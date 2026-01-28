"""
Bowtie1-only off-target counting for EXCAVATE-HT

What this gives you:
- Exact-match counts across a Bowtie1-indexed *genome* FASTA (guide+PAM variants, both strands)
- Exactly-1-bp-mismatch-in-guide counts on a Bowtie1-indexed *chromosome* FASTA (PAM exact, both strands)
- Shared SAM parsing
- Appends two columns to your guides df:
    - exact_matches_genome
    - mm1_matches_chr   (or a custom column name)

Assumptions:
- df has column 'gRNA'
- cas_obj has fields: pam_three_prime, pam_five_prime, pam_length, is_three_prime, is_five_prime, name
- Bowtie1 is installed (binaries: bowtie, bowtie-build)
"""

from __future__ import annotations

import os
import re
import shutil
import subprocess
from collections import defaultdict
from itertools import product
from pathlib import Path
from typing import Dict, List, Tuple, Optional
import urllib.request
import zipfile

from Bio.Seq import Seq

# get chromosome name for column and file naming
def chrom_name_from_fasta(fasta_path: str) -> str:
    """
    Extract chromosome/contig name from first FASTA header.

    Example header:
        >chr1
        >chr1 dna:chromosome chromosome:GRCh38:1:1:248956422:1 REF

    Returns:
        chr1
    """
    with open(fasta_path, "r") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                # take everything after '>' up to first whitespace
                return line[1:].split()[0]

    raise ValueError(f"No FASTA header found in {fasta_path}")
# ----------------------------
# PAM expansion + target generation
# ----------------------------

def _strip_named_groups(pat: str) -> str:
    # "(?P<p0>...)" -> "..."
    return re.sub(r"\(\?P<[^>]+>([^)]+)\)", r"\1", pat)


def expand_pam_regex(pam_regex: str) -> List[str]:
    """
    Expand regex-like PAM strings into concrete PAMs.
    Supports:
      - bracket classes [ACGT], [AG], [TGC], etc.
      - literals
      - top-level | alternatives
      - named groups (?P<p0>...)
    """
    pam_regex = _strip_named_groups(pam_regex)
    alts = pam_regex.split("|")

    out = set()
    for alt in alts:
        alt = alt.strip()
        tokens: List[List[str]] = []
        i = 0
        while i < len(alt):
            if alt[i] == "[":
                j = alt.index("]", i)
                tokens.append(list(alt[i + 1 : j]))
                i = j + 1
            else:
                tokens.append([alt[i]])
                i += 1

        for tup in product(*tokens):
            out.add("".join(tup).upper())

    return sorted(out)


def _targets_exact_for_guide(guide: str, cas_obj, pam3: List[str], pam5: List[str]) -> Tuple[List[str], List[str]]:
    """
    Returns (fwd_targets, rev_targets) for *exact* matching, mirroring your old regex logic.
    """
    g = guide.upper()
    rc = str(Seq(g).reverse_complement()).upper()

    if cas_obj.is_three_prime:
        # forward: guide + PAM3
        fwd = [g + p for p in pam3]
        # omitting reverse targets because bowtie looks on rc strand anyway
        # reverse: PAM5 + rc(guide)
        # rev = [p + rc for p in pam5]
    elif cas_obj.is_five_prime:
        # forward: PAM5 + guide
        fwd = [p + g for p in pam5]
        # omitting reverse targets because bowtie looks on rc strand anyway
        # reverse: rc(guide) + PAM3
        # rev = [rc + p for p in pam3]
    else:
        raise ValueError("Cas PAM orientation not defined")

    return fwd


def one_mismatch_variants(seq: str) -> List[str]:
    """
    Exactly-1 substitution variants (Hamming distance == 1).
    For length L -> 3*L variants.
    """
    seq = seq.upper()
    bases = ("A", "C", "G", "T")
    out: List[str] = []
    for i, b in enumerate(seq):
        for nb in bases:
            if nb != b:
                out.append(seq[:i] + nb + seq[i + 1 :])
    return out


def _targets_mm1_for_guide(guide: str, cas_obj, pam3: List[str], pam5: List[str]) -> Tuple[List[str], List[str]]:
    """
    Returns (fwd_targets, rev_targets) for *exactly-1 mismatch in guide* (PAM exact),
    mirroring your old regex approach (which never included the exact guide).
    """
    g = guide.upper()
    rc = str(Seq(g).reverse_complement()).upper()

    g_vars = one_mismatch_variants(g)
    rc_vars = one_mismatch_variants(rc)

    if cas_obj.is_three_prime:
        # forward: (guide_variant) + PAM3
        fwd = [gv + p for gv in g_vars for p in pam3]
        # omitting reverse targets because bowtie looks on rc strand anyway
        # # reverse: PAM5 + (rc_guide_variant)
        # rev = [p + rv for rv in rc_vars for p in pam5]
    elif cas_obj.is_five_prime:
        # forward: PAM5 + (guide_variant)
        fwd = [p + gv for gv in g_vars for p in pam5]
        # omitting reverse targets because bowtie looks on rc strand anyway
        # # reverse: (rc_guide_variant) + PAM3
        # rev = [rv + p for rv in rc_vars for p in pam3]
    else:
        raise ValueError("Cas PAM orientation not defined")

    return fwd


def write_queries_fasta(
    guides: List[str],
    cas_obj,
    out_fa: str,
    mode: str = "exact",   # "exact" or "mm1"
) -> Tuple[Dict[str, int], int]:
    """
    Write a FASTA of queries (guide+PAM variants) for Bowtie.
    Returns:
      - qname_to_guide: maps query record name -> guide index
      - n_queries written
    """
    guides = [g.upper() for g in guides]
    if not guides:
        raise ValueError("No guides provided.")

    pam3 = expand_pam_regex(cas_obj.pam_three_prime)
    pam5 = expand_pam_regex(cas_obj.pam_five_prime)

    # sanity checks on PAM lengths
    if any(len(p) != cas_obj.pam_length for p in pam3):
        raise ValueError(f"{cas_obj.name}: pam_three_prime expansion not length {cas_obj.pam_length}")
    if any(len(p) != cas_obj.pam_length for p in pam5):
        raise ValueError(f"{cas_obj.name}: pam_five_prime expansion not length {cas_obj.pam_length}")

    out_path = Path(out_fa)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    qname_to_guide: Dict[str, int] = {}
    n = 0

    with open(out_path, "w") as f:
        for gi, g in enumerate(guides):
            if mode == "exact":
                fwd_targets = _targets_exact_for_guide(g, cas_obj, pam3, pam5)
                tag = "EX"
            elif mode == "mm1":
                fwd_targets = _targets_mm1_for_guide(g, cas_obj, pam3, pam5)
                tag = "MM1"
            else:
                raise ValueError("mode must be 'exact' or 'mm1'")

            for j, t in enumerate(fwd_targets):
                qn = f"g{gi}|F|{tag}|{j}"
                f.write(f">{qn}\n{t}\n")
                qname_to_guide[qn] = gi
                n += 1

    return qname_to_guide, n


# ----------------------------
# Bowtie1 index + run + parse
# ----------------------------

def get_bowtie1_version() -> tuple[int, int, int]:
    """
    Returns Bowtie version as (major, minor, patch).
    Robust to outputs like:
      - 'bowtie version 1.3.1'
      - 'bowtie-align-s version 1.3.1'
    """
    try:
        p = subprocess.run(
            ["bowtie", "--version"],
            capture_output=True,
            text=True,
            check=True,
        )
    except FileNotFoundError:
        raise RuntimeError("bowtie not found on PATH.")
    except subprocess.CalledProcessError as e:
        raise RuntimeError(f"bowtie --version failed:\n{e.stdout}\n{e.stderr}")

    text = (p.stdout or "") + "\n" + (p.stderr or "")

    # Match 'version 1.2.3' anywhere
    m = re.search(r"\bversion\s+(\d+)\.(\d+)\.(\d+)\b", text)
    if not m:
        raise RuntimeError(f"Could not parse Bowtie version from:\n{text}")

    return tuple(map(int, m.groups()))


def ensure_bowtie1_version(min_version=(1, 2, 3)):
    v = get_bowtie1_version()
    if v < min_version:
        raise RuntimeError(
            f"Bowtie >= {'.'.join(map(str, min_version))} required to use .bt2 indexes. "
            f"Detected version {'.'.join(map(str, v))}."
        )


HG38_BT2_URL = "https://genome-idx.s3.amazonaws.com/bt/GRCh38_noalt_as.zip"
HG38_PREFIX_NAME = "GRCh38_noalt_as"

import urllib.request
import zipfile
from pathlib import Path

HG38_BT2_URL = "https://genome-idx.s3.amazonaws.com/bt/GRCh38_noalt_as.zip"
HG38_PREFIX_NAME = "GRCh38_noalt_as"

def _bt2_forward_files(prefix: Path) -> list[Path]:
    return [Path(f"{prefix}.{i}.bt2") for i in range(1, 5)]

def download_hg38_bt2(outdir: str | Path) -> str:
    outdir = Path(outdir)
    index_dir = outdir / "bowtie_index"
    index_dir.mkdir(parents=True, exist_ok=True)

    # Preferred prefix location (no nested folder)
    preferred_prefix = index_dir / HG38_PREFIX_NAME

    # If already present at preferred location, done
    if all(p.exists() for p in _bt2_forward_files(preferred_prefix)):
        print(f"hg38 Bowtie indexes already present at: {preferred_prefix}")
        return str(preferred_prefix)

    zip_path = index_dir / f"{HG38_PREFIX_NAME}.zip"

    print("Downloading hg38 Bowtie indexes (~3–4 GB). This may take a few minutes...")
    urllib.request.urlretrieve(HG38_BT2_URL, zip_path)

    print("Extracting hg38 Bowtie indexes...")
    with zipfile.ZipFile(zip_path, "r") as z:
        z.extractall(index_dir)

    # Clean up zip to save space
    zip_path.unlink(missing_ok=True)

    # After extraction, indexes may be in a nested folder.
    # Find the location of "<name>.1.bt2" and infer the prefix from that.
    candidates = list(index_dir.rglob(f"{HG38_PREFIX_NAME}.1.bt2"))
    if not candidates:
        # Helpful debug listing
        sample = list(index_dir.rglob("*.bt2"))[:20]
        raise RuntimeError(
            "hg38 download completed but expected .bt2 files not found.\n"
            f"Searched under: {index_dir}\n"
            f"Found bt2 files (sample): {[str(p) for p in sample]}"
        )

    # Pick the first match and compute its prefix path by stripping ".1.bt2"
    one_bt2 = candidates[0]
    detected_prefix = Path(str(one_bt2).replace(".1.bt2", ""))

    # Validate that the full forward set exists for the detected prefix
    if not all(p.exists() for p in _bt2_forward_files(detected_prefix)):
        raise RuntimeError(
            "Found an hg38 bt2 index file but not the full expected set of forward indexes.\n"
            f"Detected prefix: {detected_prefix}\n"
            f"Missing: {[str(p) for p in _bt2_forward_files(detected_prefix) if not p.exists()]}"
        )

    print(f"hg38 Bowtie indexes ready at: {detected_prefix}")
    return str(detected_prefix)


def bowtie_index_exists(prefix: str) -> bool:
    """
    Returns True if either Bowtie1 (.ebwt) or Bowtie2 (.bt2) index exists.
    Requires forward indexes only; reverse indexes optional for exact matching.
    """
    p = Path(prefix)

    ebwt = [
        f"{p}.1.ebwt",
        f"{p}.2.ebwt",
        f"{p}.3.ebwt",
        f"{p}.4.ebwt",
    ]

    bt2 = [
        f"{p}.1.bt2",
        f"{p}.2.bt2",
        f"{p}.3.bt2",
        f"{p}.4.bt2",
    ]

    return all(Path(x).exists() for x in ebwt) or all(Path(x).exists() for x in bt2)


def ensure_bowtie1_index(fasta_path: str, index_prefix: str) -> None:
    """
    Ensure Bowtie index exists.

    - Accepts either .ebwt or .bt2.
    - If missing, builds Bowtie1 (.ebwt) index.
    - Requires Bowtie >=1.2.3 (so .bt2 can also be used if provided).
    """
    # Make sure bowtie itself is modern enough
    ensure_bowtie1_version()

    # If user already provided either ebwt or bt2, we're done
    if bowtie_index_exists(index_prefix):
        return

    # Otherwise build ebwt
    if shutil.which("bowtie-build") is None:
        raise RuntimeError("bowtie-build not found. Install Bowtie1 or load the module.")

    Path(index_prefix).parent.mkdir(parents=True, exist_ok=True)

    subprocess.run(
        ["bowtie-build", fasta_path, index_prefix],
        check=True,
    )


def run_bowtie1(
    index_prefix: str,
    queries_fa: str,
    sam_out: str,
    threads: int = 1,
    max_alignments: Optional[int] = None,  # if you want to cap (None = all)
) -> None:
    """
    Bowtie1 run for exact matching of provided query sequences.
    We always use -v 0 because:
      - for exact mode, queries are exact targets
      - for mm1 mode, queries are explicitly enumerated 1-mismatch variants
    """
    if shutil.which("bowtie") is None:
        raise RuntimeError("bowtie (Bowtie1) not found on PATH. Install Bowtie1 or load the module.")

    cmd = [
        "bowtie",
        "-x", index_prefix,
        "-v", "0",
        "-a",
        "--best", "--strata",
        "-p", str(max(1, int(threads))),
        "-f",
        "-S",
        queries_fa,
        sam_out,
    ]

    if max_alignments is not None:
        # Bowtie1: -m suppresses reads with >m alignments; not what you want usually.
        # So we DON'T use -m here. If you truly want a cap, you'll need a different strategy.
        raise ValueError("Bowtie1 does not support a clean 'cap but still report' mode with -a. Leave max_alignments=None.")

    print("Running bowtie command:", " ".join(cmd))
    p = subprocess.run(cmd, capture_output=True, text=True)
    if p.returncode != 0:
        raise RuntimeError(
            "Bowtie failed.\n"
            f"CMD: {' '.join(cmd)}\n"
            f"STDOUT:\n{p.stdout}\n"
            f"STDERR:\n{p.stderr}\n"
        )


def parse_sam_qname_counts(sam_path: str) -> Dict[str, int]:
    """
    Returns {QNAME: number_of_alignments} for mapped alignments in SAM.
    """
    counts: Dict[str, int] = defaultdict(int)
    with open(sam_path, "r") as f:
        for line in f:
            if not line or line.startswith("@"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 2:
                continue
            qname = fields[0]
            flag = int(fields[1])
            if flag & 0x4:  # unmapped
                continue
            counts[qname] += 1
    return dict(counts)


def sum_counts_per_guide(qname_to_guide: Dict[str, int], q_counts: Dict[str, int], n_guides: int) -> List[int]:
    out = [0] * n_guides
    for qname, hits in q_counts.items():
        out[qname_to_guide[qname]] += hits
    return out


# ----------------------------
# High-level pipeline function (what you asked for)
# ----------------------------

def add_bowtie_offtargets(
    df,
    cas_obj,
    genome_index_prefix: str,
    chr_index_prefix: str,
    *,
    # Only needed if you want auto-build; safe to pass None and build elsewhere.
    genome_fasta_for_autobuild: Optional[str] = None,
    chr_fasta_for_autobuild: Optional[str] = None,
    ensure_index: bool = False,
    threads_genome: int = 4,
    threads_chr: int = 2,
    tmp_dir: str = "tmp_offtargets",
    exact_col: str = "exact_matches_in_genome",
    mm1_col: str = "1bp_mismatch_matches_in_chr",
):
    """
    Appends:
      - exact_col: exact matches across genome (Bowtie1)
      - mm1_col: exactly-1-bp mismatch-in-guide matches on chromosome (Bowtie1)

    Requires df['gRNA'].
    """
    import pandas as pd  # keep module import lighter

    if "gRNA" not in df.columns:
        raise ValueError("df must have a 'gRNA' column")

    guides = [g.upper() for g in df["gRNA"].tolist()]
    n_guides = len(guides)

    tmp = Path(tmp_dir)
    tmp.mkdir(parents=True, exist_ok=True)

    # Optionally ensure indexes exist
    if ensure_index:
        if genome_fasta_for_autobuild is None or chr_fasta_for_autobuild is None:
            raise ValueError("If ensure_index=True, provide both genome_fasta_for_autobuild and chr_fasta_for_autobuild.")
        ensure_bowtie1_index(genome_fasta_for_autobuild, genome_index_prefix)
        ensure_bowtie1_index(chr_fasta_for_autobuild, chr_index_prefix)

    # ---- Exact genome-wide ----
    exact_fa = tmp / f"queries_exact_{cas_obj.name}.fa"
    exact_sam = tmp / f"hits_exact_{cas_obj.name}.sam"

    qmap_exact, n_exact = write_queries_fasta(guides, cas_obj, str(exact_fa), mode="exact")
    run_bowtie1(genome_index_prefix, str(exact_fa), str(exact_sam), threads=threads_genome)
    qcounts_exact = parse_sam_qname_counts(str(exact_sam))
    exact_per_guide = sum_counts_per_guide(qmap_exact, qcounts_exact, n_guides)

    # ---- 1-mismatch (guide-only) on chromosome ----
    mm1_fa = tmp / f"queries_mm1_{cas_obj.name}.fa"
    mm1_sam = tmp / f"hits_mm1_{cas_obj.name}.sam"

    qmap_mm1, n_mm1 = write_queries_fasta(guides, cas_obj, str(mm1_fa), mode="mm1")
    run_bowtie1(chr_index_prefix, str(mm1_fa), str(mm1_sam), threads=threads_chr)
    qcounts_mm1 = parse_sam_qname_counts(str(mm1_sam))
    mm1_per_guide = sum_counts_per_guide(qmap_mm1, qcounts_mm1, n_guides)

    out = df.copy()
    out[exact_col] = exact_per_guide
    out[mm1_col] = mm1_per_guide

    # (optional) attach debug metadata columns if you want
    # out["_n_exact_queries"] = n_exact
    # out["_n_mm1_queries"] = n_mm1

    return out


# ----------------------------
# Example usage
# ----------------------------
"""
# You already built indexes:
#   bowtie-build hg38.fa /path/to/index/hg38_bt1
#   bowtie-build ch1sequence.fasta /path/to/index/ch1_bt1

df_out = add_bowtie_offtargets(
    df_guides,
    SpCas9,
    genome_index_prefix="/path/to/index/hg38_bt1",
    chr_index_prefix="/path/to/index/ch1_bt1",
    ensure_index=False,                 # you said you already built
    threads_genome=8,
    threads_chr=4,
    tmp_dir="tmp_offtargets",
    exact_col="exact_matches_genome",
    mm1_col="mm1_matches_chr1",
)

# df_out now has both columns appended.
"""
