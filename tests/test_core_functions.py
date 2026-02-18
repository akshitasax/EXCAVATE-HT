"""
Unit tests for core EXCAVATE-HT functions using toy data.

Toy chromosome (300 nt, all A's except two deliberately placed PAMs):
  - positions 53-55  (1-based): "CGG"    → SpCas9 PAM (NGG) near SNP at pos 50
  - positions 153-158 (1-based): "AAGAGT" → SaCas9 PAM (NNGRRT) near SNP at pos 150

For SpCas9 (pam_len=3, guide_len=20, SNP at pos 50):
  guide occupies positions 33-52, PAM at 53-55
  SNP is at guide position 18 from the 5' end (3 nt from PAM)

For SaCas9 (pam_len=6, guide_len=20, SNP at pos 150):
  guide occupies positions 133-152, PAM at 153-158
  SNP is at guide position 18 from the 5' end (3 nt from PAM)

Neither PAM creates spurious guides for the other Cas species (verified by
inspecting the snpregion sequence and PAM regex patterns).
"""

import pandas as pd
import pytest
from Bio.Seq import Seq

from excavate import ap

# ─────────────────────────────────────────────────────────────────────────────
# Toy chromosome (module-level constant)
# ─────────────────────────────────────────────────────────────────────────────

def _build_chrom_seq() -> Seq:
    """300-nt all-A sequence with two deliberately placed PAMs."""
    seq = list("A" * 300)
    seq[52:55] = list("CGG")      # SpCas9 PAM at 1-based positions 53-55
    seq[152:158] = list("AAGAGT") # SaCas9 PAM at 1-based positions 153-158
    return Seq("".join(seq))

CHROM_SEQ = _build_chrom_seq()

# ─────────────────────────────────────────────────────────────────────────────
# Helper: build a gens DataFrame as create_gens() would produce
# ─────────────────────────────────────────────────────────────────────────────

def make_gens_df(pos: int, ref: str, alt: str, af: float, rsid: str = "rs_test") -> pd.DataFrame:
    """Single-row population-style gens DataFrame (no bcftools needed)."""
    return pd.DataFrame({
        "chrom":           ["chr1"],
        "pos":             [pos],
        "snp_id":          [rsid],
        "ref":             [ref],
        "alt":             [alt],
        "alt AF":          [af],
        "present alleles": [(f"{ref}(ref)", f"{alt}(alt)")],
    })


# ─────────────────────────────────────────────────────────────────────────────
# getaltseq
# ─────────────────────────────────────────────────────────────────────────────

def test_getaltseq_allele1_keeps_ref():
    """allele1 retains the reference base at the SNP position."""
    gens = make_gens_df(50, "A", "T", 0.3, "rs1")
    allele1 = ap.getaltseq(gens, CHROM_SEQ, "allele1")
    assert allele1[49] == "A"  # 0-based index 49 = 1-based position 50


def test_getaltseq_allele2_carries_alt():
    """allele2 substitutes the alt base at the SNP position."""
    gens = make_gens_df(50, "A", "T", 0.3, "rs1")
    allele2 = ap.getaltseq(gens, CHROM_SEQ, "allele2")
    assert allele2[49] == "T"


def test_getaltseq_does_not_mutate_elsewhere():
    """Only the SNP position changes; the rest of the sequence is unchanged."""
    gens = make_gens_df(50, "A", "T", 0.3, "rs1")
    allele2 = ap.getaltseq(gens, CHROM_SEQ, "allele2")
    # PAM at positions 53-55 (0-based 52-54) must be untouched
    assert str(allele2[52:55]) == "CGG"


# ─────────────────────────────────────────────────────────────────────────────
# SpCas9 guide-finding (SNP at pos 50, PAM CGG at pos 53-55)
# ─────────────────────────────────────────────────────────────────────────────

def test_spcas9_finds_exactly_one_guide_allele1():
    """SpCas9 finds exactly one + strand guide for the ref allele."""
    gens = make_gens_df(50, "A", "T", 0.3, "rs1")
    allele1 = ap.getaltseq(gens, CHROM_SEQ, "allele1")
    guides = ap.find_guides(gens, allele1, ap.SpCas9, 10, 20)
    assert len(guides) == 1


def test_spcas9_finds_exactly_one_guide_allele2():
    """SpCas9 finds exactly one + strand guide for the alt allele."""
    gens = make_gens_df(50, "A", "T", 0.3, "rs1")
    allele2 = ap.getaltseq(gens, CHROM_SEQ, "allele2")
    guides = ap.find_guides(gens, allele2, ap.SpCas9, 10, 20)
    assert len(guides) == 1


def test_spcas9_guide_sequence_allele1():
    """Ref-allele SpCas9 guide is 20 A's (all ref at the SNP)."""
    gens = make_gens_df(50, "A", "T", 0.3, "rs1")
    allele1 = ap.getaltseq(gens, CHROM_SEQ, "allele1")
    guides = ap.find_guides(gens, allele1, ap.SpCas9, 10, 20)
    assert guides["gRNA"].iloc[0] == "A" * 20


def test_spcas9_guide_sequence_allele2():
    """Alt-allele SpCas9 guide has T at position 18 from the 5' end."""
    gens = make_gens_df(50, "A", "T", 0.3, "rs1")
    allele2 = ap.getaltseq(gens, CHROM_SEQ, "allele2")
    guides = ap.find_guides(gens, allele2, ap.SpCas9, 10, 20)
    # Guide positions 1-17: A, position 18: T (SNP), positions 19-20: A
    expected = "A" * 17 + "T" + "A" * 2
    assert guides["gRNA"].iloc[0] == expected


def test_spcas9_guide_metadata():
    """SpCas9 guide has correct PAM, strand, genomic coordinates, and IDs."""
    gens = make_gens_df(50, "A", "T", 0.3, "rs1")
    allele1 = ap.getaltseq(gens, CHROM_SEQ, "allele1")
    guides = ap.find_guides(gens, allele1, ap.SpCas9, 10, 20)
    row = guides.iloc[0]

    assert row["PAM"]          == "CGG"
    assert row["strand"]       == "+"
    assert row["start"]        == 33       # 1-based guide start
    assert row["end"]          == 52       # 1-based guide end (inclusive)
    assert row["chrom"]        == "chr1"
    assert row["SNP position"] == 50
    assert row["present allele"] == "A (ref)"
    assert row["variation"]    == "A>T"
    assert row["guide ID"]     == "chr1_33_+_A_SpCas9_20nt"


def test_spcas9_no_guides_near_sacas9_pam():
    """SpCas9 finds no guide near the SaCas9 PAM at pos 153-158."""
    gens = make_gens_df(150, "A", "G", 0.25, "rs2")
    allele1 = ap.getaltseq(gens, CHROM_SEQ, "allele1")
    guides = ap.find_guides(gens, allele1, ap.SpCas9, 10, 20)
    assert len(guides) == 0


# ─────────────────────────────────────────────────────────────────────────────
# SaCas9 guide-finding (SNP at pos 150, PAM AAGAGT at pos 153-158)
# ─────────────────────────────────────────────────────────────────────────────

def test_sacas9_finds_exactly_one_guide():
    """SaCas9 finds exactly one + strand guide for the ref allele."""
    gens = make_gens_df(150, "A", "G", 0.25, "rs2")
    allele1 = ap.getaltseq(gens, CHROM_SEQ, "allele1")
    guides = ap.find_guides(gens, allele1, ap.SaCas9, 10, 20)
    assert len(guides) == 1


def test_sacas9_guide_metadata():
    """SaCas9 guide has correct PAM, strand, genomic coordinates, and IDs."""
    gens = make_gens_df(150, "A", "G", 0.25, "rs2")
    allele1 = ap.getaltseq(gens, CHROM_SEQ, "allele1")
    guides = ap.find_guides(gens, allele1, ap.SaCas9, 10, 20)
    row = guides.iloc[0]

    assert row["gRNA"]         == "A" * 20
    assert row["PAM"]          == "AAGAGT"
    assert row["strand"]       == "+"
    assert row["start"]        == 133
    assert row["end"]          == 152
    assert row["chrom"]        == "chr1"
    assert row["SNP position"] == 150
    assert row["present allele"] == "A (ref)"
    assert row["variation"]    == "A>G"
    assert row["guide ID"]     == "chr1_133_+_A_SaCas9_20nt"


def test_sacas9_alt_allele_guide_sequence():
    """Alt-allele SaCas9 guide has G at position 18 from the 5' end."""
    gens = make_gens_df(150, "A", "G", 0.25, "rs2")
    allele2 = ap.getaltseq(gens, CHROM_SEQ, "allele2")
    guides = ap.find_guides(gens, allele2, ap.SaCas9, 10, 20)
    expected = "A" * 17 + "G" + "A" * 2
    assert guides["gRNA"].iloc[0] == expected


def test_sacas9_no_guides_near_spcas9_pam():
    """SaCas9 finds no guide near the SpCas9 PAM at pos 53-55."""
    gens = make_gens_df(50, "A", "T", 0.3, "rs1")
    allele1 = ap.getaltseq(gens, CHROM_SEQ, "allele1")
    guides = ap.find_guides(gens, allele1, ap.SaCas9, 10, 20)
    assert len(guides) == 0


# ─────────────────────────────────────────────────────────────────────────────
# BED format output
# ─────────────────────────────────────────────────────────────────────────────

def _spcas9_allele1_guides() -> pd.DataFrame:
    gens = make_gens_df(50, "A", "T", 0.3, "rs1")
    allele1 = ap.getaltseq(gens, CHROM_SEQ, "allele1")
    return ap.find_guides(gens, allele1, ap.SpCas9, 10, 20)


def test_bed_format_columns():
    """output_bed_format returns the six standard BED columns."""
    bed = ap.output_bed_format(_spcas9_allele1_guides())
    assert list(bed.columns) == ["chrom", "start", "end", "name", "score", "strand"]


def test_bed_format_start_is_zero_based():
    """BED start is 0-based (guide 1-based start minus 1)."""
    bed = ap.output_bed_format(_spcas9_allele1_guides())
    assert bed["start"].iloc[0] == 32  # 33 - 1


def test_bed_format_end_is_1based_inclusive():
    """BED end is the 1-based inclusive guide end (unchanged from CSV)."""
    bed = ap.output_bed_format(_spcas9_allele1_guides())
    assert bed["end"].iloc[0] == 52


def test_bed_format_score_and_name():
    """BED score is 500 and name matches guide ID."""
    bed = ap.output_bed_format(_spcas9_allele1_guides())
    assert bed["score"].iloc[0] == 500
    assert bed["name"].iloc[0]  == "chr1_33_+_A_SpCas9_20nt"


def test_bed_written_without_header(tmp_path):
    """BED file written to disk has no header — first line starts with chrom."""
    bed = ap.output_bed_format(_spcas9_allele1_guides())
    out = tmp_path / "guides.bed"
    bed.to_csv(out, sep="\t", header=False, index=False)
    first_line = out.read_text().splitlines()[0]
    # Must start with the chromosome field, not a column name
    assert first_line.startswith("chr1\t")
    assert "chrom" not in first_line


def test_bed_written_tab_separated(tmp_path):
    """BED file is tab-separated (not comma-separated)."""
    bed = ap.output_bed_format(_spcas9_allele1_guides())
    out = tmp_path / "guides.bed"
    bed.to_csv(out, sep="\t", header=False, index=False)
    first_line = out.read_text().splitlines()[0]
    fields = first_line.split("\t")
    assert len(fields) == 6


# ─────────────────────────────────────────────────────────────────────────────
# all_guides_var_info — variant annotation
# ─────────────────────────────────────────────────────────────────────────────

def test_var_info_in_vcf_column_present():
    """Guides whose SNP is in the gens_dict get 'Y' in the In-VCF column."""
    gens = make_gens_df(50, "A", "T", 0.3, "rs1")
    allele1 = ap.getaltseq(gens, CHROM_SEQ, "allele1")
    allele2 = ap.getaltseq(gens, CHROM_SEQ, "allele2")
    g1 = ap.find_guides(gens, allele1, ap.SpCas9, 10, 20)
    g2 = ap.find_guides(gens, allele2, ap.SpCas9, 10, 20)
    all_guides = pd.concat([g1, g2], ignore_index=True)

    gens_dict = {"gens_df_population": gens}
    annotated = ap.all_guides_var_info(gens_dict, all_guides)

    assert "In population" in annotated.columns
    assert (annotated["In population"] == "Y").all()


def test_var_info_in_vcf_column_absent():
    """Guides whose SNP is NOT in a gens_df get 'N' in that column."""
    gens_snp1 = make_gens_df(50,  "A", "T", 0.3,  "rs1")
    gens_snp2 = make_gens_df(150, "A", "G", 0.25, "rs2")

    allele1_snp2 = ap.getaltseq(gens_snp2, CHROM_SEQ, "allele1")
    guides_snp2  = ap.find_guides(gens_snp2, allele1_snp2, ap.SaCas9, 10, 20)

    # Annotate snp2 guides using only the snp1 gens_dict → SNP position won't match
    gens_dict = {"gens_df_population": gens_snp1}
    annotated = ap.all_guides_var_info(gens_dict, guides_snp2)

    assert (annotated["In population"] == "N").all()


def test_var_info_alt_af_populated():
    """alt allele frequency column is filled from the population gens_df."""
    gens = make_gens_df(50, "A", "T", 0.3, "rs1")
    allele1 = ap.getaltseq(gens, CHROM_SEQ, "allele1")
    guides = ap.find_guides(gens, allele1, ap.SpCas9, 10, 20)

    gens_dict = {"gens_df_population": gens}
    annotated = ap.all_guides_var_info(gens_dict, guides)

    assert annotated["alt allele frequency"].iloc[0] == "0.3"


# ─────────────────────────────────────────────────────────────────────────────
# targetable_vars — summary table
# ─────────────────────────────────────────────────────────────────────────────

def _make_annotated_guides() -> pd.DataFrame:
    """Return annotated SpCas9 guides for SNP at pos 50 (both alleles)."""
    gens = make_gens_df(50, "A", "T", 0.3, "rs1")
    allele1 = ap.getaltseq(gens, CHROM_SEQ, "allele1")
    allele2 = ap.getaltseq(gens, CHROM_SEQ, "allele2")
    g1 = ap.find_guides(gens, allele1, ap.SpCas9, 10, 20)
    g2 = ap.find_guides(gens, allele2, ap.SpCas9, 10, 20)
    all_guides = pd.concat([g1, g2], ignore_index=True)
    all_guides["alt allele frequency"] = "0.3"
    return all_guides


def test_targetable_vars_columns():
    """Summary table has the four required columns."""
    guides = _make_annotated_guides()
    snplist = pd.DataFrame({
        "rsID":                ["rs1"],
        "SNP position":        [50],
        "alt allele frequency": [0.3],
    })
    summary = ap.targetable_vars(guides, snplist)
    count_col = "no. of guides found (with ref or alt allele)"
    expected = {"rsID", "SNP position", "alt allele frequency", count_col}
    assert expected.issubset(set(summary.columns))


def test_targetable_vars_snp_with_guides():
    """SNP at position 50 has two guides (one per allele)."""
    guides = _make_annotated_guides()
    snplist = pd.DataFrame({
        "rsID":                ["rs1"],
        "SNP position":        [50],
        "alt allele frequency": [0.3],
    })
    summary  = ap.targetable_vars(guides, snplist)
    count_col = "no. of guides found (with ref or alt allele)"
    count = summary[summary["SNP position"] == 50][count_col].iloc[0]
    assert count == 2  # one guide per allele (ref + alt)


def test_targetable_vars_snp_without_guides():
    """SNP at position 80 (no PAM nearby) reports 0 guides."""
    guides = _make_annotated_guides()
    snplist = pd.DataFrame({
        "rsID":                ["rs1",  "rs_none"],
        "SNP position":        [50,     80],
        "alt allele frequency": [0.3,   0.2],
    })
    summary   = ap.targetable_vars(guides, snplist)
    count_col = "no. of guides found (with ref or alt allele)"
    count = summary[summary["SNP position"] == 80][count_col].iloc[0]
    assert count == 0
