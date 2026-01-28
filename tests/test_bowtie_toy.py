import tempfile
from pathlib import Path
import subprocess
import pandas as pd
from excavate import ap

# Adjust import to wherever your offtargets module lives
from excavate import offtargets as ot


def write_toy_genome(tmpdir: Path):
    genome = """>chrTest
CTGTCCCCAGGTGGGGGTGCAGG
"""
    genome_fa = tmpdir / "genome.fa"
    genome_fa.write_text(genome)
    return genome_fa


def test_exact_matches():
    with tempfile.TemporaryDirectory() as d:
        d = Path(d)

        # 1. Toy genome
        genome_fa = write_toy_genome(d)

        # 2. Two guides: one present, one absent
        guides = pd.DataFrame(
            {
                "gRNA": [
                    "CTGTCCCCAGGTGGGGGTGC",  # exact match once
                    "TTTTTTTTTTTTTTTTTTTT",  # not present
                ]
            }
        )

        # 3. Index prefix
        index_prefix = str(d / "toy_bt")

        # 4. Dummy cas_obj (adjust if your function needs real fields)

        cas_obj = ap.SpCas9

        # 5. Run exact matching only
        out = ot.add_bowtie_offtargets(
            guides,
            cas_obj,
            genome_index_prefix=index_prefix,
            chr_index_prefix=index_prefix,
            ensure_index=True,
            genome_fasta_for_autobuild=str(genome_fa),
            chr_fasta_for_autobuild=str(genome_fa),
            threads_genome=1,
            threads_chr=1,
            tmp_dir=str(d / "tmp"),
            exact_col="exact",
            mm1_col="mm1",
        )

        print(out)

        # 6. Assertions
        assert out.loc[0, "exact"] == 1
        assert out.loc[1, "exact"] == 0

        print("\n✅ Exact matching test passed!")


if __name__ == "__main__":
    test_exact_matches()
