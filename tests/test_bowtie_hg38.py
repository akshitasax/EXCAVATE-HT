import os
import tempfile
from pathlib import Path
import pandas as pd

from excavate import offtargets as ot
from excavate import ap


def test_hg38_exact_two_guides():
    """
    Integration test against real hg38 Bowtie index.

    Requires env var:
        HG38_BT2_PREFIX=/full/path/to/hg38

    This should point to files like:
        hg38.1.bt2, hg38.2.bt2, ...
    """

    hg38_prefix='/Users/akshitasaxena/Downloads/ncbi_dataset_hg38/ncbi_dataset/data/GCF_000001405.26/GRCh38_noalt_as/GRCh38_noalt_as'
    ch1_prefix='/Users/akshitasaxena/Downloads/excavate/excavate_paper_input/012626_outputs/ch1_hg38_bt1'
    # Two guides:
    # 1) poly-A 20mer: should have MANY exact matches in hg38
    # 2) unlikely random 20mer: very likely zero matches
    guides = pd.DataFrame(
        {
            "gRNA": [
                "AAAAAAAAAAAAAAAAAAAA",  # should hit lots of places
                "GCTTAGCGTACGATCGATGC", # probably zero
                "TTGTATTTTTAGTAGATATG",  # 731
                "CAGTTGGCCAAGGGACAATG" # 1
            ]
        }
    )


    cas_obj = ap.SpCas9

    with tempfile.TemporaryDirectory() as td:
        td = Path(td)

        out = ot.add_bowtie_offtargets(
            guides,
            cas_obj,
            genome_index_prefix=hg38_prefix,
            chr_index_prefix=ch1_prefix,   # unused here but required by API
            ensure_index=False,             # hg38 already indexed
            genome_fasta_for_autobuild=None,
            chr_fasta_for_autobuild=None,
            threads_genome=2,
            threads_chr=2,
            tmp_dir=str(td / "tmp"),
            exact_col="exact",
            mm1_col="mm1",
        )

        print(out)

        # Poly-A should have >= 1 exact hit
        assert out.loc[0, "exact"] > 0

        # Random guide should usually be zero
        assert out.loc[1, "exact"] == 0

        # mfn2 guide mm1 should be exactly 7906
        assert out.loc[2, "mm1"] == 7906

        # mfn2 guide2 exact should be exactly 1
        assert out.loc[3, "exact"] == 1

        print("\n✅ hg38 exact-match test passed!")


if __name__ == "__main__":
    test_hg38_exact_two_guides()
