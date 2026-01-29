Parameters
==========

The command-line interface is installed as ``excavate-ht``.
Run the following in Terminal to print descriptions of parameters (also described below).
.. code-block:: bash

   excavate-ht --help
   excavate-ht generate --help
   excavate-ht pair --help

EXCAVATE-HT Modes
=================

EXCAVATE-HT has two modes: ``generate`` and ``pair``.

excavate-ht generate
--------------------

The ``generate`` mode creates a complete gRNA library using inputted variant data.

Required arguments
~~~~~~~~~~~~~~~~~~

``--vcf``
   Path to the VCF file. Comma-separated if more than one
   (e.g., ``cell_line.vcf.gz,1000genomes.vcf.gz``).
   If using a phased cell-line VCF, put that first.
   *Required.*

``--var_type``
   Variant type: ``cell-line`` or ``population``.
   Can specify multiple, comma-separated
   (e.g., ``cell-line,population1,population2``).
   *Required.*

``--chrom_fa``
   Path to the chromosome FASTA file for your locus of interest
   (e.g., ``chr1sequence.fasta``).
   *Required.*

``--locus``
   Genomic locus in format ``chr#:start-end``.
   *Required.*

``-o``, ``--output_dir``
   Output directory for results. Folder name if you are already in the excavate folder.
   *Required.*

Cas protein and PAM configuration
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

You must specify either ``--cas`` OR both ``--pam-list`` and ``--orient``.

``--cas``
   Specify Cas protein. One of:
   ``SpCas9``, ``SpCas9_NG``, ``enAsCas12a``, or ``SaCas9``.

``--pam-list``
   PAM sequences (if not using a supported Cas protein).
   Comma-separated. Use IUPAC codes.
   Required if ``--cas`` is not specified.

``--orient``
   PAM orientation. Choices: ``3prime`` or ``5prime``.
   Required if ``--pam-list`` is specified.

Guide design parameters
~~~~~~~~~~~~~~~~~~~~~~~

``--af-threshold``
   Allele frequency threshold between 0 and 1. A buffer of 0.01 is subtracted from the AF threshold
   to account for rounded values in VCF files.
   *Default:* ``0.1``

``-g``, ``--guide-length``
   Guide length in base pairs.
   *Default:* ``20``

``-m``, ``--max_snppos_in_protospacer``
   Maximum distance (bp) of SNP from PAM sequence.
   *Default:* ``10``

Off-target analysis
~~~~~~~~~~~~~~~~~~~

``--off-targets``
   Enable off-target analysis: counts exact matches genome-wide and 1-bp mismatches in the chromosome of interest. Uses Bowtie1.
   *Flag (no value needed).*

``--download-hg38-index``
   Download prebuilt hg38 Bowtie indexes (GRCh38_noalt_as) into ``<outdir>/bowtie_index`` and use them automatically.
   *Flag (no value needed).*

``--genome-index-prefix``
   Bowtie1 index prefix for the genome FASTA (e.g., ``/path/to/index/hg38_bt1``).
   If missing or index files are not found, EXCAVATE-HT will build it next to this prefix.

``--genome_fa``
   Path to the whole genome FASTA file for your organism.
   Required if building Bowtie indexes from scratch (when not using ``--download-hg38-index`` or ``--genome-index-prefix``).

``--chrom-index-prefix``
   Bowtie1 index prefix for the chromosome FASTA (e.g., ``/path/to/index/chr1_bt1``).
   If missing or index files are not found, EXCAVATE-HT will build it next to this prefix.

``--bowtie-threads``
   Number of threads for Bowtie1. Defaults to ``NSLOTS`` if set (HPC environments), otherwise 4.

Output options
~~~~~~~~~~~~~~

``--split-phased``
   Enable splitting of gRNA libraries by cell-line phasing.
   *Flag (no value needed).*

``--summary``
   Output a summary table for each gRNA library.
   *Flag (no value needed).*

``--per-vcf``
   Save single-gRNA libraries for each VCF file, split by allele.
   *Flag (no value needed).*

Pairing options
~~~~~~~~~~~~~~~

``--pairing-method``
   Enable pairing of gRNA to output a dual-guide library.
   Choices:

   - ``r``: pair all guides that target different SNPs together (default).
   - ``fp``: pair guides about a fixed point.
   - ``t``: tiled pairing - pair guides that target adjacent variants.

   *Default:* ``r``

``-f``, ``--fixed-points-list``
   One or more fixed points, comma-separated
   (genomic coordinate without chr#, e.g., ``11989251,12002042,...``)
   in your locus to pair guides around.
   Required when using ``--pairing-method fp``.

excavate-ht pair
----------------

The ``pair`` mode pairs an existing single-gRNA library to create a dual-guide library.

Required arguments
~~~~~~~~~~~~~~~~~~

``-i``, ``--input-library``
   Path to input single-gRNA library file.
   *Required.*

``-o``, ``--output_dir``
   Output directory for results. Folder name if you are already in the excavate folder.
   *Required.*

Pairing options
~~~~~~~~~~~~~~~

``--pairing-method``
   Enable pairing of gRNA to output a dual-guide library.
   Choices:

   - ``r``: pair all guides that target different SNPs together (default).
   - ``fp``: pair guides about a fixed point.
   - ``t``: tiled pairing - pair guides that target adjacent variants.

   *Default:* ``r``

``-f``, ``--fixed-points-list``
   One or more fixed points, comma-separated
   (genomic coordinate without chr#, e.g., ``11989251,12002042,...``)
   in your locus to pair guides around.
   Required when using ``--pairing-method fp``.
