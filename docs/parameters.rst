Usage
=====

The command-line interface is installed as ``excavate-ht``.
Run the following in Terminal to print descriptions of parameters (also described below).
.. code-block:: bash

   excavate-ht --help
   excavate-ht generate --help
   excavate-ht pair --help

EXCAVATE-HT Modes
-----------------

EXCAVATE-HT has two modes, generate and pair.

excavate-ht generate
--------------------

Positional arguments
~~~~~~~~~~~~~~~~~~~~

The following arguments should be inputted in the following order:

**vcf_file**
   Path to the VCF file. Comma-separated if more than one
   (e.g., ``cell_line.vcf.gz,1000genomes.vcf.gz``).
   If using a phased cell-line VCF, put that first.

**var_type**
   Variant type: ``cell-line`` or ``population``.
   Can specify multiple, comma-separated
   (e.g., ``cell-line,population1,population2``).

**fa_file**
   Path to the chromosome FASTA file for your locus of interest
   (e.g., ``ch1sequence.fasta``).

**genome_fa**
   Path to the whole genome FASTA file for your organism.

**locus**
   Genomic locus in format ``chr#:start-end``.

Keyword arguments
~~~~~~~~~~~~~~~~~

The following arguments can be input in any order:

**--af-threshold**
   Allele frequency threshold between 0 and 1. A buffer of 0.01 is subtracted from the AF threshold 
   to account for rounded values in VCF files.
   *Default:* ``0.1``

**--cas**
   Required if no pam-list input given
   Specify Cas protein. One of:
   ``SpCas9``, ``SpCas9_NG``, ``enAsCas12a``, or ``SaCas9``.

**--pam-list**
   Required if no cas input given
   PAM sequences (if not using a supported Cas protein).
   Comma-separated. Use IUPAC codes.

**--orient**
   PAM orientation. Choices: ``3prime`` or ``5prime``.

**-g, --guide-length**
   Guide length in base pairs.
   *Default:* ``20``

**-m, --max_snppos_in_protospacer**
   Maximum distance (bp) of SNP from PAM sequence.
   *Default:* ``10``

**--off-targets**
   Enable off-target analysis:
   counts exact matches in the genome and
   1-bp mismatches in the chromosome of interest.
   *Flag (no default, only used if set).*

**--split-phased**
   Enable splitting of gRNA libraries by cell-line phasing.
   *Flag (no default, only used if set).*

**--summary**
   Output a summary table for each gRNA library.
   *Flag (no default, only used if set).*

**--per-vcf**
   Save single-gRNA libraries for each VCF file, split by allele.
   *Flag (no default, only used if set).*

**--pairing-method**
   Enable pairing of gRNA to output a dual-guide library.  
   Choices:  
   - ``r``: pair all guides that target different SNPs together (default).  
   - ``fp``: pair guides about a fixed point.  
   - ``t``: tiled pairing: pairing guides that target adjacent variants.  
   *Default:* ``r``

**-f, --fixed-points-list**
   One or more fixed points, comma-separated
   (genomic coordinate without chr#, e.g.,
   ``11989251,12002042,...``) in your locus
   to pair guides around.

**-o, --output_dir**
   Output directory for results. Folder to be created within your working directory.  
   *Required.*

excavate-ht pair
----------------

**-i, --input-library**
   Path to input single-gRNA library file.  
   *Required.*

**--pairing-method**
   Enable pairing of gRNA to output a dual-guide library.  
   Choices:  
   - ``r``: pair all guides that target different SNPs together (default).  
   - ``fp``: pair guides about a fixed point.  
   - ``t``: tiled pairing.  
   *Default:* ``r``

**-f, --fixed-points-list**
   One or more fixed points, comma-separated
   (genomic coordinate without chr#, e.g.,
   ``11989251,12002042,...``) in your locus
   to pair guides around.

**-o, --output_dir**
   Output directory for results.  
   Provide folder name if you are already in the ``excavate`` folder.  
   *Required.*
