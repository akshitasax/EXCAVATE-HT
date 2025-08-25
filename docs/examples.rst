EXCAVATE-HT Examples
--------------------

Example run: Generating a paired dual-gRNA library 
==================================================

Generating paired dual-gRNAs targeting heterozygous SNVs in an individual cell-line in the MFN2 locus, for the SpCas9 enzyme, 
paired via "random" method (all gRNA targeting different SNPs paired together) in the locus. 

.. code-block:: bash

   excavate-ht generate \
     cellline.vcf.gz \
     cell-line \
     chrom1.fasta \
     hg38.fasta \
     chr1:11975541-12019490 \
     --cas SpCas9 \
     --pairing-method r \
     --off-targets \
     --summary \
     -o paired_gRNA_output_folder


Example run: Generating a single-gRNA library, editing it, then pairing to create a dual-gRNA library 
=====================================================================================================
Generaing single-gRNAs targeting heteroyzgous SNVs in an individual cell-line and common SNPs present
in the 1000 genomes project, in the MFN2 locus, for a custom 5prime PAM sequence.

.. code-block:: bash

   excavate-ht generate \
     cellline.vcf.gz,1000genomes_chrom1.vcf.gz cell-line,population \
     chrom1.fasta \
     hg38.fasta \
     chr1:11975541-12019490 \
     --pam-list TRTV \
     --orient 5prime \
     --off-targets \
     --summary \
     -o single_gRNA_output_folder

This sgRNA library can now be manually edited. gRNA with many off-target sites, as revealed by the --off-targets analysis option, 
can be filtered out, and screening controls such as non-targeting gRNA or experimentally verified gRNA can be added in.

After you have the final library of all single-gRNA you want included in your library/screen, you can pair up this user-input 
sgRNA library to create a library of paired dual-gRNA elements, using excavate-ht pair. We will pair up gRNA about a fixed-point in the MFN2 locus.  

.. code-block:: bash

   excavate-ht pair \
     --input-library final_sgRNA_MFN2_library.csv
     --pairing-method fp\
     --fixed-points-list 11976861 \
     -o paired_MFN2_library
