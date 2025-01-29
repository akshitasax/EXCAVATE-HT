# EXCAVATE

EXtracting Common Allelic VAriants for Targeted Editing

A python-based tool to mine population and personalised variant data to generate libraries of genomic loci commonly targetable for allele-specific editing.

## Getting Started

Download:

1. **Miniconda** (https://www.anaconda.com/download/success)

- Create new conda environment with all dependencies needed to run Excavate:

    - Download ap_environment.yml
    - Run the following on Terminal:
cd <path/to/ap_environment.yml>
conda env create -f ap_environment.yml
- Verify environment was installed successfully by running:
conda env list
2. **ap.py**
3. **gen_dual_guide_lib.py**

- keep ap.py and gen_dual_guide_lib.py in the same directory

4. **Files needed to create your guide libraries**
  - Variant data (vcf.gz and vcf.gz.tbi files) for chromosome or locus of interest from a population database of your choice
  - Variant data (vcf file) for your cell-line of interest
  - Fasta file for the reference sequence of the chromosome your gene of interest falls in
  - Whole genome fasta files for off-target analysis

## Usage

gen_dual_guide_lib.py vcf_file fa_file locus var_type af_threshold out

    <vcf_file>            path to population or cell-line vcf, comma-separated if more than one (WTC.vcf, 1000genomes.vcf)
    <fa_file>             path-to-chromosome-fasta
    <hg38>                path-to-whole-genome-fasta-files
    <locus>               chr#:start-end
    <var_type>            'cell-line' or 'population', comma-separated if more than one, corresponding to each vcf file ('cell-line', 'population')
    <af_threshold>        floating point allele frequency threshold
    <out>                 directory in which to save outputted files
    <cas_types>           Cas species, comma-separated for multiple #coming soon
    <guide_length>        comma-separated guide lengths if needed for different cas enzymes #coming soon

example code:

./gen_dual_guide_lib.py Downloads/KOLF2-ARID2-A02.vcf.gz,Downloads/ALL.chr16.shapeit2_integrated_snvindels_v2a_27022019.GRCh38.phased.vcf.gz Downloads/ch16sequence.fasta chr16:31150068-31214118,16:31150068-31214118 cell-line,population 0.3 Downloads/output_folder
