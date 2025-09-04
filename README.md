# EXCAVATE-HT

EXtracting Common Allelic VAriants for Targeted Editing in High-Throughput

A python-based tool to mine population and personalised variant data to generate libraries of genomic loci commonly targetable for allele-specific editing.

Find detailed documentation here: https://excavate-ht.readthedocs.io/en/latest/
**SET UP**

**For MacOS/Unix:**

***Install miniconda*** 
1. Go to miniconda installation website (https://www.anaconda.com/download)  
2. Press skip registration  
3. Go to miniconda installers \- download for your device  
4. Now you have an installation file, open it. Follow install instructions. I recommend installation on your own user disk and not for all users of the computer.  
5. Open terminal, make sure it says (base). This means the base conda environment is active and installation was successful.  
6. Continue with downloading the genome fasta and other files needed to run excavate.

**For Windows:**

Some of the packages needed to run excavate (bcftools and bedtools) cannot be easily installed via conda on Windows (non Unix systems). Hence, it is recommended to first install and enable WSL (Windows Subsystem for Linux). This will allow you to use a Linux environment to run conda and excavate. To do this,

1. Open PowerShell as Administrator (right click on PowerShell \> Run as Administrator)  
2. Run:
   
   `wsl –install`

4. Once it's installed, restart your computer if prompted  
5. Open the Ubuntu app from Start Menu (Ubuntu is a Linux distribution downloaded when WSL was installed).
   
***Install Miniconda:***  
1. In Ubuntu, run the following:

   

   `wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh`

   `bash Miniconda3-latest-Linux-x86_64.sh`

   

2. Close and open a new window for Ubuntu. If miniconda was successfully installed, you will see (base) in the beginning of your terminal.

   

3. **Important note:** to change directory to a desired location on your Windows filesystem via Ubuntu, the path notations are a little different. For example, to change directory to your Downloads folder (located at C:\\Users\\YourUsername\\Downloads), instead of `cd Downloads`, you must run `cd /mnt/c/Users/YourUsername/Downloads. /mnt/…` allows you to access your Windows filesystem.

    

4. Continue with downloading the genome fasta and other files needed to run excavate.

**For all users:**

- Download the FASTA file for your whole genome of interest. Find hg38 here: https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000001405.26/

- Download your chromosome of interest:  
	1. Go to [https://www.ncbi.nlm.nih.gov/datasets/genome/GCF\_000001405.26/](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000001405.26/)  
	2. Go to chromosome of interest and click either the RefSeq or GenBank blue link  
	3. Click on FASTA under the title of the page  
	4. Click on Send to: at the top right corner of the page \> Complete record \> File \> make sure format is FASTA \> Create file.  
	5. This will download a sequence.fasta. I suggest changing its name to identify the chromosome downloaded, for example, change to ‘chr1sequence.fasta’

- Download the VCF files (.vcf.gz and .vcf.gz.tbi files) you will need for your cell-line of interest.

- Download the VCF files (.vcf.gz and .vcf.gz.tbi files) you will need for your population of interest. (For example, find files with data from the 1000 genomes here: https://www.internationalgenome.org/data-portal/data-collection)

- Download all scripts required to run excavate.  
    
- Make a folder in Downloads called ‘excavate’
	1. Put all code files and genome, vcf, and chromosome fasta files for your run in that folder

- Before you run EXCAVATE, your excavate folder should have:  
1. **environment.yml file**  
2. **ap.py**  
3. **main.py**  
4. **Whole genome fasta file**  
5. **Chromosome fasta file**  
6. **VCF files (cell-line or population or multiple). Have both the .gz and the .gz.tbi files in your folder.**

   

- Create and activate the excavate environment: Miniconda is a package installer and manager. It allows you to create different “environments” with different software packages that you may need for specific tasks. We will create one for excavate. To do this:  
1. Open terminal on Mac or Ubuntu on Windows.  
2. Change directory to the excavate folder in your Downloads by typing and entering this on mac:

   `cd Downloads/excavate`

   Or this on windows:

   	`cd /mnt/c/Users/YourUsername/Downloads/excavate`

   

3. Use the environment.yml file to create a new conda environment with all required packages and dependencies to run EXCAVATE. In terminal, enter:

   `conda env create -f environment.yml`

	

4. Once all packages are downloaded successfully, try activating this environment by entering:

   `conda activate excavate`

5. Now your terminal should say (excavate). The excavate environment is active\!

**RUNNING EXCAVATE**

The file ap.py contains all the functions required for excavate to run. The file main.py facilitates taking user input and running all functions in the correct order. 

The first time you run EXCAVATE, you must make the main.py script executable. To do this, enter:  
`chmod +x main.py` 

Now you can run the script.

Enter this for a description of how to use the pipeline:  
`./main.py --help`

Enter this for a description of how to use the 'generate' mode of the pipeline:  
`./main.py generate --help`

Enter this for a description of how to use the 'pair' mode of the pipeline:  
`./main.py pair --help`

Here is an example of how you can use the pipeline:

./main.py generate path/to/cellline\_vcf\_file cell-line path/to/chromosome\_fasta\_file path/to/genome\_fasta\_file chr1:11975541-12019490 \--cas SpCas9 \--pair fp \-f 11976861 \--split-phased \--summary \-o cn8\_kgp\_gRNA\_files

- Script \- make sure to add ./ first\*  
- Files\*  
  - VCF File (.gz) \+ ‘cell-line’  
  - Chromosome FASTA File  
  - Whole Genome FASTA File  
- Genomic Coordinates\*  
- Cas Information\* \<- Here you will select your Cas species by either entering –cas or –pam-list, only one is MANDATORY  
- Pairing Method \<- Today we will demo fixed point pairing  
- Fixed Point \<- Provide 1 genomic position, ideally at the midpoint of your locus of interest  
- Phasing \<- Pairs guides based on phasing   
- Output Information  
  - Summary Table \<- output stats on your gRNA library  
  - Output Folder\* \<- This folder will be placed in your excavate volder, please give it a name.
