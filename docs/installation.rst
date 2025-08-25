Installation
============

Short version:

You can install EXCAVATE-HT directly from source:

.. code-block:: bash

   git clone https://github.com/akshitasax/EXCAVATE-HT.git
   cd EXCAVATE-HT
   pip install -e .

External Tools
--------------

EXCAVATE-HT also requires external bioinformatics tools to be installed:

- bedtools
- bcftools

We recommend installing them with conda:

.. code-block:: bash

   conda install -c bioconda bedtools bcftools

Requirements
------------

- Python >= 3.9, < 3.13
- numpy
- pandas
- biopython
- regex
- pyfaidx

Long version:

Set Up
======

For MacOS/Unix
--------------

Install miniconda
~~~~~~~~~~~~~~~~

1. Go to miniconda installation website (https://www.anaconda.com/download)
2. Press skip registration
3. Go to miniconda installers - download for your device
4. Now you have an installation file, open it. Follow install instructions. I recommend installation on your own user disk and not for all users of the computer.
5. Open terminal, make sure it says (base). This means the base conda environment is active and installation was successful.
6. Continue with downloading the genome fasta and other files needed to run excavate.

For Windows
----------

Some of the packages needed to run excavate (bcftools and bedtools) cannot be easily installed via conda on Windows (non-Unix systems). Hence, it is recommended to first install and enable WSL (Windows Subsystem for Linux). This will allow you to use a Linux environment to run conda and excavate. To do this:

1. Open PowerShell as Administrator (right click on PowerShell > Run as Administrator)
2. Run:

   .. code-block:: bash

      wsl --install

3. Once it's installed, restart your computer if prompted
4. Open the Ubuntu app from Start Menu (Ubuntu is a Linux distribution downloaded when WSL was installed).

   

   .. code-block:: bash

      wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
      bash Miniconda3-latest-Linux-x86_64.sh

2. Close and open a new window for Ubuntu. If miniconda was successfully installed, you will see (base) in the beginning of your terminal.

3. **Important note:** to change directory to a desired location on your Windows filesystem via Ubuntu, the path notations are a little different. For example, to change directory to your Downloads folder (located at C:\\Users\\YourUsername\\Downloads), instead of ``cd Downloads``, you must run ``cd /mnt/c/Users/YourUsername/Downloads``. ``/mnt/...`` allows you to access your Windows filesystem.

4. Continue below:

For all users
------------

5. Create a working directory to save all files needed for your excavate-ht run. For example, you can create a folder called "excavate-ht" in your Downloads.

Create and activate the excavate environment
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Miniconda is a package installer and manager. It allows you to create different "environments" with different software packages that you may need for specific tasks. We will create one for excavate. To do this:

1. Download the environment.yml file into your excavate-ht working directory or folder
2. Open terminal on Mac or Ubuntu on Windows.
3. Change directory to the excavate-ht folder in your Downloads by typing and entering this on mac:

   .. code-block:: bash

      cd Downloads/excavate-ht

   Or this on windows:

   .. code-block:: bash

      cd /mnt/c/Users/YourUsername/Downloads/excavate-ht

4. Use the environment.yml file to create a new conda environment with all required packages and dependencies to run EXCAVATE. In terminal, enter:

   .. code-block:: bash

      conda env create -f environment.yml

5. Once all packages are downloaded successfully, try activating this environment by entering:

   .. code-block:: bash

      conda activate excavate

6. Now your terminal should say (excavate). The excavate environment is active!

Install excavate-ht
~~~~~~~~~~~~~~~~~~

You can install EXCAVATE-HT directly from source:

In terminal, within the excavate environment, type:

.. code-block:: bash

   git clone https://github.com/akshitasax/EXCAVATE-HT.git
   cd EXCAVATE-HT
   pip install -e .

EXCAVATE-HT is installed and ready to use. It will appear as a folder in your working directory.

Input data files required for your EXCAVATE-HT run
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

You can create a subfolder in your excavate-ht working directory to store input files:

* Download the FASTA file for your whole genome of interest. Find hg38 here: https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000001405.26/

* Download your chromosome of interest:

  1. Go to `https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000001405.26/ <https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000001405.26/>`_
  2. Go to chromosome of interest and click either the RefSeq or GenBank blue link
  3. Click on FASTA under the title of the page
  4. Click on Send to: at the top right corner of the page > Complete record > File > make sure format is FASTA > Create file.
  5. This will download a sequence.fasta. I suggest changing its name to identify the chromosome downloaded, for example, change to 'chr1sequence.fasta'

* Download the VCF files (.vcf.gz and .vcf.gz.tbi files) you will need for your cell-line of interest.

* Download the VCF files (.vcf.gz and .vcf.gz.tbi files) you will need for your population of interest. (For example, find files with data from the 1000 genomes project here: https://www.internationalgenome.org/data-portal/data-collection)

At this point, your working directory file structure may look like this:

.. code-block:: text

   ~/Downloads/excavate-ht/
   ├── EXCAVATE-HT/
   ├── environment.yml
   ├── input_data/
   │   ├── whole_genome.fa
   │   ├── chromosome.fa
   │   ├── cell-line.vcf.gz
   │   ├── cell-line.vcf.gz.tbi
   │   ├── population.vcf.gz
   │   └── population.vcf.gz.tbi
   └── output_results/ 

Run the following to ensure excavate-ht was installed successfully:

.. code-block:: bash

    excavate-ht --help
    
EXCAVATE-HT can be run in two modes: 'generate' and 'pair':

Run the following for a description on how to use each mode. Find descriptions of parameters and usage examples in the docs.

.. code-block:: bash

    excavate-ht generate --help

.. code-block:: bash

    excavate-ht pair --help