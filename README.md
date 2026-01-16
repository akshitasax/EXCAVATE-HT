# EXCAVATE-HT

EXtracting Common Allelic VAriants for Targeted Editing in High-Throughput

A python-based tool to mine population and personalised variant data to generate libraries of genomic loci commonly targetable for allele-specific editing.

Find detailed documentation here: https://excavate-ht.readthedocs.io/en/latest/

# Installation

## Prerequisites

Before you begin, ensure you have:
- At least 10 GB of free disk space (for miniconda, EXCAVATE-HT, and genome files)
- A stable internet connection for downloading large genome files

## Set Up

### Installing miniconda on MacOS/Unix

1. Go to the miniconda installation website: https://www.anaconda.com/download
2. Click "Skip registration"
3. Navigate to miniconda installers and download the installer for your device
4. Open the downloaded installation file and follow the install instructions. We recommend installing to your own user directory rather than for all users of the computer.
5. Open a new terminal window and verify the installation was successful by checking that `(base)` appears at the beginning of your command prompt. This indicates the base conda environment is active.

### Installing miniconda on Windows

Some packages needed to run EXCAVATE-HT (bcftools and bedtools) cannot be easily installed via conda on Windows. Therefore, you must first install WSL (Windows Subsystem for Linux) to use a Linux environment.

**Step 1: Install WSL**

1. Open PowerShell as Administrator (right-click on PowerShell > Run as Administrator)
2. Run:

   ```bash
   wsl --install
   ```

3. Restart your computer when prompted
4. After restart, open the Ubuntu app from the Start Menu (Ubuntu is automatically installed with WSL)
5. Follow the prompts to create a username and password for your Ubuntu environment

**Step 2: Install miniconda in Ubuntu**

1. In the Ubuntu terminal, download and install miniconda:

   ```bash
   wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
   bash Miniconda3-latest-Linux-x86_64.sh
   ```

2. Follow the installation prompts (press Enter to review the license, type "yes" to accept, and accept the default installation location)
3. Close and reopen the Ubuntu terminal
4. Verify installation by checking that `(base)` appears at the beginning of your command prompt

**Important note for Windows users:** To access your Windows filesystem from Ubuntu, use the `/mnt/` prefix. For example:
- Windows path: `C:\Users\YourUsername\Downloads`
- Ubuntu path: `/mnt/c/Users/YourUsername/Downloads`

To find your Windows username, open Command Prompt on Windows and run `echo %USERNAME%`.

### Create a working directory

Create a directory to store all files for your EXCAVATE-HT analysis. For example:

**On Mac/Linux:**
```bash
mkdir -p ~/Downloads/excavate-ht
cd ~/Downloads/excavate-ht
```

**On Windows (in Ubuntu terminal):**
```bash
mkdir -p /mnt/c/Users/YourUsername/Downloads/excavate-ht
cd /mnt/c/Users/YourUsername/Downloads/excavate-ht
```

Replace `YourUsername` with your actual Windows username.

### Clone the EXCAVATE-HT repository

Clone the repository to get the source code and environment configuration file:

```bash
git clone https://github.com/akshitasax/EXCAVATE-HT.git
cd EXCAVATE-HT
```

### Create and activate the excavate conda environment

The repository includes an `environment.yml` file that specifies all required dependencies. Use it to create a dedicated conda environment:

1. Create the environment (this may take several minutes):

   ```bash
   conda env create -f environment.yml
   ```

2. Activate the environment:

   ```bash
   conda activate excavate
   ```

3. Your terminal prompt should now show `(excavate)` instead of `(base)`, indicating the environment is active.

**The environment includes:**

Python dependencies:
- Python >= 3.9, < 3.13
- numpy
- pandas
- biopython
- regex
- pyfaidx

External tools:
- bedtools
- bcftools

### Install EXCAVATE-HT

With the excavate environment activated, install EXCAVATE-HT in editable mode:

```bash
pip install -e .
```

Verify the installation was successful:

```bash
excavate-ht --help
```

If you see the help message, installation was successful!

## Prepare input data files

### Set up input directory

Create a subdirectory in your working folder to store input files:

```bash
cd ~/Downloads/excavate-ht  # or your working directory
mkdir -p input_data
cd input_data
```

### Download required files

You will need the following files for your analysis:

**1. Reference genome FASTA file**

Download the human reference genome (hg38): https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000001405.26/

**2. Chromosome-specific FASTA file**

To download an individual chromosome:

1. Go to https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000001405.26/
2. Click on your chromosome of interest
3. Click either the RefSeq or GenBank blue link
4. Click "FASTA" under the page title
5. Click "Send to:" at the top right > Complete record > File > ensure format is FASTA > Create file
6. Rename the downloaded file to something descriptive (e.g., `chr1_sequence.fasta`)

**3. Cell line VCF files**

Download both the compressed VCF file and its index:
- `cell-line.vcf.gz`
- `cell-line.vcf.gz.tbi`

**4. Population VCF files**

Download population variant data (e.g., from the 1000 Genomes Project): https://www.internationalgenome.org/data-portal/data-collection

You'll need:
- `population.vcf.gz`
- `population.vcf.gz.tbi`

**Note:** Genome files can be large (several GB). Download times will vary based on your internet connection.

### Expected directory structure

After setup, your directory structure should look like this:

```text
~/Downloads/excavate-ht/
├── EXCAVATE-HT/
│   ├── environment.yml
│   ├── (other repository files)
├── input_data/
│   ├── whole_genome.fa
│   ├── chr1_sequence.fasta
│   ├── cell-line.vcf.gz
│   ├── cell-line.vcf.gz.tbi
│   ├── population.vcf.gz
│   └── population.vcf.gz.tbi
└── output_results/
```

## Using EXCAVATE-HT

EXCAVATE-HT can be run in two modes: `generate` and `pair`.

To see detailed help for each mode:

```bash
excavate-ht generate --help
```

```bash
excavate-ht pair --help
```

For parameter descriptions and usage examples, see the full documentation.

## Troubleshooting

**Issue: `(base)` doesn't appear after installing miniconda**
- Close and reopen your terminal
- If still not showing, run: `conda init` and restart your terminal

**Issue: `conda: command not found`**
- The conda installation directory may not be in your PATH. Run the installer again and ensure you answer "yes" when asked to initialize conda.

**Issue: Environment creation fails**
- Ensure you have a stable internet connection
- Try: `conda clean --all` then retry `conda env create -f environment.yml`

**Issue: Permission denied errors on Windows/WSL**
- Make sure you're working in a directory you have write access to (like `/mnt/c/Users/YourUsername/`)

**Issue: Cannot find environment.yml**
- Ensure you're in the EXCAVATE-HT directory: `cd EXCAVATE-HT`
- Verify the file exists: `ls environment.yml`

For additional help, please open an issue on the GitHub repository.