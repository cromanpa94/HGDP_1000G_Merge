# Merging HGDP and 1000G

This repo will take you through the steps to merge the HGDP and 1000G reference files into a single plink binary file.
Once the repo has been downloaded make sure that you meet all of the requirements, and download the necessary files in their respective folders.
To do that, first go to the `DataBases` folder, and read the README files indicating what needs to be downloaded on each (both the `HGDP` and the `1000G` folders).
After everything has been downloaded, you can start running the [script](https://nbviewer.jupyter.org/github/tomszar/HGDP_1000G_Merge/blob/master/Code/2018-05-MergeGenotypes.ipynb), located in the `Code` folder. 
This is a python notebook, so you can interactively run it, and modify it to your needs.
In summary the script will follow these steps:
- Transform the HGDP into plink files
- LifOver the HGDP from hg18 to hg19
- Extract only the SNPs found in the HGDP from the 1000G vcf files
- Concatenate the different chromosomes and export to plink files
- Merge the HGDP and 1000G

## Requirements

This script was ran on a Linux machine, using Ubuntu 18.04.
You will need the following programs:
- [Python 3.x](https://www.python.org/downloads/): I recommend installing python 3.x using [Anaconda](https://www.anaconda.com/download/). 

For the following programs, you can use the [bioconda](https://bioconda.github.io/) channel to install them through Anaconda.
To do that, once you've installed Anaconda follow the instructions in [here](https://bioconda.github.io/).
The script will assume that all of the following programs are in your path.
- [Plink](https://www.cog-genomics.org/plink2): to install it using [bioconda](https://bioconda.github.io/recipes/plink/README.html) use the following command `conda install plink`
- [Vcftools](https://vcftools.github.io/index.html): to install it using [bioconda](https://bioconda.github.io/recipes/vcftools/README.html) use the following command `conda install vcftools`
- [Bcftools](https://samtools.github.io/bcftools/bcftools.html): to install it using [bioconda](https://bioconda.github.io/recipes/bcftools/README.html) use the following command `conda install bcftools`
- [USCS liftOver](http://genome.ucsc.edu/cgi-bin/hgLiftOver): to install it using [bioconda](https://bioconda.github.io/recipes/ucsc-liftover/README.html) use the following command `conda install ucsc-liftover`

To ease the process, there is a conda environment file in `Code/mergeref.yml`. 
With anaconda already installed you can create the same environment used to run the script:

```
conda env create -f mergeref.yml
source activate mergeref
```

## Files to download

In the `DataBases` folder you'll need to download the respective files.
In each folder (`HGDP` and `1000G`) there is a README file with the same information.

### HGDP

Download the following files and paste them in the `DataBases/HGDP` folder.
The HGDP Stanford files can be downloaded from [here](http://hagsc.org/hgdp/files.html).
You will also need to download the Sample Information from [here](https://web.stanford.edu/group/rosenberglab/data/rosenberg2006ahg/SampleInformation.txt).
Finally, you'll need to download the chain file that tells liftOver how to convert between hg18 to hg19 from [here](http://hgdownload.cse.ucsc.edu/goldenPath/hg18/liftOver/hg18ToHg19.over.chain.gz).

### 1000G

Download the following files and and paste them in the `DataBases/1000G` folder.
The 1000G Phase 3 files can be downloaded from [here](ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/).

## Steps

[Here](https://nbviewer.jupyter.org/github/tomszar/HGDP_1000G_Merge/blob/master/Code/2018-05-MergeGenotypes.ipynb) you can see the steps for merging the HGDP and 1000G databases.





# Merging HGDP and 1000 Genomes

This script will walk you through the merging of the HGDP and 1000 Genomes genotypes.
We will follow these steps:
- Download HGDP files, and transform them into a plink file
- Download 1000G files and keep only SNPs found in the HGDP files
- Merge the HGDP and 1000G files

## Needed packages

To run this script you'll need [plink](https://www.cog-genomics.org/plink2), [vcftools](https://vcftools.github.io/index.html), [bcftools](https://samtools.github.io/bcftools/bcftools.html).
All of these can be found in [bioconda](https://bioconda.github.io/)

## Preliminaries

First let's import modules and set up paths. Open Python.

```{python}
import glob, os, shutil, subprocess, csv, time
import pandas as pd
import numpy as np
```

```{python}
projpath = os.path.realpath("../HGDP_1000G_Merge")
pathhgdp = os.path.join(projpath, "DataBases", "HGDP")
path1000 = os.path.join(projpath, "DataBases", "1000G")
print(projpath)
```

## HGDP

The HGDP files (Stanford) were downloaded from [here](http://hagsc.org/hgdp/files.html), and the sample list file from [here](http://www.stanford.edu/group/rosenberglab/data/rosenberg2006ahg/SampleInformation.txt).
The script to transform the HGDP data to plink format is called HGDPtoPlink.sh and was modified from [here](http://www.harappadna.org/2011/02/hgdp-to-ped-conversion/).
The HGDP data uses coordinates from build 36.1 (a list of assemblies can be found [here](https://genome.ucsc.edu/FAQ/FAQreleases.html))

```{python}
#Run the HGDPtoPlink script
#It can take a while
os.chdir(os.path.join(projpath, "Code"))
subprocess.run(["bash", "HGDPtoPlink.sh"])
```

Because the 1000G uses the GRCh37 assembly (fasta file can be found [here](http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/) as human_g1k_v37.fasta.gz) we'll need to liftover the HGDP coordinates.
To do that we'll use UCSC [liftOver](http://genome.ucsc.edu/cgi-bin/hgLiftOver) already installed using bioconda, and [liftOverPlink](https://github.com/sritchie73/liftOverPlink) as a wrapper to work with plink files (`ped` and `map` formats).
The chain file that tells liftOver how to convert between hg18 and hg19 can be downloaded [here](http://hgdownload.cse.ucsc.edu/goldenPath/hg18/liftOver/hg18ToHg19.over.chain.gz).

This part is run in bash. Exit python.
```{bash}
#Using liftover
Code/liftoverplink --map Databases/HGDP/hgdp940.map --out Databases/HGDP/lifted --chain Databases/HGDP/hg18ToHg19.over.chain.gz 
Code/rmBadLifts --map Databases/HGDP/lifted.map --out Databases/HGDP/good_lifted.map --log Databaes/HGDP/bad_lifted.dat
```

Back to python.
```{python}
os.chdir(pathhgdp)

#Creating a list of snps to include in lifted version
snps = pd.read_csv("good_lifted.map", sep = "\t", header = None)
snps.iloc[:,1].to_csv("snplist.txt", index = False)
#Excluding snps and creating binary file
subprocess.run(["plink", "--file", "hgdp940", "--recode", "--out", "lifted", "--extract", "snplist.txt" ])
subprocess.run(["plink", "--file", "--ped", "lifted.ped", "--map", "good_lifted.map", "--make-bed", "--out", "hgdp940hg19"])

#Removing some files
for file in glob.glob("*.ped"):
    os.remove(file)
    
for file in glob.glob("*.map"):
    os.remove(file)
```

```{python}
#Read hgdp940hg19 log file
with open("hgdp940hg19.log") as myfile:
    for num, line in enumerate(myfile, 1):
        if "variants" in line:
            print(line, end='')
    print("Finished file... \n")
```

```{python}
dest_dir = os.path.join(projpath, "Results")
for filename in glob.glob("hgdp940hg19.*"):
    shutil.copy(filename, dest_dir)
```

## 1000G

The 1000G Phase 3 files were downloaded from [here](ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/).
First, we will extract the list of snps (`snplist.txt`) from the HGDP dataset for each chromosome of the 1000G samples using vcftools.
Then we will concatenate the different autosomal chromosomes in one file and convert it into a plink binary file using bcftools and plink.

```{python}
#Extracting snps found in the HGDP to the 1000G files
#Takes time
os.chdir(path1000)
shutil.copy(os.path.join(pathhgdp, "snplist.txt"), path1000)
for file in glob.glob("*chr[0-9]*.gz"):
    outname = file.split(".")[1] + "_extracted"
    subprocess.run(["vcftools", "--gzvcf", file, "--snps", "snplist.txt", "--recode", "--out", outname])
```

```{python}
concatfiles = glob.glob("chr[0-9]*.recode.vcf")
function = ["bcftools", "concat", "-o", "1000g.vcf.gz", "-Oz"]
function.extend(concatfiles)
subprocess.run(function)
```

```{python}
#Convert to binary plink
subprocess.run(["plink", "--vcf", "1000g.vcf.gz", "--make-bed", "--out", "1000Ghg19" ])
#Updating fam file
allfam = pd.read_csv("integrated_call_samples_v2.20130502.ALL.ped", header = None, skiprows = 1, sep = "\t")
oldfam = pd.read_csv("1000Ghg19.fam", header = None, sep = " ")
updatedfam = pd.merge(oldfam, allfam, how = "inner", left_on = 1, right_on = 1)
updatedfam.iloc[:,[6,1,7,8,9,5]].to_csv("1000Ghg19.fam", sep = " ", header = False, index = False)
os.remove("1000g.vcf.gz")
```

```{python}
for file in glob.glob("chr*.recode.vcf"):
    os.remove(file)
    
dest_dir = os.path.join(projpath, "Results")
for filename in glob.glob("1000Ghg19.*"):
    shutil.copy(filename, dest_dir)
```

## Merge reference samples

Now we will merge the 1000G and HGDP databases, both using the hg19 coordinates and with related people removed.
We will attempt an initial merge, and then flip strands for the problematic SNPs. 
Then we will try a second attempt, and remove the still problematic snps from both databases, to then attempt the final merge.

```{python}
os.chdir(os.path.join(projpath, "Results"))
#Preliminary merge
subprocess.run(["plink", "--bfile", "1000Ghg19", "--bmerge", "hgdp940hg19", "--make-bed", "--out", "hgdp1000ghg19"])
#Flip snps
subprocess.run(["plink", "--bfile", "1000Ghg19", "--flip", "hgdp1000ghg19-merge.missnp", "--make-bed", "--out", "1000Ghg19_flip"])
subprocess.run(["plink", "--bfile", "1000Ghg19_flip", "--bmerge", "hgdp940hg19", "--make-bed", "--out", "hgdp1000ghg19"])
#Removing snps
for file in glob.glob("*.bed"):
    outname = file.split(".")[0] + "_temp"
    subprocess.run(["plink", "--bfile", file.split(".")[0], "--exclude", "hgdp1000ghg19-merge.missnp", "--make-bed", "--out", outname])

#Removing temp files
subprocess.run(["plink", "--bfile", "hgdp940hg19_temp", "--bmerge", "1000Ghg19_flip_temp", "--make-bed", "--out", "hgdp1000ghg19"])
for file in glob.glob("*_temp*"):
    os.remove(file)
    
for file in glob.glob("*_flip*"):
    os.remove(file)

with open("hgdp1000ghg19.log", 'r') as fin:
    file_contents = fin.read()
    print(file_contents)
```

