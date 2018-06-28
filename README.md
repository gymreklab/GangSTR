# GangSTR

GangSTR is a tool for genome-wide profiling tandem repeats from short reads. A key advantage of GangSTR over existing tools (e.g. [lobSTR](https://github.com/mgymrek/lobstr-code) or [hipSTR](https://github.com/tfwillems/HipSTR)) is that it can handle repeats that are longer than the read length.

GangSTR takes aligned reads (BAM) and a set of repeats in the reference genome as input and outputs a VCF file containing genotypes for each locus.

For questions on installation or usage, please open an issue, submit a pull request, or contact Nima Mousavi (mousavi@ucsd.edu).

[Download](#download)

[Install](#install)

[Usage](#usage)

[File formats](#formats)

[References](#references)

[Citation](#citation)


<a name="download"></a>
## Download

The latest GangSTR release is [GangSTR-1.2](https://github.com/gymreklab/GangSTR/releases/download/v1.2/GangSTR-1.2.tar.gz).

For a list of previous releases and changelog see the [releases page](https://github.com/gymreklab/GangSTR/releases).

For a list of TR references available, see [references](#references) below. 

<a name="install"></a>
## Install

### Basic installation
GangSTR requires third party packages [nlopt](https://nlopt.readthedocs.io/en/latest/), [gsl](https://www.gnu.org/software/gsl/doc/html/index.html), and [htslib](http://www.htslib.org//). See below for details on [installing dependencies](#dependencies) or for [installing without root access](#noroot). Once the dependencies are installed, you can download and compile GangSTR using the following steps:

```
# 1. Download and unpack the tarball
wget https://github.com/gymreklab/GangSTR/releases/download/v1.2/GangSTR-1.2.tar.gz
tar -xzvf GangSTR-1.2.tar.gz
cd GangSTR-1.2

# 2. Compile
./configure
make
make install
```

This will install the binary `GangSTR` to the default `PATH` location (e.g. `/usr/bin/GangSTR`).

<a name="dependencies"></a>
### Dependencies
GangSTR relies on three external libraries: [nlopt](https://nlopt.readthedocs.io/en/latest/), [gsl](https://www.gnu.org/software/gsl/doc/html/index.html), and [htslib](http://www.htslib.org//). These can be installed using the commands below. Additional detailed instructions for installing each package are available on the project sites. 

```
# Install GSL 
wget ftp://ftp.gnu.org/gnu/gsl/gsl-2.5.tar.gz
tar -xzvf gsl-2.5.tar.gz
cd gsl-2.5/
./configure
make
make install # may require sudo access

# Install NLOPT 
wget http://ab-initio.mit.edu/nlopt/nlopt-2.4.2.tar.gz
tar -xzvf nlopt-2.4.2.tar.gz
cd nlopt-2.4.2
./configure
make
make install # may require sudo access

# Install HTSLIB 
wget https://github.com/samtools/htslib/releases/download/1.8/htslib-1.8.tar.bz2
tar -xjvf htslib-1.8.tar.bz2
cd htslib-1.8/
./configure
make
make install # may require sudo access
```

<a name="noroot"></a>
### Installing locally (e.g. without root access)

GangSTR and all dependencies can be compiled and installed without root access by changing the installation prefix to your home directory (or another place you have permissions to write). You may also need to edit the environment variable `PKG_CONFIG_PATH`. The following give complete instructions for compiling locally without root access.

```
# Install GSL 
wget ftp://ftp.gnu.org/gnu/gsl/gsl-2.5.tar.gz
tar -xzvf gsl-2.5.tar.gz
cd gsl-2.5/
./configure --prefix=/home/username
make
make install 

# Install NLOPT 
wget http://ab-initio.mit.edu/nlopt/nlopt-2.4.2.tar.gz
tar -xzvf nlopt-2.4.2.tar.gz
cd nlopt-2.4.2
./configure --prefix=/home/username
make
make install 

# Install HTSLIB 
wget https://github.com/samtools/htslib/releases/download/1.8/htslib-1.8.tar.bz2
tar -xjvf htslib-1.8.tar.bz2
cd htslib-1.8/
./configure --prefix=/home/username
make
make install 

# Update PKG_CONFIG_PATH just in case
export PKG_CONFIG_PATH=/home/username/lib/pkgconfig

# Install GangSTR
wget https://github.com/gymreklab/GangSTR/releases/download/v1.2/GangSTR-1.2.tar.gz
tar -xzvf GangSTR-1.2.tar.gz
cd GangSTR-1.2
./configure --prefix=/home/username
make
make install 
```

<a name="usage"></a>
## Usage

<a name="formats"></a>
## File formats

GangSTR takes as input a BAM file of short read alignments, a reference set of TRs, and a reference genome, and outputs genotypes in a VCF file. Each of these formats is described below.

### BAM (`--bam`)
GangSTR requires a [BAM](https://samtools.github.io/hts-specs/SAMv1.pdf) file produced by an indel-sensitive aligner. The BAM file must be sorted and indexed e.g. by using `samtools sort` and `samtools index`. GangSTR currently only processes a single sample at a time.

### FASTA Reference genome (`--ref`)
You must input a reference genome in FASTA format. This must be the same reference build used to align the sequences in the BAM file.

### TR regions (`--regions`)
GangSTR requires a reference set of regions to genotype. This is a BED-like file with the following required columns:

1. The name of the chromosome on which the STR is located
2. The start position of the STR on its chromosome
3. The end position of the STR on its chromosome
4. The motif length
5. The repeat motif
6. The sequence of the repeat region in the reference genome

Lines beginning with `*` in the reference are ignored.

Below is an example file which contains 5 TR loci. Standard references for hg19 and GRCh38 can be obtained [below](#references).
**NOTE: The table header is for descriptive purposes. The BED file should not have a header**

| **CHROM** | **START** | **END** | **MOTIF_LEN** | **MOTIF** | **REFSEQ** |
|-----------|-----------|---------|----------------|----------|------------|
|chr1	|10689	|10700|	5	|CGCGC|	CGCGCCGCGCCG|
| chr1  |  28589  | 28603  | 1 |      T    |   TTTTTTTTTTTTTTT|
|chr4  |  11173|   11194  | 11   |   CGCCGGCGCGG |    CGCCGGCGCGGCGCCGGGGCGG|
|chr4   | 150889 | 150909 | 2    |   TG      |TGTGTGTGTGTGTGTGTGTGT|
|chr20  | 72993 |  73011|   5   |    TACTA  | TACTACAATATACTATACT|

### VCF (output)
For more information on VCF file format, see the [VCF spec](http://samtools.github.io/hts-specs/VCFv4.2.pdf). In addition to standard VCF fields, GangSTR adds custom fields described below.

#### INFO fields

INFO fields contain aggregated statistics about each TR. The following custom fields are added:

| **FIELD** | **DESCRIPTION** |
|-----------|------------------|
| END | End position of the TR |
| RU| Repeat motif | 
| REF| Reference copy number (number of repeat units| 

#### FORMAT fields
FORMAT fields contain information specific to each genotype call. The following custom fields are added:

| **FIELD** | **DESCRIPTION** |
|-----------|------------------|
| GB | Base pair length differences of genotype from reference for each allele |
| CI| 95% confidence intervals for each allele | 
| RC| Number of reads in each class (enclosing, spanning, FRR, flanking)| 
| Q| Minimum negative likelihood| 
| INS| Insert size mean and stddev at the locus| 


<a name="references"></a>
## GangSTR reference files

The following lists available references created using Tandem Repeats Finder. We update the reference periodically with additional loci or annotation changes. Please see the Changelog file for details. Note references must be unzipped before using with GangSTR.

| **Reference build** | **Version** | **Link** |
| --------------------| ------------|----------|
| hg19 | ver9 | [hg19_ver9.bed.gz](https://s3.amazonaws.com/gangstr/hg19_ver9.bed.gz) |
| hs38 | ver8 | [hg37_ver8.bed.gz](https://s3.amazonaws.com/gangstr/hs37_ver8.bed.gz) |


<a name="citation"></a>
## Citing GangSTR
