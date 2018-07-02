# **UNDER CONSTRUCTION** OFFICIAL GANGSTR RELEASE COMING SOON

# GangSTR

GangSTR is a tool for genome-wide profiling tandem repeats from short reads. A key advantage of GangSTR over existing tools (e.g. [lobSTR](https://github.com/mgymrek/lobstr-code) or [hipSTR](https://github.com/tfwillems/HipSTR)) is that it can handle repeats that are longer than the read length.

GangSTR takes aligned reads (BAM) and a set of repeats in the reference genome as input and outputs a VCF file containing genotypes for each locus.

For questions on installation or usage, please open an issue, submit a pull request, or contact Nima Mousavi (mousavi@ucsd.edu).

[Download](#download) | [Install](#install) | [Usage](#usage) | [File formats](#formats) | [References](#references)

<a name="download"></a>
## Download

The latest GangSTR release is GangSTR-1.3. Download [GangSTR-1.3.4-a49a.tar.gz](https://github.com/gymreklab/GangSTR/releases/download/test/GangSTR-1.3.4-a49a.tar.gz).

For a list of previous releases and changelog see the [releases page](https://github.com/gymreklab/GangSTR/releases).

For a list of TR references available, see [references](#references) below. 

<a name="install"></a>

## Basic Install

GangSTR requires third party packages [nlopt](https://nlopt.readthedocs.io/en/latest/), [gsl](https://www.gnu.org/software/gsl/doc/html/index.html), and [htslib](http://www.htslib.org//). The built-in script `install-gangstr.sh` installs these for you before compiling and installing GangSTR. Both UNIX and Mac OSX are supported.

If you are running as root:
```
tar -xzvf GangSTR-X.X.tar.gz
cd GangSTR-X.X
sudo ./install-gangstr.sh
```

If you are installing locally (e.g. on a cluster where you don't have root access):
```
tar -xzvf GangSTR-X.X.tar.gz
cd GangSTR-X.X
./install-gangstr.sh PREFIX
```

where `PREFIX` is a place you have write permissions. In most cases this will be your home directory, e.g. `$HOME`. If you install locally, make sure `$PREFIX/bin` is on your `PATH`.

Typing `GangSTR --help` should show a help message if GangSTR was successfully installed.

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
GangSTR requires a reference set of regions to genotype. This is a BED-like file with the following columns:

1. The name of the chromosome on which the STR is located
2. The start position of the STR on its chromosome
3. The end position of the STR on its chromosome
4. The motif length
5. The repeat motif

Below is an example file which contains 5 TR loci. Standard references for hg19 and GRCh38 can be obtained [below](#references).
**NOTE: The table header is for descriptive purposes. The BED file should not have a header**

| **CHROM** | **START** | **END** | **MOTIF_LEN** | **MOTIF** | **OFFTARGET** |
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
| hg19 | ver8 | [hg19_ver8.bed.gz](https://s3.amazonaws.com/gangstr/hg19_ver8.bed.gz) | 
| hs37 | ver8 | [hs37_ver8.bed.gz](https://s3.amazonaws.com/gangstr/hs37_ver8.bed.gz) |
| hg38 | ver5 | [hg38_ver5.bed.gz](https://s3.amazonaws.com/gangstr/hg38_ver5.bed.gz) |

