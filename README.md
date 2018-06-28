# GangSTR

GangSTR is a tool for genome-wide profiling STRs and short VNTRs from short reads. A key advantage of GangSTR over existing tools (e.g. [lobSTR](https://github.com/mgymrek/lobstr-code) or [hipSTR](https://github.com/tfwillems/HipSTR)) is that it can handle repeats that are longer than the read length.

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

<a name="install"></a>
## Install

GangSTR requires third party packages [nlopt](https://nlopt.readthedocs.io/en/latest/), [gsl](https://www.gnu.org/software/gsl/doc/html/index.html), and [htslib](http://www.htslib.org//). Before proceeding these should be installed. See below for details on installing dependencies. Once the dependencies are installed, you can download and compile GangSTR using the following steps:

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

### Dependencies

### Installing locally (e.g. without root access)

<a name="usage"></a>
## Usage

<a name="formats"></a>
## File formats

<a name="references"></a>
## GangSTR reference files

<a name="citation"></a>
## Citing GangSTR
