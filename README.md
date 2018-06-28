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

<a name="references"></a>
## GangSTR reference files

<a name="citation"></a>
## Citing GangSTR
