#!/bin/bash
# Copyright (C) 2017 Melissa Gymrek <mgymrek@ucsd.edu>
# and Nima Mousavi (mousavi@ucsd.edu)
#
# This file is part of GangSTR.
#
# GangSTR is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# GangSTR is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with GangSTR.  If not, see <http://www.gnu.org/licenses/>.

die()
{
    BASE=$(basename "$0")
    echo "$BASE error: $1" >&2
    exit 1
}

check_program()
{
	PROG="$1"
	# If the program is found in the $PATH, just return
	which "$PROG" >/dev/null && return

	# Otherwise, show an error with some helpful information
	echo "
--- GangSTR Compilation Error ---
You are trying to install GangSTR.
A required program '$PROG' was not found.
To run this script, the following programs are needed:
   wget
   A C++ compiler
   pkg-config
You additionally must have an internet connection
---
"
	exit 1
}

usage()
{
    BASE=$(basename -- "$0")
    echo "Install GangSTR and dependencies
Usage:
   ./$BASE [PREFIX]

If a PREFIX is provided, will install GangSTR and all dependencies to that directory. Otherwise install do default directory.

IF YOU ARE INSTALLING LOCALLY, it is recommended to use your home directory as the prefix. e.g.:
     ./$BASE /home/<username>

IF YOU ARE INSTALLING AS ROOT, you don't need to pass a prefix, but probably need to use sudo, e,g.:
     sudo ./$BASE
"
    exit 0
}

for PROG in wget make pkg-config;
do
  check_program $PROG
done

PREFIX=$1

if [ "x$1" == "x-h" -o "x$1" == "x--help" ]
then
    usage
fi

# If PREFIX not set, use default.
if [ "x$PREFIX" == "x" ]
then
    PREFIX=/usr/local
else
    export PKG_CONFIG_PATH=$PREFIX/lib/pkgconfig
fi

# Check that prefix exists
if [ ! -d "$PREFIX" ]; then
    die "Directory $PREFIX does not exist. Try running: ./install-gangstr.sh $HOME"
fi

echo "[install-deps.sh] Compiling GangSTR+dependencies with prefix=$PREFIX"

mkdir -p dependences || die "Could not make dependencies directory"

# Install GSL 
echo "[install-deps.sh] Compiling GSL..."
cd dependences
wget ftp://ftp.gnu.org/gnu/gsl/gsl-2.5.tar.gz || die "Error downloading GSL"
tar -xzvf gsl-2.5.tar.gz || die "Error unzipping GSL"
cd gsl-2.5/ || die "Error navigating to GSL directory"
./configure --prefix=$PREFIX || die "Error configuring GSL"
make || die "Error compiling GSL"
make install || die "Error installing GSL"

# Install NLOPT 
echo "[install-deps.sh] Compiling GSL..."
cd ../ || die "Error navigating to dependencies"
wget http://ab-initio.mit.edu/nlopt/nlopt-2.4.2.tar.gz || die "Error downloading NLOPT"
tar -xzvf nlopt-2.4.2.tar.gz || die "Error unzipping NLOPT"
cd nlopt-2.4.2 || die "Error navigating to NLOPT directory"
./configure --prefix=$PREFIX || die "Error configuring NLOPT"
make || die "Error compiling NLOPT"
make install || die "Error installing NLOPT"

# Install HTSLIB 
cd ../ || die "Error navigating to dependencies"
wget https://github.com/samtools/htslib/releases/download/1.8/htslib-1.8.tar.bz2 || die "Error downloading HTSLIB"
tar -xjvf htslib-1.8.tar.bz2 || die "Error unzipping HTSLIB"
cd htslib-1.8/ || die "Error navigating to HTSLIB"
./configure --disable-lzma --disable-bz2 --prefix=$PREFIX || die "Error configuring HTSLIB"
make || die "Error compiling HTSLIB"
make install || die "Error installing HTSLIB"
cd ../../

# Install GangSTR
./configure --prefix=$PREFIX || die "Error configuring GangSTR"
make || die "Error compiling GangSTR"
make install || die "Error installing GangSTR"

echo "[install-deps.sh] Success! Type GangSTR --help to make sure GangSTR is installed"

exit 0
