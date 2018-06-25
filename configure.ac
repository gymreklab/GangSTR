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

AC_INIT([GangSTR], [1.0], [mgymrek@ucsd.edu])
AM_INIT_AUTOMAKE([foreign -Wall -Werror subdir-objects])
AC_PROG_CC
AC_PROG_CXX
AC_PROG_RANLIB
AM_PROG_AR
AM_PROG_CC_C_O 

# Set variables using shell commands
GIT_VERSION=${PACKAGE_VERSION}
AC_SUBST(GIT_VERSION)
MACHTYPE=`uname -m`
AC_SUBST(MACHTYPE)

# Check for installed packages
#PKG_CHECK_MODULES([CPPUNIT],[cppunit])
PKG_CHECK_MODULES([GSL],[gsl])
AM_PATH_CPPUNIT(1.9.6)

# Custom Flags to add to CFLAGS/CXXFLAGS
GANGSTR_CFLAGS="-D_GIT_VERSION=\"\\\"${GIT_VERSION}\\\"\" -D_MACHTYPE=\"\\\"${MACHTYPE}\\\"\""

# Set the final value for CFLAGS/CXXFLAGS
CFLAGS="$CFLAGS $GANGSTR_CFLAGS"
CXXFLAGS="$CXXFLAGS $GANGSTR_CFLAGS"
#CFLAGS="$GANGSTR_CFLAGS"
#CXXFLAGS="$GANGSTR_CFLAGS"

AC_CONFIG_HEADERS([config.h])
AC_CONFIG_FILES([Makefile src/Makefile])
AC_OUTPUT
AC_MSG_RESULT([    CFLAGS:  $CFLAGS])
AC_MSG_RESULT([    CXXFLAGS:  $CXXFLAGS])
AC_MSG_RESULT([    GANGSTR_CFLAGS:  $GANGSTR_CFLAGS])
AC_MSG_RESULT([    GSL_LIBS:  $GSL_LIBS])
AC_MSG_RESULT([    GSL_CFLAGS:  $GSL_CFLAGS])
