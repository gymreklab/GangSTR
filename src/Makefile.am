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

HTS_LIBS = -Ilib/htslib -lhts

AM_CPPFLAGS = -I$(top_srcdir) $(GSL_CFLAGS)
AM_LDFLAGS = $(GSL_LIBS)

bin_PROGRAMS = GangSTR

noinst_LIBRARIES = libgangstr.a

libgangstr_a_CPPFLAGS = $(AM_CPPFLAGS)
libgangstr_a_SOURCES = \
	common.h common.cpp \
	options.h options.cpp \
	locus.h locus.cpp \
	region_reader.h region_reader.cpp \
	ref_genome.h ref_genome.cpp \
	genotyper.h genotyper.cpp \
	read_class.h read_class.cpp \
	frr_class.h frr_class.cpp \
	enclosing_class.h enclosing_class.cpp \
	spanning_class.h spanning_class.cpp \
	likelihood_maximizer.h likelihood_maximizer.cpp \
	mathops.h mathops.cpp \
	fastonebigheader.h \
	read_extractor.h read_extractor.cpp \
	bam_io.h bam_io.cpp \
	stringops.h stringops.cpp \
	read_pair.h read_pair.cpp \
	realignment.h realignment.cpp

GangSTR_LDFLAGS = $(AM_LDFLAGS) $(HTS_LIBS)
GangSTR_SOURCES = main_gangstr.cpp
GangSTR_CPPFLAGS = $(AM_CPPFLAGS)
GangSTR_CFLAGS = $(AM_CFLAGS)
GangSTR_CXXFLAGS = $(AM_CXXFLAGS)
GangSTR_LDADD = libgangstr.a $(GSL_LIBS)

TESTS = GangSTRTests
AUTOMAKE_OPTIONS = serial-tests
check_PROGRAMS=$(TESTS)
GangSTRTests_SOURCES = \
	read_class.h read_class.cpp \
	enclosing_class.h enclosing_class.cpp \
	frr_class.h frr_class.cpp \
	spanning_class.h spanning_class.cpp \
	mathops.h mathops.cpp \
	fastonebigheader.h \
	likelihood_maximizer.h likelihood_maximizer.cpp \
	tests/ReadClass_test.h tests/ReadClass_test.cpp \
	tests/LiketestslihoodMaximizer_test.sh tests/LikelihoodMaximizer_test.cpp \
	tests/main_test.cpp
GangSTRTests_CXXFLAGS = $(CPPUNIT_CFLAGS) $(GSL_CFLAGS)
GangSTRTests_LDFLAGS = $(CPPUNIT_LIBS) -ldl $(GSL_LIBS)
