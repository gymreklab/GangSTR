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

AM_CPPFLAGS = -I$(top_srcdir) $(HTS_CPPFLAGS) ${BAM_CPPFLAGS}
AM_LDFLAGS = ${BAM_LDFLAGS} ${HTS_LDFLAGS}

bin_PROGRAMS = GangSTR

noinst_LIBRARIES = libgangstr.a

libgangstr_a_SOURCES = \
	common.h common.cpp \
	options.h options.cpp \
	locus.h locus.cpp \
	region_reader.h region_reader.cpp \
	bam_reader.h bam_reader.cpp \
	ref_genome.h ref_genome.cpp \
	genotyper.h genotyper.cpp \
	read_class.h read_class.cpp \
	frr_class.h frr_class.cpp \
	enclosing_class.h enclosing_class.cpp \
	spanning_class.h spanning_class.cpp \
	likelihood_maximizer.h likelihood_maximizer.cpp \
	mathops.h mathops.cpp \
	fastonebigheader.h

GangSTR_LDFLAGS = ${AM_LDFLAGS}
GangSTR_SOURCES = main_gangstr.cpp
GangSTR_LDADD = libgangstr.a

TESTS = GangSTRTests
check_PROGRAMS=$(TESTS)
GangSTRTests_SOURCES = \
	tests/ReadClass_test.h tests/ReadClass_test.cpp \
	tests/main_test.cpp
GangSTRTests_CXXFLAGS = $(CPPUNIT_CFLAGS)
GangSTRTests_LDFLAGS = $(CPPUNIT_LIBS) -ldl
