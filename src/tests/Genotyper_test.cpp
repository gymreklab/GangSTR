/*
Copyright (C) 2017 Melissa Gymrek <mgymrek@ucsd.edu>
and Nima Mousavi (mousavi@ucsd.edu)

This file is part of GangSTR.

GangSTR is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

GangSTR is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with GangSTR.  If not, see <http://www.gnu.org/licenses/>.
*/

#include "src/tests/Genotyper_test.h"

// Registers the fixture into the 'registry'
CPPUNIT_TEST_SUITE_REGISTRATION(GenotyperTest);

void GenotyperTest::setUp() {
  test_dir = getenv("GANGSTR_TEST_DIR");
  locus.chrom = "3";
  locus.start = 201;
  locus.end = 230;
  locus.motif = "CAG";
}

void GenotyperTest::tearDown() {}

void GenotyperTest::test_SetFlanks() {
  Options options;
  options.flanklen = 99;
  std::string fastafile = test_dir + "/test.fa";
  RefGenome refgenome(fastafile);
  Genotyper genotyper(refgenome, options);
  std::string pre_flank = "GGAGCGGAAAGAATGTCGGAGCGGGCCGCGGATGACGTCAGGGGGGAGCCGCGCCGCGCGGCGGCGGCGGCGGGCGGAGCAGCGGCCGCGGCCGCCCGG";
  std::string post_flank = "CCGCCGCCTCCGCAGCCCCAGCGGCAGCAGCACCCGCCACCGCCGCCACGGCGCACACGGCCGGAGGACGGCGGGCCCGGCGCCGCCTCCACCTCGGCC";
  if (!genotyper.SetFlanks(&locus)) {
    CPPUNIT_FAIL("SetFlanks returned false unexpectedly");
  }
  CPPUNIT_ASSERT_EQUAL(pre_flank, locus.pre_flank);
  CPPUNIT_ASSERT_EQUAL(post_flank, locus.post_flank);
}

void GenotyperTest::test_ProcessLocus() {
  CPPUNIT_FAIL("test_ProcessLocus() not implemented");
}
