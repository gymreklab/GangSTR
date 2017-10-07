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
#include "src/bam_io.h"

// Registers the fixture into the 'registry'
CPPUNIT_TEST_SUITE_REGISTRATION(GenotyperTest);

void GenotyperTest::setUp() {
  test_dir = getenv("GANGSTR_TEST_DIR");
  // locus.chrom = "3";
  // locus.start = 201;
  // locus.end = 230;
  // locus.motif = "CAG";
  locus.chrom = "19";
  // locus.start = 13318672;
  // locus.end = 13318710;
  locus.start = 5000;
  locus.end = 5039;
  locus.motif = "CTG";
  locus.period = 3;
}

void GenotyperTest::tearDown() {}

void GenotyperTest::test_SetFlanks() {
  Options options;
  options.realignment_flanklen = 99;
  options.dist_mean = 500;
  options.dist_sdev = 50;
  options.flanklen = 3000;
  options.frr_weight = 0.5;
  options.enclosing_weight = 1.0;
  options.spanning_weight = 1.0;
  options.flanking_weight = 1.0;
  options.read_len = 100;
  options.dist_max = 1000;

  locus.chrom = "3";
  locus.start = 201;
  locus.end = 230;
  locus.motif = "cag";
  locus.period = 3;
  std::string fastafile = test_dir + "/test.fa";
  RefGenome refgenome(fastafile);
  Genotyper genotyper(refgenome, options);
  std::string pre_flank = "ggagcggaaagaatgtcggagcgggccgcggatgacgtcaggggggagccgcgccgcgcggcggcggcggcgggcggagcagcggccgcggccgcccgg";
  std::string post_flank = "ccgccgcctccgcagccccagcggcagcagcacccgccaccgccgccacggcgcacacggccggaggacggcgggcccggcgccgcctccacctcggcc";
  if (!genotyper.SetFlanks(&locus)) {
    CPPUNIT_FAIL("SetFlanks returned false unexpectedly");
  }
  CPPUNIT_ASSERT_EQUAL(pre_flank, locus.pre_flank);
  CPPUNIT_ASSERT_EQUAL(post_flank, locus.post_flank);
}

void GenotyperTest::test_ProcessLocus() {
  Options options;
  options.realignment_flanklen = 99;
  options.dist_mean = 500;
  options.dist_sdev = 50;
  options.flanklen = 3000;
  options.frr_weight = 0.5;
  options.enclosing_weight = 1.0;
  options.spanning_weight = 1.0;
  options.flanking_weight = 1.0;
  options.read_len = 100;
  options.dist_max = 1000;
  std::string fastafile = test_dir + "/CACNA1A_5k_region.fa";
  // std::string fastafile = "/storage/reslsources/dbase/human/hs37d5/hs37d5.fa";
  RefGenome refgenome(fastafile);
  Genotyper genotyper(refgenome, options);

  std::string bam_file = test_dir + "/54_nc_40.sorted.bam";
  std::vector<std::string> files(0);
  files.push_back(bam_file);
  BamCramMultiReader* bamreader = new BamCramMultiReader(files, fastafile);
  CPPUNIT_ASSERT_EQUAL(genotyper.ProcessLocus(bamreader, &locus), true);
}
