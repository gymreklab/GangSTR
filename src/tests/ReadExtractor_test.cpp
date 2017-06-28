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

#include "src/bam_io.h"
#include "src/stringops.h"

#include "src/tests/ReadExtractor_test.h"

#include <vector>

// Registers the fixture into the 'registry'
CPPUNIT_TEST_SUITE_REGISTRATION(ReadExtractorTest);

void ReadExtractorTest::setUp() {
  // This gets set during "make check"
  // env variable set in ./src/Makefile.am
  test_dir = getenv("GANGSTR_TEST_DIR");
  // Set up locus
  locus.chrom = "3";
  locus.start = 63898361;
  locus.end = 63898392;
}

void ReadExtractorTest::tearDown() {}

void ReadExtractorTest::test_ExtractReads() {
  CPPUNIT_FAIL("test_ExtractReads not implemented");
}

void ReadExtractorTest::test_FindDiscardedRead() {
  std::string spanning_file = test_dir + "/test.spanning.bam";
  std::string enclosing_file = test_dir + "/test.enclosing.bam";
  std::string frr_file = test_dir + "/test.frr.bam";
  std::string discarded_file = test_dir + "/test.discarded.bam";
  std::vector<std::string> files(0);
  files.push_back(spanning_file);
  files.push_back(enclosing_file);
  files.push_back(frr_file);
  BamCramMultiReader bamreader(files);
  int32_t BUFFER = 50000; // want to capture all reads
  int32_t chrom_ref_id = bamreader.bam_header()->ref_id(locus.chrom);
  bamreader.SetRegion(locus.chrom, locus.start - BUFFER, locus.end + BUFFER);
  // Loop through alignments
  BamAlignment aln;
  int32_t insert_size;
  bool truth_is_discarded, test_is_discarded;
  while (bamreader.GetNextAlignment(aln)) {
    std::string fname = aln.Filename();
    truth_is_discarded = string_ends_with(fname, "test.discarded.bam");
    test_is_discarded = read_extractor.FindDiscardedRead(aln, chrom_ref_id, locus);
    std::stringstream msg;
    msg << "Misclassified discarded read " << aln.Name();
    CPPUNIT_ASSERT_EQUAL_MESSAGE(msg.str(), truth_is_discarded, test_is_discarded);
  }
}

void ReadExtractorTest::test_FindSpanningRead() {
  std::string spanning_file = test_dir + "/test.spanning.bam";
  std::string enclosing_file = test_dir + "/test.enclosing.bam";
  std::string frr_file = test_dir + "/test.frr.bam";
  std::string discarded_file = test_dir + "/test.discarded.bam";
  std::vector<std::string> files(0);
  files.push_back(spanning_file);
  files.push_back(enclosing_file);
  files.push_back(frr_file);
  BamCramMultiReader bamreader(files);
  int32_t BUFFER = 50000; // want to capture all reads
  int32_t chrom_ref_id = bamreader.bam_header()->ref_id(locus.chrom);
  bamreader.SetRegion(locus.chrom, locus.start - BUFFER, locus.end + BUFFER);
  // Loop through alignments
  BamAlignment aln;
  int32_t insert_size;
  bool truth_is_spanning, test_is_spanning;
  while (bamreader.GetNextAlignment(aln)) {
    std::string fname = aln.Filename();
    truth_is_spanning = string_ends_with(fname, "test.spanning.bam");
    test_is_spanning = read_extractor.FindSpanningRead(aln, chrom_ref_id, locus, &insert_size);
    std::stringstream msg;
    msg << "Misclassified spanning read " << aln.Name();
    CPPUNIT_ASSERT_EQUAL_MESSAGE(msg.str(), truth_is_spanning, test_is_spanning);
  }
}

void ReadExtractorTest::test_ProcessSingleRead() {
  CPPUNIT_FAIL("test_ProcessSingleRead not implemented");
}

// TODO add case where mate not actually in file, should return false
void ReadExtractorTest::test_RescueMate() {
  std::string frr_file = test_dir + "/test.frr.bam";
  std::vector<std::string> files(0);
  files.push_back(frr_file);
  BamCramMultiReader bamreader(files);
  int32_t BUFFER = 50000; // want to capture all reads
  bamreader.SetRegion(locus.chrom, locus.start - BUFFER, locus.end + BUFFER);
  // Loop through alignments, do several example cases
  BamAlignment aln;
  BamAlignment matealn;
  while (bamreader.GetNextAlignment(aln)) {
    read_extractor.RescueMate(&bamreader, aln, &matealn);
    std::stringstream msg;
    msg << "RescueMate failed for " << aln.Name();
    if (aln.Name() == "ATXN7_52_cov60_dist500_DIP_const70_70_constAllele_3065_3556_0:0:0_0:0:0_136") {
      CPPUNIT_ASSERT_EQUAL_MESSAGE(msg.str(), matealn.RefID(), aln.MateRefID());
      CPPUNIT_ASSERT_EQUAL_MESSAGE(msg.str(), matealn.Position(), 53253385);
    }
    if (aln.Name() == "ATXN7_52_cov60_dist500_DIP_const70_70_constAllele_2667_3154_0:0:0_0:0:0_a8") {
      CPPUNIT_ASSERT_EQUAL_MESSAGE(msg.str(), matealn.RefID(), aln.MateRefID());
      CPPUNIT_ASSERT_EQUAL_MESSAGE(msg.str(), matealn.Position(), 53253384);
    }
  }
}


