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
#include "src/likelihood_maximizer.h"

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
  locus.start = 63898362;
  locus.end = 63898391;
  locus.pre_flank = "taggagcggaaagaatgtcggagcgggccgcggatgacgtcaggggggagccgcgccgcgcggcggcggcggcgggcggagcagcggccgcggccgcccgg";
  locus.post_flank = "ccgccgcctccgcagccccagcggcagcagcacccgccaccgccgccacggcgcacacggccggaggacggcgggcccggcgccgcctccacctcggccgc";
  locus.motif = "cag";
  locus.period = 3;
  std::string answers_file = test_dir + "/test_pair_answers.tab";
  LoadAnswers(answers_file, &read_type_answers, &data_answers);
  regionsize = 5000;
  min_match = 0;
  options.dist_mean = 500;
  options.dist_sdev = 50;
  options.flanklen = 3000;
  options.read_len = 100;
  options.dist_max = 1000;
  read_extractor_ = new ReadExtractor(options);
}

void ReadExtractorTest::tearDown() {}

void ReadExtractorTest::test_ExtractReads() {
  std::string bam_file = test_dir + "/test.sorted.bam";
  std::vector<std::string> files(0);
  files.push_back(bam_file);
  BamCramMultiReader bamreader(files);
  LikelihoodMaximizer* likmax = new LikelihoodMaximizer(options);
  if (!read_extractor_->ExtractReads(&bamreader, locus, regionsize, min_match, likmax)) {
    CPPUNIT_FAIL("ExtractReads unexpectedly returned false");
  }
  // Check each class has at least as many as the python code found
  CPPUNIT_ASSERT(likmax->GetEnclosingDataSize()>=10);
  CPPUNIT_ASSERT(likmax->GetSpanningDataSize()>=40);
  CPPUNIT_ASSERT(likmax->GetFRRDataSize()>=30);
  CPPUNIT_ASSERT(likmax->GetFlankingDataSize()>=55);
}

void ReadExtractorTest::LoadAnswers(const std::string& answers_file,
				    std::map<std::string, ReadType>* read_type_answers,
				    std::map<std::string, int32_t>* data_answers) {
  std::ifstream* freader = new std::ifstream(answers_file.c_str());
  while (true) {
    std::string line;
    if (!std::getline(*freader, line)) {
      break;
    }
    std::vector<std::string> items;
    split_by_delim(line, '\t', items);
    std::string read_name = "1_" + items[3];
    ReadType rt = (ReadType)atoi(items[4].c_str());
    int32_t data = (int32_t)atoi(items[5].c_str());
    read_type_answers->insert(std::pair<std::string, ReadType>(read_name, rt));
    data_answers->insert(std::pair<std::string, int32_t>(read_name, data));
  }
}

void ReadExtractorTest::test_ProcessReadPairs() {
  // Load Bam file
  std::string bam_file = test_dir + "/test.sorted.bam";
  std::vector<std::string> files(0);
  files.push_back(bam_file);
  BamCramMultiReader bamreader(files);
  std::map<std::string, ReadPair> read_pairs;
  // Test each pair
  if (!read_extractor_->ProcessReadPairs(&bamreader, locus, regionsize, min_match, &read_pairs)) {
    CPPUNIT_FAIL("ProcessReadPairs returned false unexpectedly.");
  }
  for (std::map<std::string, ReadPair>::const_iterator iter = read_pairs.begin();
       iter != read_pairs.end(); iter++) {
    std::stringstream msg;
    msg << "Misclassified " << iter->first
	<< " Found read type " << iter->second.read_type
	<< " Found data value " << iter->second.data_value;
    if (read_type_answers.find(iter->first) != read_type_answers.end()) {
      CPPUNIT_ASSERT_EQUAL_MESSAGE(msg.str(), read_type_answers[iter->first], iter->second.read_type);
      if (read_type_answers[iter->first] == RC_SPAN) {
	bool correct = abs(data_answers[iter->first]-iter->second.data_value) <= 2; // wiggle room for start coords?
	CPPUNIT_ASSERT_MESSAGE(msg.str(), correct);
      } else if (read_type_answers[iter->first] == RC_ENCL) {
	// for enclose, must be exactly right
	CPPUNIT_ASSERT_EQUAL_MESSAGE(msg.str(), data_answers[iter->first], iter->second.data_value);
      } else if (read_type_answers[iter->first] == RC_FRR) {
	bool correct = abs(data_answers[iter->first]-iter->second.data_value) <= 102; // account for read length bug + wiggle room
	CPPUNIT_ASSERT_MESSAGE(msg.str(), correct);
      } else {
	CPPUNIT_FAIL("Shouldn't get here");
      }
    }
  }
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
  files.push_back(discarded_file);
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
    test_is_discarded = read_extractor_->FindDiscardedRead(aln, chrom_ref_id, locus);
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
    test_is_spanning = read_extractor_->FindSpanningRead(aln, chrom_ref_id, locus, &insert_size);
    std::stringstream msg;
    msg << "Misclassified spanning read " << aln.Name();
    CPPUNIT_ASSERT_EQUAL_MESSAGE(msg.str(), truth_is_spanning, test_is_spanning);
  }
}

void ReadExtractorTest::test_ProcessSingleRead() {
  std::string spanning_file = test_dir + "/test.spanning.single.bam";
  std::string enclosing_file = test_dir + "/test.enclosing.single.bam";
  std::string frr_file = test_dir + "/test.frr.single.bam";
  std::vector<std::string> files(0);
  files.push_back(spanning_file);
  files.push_back(enclosing_file);
  files.push_back(frr_file);
  BamCramMultiReader bamreader(files);
  int32_t BUFFER = 50000; // want to capture all reads
  int32_t chrom_ref_id = bamreader.bam_header()->ref_id(locus.chrom);
  bamreader.SetRegion(locus.chrom, locus.start - BUFFER, locus.end + BUFFER);
  BamAlignment aln;
  ReadType read_type;
  SingleReadType srt;
  int32_t score_value, data_value, nCopy_value;
  bool truth_is_enclosing;
  while (bamreader.GetNextAlignment(aln)) {
    std::string fname = aln.Filename();
    truth_is_enclosing = string_ends_with(fname, "test.enclosing.single.bam");
    if (!read_extractor_->ProcessSingleRead(aln, chrom_ref_id, locus, min_match,
      &data_value, &nCopy_value, &score_value, &read_type, &srt)) {
      CPPUNIT_FAIL("ProcessSingleRead returned false unexpectedly");
    }
    std::stringstream msg;
    msg << "Misclassified discarded read " << aln.Name();
    if (truth_is_enclosing) {
      int64_t true_value;
      std::string tag_id = "nc";
      aln.GetIntTag(tag_id.c_str(), true_value);
      CPPUNIT_ASSERT_EQUAL_MESSAGE(msg.str(), RC_ENCL, read_type);
      CPPUNIT_ASSERT_EQUAL_MESSAGE(msg.str(), (int32_t)true_value, data_value);
    } else if (aln.Name() == "ATXN7_52_cov60_dist500_DIP_const70_20_altAllele_2920_3331_0:0:0_0:0:0_219" ||
	       aln.Name() == "ATXN7_52_cov60_dist500_DIP_const70_70_constAllele_2838_3270_0:0:0_0:0:0_3b") {
      int64_t true_is;
      std::string tag_is = "is";
      aln.GetIntTag(tag_is.c_str(), true_is);
      bool correct = abs((int32_t)true_is - data_value) <= 2; // allow some wiggle room for now
      CPPUNIT_ASSERT_EQUAL_MESSAGE(msg.str(), RC_SPAN, read_type);
      CPPUNIT_ASSERT_MESSAGE(msg.str(), correct);
    } 
  }
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
    read_extractor_->RescueMate(&bamreader, aln, &matealn);
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


