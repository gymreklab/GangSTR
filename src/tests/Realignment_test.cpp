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

#include "src/tests/Realignment_test.h"

#include <sstream>

// Registers the fixture into the 'registry'
CPPUNIT_TEST_SUITE_REGISTRATION(RealignmentTest);

void RealignmentTest::setUp() {}

void RealignmentTest::tearDown() {}

std::string RealignmentTest::ConstructSeq(const std::string& prefix,
			 const std::string& suffix,
			 const std::string& motif,
			 const int32_t& nCopy) {
  std::stringstream ss;
  ss << prefix;
  for (int i=0; i<nCopy; i++) {
    ss << motif;
  }
  ss << suffix;
  return ss.str();
}

void RealignmentTest::test_ExpansionAwareRealign() {
  std::string pre_flank = "ACTAGCTACTCATCCA";
  std::string post_flank = "ATCATCGACTACGACT";
  std::string motif = "CAG";
  std::string seq;
  int32_t nCopy, pos, score;
  // Case 1 - enclosing
  for (int32_t i=0; i<20; i++) {
    seq = ConstructSeq(pre_flank, post_flank, motif, i);
    if (!expansion_aware_realign(seq, pre_flank, post_flank, motif,
				 &nCopy, &pos, &score)) {
      CPPUNIT_FAIL("expansion_aware_realign returned false unexpectedly");
    }
    CPPUNIT_ASSERT_EQUAL(i, nCopy);
    CPPUNIT_ASSERT_EQUAL(0, pos);
    CPPUNIT_ASSERT_EQUAL((int32_t)seq.size()*MATCH_SCORE, score);
  }
  // Case 2 - preflank
  for (int32_t i=0; i<50; i++) {
    seq = ConstructSeq(pre_flank, "", motif, i);
    if (!expansion_aware_realign(seq, pre_flank, post_flank, motif,
				 &nCopy, &pos, &score)) {
      CPPUNIT_FAIL("expansion_aware_realign returned false unexpectedly");
    }
    CPPUNIT_ASSERT_EQUAL(i, nCopy);
    CPPUNIT_ASSERT_EQUAL(0, pos);
    CPPUNIT_ASSERT_EQUAL((int32_t)seq.size()*MATCH_SCORE, score);
  }
  // Case 3 - postflank
  for (int32_t i=0; i<50; i++) {
    seq = ConstructSeq("", post_flank, motif, i);
    if (!expansion_aware_realign(seq, pre_flank, post_flank, motif,
				 &nCopy, &pos, &score)) {
      CPPUNIT_FAIL("expansion_aware_realign returned false unexpectedly");
    }
    CPPUNIT_ASSERT_EQUAL(i, nCopy);
  }
  // Case 3 - IRR
  for (int32_t i=0; i<50; i++) {
    seq = ConstructSeq("", "", motif, i);
    if (!expansion_aware_realign(seq, pre_flank, post_flank, motif,
				 &nCopy, &pos, &score)) {
      CPPUNIT_FAIL("expansion_aware_realign returned false unexpectedly");
    }
    CPPUNIT_ASSERT_EQUAL(i, nCopy);
  }
  // Case 4 - random sequence
  seq ="NNNNNNNNNNN";
  if (!expansion_aware_realign(seq, pre_flank, post_flank, motif,
			       &nCopy, &pos, &score)) {
    CPPUNIT_FAIL("expansion_aware_realign returned false unexpectedly");
  }
  CPPUNIT_ASSERT_EQUAL(0, score);
}

void RealignmentTest::test_SmithWaterman() {
  std::string seq1, seq2;
  // Case 1
  seq1 = "ACGT";
  seq2 = "ACGT";
  int32_t pos, score;
  if (!smith_waterman(seq1, seq2, &pos, &score)) {
    CPPUNIT_FAIL("smith_waterman returned false unexpectedly");
  }
  CPPUNIT_ASSERT_EQUAL(MATCH_SCORE*4, score);
  CPPUNIT_ASSERT_EQUAL(0, pos);
  // Case 2
  seq1 = "ATCACGT";
  seq2 = "ATCGCGT";
  if (!smith_waterman(seq1, seq2, &pos, &score)) {
    CPPUNIT_FAIL("smith_waterman returned false unexpectedly");
  }
  CPPUNIT_ASSERT_EQUAL(MATCH_SCORE*6+MISMATCH_SCORE, score);
  CPPUNIT_ASSERT_EQUAL(0, pos);
  // Case 3
  seq1 = "ATCACGT";
  seq2 = "CACGT";
  if (!smith_waterman(seq1, seq2, &pos, &score)) {
    CPPUNIT_FAIL("smith_waterman returned false unexpectedly");
  }
  CPPUNIT_ASSERT_EQUAL(MATCH_SCORE*5, score);
  CPPUNIT_ASSERT_EQUAL(2, pos);
}

void RealignmentTest::test_CreateScoreMatrix() {
  int32_t current_score;
  int32_t start_pos;
  // Case 1 - trivial matching
  std::string seq1 = "ACT";
  std::string seq2 = "ACT";
  int32_t rows = (int32_t)seq1.size() + 1;
  int32_t cols = (int32_t)seq2.size() + 1;
  std::vector<std::vector<int32_t> > score_matrix;
  score_matrix.resize(rows, std::vector<int32_t>(cols, 0));
  if (!create_score_matrix(rows, cols, seq1, seq2,
			   &score_matrix, &start_pos, &current_score)) {
    CPPUNIT_FAIL("create_score_matrix failed unexpectedly");
  }
  CPPUNIT_ASSERT_EQUAL(MATCH_SCORE*3, current_score);
  CPPUNIT_ASSERT_EQUAL(3, start_pos);

  // Case 2 - almost trivial matching
  seq1 = "ACTG";
  seq2 = "ACT";
  rows = (int32_t)seq1.size() + 1;
  cols = (int32_t)seq2.size() + 1;
  std::vector<std::vector<int32_t> > score_matrix2;
  score_matrix2.resize(rows, std::vector<int32_t>(cols, 0));
  if (!create_score_matrix(rows, cols, seq1, seq2,
			   &score_matrix2, &start_pos, &current_score)) {
    CPPUNIT_FAIL("create_score_matrix failed unexpectedly");
  }
  CPPUNIT_ASSERT_EQUAL(MATCH_SCORE*3, current_score);
  CPPUNIT_ASSERT_EQUAL(3, start_pos);
  // Case 3
  seq1 = "ACCTGA";
  seq2 = "CCTG";
  rows = (int32_t)seq1.size() + 1;
  cols = (int32_t)seq2.size() + 1;
  std::vector<std::vector<int32_t> > score_matrix3;
  score_matrix3.resize(rows, std::vector<int32_t>(cols, 0));
  if (!create_score_matrix(rows, cols, seq1, seq2,
			   &score_matrix3, &start_pos, &current_score)) {
    CPPUNIT_FAIL("create_score_matrix failed unexpectedly");
  }
  CPPUNIT_ASSERT_EQUAL(MATCH_SCORE*4, current_score);
  CPPUNIT_ASSERT_EQUAL(5, start_pos);
  // Case 4
  seq1 = "ACCTGAT";
  seq2 = "ACTGAT";
  rows = (int32_t)seq1.size() + 1;
  cols = (int32_t)seq2.size() + 1;
  std::vector<std::vector<int32_t> > score_matrix4;
  score_matrix4.resize(rows, std::vector<int32_t>(cols, 0));
  if (!create_score_matrix(rows, cols, seq1, seq2,
			   &score_matrix4, &start_pos, &current_score)) {
    CPPUNIT_FAIL("create_score_matrix failed unexpectedly");
  }
  CPPUNIT_ASSERT_EQUAL(MATCH_SCORE*6+GAP_SCORE, current_score);
  CPPUNIT_ASSERT_EQUAL(7, start_pos);
}

void RealignmentTest::test_CalcScore() {
  std::string seq1 = "ATCT";
  std::string seq2 = "ACT";
  int32_t rows = (int32_t)seq1.size() + 1;
  int32_t cols = (int32_t)seq2.size() + 1;
  std::vector<std::vector<int32_t> > score_matrix;
  score_matrix.resize(rows, std::vector<int32_t>(cols, 0));
  // Walk through a couple steps of filling in the matrix
  int32_t i=1;
  int32_t j=1;
  if (!calc_score(i, j, seq1, seq2, &score_matrix)) {
    CPPUNIT_FAIL("calc_score returned false unexpectedly");
  }
  CPPUNIT_ASSERT_EQUAL(MATCH_SCORE, score_matrix.at(i).at(j));
  i=2;
  j=1;
  if (!calc_score(i, j, seq1, seq2, &score_matrix)) {
    CPPUNIT_FAIL("calc_score returned false unexpectedly");
  }
  CPPUNIT_ASSERT_EQUAL(0, score_matrix.at(i).at(j));
  i=2;
  j=2;
  if (!calc_score(i, j, seq1, seq2, &score_matrix)) {
    CPPUNIT_FAIL("calc_score returned false unexpectedly");
  }
  CPPUNIT_ASSERT_EQUAL(2, score_matrix.at(i).at(j));
  i=3;
  j=2;
  if (!calc_score(i, j, seq1, seq2, &score_matrix)) {
    CPPUNIT_FAIL("calc_score returned false unexpectedly");
  }
  CPPUNIT_ASSERT_EQUAL(3, score_matrix.at(i).at(j));
}

void RealignmentTest::test_ClassifyRealignedRead() {
  std::string pre_flank = "ACTAGCTACTCATCCA";
  std::string post_flank = "ATCATCGACTACGACT";
  std::string motif = "CAG";
  std::string seq;
  int32_t nCopy, pos, score;
  SingleReadType src;
  // Case 1 - enclosing
  for (int32_t i=1; i<20; i++) {
    seq = ConstructSeq(pre_flank, post_flank, motif, i);
    expansion_aware_realign(seq, pre_flank, post_flank, motif,
			    &nCopy, &pos, &score);
    classify_realigned_read(seq, motif, pos, nCopy, score, (int32_t)pre_flank.size(),
			    &src);
    CPPUNIT_ASSERT_EQUAL(SR_ENCLOSING, src);
  }
  // Case 2 - preflank
  for (int32_t i=1; i<20; i++) {
    seq = ConstructSeq(pre_flank, "", motif, i);
    expansion_aware_realign(seq, pre_flank, post_flank, motif,
			    &nCopy, &pos, &score);
    classify_realigned_read(seq, motif, pos, nCopy, score, (int32_t)pre_flank.size(),
			    &src);
    CPPUNIT_ASSERT_EQUAL(SR_PREFLANK, src);
  }
  // Case 3 - postflank
  for (int32_t i=1; i<20; i++) {
    seq = ConstructSeq("", post_flank, motif, i);
    expansion_aware_realign(seq, pre_flank, post_flank, motif,
			    &nCopy, &pos, &score);
    classify_realigned_read(seq, motif, pos, nCopy, score, (int32_t)pre_flank.size(),
			    &src);
    CPPUNIT_ASSERT_EQUAL(SR_POSTFLANK, src);
  }
  // Case 4 - IRR
  for (int32_t i=1; i<20; i++) {
    seq = ConstructSeq("", "", motif, i);
    expansion_aware_realign(seq, pre_flank, post_flank, motif,
			    &nCopy, &pos, &score);
    classify_realigned_read(seq, motif, pos, nCopy, score, (int32_t)pre_flank.size(),
			    &src);
    CPPUNIT_ASSERT_EQUAL(SR_IRR, src);
  }
  // Case 5 - unknown
  seq = "NNNNNNNNNNN";
  expansion_aware_realign(seq, pre_flank, post_flank, motif,
			  &nCopy, &pos, &score);
  classify_realigned_read(seq, motif, pos, nCopy, score, (int32_t)pre_flank.size(),
			  &src);
  CPPUNIT_ASSERT_EQUAL(SR_UNKNOWN, src);
  seq = "ATACGTACGATCTACAG";
  expansion_aware_realign(seq, pre_flank, post_flank, motif,
			  &nCopy, &pos, &score);
  classify_realigned_read(seq, motif, pos, nCopy, score, (int32_t)pre_flank.size(),
			  &src);
  CPPUNIT_ASSERT_EQUAL(SR_UNKNOWN, src);
}


