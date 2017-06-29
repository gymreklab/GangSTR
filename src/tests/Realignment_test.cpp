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

// Registers the fixture into the 'registry'
CPPUNIT_TEST_SUITE_REGISTRATION(RealignmentTest);

void RealignmentTest::setUp() {}

void RealignmentTest::tearDown() {}

void RealignmentTest::test_ExpansionAwareRealign() {
  CPPUNIT_FAIL("test_ExpansionAwareRealign() not implemented");
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
  CPPUNIT_FAIL("test_ClassifyRealignedRead() not implemented");
}


