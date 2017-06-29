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

#ifndef SRC_TESTS_READEXTRACTOR_H__
#define SRC_TESTS_READEXTRACTOR_H__

#include <cppunit/extensions/HelperMacros.h>

#include "src/locus.h"
#include "src/read_extractor.h"

#include <string>

class ReadExtractorTest: public CppUnit::TestFixture {
  CPPUNIT_TEST_SUITE(ReadExtractorTest);
  CPPUNIT_TEST(test_ExtractReads);
  CPPUNIT_TEST(test_ProcessReadPairs);
  CPPUNIT_TEST(test_FindDiscardedRead);
  CPPUNIT_TEST(test_FindSpanningRead);
  CPPUNIT_TEST(test_ProcessSingleRead);
  CPPUNIT_TEST(test_RescueMate);
  CPPUNIT_TEST_SUITE_END();

 public:
  void setUp();
  void tearDown();
 private:
  void test_ExtractReads();
  void test_ProcessReadPairs();
  void test_FindDiscardedRead();
  void test_FindSpanningRead();
  void test_ProcessSingleRead();
  void test_RescueMate();
  void LoadAnswers(const std::string& answers_file,
		   std::map<std::string, ReadType>* read_type_answers,
		   std::map<std::string, int32_t>* data_answers);
  std::string test_dir;
  ReadExtractor read_extractor;
  Locus locus;
};

#endif //  SRC_TESTS_READEXTRACTOR_H_
