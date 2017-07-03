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

#ifndef SRC_TESTS_READCLASS_H__
#define SRC_TESTS_READCLASS_H__

#include <cppunit/extensions/HelperMacros.h>

#include "src/enclosing_class.h"
#include "src/spanning_class.h"
#include "src/frr_class.h"

class ReadClassTest: public CppUnit::TestFixture {
  CPPUNIT_TEST_SUITE(ReadClassTest);
  CPPUNIT_TEST(test_AddData);
  CPPUNIT_TEST(test_Reset);
  CPPUNIT_TEST(test_SpanClassProb);
  CPPUNIT_TEST(test_SpanReadProb);
  CPPUNIT_TEST(test_FRRClassProb);
  CPPUNIT_TEST(test_FRRReadProb);
  CPPUNIT_TEST(test_EnclosingClassProb);
  CPPUNIT_TEST(test_EnclosingReadProb);
  CPPUNIT_TEST(test_GetClassLogLikelihood);
  CPPUNIT_TEST(test_GetAlleleLogLikelihood);
  CPPUNIT_TEST_SUITE_END();

 public:
  void setUp();
  void tearDown();
  void test_AddData();
  void test_Reset();
  void test_SpanClassProb();
  void test_SpanReadProb();
  void test_FRRClassProb();
  void test_FRRReadProb();
  void test_EnclosingClassProb();
  void test_EnclosingReadProb();
  void test_GetClassLogLikelihood();
  void test_GetAlleleLogLikelihood();
 private:
  EnclosingClass encl_class_;
  SpanningClass span_class_;
  FRRClass frr_class_;
  int32_t read_len;
  int32_t motif_len;
  int32_t ref_count;
};

#endif //  SRC_TESTS_READCLASS_H_
