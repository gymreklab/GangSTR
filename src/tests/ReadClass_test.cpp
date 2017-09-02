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

#include "src/tests/ReadClass_test.h"
#include <math.h>
#include <iostream>
using namespace std;
// Registers the fixture into the 'registry'
CPPUNIT_TEST_SUITE_REGISTRATION(ReadClassTest);

void ReadClassTest::setUp() {
  Options options;

  options.dist_mean = 400;
  options.dist_sdev = 50;
  options.stutter_up = 0.01;
  options.stutter_down = 0.02;
  options.stutter_p = 0.95;
  options.flanklen = 2000;
  options.realignment_flanklen = 100;
  options.frr_weight = 0.8;
  options.enclosing_weight = 1.0;
  options.spanning_weight = 1.0;
  options.verbose = false;

  encl_class_.SetOptions(options);
  span_class_.SetOptions(options);
  frr_class_.SetOptions(options);
  encl_class_.Reset();
  span_class_.Reset();
  frr_class_.Reset();
  read_len = 100;
  motif_len = 3;
  ref_count = 10;
  ploidy = 2;
}

void ReadClassTest::tearDown() {}

void ReadClassTest::test_AddData() {
  int32_t test_data1 = 10;
  int32_t test_data2 = 20;
  int32_t test_data3 = 30;
  encl_class_.AddData(test_data1);
  encl_class_.AddData(test_data2);
  span_class_.AddData(test_data1);
  span_class_.AddData(test_data2);
  span_class_.AddData(test_data3);
  frr_class_.AddData(test_data1);
  CPPUNIT_ASSERT_EQUAL((int)encl_class_.GetDataSize(), 2);
  CPPUNIT_ASSERT_EQUAL((int)span_class_.GetDataSize(), 3);
  CPPUNIT_ASSERT_EQUAL((int)frr_class_.GetDataSize(), 1);
}

void ReadClassTest::test_Reset() {
  int32_t test_data1 = 10;
  int32_t test_data2 = 20;
  int32_t test_data3 = 30;
  encl_class_.AddData(test_data1);
  encl_class_.AddData(test_data2);
  span_class_.AddData(test_data1);
  span_class_.AddData(test_data2);
  span_class_.AddData(test_data3);
  frr_class_.AddData(test_data1);
  encl_class_.Reset();
  span_class_.Reset();
  frr_class_.Reset();
  CPPUNIT_ASSERT_EQUAL((int)encl_class_.GetDataSize(), 0);
  CPPUNIT_ASSERT_EQUAL((int)span_class_.GetDataSize(), 0);
  CPPUNIT_ASSERT_EQUAL((int)frr_class_.GetDataSize(), 0);
}
// NOTE:
// exp: ATXN7_18_class2_cov50_dist400
void ReadClassTest::test_SpanClassProb() {
  // Example - TODO change
  int32_t allele = 25;
  double log_class_prob = 0.5;
  span_class_.GetLogClassProb(allele, read_len, motif_len, &log_class_prob);
  CPPUNIT_ASSERT_EQUAL(roundf(log_class_prob*pow(10,13))/pow(10,13), roundf(log(0.0838726946383)*pow(10,13))/pow(10,13));
  // CPPUNIT_FAIL("test_SpanClassProb not implemented");
}

void ReadClassTest::test_SpanReadProb() {
  // Example - TODO change
  int32_t allele = 25;
  double log_allele_prob = 0.5;
  int32_t data = 450;
  span_class_.GetLogReadProb(allele, data, read_len, motif_len, ref_count, &log_allele_prob);
  CPPUNIT_ASSERT_EQUAL(roundf(log_allele_prob*pow(10,13))/pow(10,13), roundf(log(0.00131231629549)*pow(10,13))/pow(10,13));
  // CPPUNIT_FAIL("test_SpanClassProb not implemented");
}

void ReadClassTest::test_FRRClassProb() {
  // Example - TODO change
  int32_t allele = 45;
  double log_class_prob = 0.5;
  frr_class_.GetLogClassProb(allele, read_len, motif_len, &log_class_prob);
  CPPUNIT_ASSERT_EQUAL(roundf(log_class_prob*pow(10,13))/pow(10,13), roundf(log(0.0177896348168/2)*pow(10,13))/pow(10,13));
  // CPPUNIT_FAIL("test_FRRClassProb not implemented");
}

void ReadClassTest::test_FRRReadProb() {
  // Example - TODO change
  int32_t allele = 45;
  double log_allele_prob = 0.5;
  int32_t data = 80;
  frr_class_.GetLogReadProb(allele, data, read_len, motif_len, ref_count, &log_allele_prob);
  CPPUNIT_ASSERT_EQUAL(roundf(log_allele_prob*pow(10,13))/pow(10,13), roundf(log(0.0363690786878)*pow(10,13))/pow(10,13));
  // CPPUNIT_FAIL("test_FRRClassProb not implemented");
}

void ReadClassTest::test_EnclosingClassProb() {
  // Example - TODO change
  int32_t allele = 25;
  double log_class_prob = 0.5;
  encl_class_.GetLogClassProb(allele, read_len, motif_len, &log_class_prob);
  CPPUNIT_ASSERT_EQUAL(roundf(log_class_prob*pow(10,13))/pow(10,13), roundf(log(0.0129032258065/2)*pow(10,13))/pow(10,13));
  // CPPUNIT_FAIL("test_EnclosingClassProb not implemented");
}

void ReadClassTest::test_EnclosingReadProb() {
  // Example - TODO change
  int32_t allele = 25;
  double log_allele_prob = 0.5;
  int32_t data = 25;
  encl_class_.GetLogReadProb(allele, data, read_len, motif_len, ref_count, &log_allele_prob);
  CPPUNIT_ASSERT_EQUAL(roundf(log_allele_prob*pow(10,13))/pow(10,13), roundf(log(0.97)*pow(10,13))/pow(10,13));
  // CPPUNIT_FAIL("test_EnclosingClassProb not implemented");
}

void ReadClassTest::test_GetClassLogLikelihood() {
  int32_t test_data1 = 10;
  int32_t test_data2 = 20;
  int32_t test_data3 = 30;
  int32_t test_data4 = 40;
  int32_t allele1 = 20;
  int32_t allele2 = 50;
  double class_ll;
  encl_class_.Reset();
  encl_class_.AddData(test_data1);
  encl_class_.AddData(test_data2);
  encl_class_.AddData(test_data3);
  encl_class_.AddData(test_data4);
  encl_class_.GetClassLogLikelihood(allele1, allele2, read_len, motif_len, ref_count, ploidy, &class_ll);
  CPPUNIT_ASSERT_EQUAL(roundf(class_ll*pow(10,11))/pow(10,11), roundf(-145.199557345*pow(10,11))/pow(10,11));

  span_class_.Reset();
  span_class_.AddData(test_data2);
  span_class_.AddData(test_data3);
  span_class_.AddData(test_data4);
  // double allele_ll;
  // cout<<frr_class_.GetLogReadProb(allele1, test_data1, &allele_ll);
  // cout<<endl<<allele_ll<<endl;

  span_class_.GetClassLogLikelihood(allele1, allele2, read_len, motif_len, ref_count, ploidy, &class_ll);
  CPPUNIT_ASSERT_EQUAL(roundf(class_ll*pow(10,4))/pow(10,4), roundf(-62.3921692604*pow(10,4))/pow(10,4));

  frr_class_.Reset();
  frr_class_.AddData(test_data1);
  frr_class_.AddData(test_data2);
  frr_class_.AddData(test_data3);
  frr_class_.AddData(test_data4);
  frr_class_.GetClassLogLikelihood(allele1, allele2, read_len, motif_len, ref_count, ploidy, &class_ll);
  CPPUNIT_ASSERT_EQUAL(roundf(class_ll*pow(10,11))/pow(10,11), roundf(-40.8239143859*pow(10,11))/pow(10,11));

  // std::cout<<std::endl<<class_ll<<std::endl;
  // CPPUNIT_FAIL("test_GetClassLogLikelihood not implemented");
}

void ReadClassTest::test_GetAlleleLogLikelihood() {
  // Example - TODO change
  int32_t allele = 25;
  double allele_ll = 0.5;
  int32_t data = 24;
  encl_class_.GetAlleleLogLikelihood(allele, data, read_len, motif_len, ref_count, &allele_ll);
  CPPUNIT_ASSERT_EQUAL(roundf(allele_ll*pow(10,11))/pow(10,11), roundf(-9.00674141673*pow(10,11))/pow(10,11));
  allele = 55;
  data = 430;
  span_class_.GetAlleleLogLikelihood(allele, data, read_len, motif_len, ref_count, &allele_ll);
  CPPUNIT_ASSERT_EQUAL(roundf(allele_ll*pow(10,11))/pow(10,11), roundf(-13.1016086835*pow(10,11))/pow(10,11));
  allele = 55;
  data = 80;
  frr_class_.GetAlleleLogLikelihood(allele, data, read_len, motif_len, ref_count, &allele_ll);
  CPPUNIT_ASSERT_EQUAL(roundf(allele_ll*pow(10,11))/pow(10,11), roundf(-6.17069762893*pow(10,11))/pow(10,11));
  // CPPUNIT_FAIL("test_GetAlleleLogLikelihood not implemented");
}

