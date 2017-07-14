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

#include "src/tests/LikelihoodMaximizer_test.h"

#include <math.h>
// Registers the fixture into the 'registry'
CPPUNIT_TEST_SUITE_REGISTRATION(LikelihoodMaximizerTest);

void LikelihoodMaximizerTest::setUp() {
  Options options;
  options.dist_mean = 400;
  options.dist_sdev = 50;
  options.stutter_up = 0.01;
  options.stutter_down = 0.02;
  options.stutter_p = 0.95;
  options.flanklen = 2000;
  options.realignment_flanklen = 100;
  options.frr_weight = 1.0;
  options.enclosing_weight = 1.0;
  options.spanning_weight = 1.0;
  options.verbose = false;
  likelihood_maximizer_ = new LikelihoodMaximizer(options);
  likelihood_maximizer_->Reset();
  read_len = 100;
  motif_len = 3;
  ref_count = 10;
}

void LikelihoodMaximizerTest::tearDown() {}

void LikelihoodMaximizerTest::test_Reset() {
  int32_t test_data1 = 10;
  int32_t test_data2 = 20;
  int32_t test_data3 = 30;
  likelihood_maximizer_->Reset();
  likelihood_maximizer_->AddEnclosingData(test_data1);
  likelihood_maximizer_->AddSpanningData(test_data2);
  likelihood_maximizer_->AddFRRData(test_data3);
  likelihood_maximizer_->Reset();
  CPPUNIT_ASSERT_EQUAL((int)likelihood_maximizer_->GetEnclosingDataSize(), 0);
  CPPUNIT_ASSERT_EQUAL((int)likelihood_maximizer_->GetSpanningDataSize(), 0);
  CPPUNIT_ASSERT_EQUAL((int)likelihood_maximizer_->GetFRRDataSize(), 0);
}

void LikelihoodMaximizerTest::test_AddEnclosingData() {
  int32_t test_data1 = 10;
  int32_t test_data2 = 20;
  int32_t test_data3 = 30;
  likelihood_maximizer_->Reset();
  likelihood_maximizer_->AddEnclosingData(test_data1);
  likelihood_maximizer_->AddEnclosingData(test_data2);
  likelihood_maximizer_->AddEnclosingData(test_data3);
  CPPUNIT_ASSERT_EQUAL((int)likelihood_maximizer_->GetEnclosingDataSize(), 3);
  CPPUNIT_ASSERT_EQUAL((int)likelihood_maximizer_->GetSpanningDataSize(), 0);
  CPPUNIT_ASSERT_EQUAL((int)likelihood_maximizer_->GetFRRDataSize(), 0);
}

void LikelihoodMaximizerTest::test_AddSpanningData() {
  int32_t test_data1 = 10;
  likelihood_maximizer_->Reset();
  likelihood_maximizer_->AddSpanningData(test_data1);
  CPPUNIT_ASSERT_EQUAL((int)likelihood_maximizer_->GetSpanningDataSize(), 1);
  CPPUNIT_ASSERT_EQUAL((int)likelihood_maximizer_->GetEnclosingDataSize(), 0);
  CPPUNIT_ASSERT_EQUAL((int)likelihood_maximizer_->GetFRRDataSize(), 0);
}

void LikelihoodMaximizerTest::test_AddFRRData() {
  int32_t test_data1 = 10;
  int32_t test_data2 = 20;
  likelihood_maximizer_->Reset();
  likelihood_maximizer_->AddFRRData(test_data1);
  likelihood_maximizer_->AddFRRData(test_data2);
  CPPUNIT_ASSERT_EQUAL((int)likelihood_maximizer_->GetFRRDataSize(), 2);
  CPPUNIT_ASSERT_EQUAL((int)likelihood_maximizer_->GetEnclosingDataSize(), 0);
  CPPUNIT_ASSERT_EQUAL((int)likelihood_maximizer_->GetSpanningDataSize(), 0);
}

void LikelihoodMaximizerTest::test_GetGenotypeNegLogLikelihood() {
  int32_t allele1 = 10, allele2 = 30;
  int32_t test_data1 = 10;
  int32_t test_data2 = 380;
  int32_t test_data3 = 80;
  likelihood_maximizer_->Reset();
  likelihood_maximizer_->AddEnclosingData(test_data1);
  likelihood_maximizer_->AddSpanningData(test_data2);
  // likelihood_maximizer_->AddFRRData(test_data3);   // TODO add an example with all reads

  double gt_ll;
  if (!likelihood_maximizer_->GetGenotypeNegLogLikelihood(allele1, allele2,
          read_len, motif_len, ref_count, &gt_ll)){
    CPPUNIT_FAIL( "Running GetGenotypeNegLogLikelihood failed." );
  }

  // double gt_ll2;
  // allele1 = 45;
  // allele2 = 65;
  // likelihood_maximizer_->Reset();
  // likelihood_maximizer_->AddSpanningData(test_data2);
  // likelihood_maximizer_->AddFRRData(test_data3);
  // if (!likelihood_maximizer_->GetGenotypeNegLogLikelihood(allele1, allele2,
  //         read_len, motif_len, ref_count, &gt_ll2)){
  //   CPPUNIT_FAIL( "Running GetGenotypeNegLogLikelihood failed." );
  // }
  
  
  CPPUNIT_ASSERT_EQUAL(roundf(gt_ll * 1000)/1000, roundf(12.1668285343*1000)/1000);
  // CPPUNIT_ASSERT_EQUAL(roundf(gt_ll2 * 100)/100, roundf(15.2104894117*100)/100);

}


