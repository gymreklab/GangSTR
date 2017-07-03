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

// Registers the fixture into the 'registry'
CPPUNIT_TEST_SUITE_REGISTRATION(LikelihoodMaximizerTest);

void LikelihoodMaximizerTest::setUp() {
  Options options;
  likelihood_maximizer_ = new LikelihoodMaximizer(options);
  likelihood_maximizer_->Reset();
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
  CPPUNIT_FAIL( "test_GetGenotypeNegLogLikelihood not implemented" );
}


