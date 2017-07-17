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

#ifndef SRC_TESTS_LIKELIHOODMAXIMIZER_H__
#define SRC_TESTS_LIKELIHOODMAXIMIZER_H__

#include <cppunit/extensions/HelperMacros.h>

#include "src/locus.h"
#include "src/options.h"
#include "src/likelihood_maximizer.h"
#include "src/read_extractor.h"
#include "src/ref_genome.h"
#include "src/genotyper.h"

class LikelihoodMaximizerTest: public CppUnit::TestFixture {
  CPPUNIT_TEST_SUITE(LikelihoodMaximizerTest);
  CPPUNIT_TEST(test_Reset);
  CPPUNIT_TEST(test_AddEnclosingData);
  CPPUNIT_TEST(test_AddSpanningData);
  CPPUNIT_TEST(test_AddFRRData);
  CPPUNIT_TEST(test_GetGenotypeNegLogLikelihood);
  CPPUNIT_TEST(test_OptimizeLikelihood);
  CPPUNIT_TEST_SUITE_END();

 public:
  void setUp();
  void tearDown();
  void test_Reset();
  void test_AddEnclosingData();
  void test_AddSpanningData();
  void test_AddFRRData();
  void test_GetGenotypeNegLogLikelihood();
  void test_OptimizeLikelihood();

  int read_len, motif_len, ref_count;
 private:
  LikelihoodMaximizer* likelihood_maximizer_;
  Locus locus;
  std::string test_dir;
  Options options;

};

#endif //  SRC_TESTS_LIKELIHOODMAXIMIZER_H_
