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

#include "src/bam_io.h"
#include <math.h>

#include <iostream>
using namespace std;

// Registers the fixture into the 'registry'
CPPUNIT_TEST_SUITE_REGISTRATION(LikelihoodMaximizerTest);

void LikelihoodMaximizerTest::setUp() {
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
  options.min_match = 0;
  options.read_len = 100;
  likelihood_maximizer_ = new LikelihoodMaximizer(options);
  likelihood_maximizer_->Reset();
  read_len = 100;
  motif_len = 3;
  ref_count = 10;
  resampled = false;


  test_dir = getenv("GANGSTR_TEST_DIR");
  locus.chrom = "19";
  locus.start = 5000;
  locus.end = 5039;
  locus.motif = "CTG";
  locus.period = 3;
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
  /*
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
          read_len, motif_len, ref_count, resampled, &gt_ll)){
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
  */
}

void LikelihoodMaximizerTest::test_OptimizeLikelihood() {
  options.dist_mean = 500;
  options.dist_sdev = 50;
  options.flanklen = 3000;
  options.frr_weight = 0.5;
  options.enclosing_weight = 1.0;
  options.spanning_weight = 1.0;
  options.flanking_weight = 1.0;
  options.read_len = 100;
  options.dist_max = 1000;

  LikelihoodMaximizer* likelihood_maximizer_opt = new LikelihoodMaximizer(options);
  likelihood_maximizer_opt->Reset();

  std::string fastafile = test_dir + "/CACNA1A_5k_region.fa";
  RefGenome refgenome(fastafile);

  ReadExtractor* read_extractor = new ReadExtractor(options);

  std::string bam_file = test_dir + "/54_nc_40.sorted.bam";
  std::vector<std::string> files(0);
  files.push_back(bam_file);
  BamCramMultiReader* bamreader = new BamCramMultiReader(files, fastafile);

  // Load preflank and postflank to locus
  if (!refgenome.GetSequence(locus.chrom,
            locus.start-options.realignment_flanklen-1,
            locus.start-2,
            &locus.pre_flank)) {
    CPPUNIT_FAIL( "Running OptimizeLikelihood failed." );
  }
  if (!refgenome.GetSequence(locus.chrom,
            locus.end,
            locus.end+options.realignment_flanklen-1,
            &locus.post_flank)) {
    CPPUNIT_FAIL( "Running OptimizeLikelihood failed." );
  }


  // Load all read data
  likelihood_maximizer_opt->Reset();
  if (!read_extractor->ExtractReads(bamreader, locus, options.regionsize,
            options.min_match, likelihood_maximizer_opt)) {
    CPPUNIT_FAIL( "Running OptimizeLikelihood failed." );
  }
  // Maximize the likelihood
  int32_t allele1, allele2;
  int32_t read_len = options.read_len;

  int32_t ref_count = (int32_t)((locus.end-locus.start+1)/locus.motif.size());
  double min_negLike;
  // TODO remake test case using new optimizer function
  /* 
  if (!likelihood_maximizer_opt->OptimizeLikelihood(read_len, (int32_t)(locus.motif.size()),
            ref_count, resampled,
            &allele1, &allele2, &min_negLike)) {
    CPPUNIT_FAIL( "Running OptimizeLikelihood failed." );
    }*/
  //  CPPUNIT_ASSERT_EQUAL(allele1, 31);
  //CPPUNIT_ASSERT_EQUAL(allele2, 60);
  //CPPUNIT_ASSERT_EQUAL(roundf(min_negLike * 100)/100, roundf(1725.53*100)/100); 
}


