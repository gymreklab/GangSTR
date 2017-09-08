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

#include <iostream>

#include "src/genotyper.h"
#include "src/mathops.h"

using namespace std;

Genotyper::Genotyper(RefGenome _refgenome,
		    Options& _options) {
  refgenome = &_refgenome;
  options = &_options;
  read_extractor = new ReadExtractor();
  likelihood_maximizer = new LikelihoodMaximizer(_options);
}

bool Genotyper::SetFlanks(Locus* locus) {
  if (!refgenome->GetSequence(locus->chrom,
			      locus->start-options->realignment_flanklen-1,
			      locus->start-2,
			      &locus->pre_flank)) {
    return false;
  }
  if (!refgenome->GetSequence(locus->chrom,
			      locus->end,
			      locus->end+options->realignment_flanklen-1,
			      &locus->post_flank)) {
    return false;
  }
  return true;
}

bool Genotyper::ProcessLocus(BamCramMultiReader* bamreader, Locus* locus) {
  // Load preflank and postflank to locus
  if (options->verbose) {
    PrintMessageDieOnError("\tSetting flanking regions", M_PROGRESS);
  }
  if (!SetFlanks(locus)) {
    return false;
  }

  // Compute insert size distribution
  if (options->verbose) {
    PrintMessageDieOnError("\tComputing insert size distribution", M_PROGRESS);
  }
  double mean, std_dev;
  likelihood_maximizer->Reset();
  if (!read_extractor->ComputeInsertSizeDistribution(bamreader, *locus,
            &mean, &std_dev)) {
    return false;
  }
  // TODO upgrade ComputeInsertSizeDistribution to bwa mem edition
  options->dist_mean = mean;
  options->dist_sdev = std_dev;
  likelihood_maximizer->UpdateOptions();
  if (options->verbose) {
    stringstream ss;
    ss << "\t\tMean=" << mean << " SD=" << std_dev;
    PrintMessageDieOnError(ss.str(), M_PROGRESS);
  }

  // Load all read data
  if (options->verbose) {
    PrintMessageDieOnError("\tLoading read data", M_PROGRESS);
  }
  if (!read_extractor->ExtractReads(bamreader, *locus, likelihood_maximizer->options->regionsize,
				    likelihood_maximizer)) {
    return false;
  }
  
  // Maximize the likelihood
  if (options->verbose) {
    PrintMessageDieOnError("\tMaximizing likelihood", M_PROGRESS);
  }
  int32_t allele1, allele2;
  int32_t read_len = read_extractor->guessed_read_length;
  int32_t ref_count = (int32_t)((locus->end-locus->start+1)/locus->motif.size());
  double min_negLike, lob1, lob2, hib1, hib2;
  bool resampled = false;

  if (!likelihood_maximizer->OptimizeLikelihood(read_len, (int32_t)(locus->motif.size()),
						ref_count, resampled,
						&allele1, &allele2, &min_negLike)) {
    return false;
  }
  cout<<">>Genotyper Results:\t"<<allele1<<", "<<allele2<<"\tlikelihood = "<<min_negLike<<"\n";

  // for (int jj = 0; jj < 10; jj++){
  //   likelihood_maximizer->ResampleReadPool();
  //   resampled = true;
  //   if (!likelihood_maximizer->OptimizeLikelihood(read_len, (int32_t)(locus->motif.size()),
  //             ref_count, resampled,
  //             &allele1, &allele2, &min_negLike)) {
  //     return false;
  //   }
  //   cout<<">>Resampled Results:\t"<<allele1<<", "<<allele2<<"\tlikelihood = "<<min_negLike<<"\n";
  // }
  if (!likelihood_maximizer->GetConfidenceInterval(read_len, (int32_t)(locus->motif.size()),
            ref_count, allele1, allele2,
            &lob1, &hib1, &lob2, &hib2)) {
    return false;
  }
  cout<<"@@Small Allele Bound:\t["<<allele1 - lob1<<", "<<allele1 - hib1<<"]\n";
  cout<<"@@Large Allele Bound:\t["<<allele2 - lob2<<", "<<allele2 - hib2<<"]\n";

  return true;
}

void Genotyper::Debug(BamCramMultiReader* bamreader) {
  cerr << "testing refgenome" << endl;
  std::string seq;
  refgenome->GetSequence("3", 63898261, 63898360, &seq);
  cerr << seq << endl;
  cerr << "testing bam" << endl;
  bamreader->SetRegion("1", 0, 10000);
  BamAlignment aln;
  if (bamreader->GetNextAlignment(aln)) { // Requires SetRegion was called
    std::string testread = aln.QueryBases();
    cerr << testread << endl;
  } else {
    cerr << "testing bam failed" << endl;
  }
  cerr << "testing GSL" << endl;
  double x = TestGSL();
  cerr << "gsl_ran_gaussian_pdf(0, 1) " << x << endl;
  double y = TestNLOPT();
}

Genotyper::~Genotyper() {
  delete read_extractor;
  delete likelihood_maximizer;
}
