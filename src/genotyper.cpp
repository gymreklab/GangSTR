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

#include <set>
using namespace std;

Genotyper::Genotyper(RefGenome& _refgenome,
		     Options& _options,
		     SampleInfo& _sample_info) {
  refgenome = &_refgenome;
  options = &_options;
  sample_info = &_sample_info;
  read_extractor = new ReadExtractor(_options, *sample_info);
  std::set<std::string> rg_samples = sample_info->GetSamples();
  for (std::set<std::string>::iterator it=rg_samples.begin();
       it != rg_samples.end(); it++) {
    sample_likelihood_maximizers[*it] = new LikelihoodMaximizer(_options, *sample_info, *it);
  }
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
  int32_t read_len = sample_info->GetReadLength();

  // Load preflank and postflank to locus
  if (options->verbose) {
    PrintMessageDieOnError("\tSetting flanking regions", M_PROGRESS);
  }
  if (!SetFlanks(locus)) {
    return false;
  }

  for (std::map<std::string, LikelihoodMaximizer*>::iterator it = sample_likelihood_maximizers.begin();
       it != sample_likelihood_maximizers.end(); it++) {
    it->second->Reset();
  }

  // Infer GC bin
  int gcbin = -1;
  if (options->model_gc_cov) {
    std::string seq;
    if (refgenome->GetSequence(locus->chrom,
			       int((locus->start+locus->end)/2-options->gc_region_len/2),
			       int((locus->start+locus->end)/2+options->gc_region_len/2),
			       &seq)) {
      float gc = GetGC(seq);
      gcbin = int(floor(gc/options->gc_bin_size));
    }
  }
  // Load all read data
  if (options->verbose) {
    PrintMessageDieOnError("\tLoading read data", M_PROGRESS);
  }
  if (!read_extractor->ExtractReads(bamreader, *locus, options->regionsize,
				    options->min_match, sample_likelihood_maximizers)) {
    return false;
  }


  // Maximize the likelihood
  if (options->verbose) {
    PrintMessageDieOnError("\tMaximizing likelihood", M_PROGRESS);
  }
  int32_t allele1, allele2;
  int32_t ref_count = (int32_t)((locus->end-locus->start+1)/locus->motif.size());
  double min_negLike, lob1, lob2, hib1, hib2;
  double a1_se, a2_se;
  bool resampled = false;
  std::set<std::string> rg_samples = sample_info->GetSamples();
  // First set grid size
  int32_t sample_min_allele, sample_max_allele;
  int32_t min_allele = 100000;
  int32_t max_allele = 0;
  for (std::set<std::string>::iterator it = rg_samples.begin();
       it != rg_samples.end(); it++) {
    const std::string samp = *it;
    if (!sample_likelihood_maximizers[samp]->InferGridSize(read_len, (int32_t)(locus->motif.size()))) {
      PrintMessageDieOnError("Error inferring grid size", M_PROGRESS);
      return false;
    }
    sample_likelihood_maximizers[samp]->GetGridSize(&sample_min_allele, &sample_max_allele);
    if (sample_max_allele > max_allele) max_allele = sample_max_allele;
    if (sample_min_allele < min_allele) min_allele = sample_min_allele;
  }
  locus->grid_min_allele = min_allele;
  locus->grid_max_allele = max_allele;
  // Then max likelihood for each sample
  for (std::set<std::string>::iterator it = rg_samples.begin();
       it != rg_samples.end(); it++) {
    const std::string samp = *it;
    sample_likelihood_maximizers[samp]->SetGridSize(min_allele, max_allele);
    if (gcbin != -1) {
      sample_likelihood_maximizers[samp]->SetCoverage(sample_info->GetGCCoverage(samp)[gcbin]);
    }
    try {
      if (!sample_likelihood_maximizers[samp]->OptimizeLikelihood(read_len, 
								  (int32_t)(locus->motif.size()),
								  ref_count, 
								  resampled, 
								  options->ploidy, 
								  0,
								  locus->offtarget_share,
								  &allele1, 
								  &allele2, 
								  &min_negLike)) {
	locus->called[samp] = false;
      }
      locus->allele1[samp] = allele1;
      locus->allele2[samp] = allele2;
      locus->min_neg_lik[samp] = min_negLike;
      locus->enclosing_reads[samp] = sample_likelihood_maximizers[samp]->GetEnclosingDataSize();
      locus->spanning_reads[samp] = sample_likelihood_maximizers[samp]->GetSpanningDataSize();
      locus->frr_reads[samp] = sample_likelihood_maximizers[samp]->GetFRRDataSize();
      locus->flanking_reads[samp] = sample_likelihood_maximizers[samp]->GetFlankingDataSize();
      locus->depth[samp] = sample_likelihood_maximizers[samp]->GetReadPoolSize();
      
      if (options->num_boot_samp > 0){
	if (options->verbose) {
	  PrintMessageDieOnError("\tGetting confidence intervals", M_PROGRESS);
	}
	try{
	  if (!sample_likelihood_maximizers[samp]->GetConfidenceInterval(read_len, (int32_t)(locus->motif.size()),
									 ref_count, allele1, allele2, *locus,
									 &lob1, &hib1, &lob2, &hib2, &a1_se, &a2_se)) {
	    locus->called[samp] = false;
	  }
	  locus->lob1[samp] = lob1;
	  locus->lob2[samp] = lob2;
	  locus->hib1[samp] = hib1;
	  locus->hib2[samp] = hib2;
	  locus->a1_se[samp] = a1_se;
	  locus->a2_se[samp] = a2_se;
	  
	  stringstream msg;
	  msg<<"\tGenotyper Results:  "<<allele1<<", "<<allele2<<"\tlikelihood = "<<min_negLike;
	  PrintMessageDieOnError(msg.str(), M_PROGRESS);
	  if (options->verbose) {
	    msg.clear();
	    msg.str(std::string());
	    msg<<"\tSmall Allele Bound: ["<<lob1<<", "<<hib1<<"]";
	    PrintMessageDieOnError(msg.str(), M_PROGRESS);
	    msg.clear();
	    msg.str(std::string());
	    msg<<"\tLarge Allele Bound: ["<<lob2<<", "<<hib2<<"]";
	    PrintMessageDieOnError(msg.str(), M_PROGRESS);
	  }
	}
	catch (std::exception &exc){
	  if (options->verbose) {
	    stringstream msg;
	    msg<<"\tEncountered error("<< exc.what() <<") in likelihood maximization. Skipping locus";
	    PrintMessageDieOnError(msg.str(), M_PROGRESS);
	  }
	  locus->called[samp] = false;
	}
      }
    }
    catch (std::exception &exc){
      if (options->verbose) {
	stringstream msg;
	msg<<"\tEncountered error("<< exc.what() <<") in likelihood maximization. Skipping locus";
	PrintMessageDieOnError(msg.str(), M_PROGRESS);
      }
      locus->called[samp] = false;
      
    }
    locus->called[samp] = true;
  }
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
  //  double y = TestNLOPT();
}

Genotyper::~Genotyper() {
  delete read_extractor;
  for (std::map<std::string, LikelihoodMaximizer*>::iterator it = sample_likelihood_maximizers.begin();
       it != sample_likelihood_maximizers.end(); it++) {
    delete it->second;
  }
}
