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
		     SampleInfo& _sample_info,
		     STRInfo& _str_info) {
  refgenome = &_refgenome;
  options = &_options;
  sample_info = &_sample_info;
  str_info = &_str_info;
  read_extractor = new ReadExtractor(_options, *sample_info);
  std::set<std::string> rg_samples = sample_info->GetSamples();
  for (std::set<std::string>::iterator it=rg_samples.begin();
       it != rg_samples.end(); it++) {
    SampleProfile sp;
    if (sample_info->GetSampleProfile(*it, &sp)) {
      sample_likelihood_maximizers[*it] = new LikelihoodMaximizer(_options, sp, sample_info->GetReadLength());
    } else {
      PrintMessageDieOnError("Could not find sample profile for " + *it, M_ERROR, false);
    }
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

bool Genotyper::SetGGL(Locus& locus, const std::string& samp) {
  std::map<std::pair<int32_t, int32_t>, double> gridlik;
  for (int32_t a1=locus.grid_min_allele; a1<=locus.grid_max_allele; a1++) {
    for (int32_t a2=a1; a2<=locus.grid_max_allele; a2++) {
      double gt_ll;
      if (!sample_likelihood_maximizers[samp]->GetGenotypeNegLogLikelihood(a1, a2, false, &gt_ll)) {
	return false;
      }
      std::pair<int32_t, int32_t> gt(a1, a2);
      gridlik[gt] = gt_ll/log(10)*-1;
    }
  }
  locus.grid_likelihoods[samp] = gridlik;
  return true;
}

bool Genotyper::ProcessLocus(BamCramMultiReader* bamreader, Locus* locus) {
  int32_t read_len = sample_info->GetReadLength();

  // Load preflank and postflank to locus
  if (options->verbose) {
    PrintMessageDieOnError("\tSetting flanking regions", M_PROGRESS, options->quiet);
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
    PrintMessageDieOnError("\tLoading read data", M_PROGRESS, options->quiet);
  }
  if (!read_extractor->ExtractReads(bamreader, *locus, options->regionsize,
				    options->min_match, sample_likelihood_maximizers)) {
    return false;
  }


  std::set<std::string> rg_samples = sample_info->GetSamples();
  // First set grid size
  int32_t sample_min_allele, sample_max_allele;
  int32_t min_allele = 100000;
  int32_t max_allele = 0;
  int32_t ref_count = (int32_t)((locus->end-locus->start+1)/locus->motif.size());
  for (std::set<std::string>::iterator it = rg_samples.begin();
       it != rg_samples.end(); it++) {
    const std::string samp = *it;
    if (gcbin != -1) {
      sample_likelihood_maximizers[samp]->SetLocusParams(str_info->GetSTRInfo(locus->chrom, locus->start),
							 sample_info->GetGCCoverage(samp, gcbin),
							 sample_info->GetReadLength(), (int32_t)(locus->motif.size()),
							 ref_count);
    } else {
      sample_likelihood_maximizers[samp]->SetLocusParams(str_info->GetSTRInfo(locus->chrom, locus->start),
							 sample_info->GetCoverage(samp),
							 sample_info->GetReadLength(), (int32_t)(locus->motif.size()),
							 ref_count);
    }
    if (!sample_likelihood_maximizers[samp]->InferGridSize() ) {
      PrintMessageDieOnError("Error inferring grid size", M_PROGRESS, options->quiet);
      return false;
    }
    sample_likelihood_maximizers[samp]->GetGridSize(&sample_min_allele, &sample_max_allele);
    if (sample_max_allele > max_allele) max_allele = sample_max_allele;
    if (sample_min_allele < min_allele) min_allele = sample_min_allele;
  }
  // Set locus info to output to VCF
  locus->grid_min_allele = min_allele;
  locus->grid_max_allele = max_allele;
  locus->expansion_threshold = str_info->GetExpansionThreshold(locus->chrom, locus->start);
  locus->stutter_up = str_info->GetStutterUp(locus->chrom, locus->start);
  locus->stutter_down = str_info->GetStutterDown(locus->chrom, locus->start);
  locus->stutter_p = str_info->GetStutterP(locus->chrom, locus->start);
  // Maximize the likelihood
  if (options->verbose) {
    PrintMessageDieOnError("\tMaximizing likelihood", M_PROGRESS, options->quiet);
  }
  int32_t allele1, allele2;
  double min_negLike, lob1, lob2, hib1, hib2, q_score;
  double a1_se, a2_se;
  bool resampled = false;
  for (std::set<std::string>::iterator it = rg_samples.begin();
       it != rg_samples.end(); it++) {
    const std::string samp = *it;
    locus->called[samp] = false;
    sample_likelihood_maximizers[samp]->SetGridSize(min_allele, max_allele);
    try {
      if (!sample_likelihood_maximizers[samp]->OptimizeLikelihood(resampled, options->ploidy,
								  0,
								  locus->offtarget_share,
								  &allele1, 
								  &allele2, 
								  &min_negLike)) {
	continue;
      }
      std::vector<double> sample_prob_vec;
      if (!sample_likelihood_maximizers[samp]->GetExpansionProb(&sample_prob_vec, locus->expansion_threshold)) {
	sample_prob_vec.clear();
	sample_prob_vec.push_back(-1.0);
	sample_prob_vec.push_back(-1.0);
	sample_prob_vec.push_back(-1.0);
      }
      if (!options->skip_qscore){
	if (!sample_likelihood_maximizers[samp]->GetQScore(allele1, allele2, &q_score)){ 
	  q_score = -1;
	  PrintMessageDieOnError("\tProblem setting quality scores", M_WARNING, options->quiet);
	}
      }
      else{
	q_score = -1;
      }
      locus->q_scores[samp] = q_score;
      locus->expansion_probs[samp] = sample_prob_vec;
      locus->allele1[samp] = allele1;
      locus->allele2[samp] = allele2;
      locus->min_neg_lik[samp] = min_negLike;
      locus->enclosing_reads[samp] = sample_likelihood_maximizers[samp]->GetEnclosingDataSize();
      locus->spanning_reads[samp] = sample_likelihood_maximizers[samp]->GetSpanningDataSize();
      locus->frr_reads[samp] = sample_likelihood_maximizers[samp]->GetFRRDataSize();
      locus->flanking_reads[samp] = sample_likelihood_maximizers[samp]->GetFlankingDataSize();
      locus->enclosing_reads_dict[samp] = sample_likelihood_maximizers[samp]->GetEnclosingReadDictStr();
      locus->flanking_reads_dict[samp] = sample_likelihood_maximizers[samp]->GetFlankingReadDictStr();


      locus->depth[samp] = sample_likelihood_maximizers[samp]->GetReadPoolSize();
      locus->called[samp] = true;
      if (allele1 <= 0 and allele2 <= 0){
	PrintMessageDieOnError("\tProblem maximizing likelihood. Skipping locus", M_WARNING, options->quiet);
	locus->called[samp] = false;
	continue;
      }
      if (options->include_ggl && !SetGGL(*locus, samp)) {
	PrintMessageDieOnError("\tProblem setting genotype likelihoods", M_WARNING, options->quiet);
      }
      if (options->num_boot_samp > 0){
	if (options->verbose) {
	  PrintMessageDieOnError("\tGetting confidence intervals", M_PROGRESS, options->quiet);
	}
	try{
	  if (!sample_likelihood_maximizers[samp]->GetConfidenceInterval(allele1, allele2, *locus,
									 &lob1, &hib1, &lob2, &hib2, &a1_se, &a2_se)) {
	    locus->called[samp] = false;
	    continue;
	  }
	  locus->lob1[samp] = lob1;
	  locus->lob2[samp] = lob2;
	  locus->hib1[samp] = hib1;
	  locus->hib2[samp] = hib2;
	  locus->a1_se[samp] = a1_se;
	  locus->a2_se[samp] = a2_se;
	  
	  stringstream msg;
	  msg<<"\tGenotyper Results:  "<<allele1<<", "<<allele2<<"\tlikelihood = "<<min_negLike;
	  PrintMessageDieOnError(msg.str(), M_PROGRESS, options->quiet);
	  if (options->verbose) {
	    msg.clear();
	    msg.str(std::string());
	    msg<<"\tSmall Allele Bound: ["<<lob1<<", "<<hib1<<"]";
	    PrintMessageDieOnError(msg.str(), M_PROGRESS, options->quiet);
	    msg.clear();
	    msg.str(std::string());
	    msg<<"\tLarge Allele Bound: ["<<lob2<<", "<<hib2<<"]";
	    PrintMessageDieOnError(msg.str(), M_PROGRESS, options->quiet);
	  }
	}
	catch (std::exception &exc){
	  if (options->verbose) {
	    stringstream msg;
	    msg<<"\tEncountered error("<< exc.what() <<") in likelihood maximization for confidence interval. Skipping locus";
	    PrintMessageDieOnError(msg.str(), M_PROGRESS, options->quiet);
	  }
	  locus->called[samp] = false;
	}
      }
    }
    catch (std::exception &exc){
      if (options->verbose) {
	stringstream msg;
	msg<<"\tEncountered error("<< exc.what() <<") in likelihood maximization. Skipping locus";
	PrintMessageDieOnError(msg.str(), M_PROGRESS, options->quiet);
      }
      locus->called[samp] = false;
    }
  }
  return true;
}

void Genotyper::Debug(BamCramMultiReader* bamreader) {
  cerr << "testing refgenome" << endl;
  std::string seq;
  refgenome->GetSequence("3", 63898261, 63898360, &seq);
  cerr << seq << endl;
  cerr << "testing bam" << endl;
  bamreader->SetRegion("1", 0, 10000, options->trim_to_readlen);
  BamAlignment aln;
  if (bamreader->GetNextAlignment(aln, options->trim_to_readlen)) { // Requires SetRegion was called
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
