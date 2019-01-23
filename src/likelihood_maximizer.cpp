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
#include <nlopt.hpp>
// #include <nlopt.h>

#include <gsl/gsl_multimin.h>
#include <gsl/gsl_siman.h>
#include "src/likelihood_maximizer.h"
#include "src/mathops.h"
#include "src/realignment.h" // for MARGIN
#include <iostream>
#include <algorithm>
using namespace std;


LikelihoodMaximizer::LikelihoodMaximizer(const Options& _options, const SampleProfile& sp,
					 const int32_t& read_len) {
  options = &_options; // TODO remove options

  enclosing_class_.SetGlobalParams(sp, options->flanklen, options->read_prob_mode, options->hist_mode);
  frr_class_.SetGlobalParams(sp, options->flanklen, options->read_prob_mode, options->hist_mode);
  spanning_class_.SetGlobalParams(sp, options->flanklen, options->read_prob_mode, options->hist_mode);
  flanking_class_.SetGlobalParams(sp, options->flanklen, options->read_prob_mode, options->hist_mode);
  resampled_enclosing_class_.SetGlobalParams(sp, options->flanklen, options->read_prob_mode, options->hist_mode);
  resampled_frr_class_.SetGlobalParams(sp, options->flanklen, options->read_prob_mode, options->hist_mode);
  resampled_spanning_class_.SetGlobalParams(sp, options->flanklen, options->read_prob_mode, options->hist_mode);
  resampled_flanking_class_.SetGlobalParams(sp, options->flanklen, options->read_prob_mode, options->hist_mode);
    
  obj_cov = sp.coverage;
  // Set up output file
  if (options->output_bootstrap) {
    bsfile_.open((options->outprefix + ".bootstrap.tab").c_str());
  }

  // Setup random number generator
  const gsl_rng_type * T;
  gsl_rng_env_setup();
  T = gsl_rng_default;
  r = gsl_rng_alloc (T);

  gsl_rng_set(r, options->seed);
  
  //offtarget_share = 0.0;

  // Set up default grid size
  lower_bound = 1;
  upper_bound = 600;
  grid_set = false;
  grid_buffer = 3;
  grid_opt_threshold = options->grid_threshold;

  // Locus params
  locus_params_set = false;
}

void LikelihoodMaximizer::Reset() {
  enclosing_class_.Reset();
  frr_class_.Reset();
  spanning_class_.Reset();
  flanking_class_.Reset();
  offtarget_class_.Reset();
  resampled_enclosing_class_.Reset();
  resampled_frr_class_.Reset();
  resampled_spanning_class_.Reset();
  resampled_flanking_class_.Reset();
  read_pool.clear();

  locus_params_set = false;
}

void LikelihoodMaximizer::AddEnclosingData(const int32_t& data) {
  enclosing_class_.AddData(data);
  ReadRecord rec;
  rec.read_type = RC_ENCL;
  rec.data = data;
  read_pool.push_back(rec);
}
void LikelihoodMaximizer::AddSpanningData(const int32_t& data) {
  spanning_class_.AddData(data);
  ReadRecord rec;
  rec.read_type = RC_SPAN;
  rec.data = data;
  read_pool.push_back(rec);
}
void LikelihoodMaximizer::AddFRRData(const int32_t& data) {
  frr_class_.AddData(data);
  ReadRecord rec;
  rec.read_type = RC_FRR;
  rec.data = data;
  read_pool.push_back(rec);
}
void LikelihoodMaximizer::AddFlankingData(const int32_t& data) {
  flanking_class_.AddData(data);
  ReadRecord rec;
  rec.read_type = RC_BOUND;
  rec.data = data;
  read_pool.push_back(rec);
}
void LikelihoodMaximizer::AddOffTargetData(const int32_t& data) {
  offtarget_class_.AddData(data);
  ReadRecord rec;
  rec.read_type = RC_OFFT;
  rec.data = data;
  read_pool.push_back(rec);
}

void LikelihoodMaximizer::SetLocusParams(const STRLocusInfo& sli, const double& cov,
					 const int32_t& _read_len, const int32_t _motif_len,
					 const int32_t& _ref_count) {
  enclosing_class_.SetLocusParams(sli);
  frr_class_.SetLocusParams(sli);
  spanning_class_.SetLocusParams(sli);
  flanking_class_.SetLocusParams(sli);
  resampled_enclosing_class_.SetLocusParams(sli);
  resampled_frr_class_.SetLocusParams(sli);
  resampled_spanning_class_.SetLocusParams(sli);
  resampled_flanking_class_.SetLocusParams(sli);

  enclosing_class_.SetCoverage(cov);
  frr_class_.SetCoverage(cov);
  spanning_class_.SetCoverage(cov);
  flanking_class_.SetCoverage(cov);
  resampled_enclosing_class_.SetCoverage(cov);
  resampled_frr_class_.SetCoverage(cov);
  resampled_spanning_class_.SetCoverage(cov);
  resampled_flanking_class_.SetCoverage(cov);

  obj_cov = cov;
  read_len = _read_len;
  motif_len = _motif_len;
  ref_count = _ref_count;
  locus_params_set = true;
}

void LikelihoodMaximizer::PrintReadPool(){
  if (resampled_pool.size() == read_pool.size()){
    vector<ReadRecord>::iterator resamp_rec = resampled_pool.begin();
    for (vector<ReadRecord>::iterator rec = read_pool.begin();
          rec != read_pool.end(); rec++){
      cerr<<rec->read_type<<"\t"<<rec->data<<"\t|\t"
          <<resamp_rec->read_type<<"\t"<<resamp_rec->data<<endl;
      resamp_rec++;
    }
  }
  else{
    for (vector<ReadRecord>::iterator rec = read_pool.begin();
          rec != read_pool.end(); rec++){
      cerr<<rec->read_type<<"\t"<<rec->data<<endl;
    } 
  }
}

void LikelihoodMaximizer::ResampleReadPool(){
  //gsl_rng_set(r, options->seed);   // Seed reset! ~~
  int32_t pool_size = read_pool.size();
  if (!resampled_pool.empty()){
    resampled_pool.clear();
  }
  resampled_pool.resize(pool_size);

  gsl_ran_sample(r, &resampled_pool[0], pool_size, &read_pool[0], pool_size, sizeof(ReadRecord));
  
  resampled_enclosing_class_.Reset();
  resampled_frr_class_.Reset();
  resampled_spanning_class_.Reset();
  resampled_flanking_class_.Reset();
  for (vector<ReadRecord>::iterator rec = resampled_pool.begin();
        rec != resampled_pool.end(); rec ++){
    if (rec->read_type == RC_ENCL){
      resampled_enclosing_class_.AddData(rec->data);
    }
    else if (rec->read_type == RC_FRR){
      resampled_frr_class_.AddData(rec->data);
    }
    else if (rec->read_type == RC_SPAN){
      resampled_spanning_class_.AddData(rec->data);
    }
    else if (rec->read_type == RC_BOUND){
      resampled_flanking_class_.AddData(rec->data);
    }
  }

  //PrintReadPool();
}

bool LikelihoodMaximizer::GetConfidenceInterval(const int32_t& all1,
						const int32_t& all2,
						const Locus& locus,
						double* lob1, double* hib1, double* lob2, double* hib2,
						double* a1_se, double* a2_se) {
  int32_t allele1, allele2;
  // TODO allow change of alpha
  double alpha = 0.05;   // Tail error on each end
  if (all1 > all2){
    allele2 = all1;
    allele1 = all2;
  }
  else{
    allele1 = all1;
    allele2 = all2;
  }
  int32_t num_boot_samp = options->num_boot_samp;
  int32_t boot_al1_1, boot_al1_2, boot_al2_1, boot_al2_2;
  int32_t boot_al1, boot_al2;
  double min_negLike;
  std::vector<int32_t> small_alleles, large_alleles;
  for (int i = 0; i < num_boot_samp + 1; i++){
    ResampleReadPool();
    if (options->ploidy == 2){
      OptimizeLikelihood(true, 1, allele1,
			 offtarget_share, 
			 &boot_al2_1, &boot_al2_2, &min_negLike);
      OptimizeLikelihood(true, 1, allele2,
			 offtarget_share, 
			 &boot_al1_1, &boot_al1_2, &min_negLike);
      if (boot_al1_1 == allele2)
	boot_al1 = boot_al1_2;
      else if (boot_al1_2 == allele2)
	boot_al1 = boot_al1_1;
      else
	cerr<< "Error running bootstrap\n";
      if (boot_al2_1 == allele1)
	boot_al2 = boot_al2_2;
      else if (boot_al2_2 == allele1)
	boot_al2 = boot_al2_1;
      else
	std::cerr<< "Error running likelihood optimization\n";
      double gt_ll1, gt_ll2;
    }
    else{ // haploid
      OptimizeLikelihood(true, 1, 0,
			 offtarget_share, 
			 &boot_al1, &boot_al2, &min_negLike);
    }

    small_alleles.push_back(boot_al1);
    large_alleles.push_back(boot_al2);
    if (options->output_bootstrap) {
      bsfile_ << locus.chrom << "\t" << locus.start << "\t" << locus.end << "\t"
	      << min(boot_al1, boot_al2) << "\t" << max(boot_al1, boot_al2) << endl;
    }
  }
  std::sort(small_alleles.begin(), small_alleles.end());
  std::sort(large_alleles.begin(), large_alleles.end());

  // Bootstrapping method from Davison and Hinkley 1997
  *lob1 = small_alleles.at(int((alpha / 2.0)  * (num_boot_samp + 1)));
  *hib1 = small_alleles.at(int((1.0 - alpha / 2.0) * (num_boot_samp + 1)));
  *lob2 = large_alleles.at(int((alpha / 2.0) * (num_boot_samp + 1)));
  *hib2 = large_alleles.at(int((1.0 - alpha / 2.0) * (num_boot_samp + 1)));
  
  
  double mean_sm_alleles = std::accumulate(small_alleles.begin(), 
					   small_alleles.end(), 0.0) / (num_boot_samp + 1);
  double mean_lg_alleles = std::accumulate(large_alleles.begin(), 
					   large_alleles.end(), 0.0) / (num_boot_samp + 1);
 
  double acum = 0;
  for (int i =0; i <= num_boot_samp; i++)
    acum += double(small_alleles[i] - mean_sm_alleles) * double(small_alleles[i] - mean_sm_alleles);
  *a1_se = std::sqrt(acum / double(num_boot_samp + 1));
  
  acum = 0;
  for (int i =0; i <= num_boot_samp; i++)
    acum += double(large_alleles[i] - mean_lg_alleles) * double(large_alleles[i] - mean_lg_alleles);
  *a2_se = std::sqrt(acum / double(num_boot_samp + 1));
  return true;  // TODO add return false cases
}

std::size_t LikelihoodMaximizer::GetEnclosingDataSize() {
  return enclosing_class_.GetDataSize();
}
std::size_t LikelihoodMaximizer::GetSpanningDataSize() {
  return spanning_class_.GetDataSize();
}
std::size_t LikelihoodMaximizer::GetFRRDataSize() {
  return frr_class_.GetDataSize();
}
std::size_t LikelihoodMaximizer::GetFlankingDataSize() {
  return flanking_class_.GetDataSize();
}
std::size_t LikelihoodMaximizer::GetOffTargetDataSize() {
  return offtarget_class_.GetDataSize();
}
std::size_t LikelihoodMaximizer::GetReadPoolSize() {
  return read_pool.size();
}

bool LikelihoodMaximizer::GetNegLikelihoodSurface(const int32_t& a_lo,
						  const int32_t& a_hi,
						  const int32_t& b_lo,
						  const int32_t& b_hi,
						  const bool& resampled,
						  double* surfaceLL) {
  double neg_inf = -100000000;

  // Checks on imputs
  if (a_hi<a_lo || b_hi<b_lo) {
    *surfaceLL = neg_inf; 
    return true;
  }

  double sum = neg_inf;
  double min = 10000.0;
  int m_a, m_b;
  double ret_val = 0.0;
  sum = neg_inf;
  for (int i = a_lo; i <= a_hi; i++){
    for (int j = b_lo; j <= b_hi; j++){
      if (!this->GetGenotypeNegLogLikelihood(i, j, resampled,
					     &ret_val)) {
	*surfaceLL = neg_inf;
	return false;
      }
      //cerr << sum << " " << ret_val << endl;
      sum = fast_log_sum_exp(sum, (-1.0) * ret_val);
      if (ret_val < min){
	min = ret_val;
	m_a = i;
	m_b = j;
      }
    }
  }
  *surfaceLL = sum;
  //  cerr << "Min: "<< min << " at: " << m_a << ", " << m_b << endl;

  return true;
}

bool LikelihoodMaximizer::GetGenotypeNegLogLikelihood(const int32_t& allele1,
						      const int32_t& allele2,
						      const bool& resampled,
						      double* gt_ll) {
  double frr_count_ll = 0.0, frr_ll, span_ll, encl_ll, flank_ll = 0.0;
  double count_weight = obj_cov;
  bool use_cov = options -> use_cov;
  int frr_count, offtarget_count = offtarget_class_.GetDataSize();

  int read_count;
  if (allele1 < 0 || allele2 < 0){
    *gt_ll = frr_class_.NEG_INF;
    return true;
  }
  
  if (!resampled){
    frr_count = frr_class_.GetDataSize();
    read_count = frr_count + enclosing_class_.GetDataSize() +
      spanning_class_.GetDataSize() + flanking_class_.GetDataSize() + 2 *  offtarget_count;

    if (read_count == 0){
      *gt_ll = frr_class_.NEG_INF;
      return true;
    }
    frr_class_.GetClassLogLikelihood(allele1, allele2, 
    				     read_len, motif_len, ref_count, 
    				     options->ploidy, &frr_ll);
    spanning_class_.GetClassLogLikelihood(allele1, allele2, 
    					  read_len, motif_len, ref_count, 
    					  options->ploidy, &span_ll);
    enclosing_class_.GetClassLogLikelihood(allele1, allele2, 
					     read_len, motif_len, ref_count, 
					   options->ploidy, &encl_ll);

    // flanking class overloads GetClassLogLikelihood function
    flanking_class_.FlankingClass::GetClassLogLikelihood(allele1, allele2, 
    							 read_len, motif_len, ref_count, 
    							 options->ploidy, &flank_ll);
    // TODO Substituting these lines changes optimization result. Find out why?!
    //if ((options->coverage > 0) && (frr_class_.GetDataSize() > 0)){
    
    if (use_cov && obj_cov > 0 && frr_count > 0){
      frr_class_.GetCountLogLikelihood(allele1, 
				       allele2,
				       read_len, 
				       motif_len, 
				       obj_cov,
				       options->ploidy, 
				       2 * offtarget_count * offtarget_share, 
				       &frr_count_ll);
      //cerr << allele1 << ","<< allele2 << "\t"<< frr_count << " " << frr_count_ll << endl; 
    }
    //    cerr << allele1 << ","<< allele2 << "\t" << flank_ll << endl; 
  }
  else {
    frr_count = resampled_frr_class_.GetDataSize();
    read_count = frr_count + resampled_enclosing_class_.GetDataSize() +
      resampled_spanning_class_.GetDataSize() + resampled_flanking_class_.GetDataSize() + 2 * offtarget_count;
    if (read_count == 0){
      *gt_ll = frr_class_.NEG_INF;
      return true;
    }
    resampled_frr_class_.GetClassLogLikelihood(allele1, allele2, 
					       read_len, motif_len, ref_count, 
					       options->ploidy, &frr_ll);
    resampled_spanning_class_.GetClassLogLikelihood(allele1, allele2, 
						    read_len, motif_len, ref_count, 
						    options->ploidy, &span_ll);
    resampled_enclosing_class_.GetClassLogLikelihood(allele1, allele2, 
						     read_len, motif_len, ref_count, 
						     options->ploidy, &encl_ll);
    // flanking class overloads GetClassLogLikelihood function
    resampled_flanking_class_.FlankingClass::GetClassLogLikelihood(allele1, allele2, 
								   read_len, motif_len, ref_count, 
								   options->ploidy, &flank_ll); 
    
    if (use_cov && obj_cov > 0 && frr_count > 0){
      resampled_frr_class_.GetCountLogLikelihood(allele1, 
				       allele2,
				       read_len, 
				       motif_len, 
				       obj_cov,
				       options->ploidy, 
				       2 * offtarget_count * offtarget_share, 
				       &frr_count_ll);
    
    }
  }
  double factor = read_count;
  factor = 1.0; // Factor is NOT applied to the count term 
  *gt_ll = -1*(1.0 / factor * (options->frr_weight * frr_ll +
			       options->spanning_weight * span_ll +
			       options->enclosing_weight * encl_ll + 
			       options->flanking_weight * flank_ll) + 
  	       count_weight * frr_count_ll);
  return true;
}

const int32_t STARTMIN = 1000000;
const int32_t STARTMAX = 0;
const int32_t DEFAULTMAX = 100;
bool LikelihoodMaximizer::InferGridSize() {
  grid_set = false;
  int32_t min_allele = STARTMIN;
  int32_t max_allele = STARTMAX;
  if (enclosing_class_.GetGridBoundaries(&min_allele, &max_allele)) {
    lower_bound = min_allele;
    upper_bound = max_allele;
    grid_set = true;
  }
  if (spanning_class_.GetGridBoundaries(&min_allele, &max_allele)) {
    lower_bound = min_allele;
    upper_bound = max_allele;
    grid_set = true;
  }
  int32_t offtarget_count = offtarget_class_.GetDataSize();
  if (frr_class_.GetGridBoundaries(&min_allele, &max_allele,
				   read_len, motif_len,
				   obj_cov, offtarget_count)) {
    lower_bound = min_allele;
    upper_bound = max_allele;
    grid_set = true;
  }
  if (flanking_class_.GetGridBoundaries(&min_allele, &max_allele)) {
    lower_bound = min_allele;
    upper_bound = max_allele;
    grid_set = true;
  }
  // Check boundaries make sense
  if (lower_bound == STARTMIN) {
    // min wasn't set. happens if no enclosing found
    lower_bound = 1;
  }
  if (upper_bound == STARTMAX) {
    // max wasn't set. happens if only spanning found
    upper_bound = DEFAULTMAX;
  }
  // Add buffer to the grid
  lower_bound = max(1, lower_bound - grid_buffer);
  upper_bound += grid_buffer;
  return grid_set;
}

void LikelihoodMaximizer::SetGridSize(const int32_t& min_allele, const int32_t max_allele) {
  lower_bound = min_allele;
  upper_bound = max_allele;
}

void LikelihoodMaximizer::GetGridSize(int32_t* min_allele, int32_t* max_allele) {
  *min_allele = lower_bound;
  *max_allele = upper_bound;
}

void LikelihoodMaximizer::InferAlleleList(std::vector<int32_t>* allele_list,
					  const int32_t& ploidy,
					  const bool& resampled, const int32_t& fix_allele) {
  allele_list->clear();
  enclosing_class_.ExtractEnclosingAlleles(allele_list);
  if (upper_bound-lower_bound <= grid_opt_threshold) {
    for (int32_t i=lower_bound; i<=upper_bound; i++) {
      allele_list->push_back(i);
    }
  } else {
    std::vector<int32_t> sublist;
    int32_t a1, a2, result;
    double minf;
    if (ploidy == 2) {
      for (std::vector<int32_t>::iterator allele_it = allele_list->begin();
	   allele_it != allele_list->end();
	   allele_it++) {	
	// 1-D optimization fixing each enclosing allele
	nlopt_1D_optimize(read_len, motif_len, ref_count, 
			  lower_bound, upper_bound, resampled, 
			  options->seed, this, *allele_it, &a1, &result, &minf);
	sublist.push_back(a1);
      }
      // 2D opt
      nlopt_2D_optimize(read_len, motif_len, ref_count, 
			lower_bound, upper_bound, resampled, 
			options->seed, this, &a1, &a2, &result, &minf);
      sublist.push_back(a1);
      sublist.push_back(a2);
      for (std::vector<int32_t>::iterator subl_it = sublist.begin();
	   subl_it != sublist.end();
	   subl_it++) {
	if(std::find(allele_list->begin(), allele_list->end(), *subl_it) == allele_list->end()) {
          allele_list->push_back(*subl_it);
	}
      }
    } else if (ploidy == 1) {
      nlopt_1D_optimize(read_len, motif_len, ref_count, 
			lower_bound, upper_bound, resampled, 
			options->seed, this, fix_allele, &a1, &result, &minf);
      allele_list->push_back(a1);
    }
  }
}
bool LikelihoodMaximizer::GetQScore(int32_t a1, int32_t a2, double* Q) {
  double surface, gt_ll;
  *Q = -1;
  if (!GetNegLikelihoodSurface(lower_bound, upper_bound,
			       lower_bound, upper_bound,
			       false,
			       &surface)) {
    return false;
  }
  if (!GetGenotypeNegLogLikelihood(a1, a2, false, &gt_ll)) {
    return false;
  }
  // adding log(2) because of uninformative prior: P(a1,a2) + P(a2,a1), unless a1 == a2
  // Posterior: Q = gt_ll * prior_a1_a2 / (surface * prior_ai_aj) 
  // if a1 == a2: prior_a1_a2 = prior_ai_aj
  // if a1 != aj: prior_a1_a2 = 2 * prior_ai_aj (to consider for both cases as described above ^^)
  if (a1 == a2){
    *Q = exp(gt_ll * -1 - surface);
  }
  else {
    *Q = exp(gt_ll * -1 - surface + log(2));
  }
  if (*Q > 1.01 or *Q < 0){
    *Q = -1;
    return false;
  }
  return true;
}
bool LikelihoodMaximizer::GetExpansionProb(std::vector<double>* prob_vec, const int32_t& exp_threshold) {
  if (exp_threshold == -1) {
    prob_vec->clear();
    prob_vec->push_back(-1);
    prob_vec->push_back(-1);
    prob_vec->push_back(-1);
    return true;
  }
  double shortshort, shortlong, longlong;

  if (!GetNegLikelihoodSurface(lower_bound, exp_threshold-1,
			       lower_bound, exp_threshold-1,
			       false,
			       &shortshort)) {
    return false;
  }
  if (!GetNegLikelihoodSurface(exp_threshold, upper_bound,
			       exp_threshold, upper_bound,
			       false,
			       &longlong)) {
    return false;
  }
  if (!GetNegLikelihoodSurface(lower_bound, exp_threshold-1,
			       exp_threshold, upper_bound,
			       false,
			       &shortlong)) {
    return false;
  }
  // Divide these by two since double counting (a,b) and (b,a)
  shortshort -= log(2); 
  longlong -= log(2); 
  //double total = log(exp(shortshort)+exp(longlong)+exp(shortlong));
  double total = fast_log_sum_exp(shortshort, fast_log_sum_exp(longlong,shortlong));
  prob_vec->clear();
  prob_vec->push_back(exp(shortshort-total));
  prob_vec->push_back(exp(shortlong-total));
  prob_vec->push_back(exp(longlong-total));
  return true;
}

bool LikelihoodMaximizer::OptimizeLikelihood(const bool& resampled, const int32_t& use_ploidy,
					     const int32_t& fix_allele,
					     const double& off_share,
					     int32_t* allele1, int32_t* allele2, double* min_negLike) {
  if (!locus_params_set) {
    PrintMessageDieOnError("Skipping locus with no params set", M_WARNING);
    return false;
  }
  if (obj_cov == -1) {
    PrintMessageDieOnError("Skipping locus with likely extreme GC content", M_WARNING);
    return false;
  }

  offtarget_share = off_share; // perc. of offtarget reads.

  // Get list of potential alleles to try
  std::vector<int32_t> allele_list;
  InferAlleleList(&allele_list, use_ploidy, resampled, fix_allele);
  /*
  if (!resampled){
    double gt_ll1;
    for (int i1 = 1; i1 < 47; i1++){
      for (int i2 = 1; i2 < 47; i2++){
	GetGenotypeNegLogLikelihood(i1, i2, resampled, &gt_ll1);
	cerr << ">> " << i1 << " " << i2 << "\t" << gt_ll1 << endl;
      }
    }
  }
  */
  if (use_ploidy == 2) {
    findBestAlleleListTuple(allele_list, use_ploidy, resampled, 0, 
			    allele1, allele2, min_negLike);
  } else if (use_ploidy == 1) {
    findBestAlleleListTuple(allele_list, use_ploidy, resampled, fix_allele, 
			    allele1, allele2, min_negLike);    
  }
  return true;
}


bool LikelihoodMaximizer::findBestAlleleListTuple(std::vector<int32_t> allele_list,
						  int32_t use_ploidy,
						  bool resampled, int32_t fix_allele,
						  int32_t* allele1, int32_t* allele2, double* min_negLike) {
  double gt_ll;
  *min_negLike = 1000000;
  int32_t best_a1 = 0, best_a2 = 0;
  if (use_ploidy == 2) {
    for (std::vector<int32_t>::iterator a1_it = allele_list.begin();
            a1_it != allele_list.end();
            a1_it++){
      // if (!resampled)
      //   cerr<<*a1_it<<endl;
      for (std::vector<int32_t>::iterator a2_it = allele_list.begin();
            a2_it != allele_list.end();
            a2_it++){
	if (*a2_it < *a1_it) continue; // want a1 the smaller allele
        GetGenotypeNegLogLikelihood(*a1_it, *a2_it, resampled, &gt_ll);
        //if (!resampled)
	//  cerr<<endl<<*a1_it<<"\t"<<*a2_it<<"\t"<<gt_ll<<endl;
	if (gt_ll < *min_negLike){
	  *min_negLike = gt_ll;
	  best_a1 = *a1_it;
	  best_a2 = *a2_it;
	}
      }
    }
  }
  else if (use_ploidy == 1) {
    best_a2 = fix_allele;
    for (std::vector<int32_t>::iterator a1_it = allele_list.begin();
            a1_it != allele_list.end();
            a1_it++){
      GetGenotypeNegLogLikelihood(*a1_it, fix_allele, resampled, &gt_ll);
      // cerr<<">> "<<fix_allele<<"\t"<<*a1_it<<"\t"<<gt_ll<<endl;
      if (gt_ll < *min_negLike){
        *min_negLike = gt_ll;
        best_a1 = *a1_it;
      }
    }
  }
  
  *allele1 = best_a1;
  *allele2 = best_a2;
  return true;    // TODO add false
}

LikelihoodMaximizer::~LikelihoodMaximizer() {
  if (options->output_bootstrap) {
    bsfile_.close();
  }
  gsl_rng_free(r);
}

double nloptNegLikelihood(unsigned n, const double *x, double *grad, void *data)
{
  if (grad) {
        cerr<< "No grad!"<<endl;
        return 0.0;
  }
  nlopt_data *d = (nlopt_data *) data;
  int read_len  = d -> read_len;
  int motif_len = d -> motif_len;
  int ref_count = d -> ref_count;
  int fix_allele = d -> fix_allele; 
  bool resampled = d -> resampled;
  LikelihoodMaximizer* lm_ptr = d -> lm_ptr;
  double gt_ll;
  if (n == 2){
    double A = x[0], B = x[1];
    if(!lm_ptr->GetGenotypeNegLogLikelihood(A, B, resampled, &gt_ll))
      return -100.0;
    else{
      return gt_ll;
    }
  }
  else{
    double A = x[0], B = fix_allele;
    if(!lm_ptr->GetGenotypeNegLogLikelihood(A, B, resampled, &gt_ll))
      return -100.0;
    else{
      return gt_ll;
    }
  
  }
}

bool nlopt_2D_optimize(const int32_t& read_len, const int32_t& motif_len,
		       const int32_t& ref_count, const int32_t& lower_bound,
		       const int32_t& upper_bound, const bool& resampled, 
		       const int& seed, LikelihoodMaximizer* lm_ptr,
		       int32_t* allele1, int32_t* allele2, int32_t* ret_result, double* minf_ret) {
  // Seed reset! ~~
  nlopt::srand(seed);
  nlopt::opt opt(nlopt::LN_COBYLA, 2);
  // opt.set_local_optimizer(nlopt::LN_COBYLA)   // TODO check nlopt::G_MLSL_LDS->multiple local
  std::vector<double> lb(2);
  lb[0] = lower_bound;
  lb[1] = lower_bound;
  opt.set_lower_bounds(lb);

  std::vector<double> ub(2);
  ub[0] = upper_bound;
  ub[1] = upper_bound;
  opt.set_upper_bounds(ub);

  nlopt_data data = nlopt_data(read_len, motif_len, ref_count, lm_ptr, 0, resampled);
  opt.set_min_objective(nloptNegLikelihood, &data);    // Change to max for maximization

  opt.set_xtol_rel(.00005);   // TODO set something appropriate
  std::vector<double> xx(2);
  double minf=100000.0, f;
  nlopt::result result;
  for (double j = 0.1; j <= 0.3; j+=0.1) {
    for (double k = 0.25; k <= 0.75; k+=0.25){
      nlopt::srand(seed);
      xx[0] = int32_t(lower_bound + j * float(upper_bound - lower_bound));
      xx[1] = int32_t(lower_bound + k * float(upper_bound - lower_bound));
      if (xx[0] > upper_bound) { xx[0] = upper_bound;}
      if (xx[1] > upper_bound) { xx[1] = upper_bound;}
      if (xx[0] < lower_bound) { xx[0] = lower_bound;}
      if (xx[1] < lower_bound) { xx[1] = lower_bound;}

      result = opt.optimize(xx, f);

      if (f < minf){
	*allele1 = int32_t(round(xx[0]));
	*allele2 = int32_t(round(xx[1]));
	*ret_result = result;
	
	minf = f;
      }
    }
  }
  
  
  *minf_ret = minf;
  return true;  // TODO add false
}

bool nlopt_1D_optimize(const int32_t& read_len, const int32_t& motif_len,
		       const int32_t& ref_count, const int32_t& lower_bound,
		       const int32_t& upper_bound, const bool& resampled, 
		       const int& seed, LikelihoodMaximizer* lm_ptr,
		       const int32_t& fix_allele, int32_t* allele1,
		       int32_t* ret_result, double* minf_ret) {
  // Seed reset! ~~
  nlopt::srand(seed);
  nlopt::opt opt(nlopt::LN_COBYLA, 1);

  std::vector<double> lb(1);
  lb[0] = lower_bound;
  opt.set_lower_bounds(lb);

  std::vector<double> ub(1);
  ub[0] = upper_bound;
  opt.set_upper_bounds(ub);

  nlopt_data data = nlopt_data(read_len, motif_len, ref_count, lm_ptr, fix_allele, resampled);
  opt.set_min_objective(nloptNegLikelihood, &data);    // Change to max for maximization

  opt.set_xtol_rel(.0005);   // TODO set something appropriate

  std::vector<double> xx(1);
  xx[0] = int32_t(lower_bound + 0.5 * float(upper_bound - lower_bound));
  if (xx[0] > upper_bound) { xx[0] = upper_bound;}
  if (xx[0] < lower_bound) { xx[0] = lower_bound;}
  double minf;
  nlopt::result result = opt.optimize(xx, minf);
  *allele1 = int32_t(round(xx[0]));
  *ret_result = result;
  *minf_ret = minf;
  return true;  // TODO add false
}
