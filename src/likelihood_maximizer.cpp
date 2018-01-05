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


LikelihoodMaximizer::LikelihoodMaximizer(Options& _options) {
  options = &_options;

  enclosing_class_.SetOptions(*options);
  frr_class_.SetOptions(*options);
  spanning_class_.SetOptions(*options);
  flanking_class_.SetOptions(*options);
  resampled_enclosing_class_.SetOptions(*options);
  resampled_frr_class_.SetOptions(*options);
  resampled_spanning_class_.SetOptions(*options);
  resampled_flanking_class_.SetOptions(*options);

  // Set up output file
  if (options->output_bootstrap) {
    bsfile_.open((options->outprefix + ".bootstrap.tab").c_str());
  }
  plotfile_.open((options->outprefix + ".plot.tab").c_str());

  // Setup random number generator
  const gsl_rng_type * T;
  gsl_rng_env_setup();
  T = gsl_rng_default;
  r = gsl_rng_alloc (T);

  gsl_rng_set(r, options->seed);
}

// // Not needed. since options are updated before creating likelihood maximizer object
// TODO delete
// void LikelihoodMaximizer::UpdateOptions(){
  
//   // std::cerr << options->dist_mean << ", " << options->dist_sdev << endl;
//   enclosing_class_.SetOptions(*options);
//   frr_class_.SetOptions(*options);
//   spanning_class_.SetOptions(*options);
//   flanking_class_.SetOptions(*options);
//   resampled_enclosing_class_.SetOptions(*options);
//   resampled_frr_class_.SetOptions(*options);
//   resampled_spanning_class_.SetOptions(*options);
//   resampled_flanking_class_.SetOptions(*options);
// }

void LikelihoodMaximizer::Reset() {
  enclosing_class_.Reset();
  frr_class_.Reset();
  spanning_class_.Reset();
  flanking_class_.Reset();
  read_pool.clear();
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

void LikelihoodMaximizer::PlotLikelihood(int32_t fix_allele,
                                          int32_t start,
                                          int32_t end,
                                          int32_t step,
                                          int32_t read_len,
                                          int32_t motif_len,
                                          int32_t ref_count){
  double gt_ll;
  for (int32_t var_allele = start; var_allele <= end; var_allele+=step){
    GetGenotypeNegLogLikelihood(fix_allele, var_allele, read_len, motif_len, ref_count, false, &gt_ll);
    plotfile_ << fix_allele << "\t" << var_allele << "\t"
        << gt_ll << endl;
  }
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

  // PrintReadPool();
}

bool LikelihoodMaximizer::GetConfidenceInterval(const int32_t& read_len, 
						const int32_t& motif_len,
						const int32_t& ref_count,
						const int32_t& all1,
						const int32_t& all2,
						const Locus& locus,
						double* lob1, double* hib1, double* lob2, double* hib2){
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
  for (int i = 0; i < num_boot_samp; i++){
    ResampleReadPool();
    if (options->ploidy == 2){
      OptimizeLikelihood(read_len, motif_len, ref_count, true, 1, allele1, &boot_al2_1, &boot_al2_2, &min_negLike);
      OptimizeLikelihood(read_len, motif_len, ref_count, true, 1, allele2, &boot_al1_1, &boot_al1_2, &min_negLike);
      if (boot_al1_1 == allele2)
	boot_al1 = boot_al1_2;
      else if (boot_al1_2 == allele2)
	boot_al1 = boot_al1_1;
      else
	cerr<< "Hell Na\n";

      if (boot_al2_1 == allele1)
	boot_al2 = boot_al2_2;
      else if (boot_al2_2 == allele1)
	boot_al2 = boot_al2_1;
      else
	cerr<< "Hell Na\n";
      double gt_ll1, gt_ll2;
      //PrintReadPool();
      GetGenotypeNegLogLikelihood(allele1, boot_al1, read_len, motif_len, ref_count, true, &gt_ll1);
      GetGenotypeNegLogLikelihood(allele2, 6, read_len, motif_len, ref_count, true, &gt_ll2);
      //cerr << allele1 << ", "<< boot_al1 << "\t" << gt_ll1  << "\n";
      //cerr << allele2 << ", "<< 6 << "\t" << gt_ll2  << "\n";
      // cerr << allele1 << "\t" << boot_al1 << endl;
      // cerr << allele2 << "\t" << boot_al2 << endl;

    }    
    // cerr<<min(boot_al1, boot_al2)<<"\t"<<max(boot_al1, boot_al2)<<endl;
    // small_alleles.push_back(min(boot_al1, boot_al2) - allele1);
    // large_alleles.push_back(max(boot_al1, boot_al2) - allele2);
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

  // TODO 0.9 or 0.1? allele1 -/+ lob1?
  // TODO allow change of 0.9 and 0.1

  // *lob1 = small_alleles.at(int(0.9 * (num_boot_samp + 1)));
  // *hib1 = small_alleles.at(int(0.1 * (num_boot_samp + 1)));
  // *lob2 = large_alleles.at(int(0.9 * (num_boot_samp + 1)));
  // *hib2 = large_alleles.at(int(0.1 * (num_boot_samp + 1)));

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
std::size_t LikelihoodMaximizer::GetReadPoolSize() {
  return read_pool.size();
}

bool LikelihoodMaximizer::GetGenotypeNegLogLikelihood(const int32_t& allele1,
						      const int32_t& allele2,
						      const int32_t& read_len,
						      const int32_t& motif_len,
						      const int32_t& ref_count,
                  const bool& resampled,
						      double* gt_ll) {
  double frr_ll, span_ll, encl_ll, flank_ll = 0.0;
  if (!resampled){
    frr_class_.GetClassLogLikelihood(allele1, allele2, read_len, motif_len, ref_count, options->ploidy, &frr_ll);
    spanning_class_.GetClassLogLikelihood(allele1, allele2, read_len, motif_len, ref_count, options->ploidy, &span_ll);
    enclosing_class_.GetClassLogLikelihood(allele1, allele2, read_len, motif_len, ref_count, options->ploidy, &encl_ll);
    // flanking class overloads GetClassLogLikelihood function
    flanking_class_.FlankingClass::GetClassLogLikelihood(allele1, allele2, read_len, motif_len, ref_count, options->ploidy, &flank_ll);
  }
  else {
    resampled_frr_class_.GetClassLogLikelihood(allele1, allele2, read_len, motif_len, ref_count, options->ploidy, &frr_ll);
    resampled_spanning_class_.GetClassLogLikelihood(allele1, allele2, read_len, motif_len, ref_count, options->ploidy, &span_ll);
    resampled_enclosing_class_.GetClassLogLikelihood(allele1, allele2, read_len, motif_len, ref_count, options->ploidy, &encl_ll);
    // flanking class overloads GetClassLogLikelihood function
    resampled_flanking_class_.FlankingClass::GetClassLogLikelihood(allele1, allele2, read_len, motif_len, ref_count, options->ploidy, &flank_ll); 
  }
  *gt_ll = -1*(options->frr_weight*frr_ll +
	       options->spanning_weight*span_ll +
	       options->enclosing_weight*encl_ll + 
         options->flanking_weight*flank_ll);
}

bool LikelihoodMaximizer::OptimizeLikelihood(const int32_t& read_len, const int32_t& motif_len,
					     const int32_t& ref_count, const bool& resampled, 
					     const int32_t& ploidy, const int32_t& fix_allele,
					     int32_t* allele1, int32_t* allele2, double* min_negLike) {

  int32_t a1, a2, result, temp;
  double minf;
  std::vector<int32_t> allele_list;
  std::vector<int32_t> sublist;
  this->enclosing_class_.ExtractEnclosingAlleles(&allele_list);
  ResampleReadPool();
  int32_t upper_bound = 500; // TODO Change 200 for number depending the parameters
  int32_t lower_bound_1d, lower_bound_2d;
  // if (allele_list.size() == 0){
  //   lower_bound_1d = 1;
  //   lower_bound_2d = 1;
  // } else if (allele_list.size() == 1){
  //   lower_bound_1d = int32_t((read_len) / motif_len);
  //   lower_bound_2d = 1;
  // }
  // else if (allele_list.size() >= 2) {
  //   lower_bound_1d = int32_t((read_len) / motif_len);
  //   lower_bound_2d = int32_t((read_len - 2 * MARGIN) / motif_len - 1);
  // }

  lower_bound_1d = 1;
  lower_bound_2d = 1;

  if (ploidy == 2){

    for (std::vector<int32_t>::iterator allele_it = allele_list.begin();
         allele_it != allele_list.end();
         allele_it++) {
      nlopt_1D_optimize(read_len, motif_len, ref_count, lower_bound_1d, upper_bound, resampled, this, *allele_it, &a1, &result, &minf);
      sublist.push_back(a1);
      // cerr<<"ER:\t"<<*allele_it<<endl;
      // cerr<<"1D:\t"<<a1<<endl;
    }
    nlopt_2D_optimize(read_len, motif_len, ref_count, lower_bound_2d, upper_bound, resampled, this, &a1, &a2, &result, &minf);
    sublist.push_back(a1);
    sublist.push_back(a2);

    // cerr<<"2D:\t"<<a1<<endl;
    // cerr<<"2D:\t"<<a2<<endl;
    for (std::vector<int32_t>::iterator subl_it = sublist.begin();
         subl_it != sublist.end();
         subl_it++) {
      if(std::find(allele_list.begin(), allele_list.end(), *subl_it) == allele_list.end()) {
          /* allele_list does not contain this sublist item */
          allele_list.push_back(*subl_it);
      }
    }
    findBestAlleleListTuple(allele_list, read_len, motif_len, ref_count, resampled, ploidy, 0,
                            allele1, allele2, min_negLike);

  }
  else if (ploidy == 1){
    nlopt_1D_optimize(read_len, motif_len, ref_count, lower_bound_1d, upper_bound, resampled, this, fix_allele, &a1, &result, &minf);
    allele_list.push_back(a1);
    findBestAlleleListTuple(allele_list, read_len, motif_len, ref_count, resampled, ploidy, fix_allele,
                            allele1, allele2, min_negLike);

  }
  if (*allele1 > *allele2){
    temp = *allele1;
    *allele1 = *allele2;
    *allele2 = temp;
  }
  // if (!resampled){
  //   double gt_ll;
  //   for (int jj = 20; jj < 100; jj += 10){
  //     GetGenotypeNegLogLikelihood(*allele1, jj, read_len, motif_len, ref_count, resampled, &gt_ll);
  //     cerr<<*allele1<<","<<jj<<"\t"<<gt_ll<<endl;
  //   }
  // }
  // cerr<<endl<<*allele1<<"\t"<<*allele2<<"\t"<<*min_negLike<<endl;
  // if (!resampled){
  //   PlotLikelihood(30, 10, 100, 5, read_len, motif_len, ref_count);
  // }
  return true;    // TODO add false
}


bool LikelihoodMaximizer::findBestAlleleListTuple(std::vector<int32_t> allele_list,
                          int32_t read_len, int32_t motif_len, int32_t ref_count, bool resampled,
			  int32_t ploidy, int32_t fix_allele,
                          int32_t* allele1, int32_t* allele2, double* min_negLike){
  double gt_ll;
  *min_negLike = 1000000;
  int32_t best_a1 = 0, best_a2 = 0;
  if (ploidy == 2){
    for (std::vector<int32_t>::iterator a1_it = allele_list.begin();
            a1_it != allele_list.end();
            a1_it++){
      // if (!resampled)
      //   cerr<<*a1_it<<endl;
      for (std::vector<int32_t>::iterator a2_it = allele_list.begin();
            a2_it != allele_list.end();
            a2_it++){
        GetGenotypeNegLogLikelihood(*a1_it, *a2_it, read_len, motif_len, ref_count, resampled, &gt_ll);
        // if (!resampled)
        //   cerr<<endl<<*a1_it<<"\t"<<*a2_it<<"\t"<<gt_ll<<endl;
          if (gt_ll < *min_negLike){
            *min_negLike = gt_ll;
            best_a1 = *a1_it;
            best_a2 = *a2_it;
          }
      }
    }
  }
  else if (ploidy == 1){
    best_a2 = fix_allele;
    for (std::vector<int32_t>::iterator a1_it = allele_list.begin();
            a1_it != allele_list.end();
            a1_it++){
      GetGenotypeNegLogLikelihood(*a1_it, fix_allele, read_len, motif_len, ref_count, resampled, &gt_ll);
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
    if(!lm_ptr->GetGenotypeNegLogLikelihood(A, B, read_len, motif_len, ref_count, resampled, &gt_ll))
      return -1.0;
    else{
      return gt_ll;
    }
  }
  else{
    double A = x[0], B = fix_allele;
    if(!lm_ptr->GetGenotypeNegLogLikelihood(A, B, read_len, motif_len, ref_count, resampled, &gt_ll))
      return -1.0;
    else{
      return gt_ll;
    }
  
  }
}

bool nlopt_2D_optimize(const int32_t& read_len, const int32_t& motif_len,
               const int32_t& ref_count, const int32_t& lower_bound,
               const int32_t& upper_bound, const bool& resampled, LikelihoodMaximizer* lm_ptr,
               int32_t* allele1, int32_t* allele2, int32_t* ret_result, double* minf_ret) {
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


  nlopt_data data[1] = nlopt_data(read_len, motif_len, ref_count, lm_ptr, 0, resampled);

  opt.set_min_objective(nloptNegLikelihood, data);    // Change to max for maximization

  opt.set_xtol_rel(1e-5);   // TODO set something appropriate

  std::vector<double> xx(2);
  xx[0] = int32_t(1.1 * (read_len / motif_len));
  xx[1] = int32_t(1.2 * (read_len / motif_len));
  double minf;
  nlopt::result result = opt.optimize(xx, minf);
  *allele1 = int32_t(xx[0]);
  *allele2 = int32_t(xx[1]);
  *ret_result = result;
  *minf_ret = minf;
  return true;  // TODO add false
}

bool nlopt_1D_optimize(const int32_t& read_len, const int32_t& motif_len,
                const int32_t& ref_count, const int32_t& lower_bound,
                const int32_t& upper_bound, const bool& resampled, LikelihoodMaximizer* lm_ptr,
                const int32_t& fix_allele, int32_t* allele1,
                int32_t* ret_result, double* minf_ret) {

  nlopt::opt opt(nlopt::LN_COBYLA, 1);

  std::vector<double> lb(1);
  lb[0] = lower_bound;
  opt.set_lower_bounds(lb);

  std::vector<double> ub(1);
  ub[0] = upper_bound;
  opt.set_upper_bounds(ub);

  nlopt_data data[1] = nlopt_data(read_len, motif_len, ref_count, lm_ptr, fix_allele, resampled);
  
  opt.set_min_objective(nloptNegLikelihood, data);    // Change to max for maximization

  opt.set_xtol_rel(1e-5);   // TODO set something appropriate

  std::vector<double> xx(1);
  xx[0] = int32_t(1.1 * (read_len / motif_len));
  double minf;
  nlopt::result result = opt.optimize(xx, minf);
  *allele1 = int32_t(xx[0]);
  *ret_result = result;
  *minf_ret = minf;
return true;  // TODO add false
}


// /// GSL siman helper functions (not complete)
// double simanEnergy(void *xp){
//   siman_data *d = (siman_data *) xp;
//   int read_len  = d -> read_len;
//   int motif_len = d -> motif_len;
//   int ref_count = d -> ref_count;
//   int fix_allele = d -> fix_allele; 
//   int A = d -> A;
//   int B = d -> B;
//   LikelihoodMaximizer* lm_ptr = d -> lm_ptr;
//   double gt_ll;
//   if(!lm_ptr->GetGenotypeNegLogLikelihood(A, B, read_len, motif_len, ref_count, &gt_ll))
//       return -1.0;
//     else{
//       return gt_ll;
//     }
// }
// double simanMetric(void *xp, void *yp){
//   siman_data *d1 = (siman_data *) xp;
//   siman_data *d2 = (siman_data *) yp;
//   int A1 = d1 -> A;
//   int B1 = d1 -> B;
//   int A2 = d2 -> A;
//   int B2 = d2 -> B;

//   return sqrt(pow(A1 - A2, 2) + pow(B1 - B2, 2))
// }
// double siman_step(const gsl_rng * r, void *xp, double step_size){
//   siman_data *old_d = (siman_data *) xp;
//   int old_A = old_A -> A;
// }
// /// GSL siman helper functions (not complete)
