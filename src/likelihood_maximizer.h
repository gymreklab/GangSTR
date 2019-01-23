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

#ifndef SRC_LIKELIHOOD_MAXIMIZER_H__
#define SRC_LIKELIHOOD_MAXIMIZER_H__

#include "src/enclosing_class.h"
#include "src/frr_class.h"
#include "src/flanking_class.h"
#include "src/options.h"
#include "src/read_class.h"
#include "src/spanning_class.h"
#include "src/read_pair.h"
#include "src/locus.h"
#include "src/sample_info.h"
#include "gsl/gsl_vector.h"
#include "gsl/gsl_rng.h"
#include "gsl/gsl_randist.h"

#include <iostream>
#include <fstream>
#include <string>
#include <numeric>

using namespace std;

// Struct for storing reads from all classes in a unified vector
struct ReadRecord{
  int32_t data;
  ReadType read_type;
};

class LikelihoodMaximizer {
 friend class Genotyper;
 public:
 LikelihoodMaximizer(const Options& _options, const SampleProfile& sp,
		     const int32_t& read_len);
 // LikelihoodMaximizer(const LikelihoodMaximizer& lm_obj); // copy constructor
  virtual ~LikelihoodMaximizer();

  // Clear data of all classes
  void Reset();
  // Add data points to each class
  void AddEnclosingData(const int32_t& data);
  void AddSpanningData(const int32_t& data);
  void AddFRRData(const int32_t& data);
  void AddFlankingData(const int32_t& data);
  void AddOffTargetData(const int32_t& data);
  // Check data size
  std::size_t GetEnclosingDataSize();
  std::size_t GetSpanningDataSize();
  std::size_t GetFRRDataSize();
  std::size_t GetFlankingDataSize();
  std::size_t GetOffTargetDataSize();
  std::size_t GetReadPoolSize();

  // Main likelihood function
  bool GetGenotypeNegLogLikelihood(const int32_t& allele1, const int32_t& allele2,
				   const bool& resampled, double* gt_ll);
  // Functions for dynamically setting grid size
  bool InferGridSize();
  void SetGridSize(const int32_t& min_allele, const int32_t max_allele);
  void GetGridSize(int32_t* min_allele, int32_t* max_allele);
  void InferAlleleList(std::vector<int32_t>* allele_list,
		       const int32_t& ploidy,
		       const bool& resampled, const int32_t& fix_allele);
  // Function for computing Q Score
  bool GetQScore(int32_t a1, int32_t a2, double* Q);
  // Functions for computing expansion probability
  bool GetExpansionProb(std::vector<double>* prob_vec, const int32_t& exp_threshold);
  bool GetNegLikelihoodSurface(const int32_t& a_lo,
			       const int32_t& a_hi,
			       const int32_t& b_lo,
			       const int32_t& b_hi,
			       const bool& resampled,
			       double* surfaceLL);
  // Main optimization function
  bool OptimizeLikelihood(const bool& resampled, const int32_t& use_ploidy,
			  const int32_t& fix_allele,
			  const double& off_share,
			  int32_t* allele1, int32_t* allele2, double* min_negLike);
  // Go over the list of the discovered alleles to find the best pair
  bool findBestAlleleListTuple(std::vector<int32_t> allele_list, int32_t use_ploidy,
			       bool resampled, int32_t fix_allele,
			       int32_t* allele1, int32_t* allele2, double* min_negLike);

  // Compute and return confidence interval with bootstrapping
  bool GetConfidenceInterval(const int32_t& allele1,
			     const int32_t& allele2,
			     const Locus& locus,
			     double* lob1, double* hib1, double* lob2, double* hib2,
			     double* a1_se, double* a2_se);

  // Set per-locus params
  void SetLocusParams(const STRLocusInfo& sli, const double& cov,
		      const int32_t& _read_len, const int32_t _motif_len,
		      const int32_t& _ref_count);

  // Print read pool
  void PrintReadPool();

  // Resample read pool with replacement
  void ResampleReadPool();

 protected:
  // Other params -> Made public for gslNegLikelihood to have access
  const Options* options;
 private:
  double obj_cov; // TODO: This is a placeholder, until we figure out how to pass coverage through sample info
  EnclosingClass enclosing_class_;
  FRRClass frr_class_;
  SpanningClass spanning_class_;
  FlankingClass flanking_class_;
  FRRClass offtarget_class_;
  std::vector<ReadRecord> read_pool;
  std::vector<ReadRecord> resampled_pool;
  EnclosingClass resampled_enclosing_class_;
  FRRClass resampled_frr_class_;
  SpanningClass resampled_spanning_class_;
  FlankingClass resampled_flanking_class_;

  // Write bootstrap samples to file
  ofstream bsfile_;
  // Random number generator
  gsl_rng * r;
  // percentage of off-target reads
  double offtarget_share;

  // Grid size
  int32_t lower_bound;
  int32_t upper_bound;
  bool grid_set;
  int32_t grid_buffer;
  int32_t grid_opt_threshold; // Above this, use nlopt rather than brute force grid search

  // Locus info
  bool locus_params_set;
  int32_t read_len;
  int32_t motif_len;
  int32_t ref_count;
};

// Helper struct for NLOPT gradient optimizer
struct nlopt_data{
  LikelihoodMaximizer* lm_ptr;
  int32_t read_len, motif_len, ref_count, n_dim, fix_allele;
  bool resampled;
  nlopt_data(int32_t Read_Len, int32_t Motif_Len, int32_t Ref_Count, 
    LikelihoodMaximizer* LM_OBJ, int32_t Fix_Allele, bool Resampled) : 
    read_len(Read_Len), motif_len(Motif_Len), ref_count(Ref_Count), 
    lm_ptr(LM_OBJ), fix_allele(Fix_Allele), resampled(Resampled) {
    }
};
// 1D gradient optimizer using NLOPT
bool nlopt_1D_optimize(const int32_t& read_len, const int32_t& motif_len,
		       const int32_t& ref_count, const int32_t& lower_bound,
		       const int32_t& upper_bound, const bool& resampled, 
		       const int& seed, LikelihoodMaximizer* lm_ptr,
		       const int32_t& fix_allele, int32_t* allele1,
		       int32_t* ret_result, double* minf_ret);
// 2D gradient optimizer using NLOPT
bool nlopt_2D_optimize(const int32_t& read_len, const int32_t& motif_len,
		       const int32_t& ref_count, const int32_t& lower_bound,
		       const int32_t& upper_bound, const bool& resampled, 
		       const int& seed, LikelihoodMaximizer* lm_ptr,
               int32_t* allele1, int32_t* allele2, int32_t* ret_result, double* minf_ret);
// Helper function for NLOPT gradient optimizer
double nloptNegLikelihood(unsigned n, const double *x, double *grad, void *data);




// // Helper struct for GSL simulated annealing
// struct siman_data{
//   LikelihoodMaximizer* lm_ptr;
//   int read_len, motif_len, ref_count, n_dim, fix_allele, A, B;
//   nlopt_data(int Read_Len, int Motif_Len, int Ref_Count, LikelihoodMaximizer* LM_OBJ, int Fix_Allele) : 
//     read_len(Read_Len), motif_len(Motif_Len), ref_count(Ref_Count), lm_ptr(LM_OBJ), fix_allele(Fix_Allele) {
//     }
// }
#endif  // SRC_LIKELIHOOD_MAXIMIZER_H__
