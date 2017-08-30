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
#include "gsl/gsl_vector.h"
#include <string>

class LikelihoodMaximizer {
 friend class Genotyper;
 public:
  LikelihoodMaximizer(Options& _options);
  // LikelihoodMaximizer(const LikelihoodMaximizer& lm_obj); // copy constructor
  virtual ~LikelihoodMaximizer();

  // Clear data of all classes
  void Reset();
  // Add data points to each class
  void AddEnclosingData(const int32_t& data);
  void AddSpanningData(const int32_t& data);
  void AddFRRData(const int32_t& data);
  void AddFlankingData(const int32_t& data);
  // Check data size
  std::size_t GetEnclosingDataSize();
  std::size_t GetSpanningDataSize();
  std::size_t GetFRRDataSize();
  // Main likelihood function
  bool GetGenotypeNegLogLikelihood(const int32_t& allele1, const int32_t& allele2,
				   const int32_t& read_len, const int32_t& motif_len,
				   const int32_t& ref_count,
				   double* gt_ll);
  // Main optimization function - TODO also return other data
  bool OptimizeLikelihood(const int32_t& read_len, const int32_t& motif_len,
			  const int32_t& ref_count,
			  int32_t* allele1, int32_t* allele2, double* min_negLike);
  // Go over the list of the discovered alleles to find the best pair
  bool findBestAlleleListTuple(std::vector<int32_t> allele_list,
                          int32_t read_len, int32_t motif_len, int32_t ref_count,
                          int32_t* allele1, int32_t* allele2, double* min_negLike);

  // Update read class options
  void UpdateOptions();
 protected:
  // Other params -> Made public for gslNegLikelihood to have access
  Options* options;

 private:
  EnclosingClass enclosing_class_;
  FRRClass frr_class_;
  SpanningClass spanning_class_;
  FlankingClass flanking_class_;
};

// Helper struct for NLOPT gradient optimizer
struct nlopt_data{
  LikelihoodMaximizer* lm_ptr;
  int read_len, motif_len, ref_count, n_dim, fix_allele;
  nlopt_data(int Read_Len, int Motif_Len, int Ref_Count, LikelihoodMaximizer* LM_OBJ, int Fix_Allele) : 
    read_len(Read_Len), motif_len(Motif_Len), ref_count(Ref_Count), lm_ptr(LM_OBJ), fix_allele(Fix_Allele) {
    }
};
// 1D gradient optimizer using NLOPT
bool nlopt_1D_optimize(const int32_t& read_len, const int32_t& motif_len,
                const int32_t& ref_count, const int32_t& lower_bound,
                const int32_t& upper_bound, LikelihoodMaximizer* lm_ptr,
                const int32_t& fix_allele, int32_t* allele1,
                int32_t* ret_result, double* minf_ret);
// 2D gradient optimizer using NLOPT
bool nlopt_2D_optimize(const int32_t& read_len, const int32_t& motif_len,
               const int32_t& ref_count, const int32_t& lower_bound,
               const int32_t& upper_bound, LikelihoodMaximizer* lm_ptr,
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
