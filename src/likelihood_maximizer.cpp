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

#include <gsl/gsl_multimin.h>
#include "src/likelihood_maximizer.h"
#include "src/mathops.h"
#include <iostream>

using namespace std;

LikelihoodMaximizer::LikelihoodMaximizer(const Options& _options) {
  options = &_options;
  enclosing_class_.SetOptions(*options);
  frr_class_.SetOptions(*options);
  spanning_class_.SetOptions(*options);
}

void LikelihoodMaximizer::Reset() {
  enclosing_class_.Reset();
  frr_class_.Reset();
  spanning_class_.Reset();
}

void LikelihoodMaximizer::AddEnclosingData(const int32_t& data) {
  enclosing_class_.AddData(data);
}
void LikelihoodMaximizer::AddSpanningData(const int32_t& data) {
  spanning_class_.AddData(data);
}
void LikelihoodMaximizer::AddFRRData(const int32_t& data) {
  frr_class_.AddData(data);
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

bool LikelihoodMaximizer::GetGenotypeNegLogLikelihood(const int32_t& allele1,
						      const int32_t& allele2,
						      const int32_t& read_len,
						      const int32_t& motif_len,
						      const int32_t& ref_count,
						      double* gt_ll) {
  double frr_ll, span_ll, encl_ll = 0.0;
  frr_class_.GetClassLogLikelihood(allele1, allele2, read_len, motif_len, ref_count, &frr_ll);
  spanning_class_.GetClassLogLikelihood(allele1, allele2, read_len, motif_len, ref_count, &span_ll);
  enclosing_class_.GetClassLogLikelihood(allele1, allele2, read_len, motif_len, ref_count, &encl_ll);
  *gt_ll = -1*(options->frr_weight*frr_ll +
	       options->spanning_weight*span_ll +
	       options->enclosing_weight*encl_ll);
}



bool LikelihoodMaximizer::OptimizeLikelihood(const int32_t& read_len, const int32_t& motif_len,
					     const int32_t& ref_count,
					     int32_t* allele1, int32_t* allele2) {
  size_t minim_dim = 2;
  const gsl_multimin_fminimizer_type * minim_type = gsl_multimin_fminimizer_nmsimplex2;
  gsl_multimin_fminimizer * minim_handle = gsl_multimin_fminimizer_alloc (minim_type, minim_dim);

  gsl_vector *ss, *x;
  gsl_multimin_function neg_log_like;

  size_t iter = 0;
  int status;
  double size;
  //LikelihoodMaximizer* par = new LikelihoodMaximizer(options); // Params = [read_len, motif_len, ref_count]
  /* Starting point */
  x = gsl_vector_alloc (2);
  gsl_vector_set (x, 0, 10.0);
  gsl_vector_set (x, 1, 30.0);

  /* Set initial step sizes to 1 */
  ss = gsl_vector_alloc (2);
  gsl_vector_set_all (ss, 4.0);

  /* Initialize method and iterate */
  neg_log_like.n = 2;
  neg_log_like.f = &(dummy_func);
  neg_log_like.params = this;

  minim_handle = gsl_multimin_fminimizer_alloc (minim_type, 2);
  gsl_multimin_fminimizer_set (minim_handle, &neg_log_like, x, ss);

  do
    {
      iter++;
      status = gsl_multimin_fminimizer_iterate(minim_handle);
      
      if (status) 
        break;

      size = gsl_multimin_fminimizer_size (minim_handle);
      status = gsl_multimin_test_size (size, 1e-2);

      if (status == GSL_SUCCESS)
        {
          printf ("converged to minimum at\n");
        }

      printf ("%5d %10.3e %10.3e f() = %7.3f size = %.3f\n", 
              iter,
              gsl_vector_get (minim_handle->x, 0), 
              gsl_vector_get (minim_handle->x, 1), 
              minim_handle->fval, size);
    }
  while (status == GSL_CONTINUE && iter < 100);
  
  gsl_vector_free(x);
  gsl_vector_free(ss);
  gsl_multimin_fminimizer_free (minim_handle);

  cout<<gsl_vector_get (minim_handle->x, 0)<<"\t"<<gsl_vector_get (minim_handle->x, 1);
  *allele1 = 2;
  // TODO
  return false;
}
LikelihoodMaximizer::~LikelihoodMaximizer() {}


double gslNegLikelihood(const gsl_vector *v, void *params)
{
  double A, B;
  LikelihoodMaximizer *p = (LikelihoodMaximizer *)params;
  double gt_ll;
  A = gsl_vector_get(v, 0);
  B = gsl_vector_get(v, 1);
 
  if(!p[0].GetGenotypeNegLogLikelihood(A, B, p[0].options->read_len, p[0].options->motif_len, p[0].options->ref_count, &gt_ll))
    return -1.0;
  else
    return gt_ll;
}