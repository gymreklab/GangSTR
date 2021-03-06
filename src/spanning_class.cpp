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

#include "src/mathops.h"
#include "src/spanning_class.h"
#include <iostream>
#include <iomanip>
#include <math.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_cdf.h>
using namespace std;

bool SpanningClass::GetLogClassProb(const int32_t& allele,
				    const int32_t& read_len, const int32_t& motif_len,
				double* log_class_prob) {
	int str_len = allele * motif_len;					// (L)
	double norm_const;
	/*
	if (!hist_mode){
	  norm_const = gsl_cdf_gaussian_P(2 * flank_len + str_len - dist_mean, dist_sdev) -
	    gsl_cdf_gaussian_P(2 * read_len - dist_mean, dist_sdev);
	}
	else{
	  norm_const = InsertSizeCDF(2 * flank_len + str_len) - InsertSizeCDF(2 * read_len);
	}
	*/
	norm_const = InsertSizeCDF(2 * flank_len + str_len) - InsertSizeCDF(2 * read_len);
	
	if (norm_const == 0 or 
	    double(2 * flank_len + str_len - 2 * read_len) == 0){
	  cerr << "SpanClassProb::Divide by Zero prevented!" << endl;
	  *log_class_prob = NEG_INF;
	  return true;
	}
	//gsl_ran_gaussian_pdf

	// int rv_dist = norm(loc = dist_mean, scale = dist_sdev)

	double coef0 = 1.0 / norm_const / double(2 * flank_len + str_len - 2 * read_len);

	double coef1 = double(dist_mean - str_len);
	double coef2 = - double(dist_sdev * dist_sdev);

	double term1, term2;
	if (2 * read_len >= str_len){
	  /*
	  if (!hist_mode){
	    term1 = gsl_cdf_gaussian_P(2 * flank_len + str_len - dist_mean, dist_sdev) - 
	      gsl_cdf_gaussian_P(2 * read_len - dist_mean, dist_sdev);
	    term2 = gsl_ran_gaussian_pdf(2 * flank_len + str_len - dist_mean, dist_sdev) -
	      gsl_ran_gaussian_pdf(2 * read_len - dist_mean, dist_sdev);
	  }
	  else{
	    term1 = InsertSizeCDF(2 * flank_len + str_len) - InsertSizeCDF(2 * read_len);
	    term2 = InsertSizePDF(2 * flank_len + str_len) - InsertSizePDF(2 * read_len);
	  }
	  */
	  term1 = InsertSizeCDF(2 * flank_len + str_len) - InsertSizeCDF(2 * read_len);
	  term2 = InsertSizePDF(2 * flank_len + str_len) - InsertSizePDF(2 * read_len);
	}
	else{
	  /*
	  if (!hist_mode){
	    term1 = gsl_cdf_gaussian_P(2 * flank_len + str_len - dist_mean, dist_sdev) - 
	      gsl_cdf_gaussian_P(str_len - dist_mean, dist_sdev);
	    term2 = gsl_ran_gaussian_pdf(2 * flank_len + str_len - dist_mean, dist_sdev) -
	      gsl_ran_gaussian_pdf(str_len - dist_mean, dist_sdev);
	  }
	  else{
	    term1 = InsertSizeCDF(2 * flank_len + str_len) - InsertSizeCDF(str_len);
	    term2 = InsertSizePDF(2 * flank_len + str_len) - InsertSizePDF(str_len);
	  }
	  */
	  term1 = InsertSizeCDF(2 * flank_len + str_len) - InsertSizeCDF(str_len);
	  term2 = InsertSizePDF(2 * flank_len + str_len) - InsertSizePDF(str_len);
	}
        
	double class_prob = coef0 * (coef1 * term1 + coef2 * term2);
	
	// cout<<endl<<class_prob<<" "<<coef0<<" "<<coef1<<" "<<coef2;
	// cout<<endl<<class_prob<<" "<<term1<<" "<<term2<<endl;

	if (class_prob > 0){
		*log_class_prob = log(class_prob);
		return true;
	}
	else if (class_prob <= 0 and class_prob > -1){ // accounting for computational errors ~0
		*log_class_prob = NEG_INF;
		return true;
	}
	else{
		return false;
	}
}

bool SpanningClass::GetLogReadProb(const int32_t& allele,
				   const int32_t& data,
				   const int32_t& read_len,
				   const int32_t& motif_len,
				   const int32_t& ref_count,
				   double* log_allele_prob) {
  
  int shift = motif_len * (allele - ref_count);
  int mean_A = dist_mean - shift;
  double allele_prob = 0.0;

  /*
  if (!hist_mode){
    allele_prob = gsl_ran_gaussian_pdf(data - mean_A, dist_sdev);
  }
  else {
    allele_prob = InsertSizePDF(data + shift);
  }
  */
  
  allele_prob = InsertSizePDF(data + shift);

  /*
  if (gsl_cdf_gaussian_P(motif_len * allele - mean_A, dist_sdev) < 1.0){
    allele_prob = 1.0 / (1.0 - gsl_cdf_gaussian_P(motif_len * allele - mean_A, dist_sdev)) * gsl_ran_gaussian_pdf(data - mean_A, dist_sdev);
  }
  else { // allele is too large to use spanning reads anyway.
    allele_prob = 0.0;
  }
  */

  //cerr << allele << "\t" << allele_prob << "\t" << InsertSizeCDF(data - shift) << endl;  
  if (allele_prob > 0){
    *log_allele_prob = log(allele_prob);
    //cerr << allele << " " << data << " " << *log_allele_prob << endl;
    return true;
  }
  else if (allele_prob == 0){
    *log_allele_prob = NEG_INF;
    //cerr << allele << " " << data << " " << *log_allele_prob << endl;
    return true;
  }
  else
    return false;
}

bool SpanningClass::GetGridBoundaries(int32_t* min_allele, int32_t* max_allele) {
  return false; // Don't use spanning reads to decide boundaries
}
