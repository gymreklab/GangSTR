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
#include "src/frr_class.h"

#include <math.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_cdf.h>
#include <iostream>

using namespace std;

bool FRRClass::GetLogClassProb(const int32_t& allele,
			       const int32_t& read_len, const int32_t& motif_len,
			       double* log_class_prob) {
	int str_len = allele * motif_len;

	if (str_len < read_len){		// condition: L > r for this read to be possible
		*log_class_prob = NEG_INF;
		return true;
	}
	// Compute normalization constant norm_const
	//	double norm_const = gsl_cdf_gaussian_P(2 * flank_len + str_len - dist_mean, dist_sdev) -
	//						gsl_cdf_gaussian_P(2 * read_len - dist_mean, dist_sdev); 
	double norm_const = InsertSizeCDF(2 * flank_len + str_len) - InsertSizeCDF(2 * read_len); 
	if (norm_const == 0 or
	    (2.0 * flank_len + str_len - 2.0 * read_len) == 0){
	  cerr << "FRRClassProb::Divide by Zero prevented!" << endl;
	  *log_class_prob = NEG_INF;
	  return true;
	}
	double coef0 = 1.0 / norm_const / (2.0 * flank_len + str_len - 2.0 * read_len);
	double coef1 = - double(dist_sdev ^ 2);
	//	double term1 = gsl_ran_gaussian_pdf(str_len - dist_mean, dist_sdev) -
	//					gsl_ran_gaussian_pdf(2 * read_len - dist_mean, dist_sdev);
	double term1 = InsertSizePDF(str_len) -	InsertSizePDF(2 * read_len);
	double coef2 = dist_mean - read_len;
	//	double term2 = gsl_cdf_gaussian_P(str_len - dist_mean, dist_sdev) - 
	//				gsl_cdf_gaussian_P(2 * read_len - dist_mean, dist_sdev);
	double term2 = InsertSizeCDF(str_len) - InsertSizeCDF(2 * read_len);
	double coef3 = str_len - read_len;
	//double term3 = gsl_cdf_gaussian_P(2 * flank_len + str_len - dist_mean, dist_sdev) - 
	//				gsl_cdf_gaussian_P(str_len - dist_mean, dist_sdev);
	double term3 = InsertSizeCDF(2 * flank_len + str_len) - InsertSizeCDF(str_len);
	double class_prob;
	if (str_len >= 2 * read_len)
		class_prob = coef0 * (coef1 * term1 + coef2 * term2 + coef3 * term3);
	else
		class_prob = coef0 * coef3 * term3;

	if (class_prob > 0){
		*log_class_prob = log(class_prob);
		return true;
	}
	else if (class_prob == 0){
		*log_class_prob = NEG_INF;
		return true;
	}
	else
		return false;
}

bool FRRClass::GetLogReadProb(const int32_t& allele,
			      const int32_t& data,
			      const int32_t& read_len,
			      const int32_t& motif_len,
			      const int32_t& ref_count,
			      double* log_allele_prob) {
	int str_len = allele * motif_len;

	if (str_len < read_len){		// redundant -> remove --> not really, because when class_prob=neg_inf this may lead to read_prob->false
		*log_allele_prob = NEG_INF;
		return true;
	}

	// Compute normalization constant norm_const
	//double norm_const = gsl_cdf_gaussian_P(2 * flank_len + str_len - dist_mean, dist_sdev) -
	//					gsl_cdf_gaussian_P(2 * read_len - dist_mean, dist_sdev); 
	double norm_const = InsertSizeCDF(2 * flank_len + str_len) - InsertSizeCDF(2 * read_len); 

	//double term1 = gsl_cdf_gaussian_P(read_len + data + str_len - dist_mean, dist_sdev) - 
	//		gsl_cdf_gaussian_P(2 * read_len + data - dist_mean, dist_sdev);
	double term1 = InsertSizeCDF(read_len + data + str_len) - InsertSizeCDF(2 * read_len + data);
	if (norm_const == 0){
	  *log_allele_prob = NEG_INF;
	  return true;
	}
	double allele_prob = 1 / norm_const * term1;
	if (allele_prob > 0){
		*log_allele_prob = log(allele_prob);
		return true;
	}
	else if (allele_prob == 0){
		*log_allele_prob = NEG_INF;
		return true;
	}
	else
		return false;
}

bool FRRClass::GetCountLogLikelihood(const int32_t& allele1,
				     const int32_t& allele2,
				     const int32_t& read_len,
				     const int32_t& motif_len,
				     const double& coverage,
				     const int32_t& ploidy,
				     const int32_t& offtarget_count,
				     double* count_ll){
  
  int32_t frr_count = GetDataSize() + offtarget_count;
  double frr_thresh = double(read_len) / double(motif_len);
  if (read_len == 0){
    cerr<< allele1 * motif_len << endl;
    cerr << "FRRCountProb::Divide by Zero prevented!" << endl;
    *count_ll = 0;
    return true;
  }
  double exp_count1 = coverage / 2.0 / double(read_len) * double(allele1 * motif_len - read_len);
  double exp_count2 = coverage / 2.0 / double(read_len) * double(allele2 * motif_len - read_len); 
  double lambda = (allele1 >= frr_thresh ? exp_count1 : 0) + 
    (allele2 >= frr_thresh ? exp_count2: 0);// Poisson parameter: Total expected number of FRRs
  
    
  double prob = 0;
  if (lambda <= 0 || allele1 <= 0 || allele2 <= 0 || frr_count <= 0){
    *count_ll = NEG_INF;
    return true;
  }
  else{
    //    *count_ll = log(pow(lambda, frr_count) * exp(-lambda) / tgamma(frr_count + 1));
    if (frr_count < 10){
      *count_ll = frr_count * log(lambda) - lambda - log(tgamma(frr_count + 1));
    }
    else{
      // Approx for log of factorial by Srinivasa Ramanujan
      double n = frr_count;
      //double log_fact = n * log(n) - n + 
      //log(n * (1 + 4 * n * (1 + 2 * n))) / 6 + 
      //  log(3.141593) / 2;

      // Approx for log of factorial: Stirling method      
      double log_fact = (n + 0.5) * log(n) - n + 0.5 * log(2 * 3.141593);
      *count_ll = frr_count * log(lambda) - lambda - log_fact;
      //std::cerr<<endl;
      //std::cerr<< ">> "<<allele1<< ": "<< exp_count1<<", "<<allele2<<": "<< exp_count2<<endl;
      //std::cerr<< ">> "<<allele1<<", "<<allele2<<"\texp:"<<lambda<<"\treal:"<<frr_count << "\tll:" << log_fact<<endl;

    }
    
    return true;
  }
  
  return true;
}
