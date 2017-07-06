#include "mathops.h"

#include <gsl/gsl_math.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

double dummy_func (const gsl_vector *v, void *params)
{
  double A;
  double *p = (double *)params;

  A = gsl_vector_get(v, 0);

  return pow(A - p[0], p[1]);
}
double TestGSL() {
  return gsl_ran_gaussian_pdf(0, 1);
}

double fast_log_sum_exp(double log_v1, double log_v2){
  if (log_v1 > log_v2){
    double diff = log_v2-log_v1;
    return diff < LOG_THRESH ? log_v1 : log_v1 + fastlog(1 + fastexp(diff));
  }
  else {
    double diff = log_v1-log_v2;
    return diff < LOG_THRESH ? log_v2 : log_v2 + fastlog(1 + fastexp(diff));
  }
}

double normal_cdf(double mean, double stdev, double x){
	
}
