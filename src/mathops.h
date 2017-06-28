#ifndef SRC_MATHOPS_H__
#define SRC_MATHOPS_H__

#include <math.h>
#include "src/fastonebigheader.h"

#include <vector>

// Testing gsl
double TestGSL();

// To accelerate logsumexp, ignore values if they're 1/1000th or less than the maximum value
const double LOG_THRESH = log(0.001);

double fast_log_sum_exp(double log_v1, double log_v2);
double normal_cdf(double mean, double stdev, double x);
#endif  // SRC_MATHOPS_H__
