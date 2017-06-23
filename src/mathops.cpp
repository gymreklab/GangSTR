#include "mathops.h"

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

