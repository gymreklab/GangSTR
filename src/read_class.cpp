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
#include "src/read_class.h"

#include <math.h>
#include <iostream>
#include <iomanip>
#include <typeinfo>
using namespace std;

ReadClass::ReadClass() {
  // Set default options
  Options default_options;
  SetOptions(default_options);
}

void ReadClass::SetOptions(const Options& options) {
  dist_mean = options.dist_mean;
  dist_sdev = options.dist_sdev;
  flank_len = options.flanklen;
  stutter_up = options.stutter_up;
  stutter_down = options.stutter_down;
  stutter_p = options.stutter_p;
  read_prob_mode = options.read_prob_mode;
}

void ReadClass::AddData(const int32_t& data) {
  read_class_data_.push_back(data);
}

/*
  Calculates log P(read_class_data_ | <allele1, allele2>) and sets class_ll

  log P(data|<allelele1, allele2>) = sum_i log P(data_i | <allele1, allele2>)
  P(data_i | <allele1, allele2> = allele1_weight*P(data_i|allele1) + allele2_weight*P(data_i|allele2)

  Return false if something goes wrong.
 */
bool ReadClass::GetClassLogLikelihood(const int32_t& allele1,
				      const int32_t& allele2,
				      const int32_t& read_len, const int32_t& motif_len,
				      const int32_t& ref_count, const int32_t& ploidy,
				      double* class_ll) {
  *class_ll = 0;
  double samp_log_likelihood, a1_ll, a2_ll;
  for (std::vector<int32_t>::iterator data_it = read_class_data_.begin();
       data_it != read_class_data_.end();
       data_it++) {
    if (!GetAlleleLogLikelihood(allele1, *data_it, read_len, motif_len, ref_count, &a1_ll)) {
      return false;
    }
    if (!GetAlleleLogLikelihood(allele2, *data_it, read_len, motif_len, ref_count, &a2_ll)) {
      return false;
    }
    // TODO delete
    // cerr<<typeid(*this).name()<<"\t";
    // cerr<<*data_it<<"\t"<<fast_log_sum_exp(log(allele1_weight_)+a1_ll, log(allele2_weight_)+a2_ll)<<endl;
    if (ploidy == 2){
      *class_ll += fast_log_sum_exp(log(allele1_weight_)+a1_ll, log(allele2_weight_)+a2_ll);
    }
    else if (ploidy == 1){
      *class_ll += log(allele1_weight_) + a1_ll;
    }
  }
  return true;
}

/*
  Calculates logP(data_i|allele) and sets allele_ll

  Return false if something goes wrong.
 */
bool ReadClass::GetAlleleLogLikelihood(const int32_t& allele,
				       const int32_t& data,
				       const int32_t& read_len, const int32_t&  motif_len,
				       const int32_t& ref_count,
				       double* allele_ll) {
  double log_class_prob, log_read_prob;
  if (!GetLogClassProb(allele, read_len, motif_len, &log_class_prob)) {
    return false;
  }
  if (!GetLogReadProb(allele, data, read_len, motif_len, ref_count, &log_read_prob)) {
    return false;
  }
  if (read_prob_mode){
    *allele_ll = log_read_prob;
  }
  else{
    *allele_ll = log_class_prob + log_read_prob;
  }
  
  return true;
}

bool ReadClass::GetLogClassProb(const int32_t& allele,
				const int32_t& read_len, const int32_t& motif_len,
				double* log_class_prob) {
  return false; // Implement in child classes
}

bool ReadClass::GetLogReadProb(const int32_t& allele,
			       const int32_t& data,
			       const int32_t& read_len, const int32_t& motif_len,
			       const int32_t& ref_count,
			       double* log_allele_prob) {
  return false; // Implement in child classes
}

void ReadClass::Reset() {
  read_class_data_.clear();
}


std::size_t ReadClass::GetDataSize() {
  return read_class_data_.size();
}

ReadClass::~ReadClass() {}
