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

using namespace std;

ReadClass::ReadClass() {}

void ReadClass::AddData(const int32_t& data) {
  read_class_data_.push_back(data);
}

bool ReadClass::GetClassLogLikelihood(const int32_t& allele1,
				   const int32_t& allele2,
				   float* class_ll) {
  *class_ll = 0;
  float samp_log_likelihood, a1_ll, a2_ll;
  for (std::vector<int32_t>::iterator data_it = read_class_data_.begin();
       data_it != read_class_data_.end();
       data_it++) {
    if (!GetAlleleLogLikelihood(allele1, *data_it, &a1_ll)) {
      return false;
    }
    if (!GetAlleleLogLikelihood(allele2, *data_it, &a2_ll)) {
      return false;
    }
    *class_ll += fast_log_sum_exp(log10(allele1_weight_)+a1_ll, log10(allele2_weight_)+a2_ll);
  }
  return true;
}

bool ReadClass::GetAlleleLogLikelihood(const int32_t& allele,
				       const int32_t& data,
				       float* allele_ll) {
  float log_class_prob, log_read_prob;
  if (!GetLogClassProb(allele, &log_class_prob)) {
    return false;
  }
  if (!GetLogReadProb(allele, data, &log_read_prob)) {
    return false;
  }
  *allele_ll = log_class_prob + log_read_prob;
  return true;
}

bool ReadClass::GetLogClassProb(const int32_t& allele,
				  float* log_class_prob) {
  return false; // Implement in child classes
}

bool ReadClass::GetLogReadProb(const int32_t& allele,
				 const int32_t& data,
				 float* log_allele_prob) {
  return false; // Implement in child classes
}

ReadClass::~ReadClass() {}
