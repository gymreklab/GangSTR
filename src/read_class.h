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

#ifndef SRC_READ_CLASS_H__
#define SRC_READ_CLASS_H__

#include <stdint.h>

#include <vector>

/*

Parent ReadClass
Individual read classes (FRR, enclosing, spanning) inherit
from this and implement their own read and class probability functions

A read class consists of:
- data (a vector of relevant values, e.g. copy number, insert size)
- a method to calculate the class log likelihood for a diploid genotype
 */
class ReadClass {
 public:
  ReadClass();
  virtual ~ReadClass();

  // Add a data point to the class data vector
  void AddData(const int32_t& data);
  // Calculate class log likelihood for diploid genotype P(data|<A,B>)
  bool GetClassLogLikelihood(const int32_t& allele1, const int32_t& allele2, double* class_ll);

 private:
  // Calculate log probability P(datapoint | allele)
  bool GetAlleleLogLikelihood(const int32_t& allele, const int32_t& data, double* allele_ll);
  // Calculate class probability for an allele - implemented in children classes
  bool GetLogClassProb(const int32_t& allele, double* log_class_prob);
  // Calculate read probability - implemented in children classes
  bool GetLogReadProb(const int32_t& allele, const int32_t& data, double* log_allele_prob);

  // Store vector of data for this class
  std::vector<int32_t> read_class_data_;
  // Allele weights. TODO: change if phasing available, would need per-read weights
  const static double allele1_weight_ = 0.5;
  const static double allele2_weight_ = 0.5;
};

#endif  // SRC_READ_CLASS_H__
