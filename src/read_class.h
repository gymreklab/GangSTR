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

class ReadClass {
 public:
  ReadClass();
  void AddData(const int32_t& data);
  bool GetClassLogLikelihood(const int32_t& allele1, const int32_t& allele2, float* class_ll);
  virtual ~ReadClass();

 private:
  bool GetAlleleLogLikelihood(const int32_t& allele, const int32_t& data, float* allele_ll);
  bool GetLogClassProb(const int32_t& allele, float* log_class_prob);
  bool GetLogReadProb(const int32_t& allele, const int32_t& data, float* log_allele_prob);

  std::vector<int32_t> read_class_data_; 
  const static float allele1_weight_ = 0.5;
  const static float allele2_weight_ = 0.5;
};

#endif  // SRC_READ_CLASS_H__
