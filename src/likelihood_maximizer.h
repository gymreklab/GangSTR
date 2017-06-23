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

#ifndef SRC_LIKELIHOOD_MAXIMIZER_H__
#define SRC_LIKELIHOOD_MAXIMIZER_H__

#include "src/enclosing_class.h"
#include "src/frr_class.h"
#include "src/read_class.h"
#include "src/spanning_class.h"

class LikelihoodMaximizer {
 public:
  LikelihoodMaximizer();
  virtual ~LikelihoodMaximizer();

  // Clear data of all classes
  void Reset();
  // Add data points to each class
  void AddEnclosingData(const int32_t& data);
  void AddSpanningData(const int32_t& data);
  void AddFRRData(const int32_t& data);
  // Check data size
  std::size_t GetEnclosingDataSize();
  std::size_t GetSpanningDataSize();
  std::size_t GetFRRDataSize();
  // Main likelihood function
  bool GetGenotypeNegLogLikelihood(const int32_t& allele1, const int32_t& allele2, double* gt_ll);

 private:
  EnclosingClass enclosing_class_;
  FRRClass frr_class_;
  SpanningClass spanning_class_;

  // Class weights - TODO?
  const static double frr_weight_ = 0.8;
  const static double enclosing_weight_ = 1.0;
  const static double spanning_weight_ = 1.0;
};

#endif  // SRC_LIKELIHOOD_MAXIMIZER_H__
