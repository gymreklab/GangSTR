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

#include "src/likelihood_maximizer.h"

LikelihoodMaximizer::LikelihoodMaximizer() {}

void LikelihoodMaximizer::Reset() {
  enclosing_class_.Reset();
  frr_class_.Reset();
  spanning_class_.Reset();
}

void LikelihoodMaximizer::AddEnclosingData(const int32_t& data) {
  enclosing_class_.AddData(data);
}
void LikelihoodMaximizer::AddSpanningData(const int32_t& data) {
  spanning_class_.AddData(data);
}
void LikelihoodMaximizer::AddFRRData(const int32_t& data) {
  frr_class_.AddData(data);
}
std::size_t LikelihoodMaximizer::GetEnclosingDataSize() {
  return enclosing_class_.GetDataSize();
}
std::size_t LikelihoodMaximizer::GetSpanningDataSize() {
  return spanning_class_.GetDataSize();
}
std::size_t LikelihoodMaximizer::GetFRRDataSize() {
  return frr_class_.GetDataSize();
}


bool LikelihoodMaximizer::GetGenotypeNegLogLikelihood(const int32_t& allele1,
						      const int32_t& allele2,
						      double* gt_ll) {
  double frr_ll, span_ll, encl_ll = 0.0;
  frr_class_.GetClassLogLikelihood(allele1, allele2, &frr_ll);
  spanning_class_.GetClassLogLikelihood(allele1, allele2, &span_ll);
  enclosing_class_.GetClassLogLikelihood(allele1, allele2, &encl_ll);
  *gt_ll = -1*(frr_weight_*frr_ll + spanning_weight_*span_ll + enclosing_weight_*encl_ll);
}

bool LikelihoodMaximizer::OptimizeLikelihood(int32_t* allele1,
					     int32_t* allele2) {
  // TODO
  return false;
}
LikelihoodMaximizer::~LikelihoodMaximizer() {}
