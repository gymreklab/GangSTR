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

#include <stdlib.h>

#include <iostream>
#include <sstream>
#include <algorithm>

#include "src/common.h"
#include "src/ref_genome.h"

using namespace std;

RefGenome::RefGenome(const std::string& _reffa) {
  // Check if file exists
  if (!file_exists(_reffa)) {
    PrintMessageDieOnError("FASTA file " + _reffa + " does not exist", M_ERROR);
  }

  // Check for index
  if (!file_exists(_reffa + ".fai")) {
    PrintMessageDieOnError("No index for FASTA file " + _reffa, M_ERROR);
  }

  // Load index
  refindex = fai_load(_reffa.c_str());
  if (refindex == NULL) {
    PrintMessageDieOnError("Failed to load FASTA index", M_ERROR);
  }
}

bool RefGenome::GetSequence(const std::string& _chrom,
			    const int32_t& _start,
			    const int32_t& _end,
			    std::string* seq) {
  int length;
  char* result = faidx_fetch_seq(refindex, _chrom.c_str(), _start, _end, &length);
  if (result == NULL) {
    stringstream ss;
    ss << "Error fetching reference sequence for " << _chrom << ":" << _start;
    PrintMessageDieOnError(ss.str(), M_ERROR);
  }
  seq->assign(result, length);
  std::transform(seq->begin(), seq->end(), seq->begin(), ::tolower);
  free((void *)result);
  return true;
}

RefGenome::~RefGenome() {}
