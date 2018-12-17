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

#ifndef SRC_REF_GENOME_H__
#define SRC_REF_GENOME_H__

#include "htslib/faidx.h"

#include <stdint.h>
#include <unistd.h>

#include <string>

class RefGenome {
 public:
  RefGenome(const std::string& _reffa);
  virtual ~RefGenome();

  bool GetSequence(const std::string& _chrom,
		   const int32_t& _start,
		   const int32_t& _end,
		   std::string* seq) const;

  const std::vector<std::string> GetChroms() const;
  const int32_t GetChromSize(const std::string& chrom) const;

 private:
  bool file_exists(std::string path) const {
    return (access(path.c_str(), F_OK) != -1);
  }

  faidx_t* refindex;
};

#endif  // SRC_REF_GENOME_H__
