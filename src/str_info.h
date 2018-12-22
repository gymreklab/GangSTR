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

#include <map>
#include <string>
#include <utility>

#include <stdint.h>

#include "src/options.h"

#ifndef STR_INFO_H_
#define STR_INFO_H_

struct STRLocusInfo {
  int32_t exp_thresh;
  double stutter_up;
  double stutter_down;
  double stutter_p;
};

class STRInfo {
 public:
  STRInfo(const Options& options);
  ~STRInfo();

  const int32_t GetExpansionThreshold(const std::string& chrom,
				      const int32_t& start);
  const double GetStutterUp(const std::string& chrom,
			    const int32_t& start);
  const double GetStutterDown(const std::string& chrom,
			      const int32_t& start);
  const double GetStutterP(const std::string& chrom,
			   const int32_t& start);

  const STRLocusInfo GetSTRInfo(const std::string& chrom, const int32_t& start);

 private:
  std::map<std::pair<std::string,int32_t>, STRLocusInfo> str_info;
  STRLocusInfo default_info;
};

#endif

