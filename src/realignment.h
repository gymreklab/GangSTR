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

#ifndef SRC_REALIGNMENT_H__
#define SRC_REALIGNMENT_H__

#include <string>

#include <stdint.h>

bool expansion_aware_realign(const std::string& seq,
			     const std::string& pre_flank,
			     const std::string& post_flank,
			     const std::string& motif,
			     int32_t* nCopy, int32_t* pos, double* score);
			     

enum SingleReadType {
  SR_PREFLANK = 0,
  SR_POSTFLANK = 1,
  SR_ENCLOSING = 2,
  SR_IRR = 3,
  SR_UNKNOWN = 4
};

bool classify_realigned_read(const std::string& seq,
			     const std::string& motif,
			     const int32_t& nCopy,
			     const double& score,
			     const int32_t read_length,
			     SingleReadType* single_read_class);

#endif  // SRC_REALIGNMENT_H__
