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

#ifndef SRC_READ_PAIR_H__
#define SRC_READ_PAIR_H__

#include "src/bam_io.h"
#include "src/realignment.h"

enum ReadType {
  RC_SPAN = 0,
  RC_ENCL = 1,
  RC_FRR = 2,
  RC_UNKNOWN = 3,
  RC_DISCARD = 4,
  RC_BOUND = 5
};

class ReadPair {
 public:
  ReadPair();
  virtual ~ReadPair();

  BamAlignment read1;
  BamAlignment read2;
  ReadType read_type;
  int32_t data_value;
  int32_t max_nCopy;    // maximum nCopy among two reads (used for flanking heuristic)
  bool found_pair;
};

#endif  // SRC_READ_PAIR_H__
