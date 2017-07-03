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

#ifndef SRC_OPTIONS_H__
#define SRC_OPTIONS_H__

#include <vector>
#include <string>

#include <stdint.h>

class Options {
 public:
  Options();
  virtual ~Options();

  // Input/output paths
  std::vector<std::string> bamfiles;
  std::string reffa;
  std::string regionsfile;
  std::string outprefix;
  // Insert sizes
  int32_t dist_mean;
  int32_t dist_sdev;
  // Stutter model - TODO later make per locus model
  double stutter_up;
  double stutter_down;
  double stutter_p;
  // Additional constants
  int32_t realignment_flanklen; // flank length to use for local realignment
  int32_t flanklen; // used in models
  int32_t regionsize; // region in bam file to search for reads around the STR
  bool verbose;
};

#endif  // SRC_OPTIONS_H__
