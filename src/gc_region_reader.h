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

#ifndef SRC_GC_REGION_READER_H__
#define SRC_GC_REGION_READER_H__

#include <fstream>
#include <string>
#include <map>
#include <vector>

#include "src/locus.h"
#include "src/ref_genome.h"

class GCRegionReader {
 public:
  GCRegionReader(const RefGenome& refgenome,
		 const float& bin_size_, const int& region_len_,
		 const int& max_regions);
  virtual ~GCRegionReader();

  bool GetGCBinLoci(std::vector<Locus>* loci,
		    const float& lb, const float& ub,
		    const int& regions_per_bin);
 private:
  float GetGC(const std::string& seq);
  
  std::vector<std::vector<Locus> > gc_bin_loci;
  float bin_size;
  int region_len;
  int max_regions;
};

#endif  // SRC_GC_REGION_READER_H__
