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
#include <algorithm>
#include <string>
#include <math.h>

#include "src/common.h"
#include "src/gc_region_reader.h"
#include "src/stringops.h"

#include <iostream>
using namespace std;

GCRegionReader::GCRegionReader(const RefGenome& refgenome,
			       const float& bin_size_, const int& region_len_,
			       const int& max_regions_) {
  bin_size = bin_size_;
  region_len = region_len_;
  max_regions = max_regions_;
  for (int i=0; i<=(1.0/bin_size); i++) {
    std::vector<Locus> loci;
    loci.clear();
    gc_bin_loci.push_back(loci);
  }
  int num_regions = 0;
  int jump_len = 50000; // jump this much between regions to sample the chrom
  const std::vector<std::string> chroms = refgenome.GetChroms();
  std::string region_seq;
  for (std::vector<std::string>::const_iterator it=chroms.begin();
       it != chroms.end(); it++) {
    const int32_t chrom_size = refgenome.GetChromSize(*it);
    for (int i=0; i<chrom_size; i+= region_len+jump_len) {
      if (!refgenome.GetSequence(*it, i, i+region_len, &region_seq)) {
	continue;
      }
      float gc = GetGC(region_seq);
      if (gc == -1) continue;
      Locus locus;
      locus.chrom = *it;
      locus.start = i;
      locus.end = i+region_len;
      gc_bin_loci[int(floor(gc/bin_size))].push_back(locus);
      num_regions += 1;
      if (num_regions >= max_regions) break;
    }
  }
}

float GCRegionReader::GetGC(const std::string& seq) {
  int total_bases = 0;
  int gc = 0;
  for (int i=0; i<seq.size(); i++) {
    if (seq[i] != 'c' && seq[i] != 'C' &&
	seq[i] != 'g' && seq[i] != 'G' &&
	seq[i] != 't' && seq[i] != 'T' &&
	seq[i] != 'a' && seq[i] != 'A') {
      return -1;
    }
    total_bases += 1;
    if (seq[i] == 'c' || seq[i] =='C' ||
	seq[i] == 'g' || seq[i] == 'G') {
      gc += 1;
    }
  }
  return float(gc)/float(total_bases);
}

bool GCRegionReader::GetGCBinLoci(std::vector<Locus>* loci,
				  const float& lb, const float& ub,
				  const int& regions_per_bin) {
  int bin = int(floor(float(lb)/bin_size));
  if (regions_per_bin < gc_bin_loci[bin].size()) {
    std::vector<Locus> l(gc_bin_loci[bin].begin(), gc_bin_loci[bin].begin()+regions_per_bin);
    *loci = l;
  } else {
    *loci = gc_bin_loci[bin];
  }
  return true;
}

GCRegionReader::~GCRegionReader() {}
