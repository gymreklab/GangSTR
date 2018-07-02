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

#include "src/common.h"
#include "src/region_reader.h"
#include "src/stringops.h"

#include <iostream>
using namespace std;

RegionReader::RegionReader(const std::string& filename) {
  freader = new std::ifstream(filename.c_str());
  if (!freader->is_open()) {
    PrintMessageDieOnError("Could not open regions file", M_ERROR);
  }
}

/*
  Read next region from regions file.
  Return true if successful, false if no more regions
 */
bool RegionReader::GetNextRegion(Locus* locus) {
  std::string line;
  std::vector<std::string> items, offtarget_regions;
  std::string offtarget_str;
  bool stat;
  if (!std::getline(*freader, line)) {
    return false;
  }
  split_by_delim(line, '\t', items);
  
  if (items.size() < 5) {
    PrintMessageDieOnError("Regions file not formatted correctly", M_ERROR);
  }
  else if (items.size() == 6){
    offtarget_str = items[5];
    split_by_delim(offtarget_str, ',', offtarget_regions);
  }
  locus->chrom = items[0];
  locus->start = atoi(items[1].c_str());
  locus->end = atoi(items[2].c_str());
  locus->period = atoi(items[3].c_str());
  locus->motif = items[4];
  std::transform(locus->motif.begin(), locus->motif.end(), locus->motif.begin(), ::tolower);
  if (offtarget_regions.size() >= 1){
    for (int off_idx = 0; off_idx < offtarget_regions.size(); off_idx++)
      {
	std::string region = offtarget_regions[off_idx];
	int pos_col = region.find(':');
	int pos_dash = region.find('-');
	if (pos_col == -1 or pos_dash == -1)
	  continue;
	locus->offtarget_set = true;
	GenomeRegion offtarget;
	offtarget.chrom = region.substr(0, pos_col);
	offtarget.start = atoi(region.substr(pos_col + 1,\
					     pos_dash - pos_col - 1).c_str());
	offtarget.end = atoi(region.substr(pos_dash + 1).c_str());
	locus->offtarget_regions.push_back(offtarget);
      }
  }
  return true;
}

/*
  Reset region reader to the top of the regions file
*/
void RegionReader::Reset(){
  freader->clear();
  freader->seekg(0, ios::beg);
}
RegionReader::~RegionReader() {
  freader->close();
  delete freader;
}
