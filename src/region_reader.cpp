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

#include "src/common.h"
#include "src/region_reader.h"

using namespace std;

RegionReader::RegionReader(const std::string& filename) {
  freader = new std::ifstream(filename.c_str());
  if (!freader->is_open()) {
    PrintMessageDieOnError("Could not open regions file", ERROR);
  }
}

/*
  Read next region from regions file.
  Return true if successful, false if no more regions
 */
bool RegionReader::GetNextRegion(Locus* locus) {
  std::string line;
  std::vector<std::string> items;
  if (!std::getline(*freader, line)) {
    return false;
  }
  split_by_delim(line, '\t', items);
  if (items.size() < 4) {
    PrintMessageDieOnError("Regions file not formatted correctly", ERROR);
  }
  locus->chrom = items[0];
  locus->start = atoi(items[1].c_str());
  locus->end = atoi(items[2].c_str());
  locus->period = atoi(items[3].c_str());
  return true;
}

RegionReader::~RegionReader() {
  freader->close();
}
