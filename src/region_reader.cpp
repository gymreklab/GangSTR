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

#include "src/region_reader.h"

using namespace std;

RegionReader::RegionReader(const std::string& filename) {
  // TODO
}

/*
  Read next region from regions file.
  Return true if successful, false if no more regions
 */
bool RegionReader::GetNextRegion(Locus* locus) {
  return false; // TODO
}

RegionReader::~RegionReader() {}
