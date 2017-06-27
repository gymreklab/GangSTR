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

#include "common.h"
#include "src/bam_reader.h"

using namespace std;

GBamReader::GBamReader(std::vector<std::string> bamfiles) {
  reader = new BamCramMultiReader();
  // Try to open the files
  /*  if (!reader->Open(bamfiles)) {
    PrintMessageDieOnError("Could not open BAM files", M_ERROR);
  }
  if (!reader->LocateIndexes()) {
    reader->CreateIndexes();
    }*/
}

std::string GBamReader::GetTestRead() {
  BamAlignment aln;
  reader->GetNextAlignment(aln);
  return aln.QueryBases;
}

GBamReader::~GBamReader() {}
