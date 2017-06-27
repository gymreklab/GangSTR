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

#include "src/read_extractor.h"
#include "src/read_pair.h"

ReadExtractor::ReadExtractor() {}

/*
  Extracts relevant reads from bamfile and 
  populates data in the likelihood_maximizer
 */
bool ReadExtractor::ExtractReads(BamCramMultiReader* bamreader,
				 const Locus& locus,
				 LikelihoodMaximizer* likelihood_maximizer) {
  // This will keep track of information for each read pair
  std::map<std::string, ReadPair> read_pairs;

  // Get bam alignments from the relevant region
  bamreader->SetRegion(locus.chrom, locus.start-REGIONSIZE, locus.end+REGIONSIZE);

  // Go through each alignment in the region
  BamAlignment alignment;
  while (bamreader->GetNextAlignment(alignment)) {
    // Check if read's mate already processed
    // TODO

    // Discard read pair if position is irrelevant
    // (e.g. both mates on same side of the STR)
    // TODO

    // Check if read is spanning
    // TODO

  }

  // TODO
  return false;
}

ReadExtractor::~ReadExtractor() {}
