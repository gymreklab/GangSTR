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

#ifndef SRC_READ_EXTRACTOR_H__
#define SRC_READ_EXTRACTOR_H__

#include "src/bam_io.h"
#include "src/locus.h"
#include "src/likelihood_maximizer.h"

class ReadExtractor {
  const static int32_t REGIONSIZE = 5000; // TODO what should this be
 public:
  ReadExtractor();
  virtual ~ReadExtractor();

  // Main function to extract reads of each class
  bool ExtractReads(BamCramMultiReader* bamreader,
		    const Locus& locus,
		    LikelihoodMaximizer* likelihood_maximizer);
 private:
  // Trim alignment read names
  std::string trim_alignment_name(const BamAlignment& aln) const;
  // Check if read should be discarded
  bool FindDiscardedRead(BamAlignment alignment,
			 const int32_t& chrom_ref_id,
			 const Locus& locus);
  // Check if read is spanning class
  bool FindSpanningRead(BamAlignment alignment,
			const int32_t& chrom_ref_id,
			const Locus& locus,
			int32_t* insert_size);
};

#endif  // SRC_READ_EXTRACTOR_H__
