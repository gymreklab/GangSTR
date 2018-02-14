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
#include "src/options.h"
#include "src/read_pair.h"

#include <iostream>
#include <fstream>

#include <math.h>

class ReadExtractor {
  friend class ReadExtractorTest;
  friend class Genotyper;
 public:
  ReadExtractor(const Options& options_);
  virtual ~ReadExtractor();
    
  bool debug = false;
// bool print_read_data = false;

  // Main function to extract reads of each class
  bool ExtractReads(BamCramMultiReader* bamreader,
		    const Locus& locus,
		    const int32_t& regionsize,
		    const int32_t& min_match, 
		    LikelihoodMaximizer* likelihood_maximizer);

 protected:
  // Trim alignment read names
  std::string trim_alignment_name(const BamAlignment& aln) const;
  
  // Process all read pairs
  bool ProcessReadPairs(BamCramMultiReader* bamreader,
			const Locus& locus, 
			const int32_t& regionsize,
			const int32_t& min_match, 
			std::map<std::string, ReadPair>* read_pairs);

  // Implemented in BamInfoExtract. TODO delete
  // // Find insert size distribution
  // bool ComputeInsertSizeDistribution(BamCramMultiReader* bamreader,
  //      const Locus& locus,
  //      double* mean, double* std_dev, int32_t* read_len);

  // Check if read should be discarded
  bool FindDiscardedRead(BamAlignment alignment,
			 const int32_t& chrom_ref_id,
			 const Locus& locus);
  // Check if read is spanning class
  bool FindSpanningRead(BamAlignment alignment,
			const int32_t& chrom_ref_id,
			const Locus& locus,
			int32_t* insert_size);
  // Check single read overlapping repeat area
  bool ProcessSingleRead(BamAlignment alignment,
			 const int32_t& chrom_ref_id,
			 const Locus& locus,
			 const int32_t &min_match,
			 int32_t* data_value,
			 int32_t* nCopy_value,
			 int32_t* score_value,
			 ReadType* read_type,
			 SingleReadType* srt);
  // Rescue mate pairs aligned elsewhere
  bool RescueMate(BamCramMultiReader* bamreader,
		  BamAlignment alignment, BamAlignment* matepair);

private:
const Options options;
ofstream readfile_;
};

#endif  // SRC_READ_EXTRACTOR_H__
