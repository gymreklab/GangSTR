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

  // Keep track of which file we're processing
  int32_t file_index = 0;
  std::string file_label = "0_";
  std::string prev_file = "";

  // Header has info about chromosome names
  const BamHeader* bam_header = bamreader->bam_header();
  const int32_t chrom_ref_id = bam_header->ref_id(locus.chrom);

  // Go through each alignment in the region
  BamAlignment alignment;
  while (bamreader->GetNextAlignment(alignment)) {
    // Check if we've moved to a different file
    if (prev_file.compare(alignment.Filename()) != 0) {
      prev_file = alignment.Filename();
      std::stringstream ss;
      ss << ++file_index << "_";
      file_label = ss.str();
    }

    // Set key to keep track of this mate pair
    std::string aln_key = file_label + trim_alignment_name(alignment);

    /*  Check if read's mate already processed */
    std::map<std::string, ReadPair>::iterator rp_iter = read_pairs.find(aln_key);
    if (rp_iter != read_pairs.end()) {
      if (rp_iter->second.read_type == RC_SPAN || 
	  rp_iter->second.read_type == RC_ENCL ||
	  rp_iter->second.read_type == RC_DISCARD) {
	// If spanning, enclosing, or discard already : nothing to do, just add read 2
	rp_iter->second.read2 = alignment;
      } else if (rp_iter->second.read_type == RC_FRR) {
	// If FRR: TODO
      } else {
	// If unknown: TODO
      }
      continue; // move on to next read
    }
  
    /* Discard read pair if position is irrelevant */
    if (FindDiscardedRead(alignment, chrom_ref_id, locus)) {
      ReadPair read_pair;
      read_pair.read_type = RC_DISCARD;
      read_pair.read1 = alignment;
      read_pairs.insert(std::pair<std::string, ReadPair>(aln_key, read_pair));
      continue;
    }

    /* Check if read is spanning */
    int32_t insert_size;
    if (FindSpanningRead(alignment, chrom_ref_id, locus, &insert_size)) {
      ReadPair read_pair;
      read_pair.read_type = RC_SPAN;
      read_pair.read1 = alignment;
      read_pair.data_value = insert_size;
      read_pairs.insert(std::pair<std::string, ReadPair>(aln_key, read_pair));
    }

    /* Check if read is enclosing */
    // TODO

    /* Check if read is FRR */
    // TODO
  }
  /*  Second pass through potential FRR reads */
  // TODO

  /* Load data into likelihood maximizer */
  for (std::map<std::string, ReadPair>::const_iterator iter = read_pairs.begin();
       iter != read_pairs.end(); iter++) {
    if (iter->second.read_type == RC_SPAN) {
      likelihood_maximizer->AddSpanningData(iter->second.data_value);
    } else if (iter->second.read_type == RC_ENCL) {
      likelihood_maximizer->AddEnclosingData(iter->second.data_value);
    } else if (iter->second.read_type == RC_FRR) {
      likelihood_maximizer->AddFRRData(iter->second.data_value);
    } else {
      continue;
    }
  }
  return false;
}

/*
  Check if read should be discarded
  Discard reads if both mates fall on same
  side of the STR.
  See 5.2_filter_spanning_only_core.py second case
  Return true if yes
*/
bool ReadExtractor::FindDiscardedRead(BamAlignment alignment,
				      const int32_t& chrom_ref_id,
				      const Locus& locus) {
  // Get read length
  int32_t read_length = (int32_t)alignment.QueryBases().size();

  // 5.2_filter_spanning_only_core.py 
  bool discard1 = alignment.RefID() == chrom_ref_id && alignment.Position() <= locus.start-read_length &&
    alignment.MateRefID() == chrom_ref_id && alignment.MatePosition() <= locus.start-read_length;
  bool discard2 = alignment.RefID() == chrom_ref_id && alignment.Position() >= locus.end &&
    alignment.MateRefID() == chrom_ref_id && alignment.MatePosition() >= locus.end;
  return discard1 || discard2;
}

/*
  Check if read is spanning
   See 5.2_filter_spanning_only_core.py
   Return true if yes
*/
bool ReadExtractor::FindSpanningRead(BamAlignment alignment,
				     const int32_t& chrom_ref_id,
				     const Locus& locus,
				     int32_t* insert_size) {
  // Get read length info
  int32_t read_length = (int32_t)alignment.QueryBases().size();

  // Similar to 5.2_filter_spanning_only_core.py:57  
  bool span1 = alignment.RefID() == chrom_ref_id && alignment.GetEndPosition() <= locus.start &&
    alignment.MateRefID() == chrom_ref_id && alignment.MatePosition() >= locus.end;
  bool span2 = alignment.RefID() == chrom_ref_id && alignment.Position() >= locus.end &&
    alignment.MateRefID() == chrom_ref_id && alignment.MatePosition() >= locus.start-read_length;
  if (span1 || span2) {
    *insert_size = abs(alignment.TemplateLength());
    return true;
  }

  // Similar to 5.2_filter_spanning_only_core.py:69
  // Allow flanking from both sides, but one side should be 
  // mapped to the correct location for this to work
  // TODO
  return false;
}

std::string ReadExtractor::trim_alignment_name(const BamAlignment& aln) const {
  std::string aln_name = aln.Name();
  if (aln_name.size() > 2){
    if (aln_name[aln_name.size()-2] == '/')
      aln_name.resize(aln_name.size()-2);
  }
  return aln_name;
}

ReadExtractor::~ReadExtractor() {}
