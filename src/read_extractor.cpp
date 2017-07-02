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

#include "src/stringops.h"
#include "src/read_extractor.h"
#include "src/realignment.h"

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
  if (!ProcessReadPairs(bamreader, locus, &read_pairs)) {
    return false;
  }

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
  return true;
}

/*
  Main function to decide what to do with each read pair
 */
bool ReadExtractor::ProcessReadPairs(BamCramMultiReader* bamreader,
				     const Locus& locus, 
				     std::map<std::string, ReadPair>* read_pairs) {
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
    if (debug) {
      std::cerr << "Processing " << alignment.Name() << std::endl;
    }
    // Check if we should skip this read
    if (alignment.IsSupplementary()) {
      continue;
    }
    
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
    std::map<std::string, ReadPair>::iterator rp_iter = read_pairs->find(aln_key);
    if (rp_iter != read_pairs->end()) {
      if (debug) {
	std::cerr << "Already found mate  " << alignment.Name()
		  << " " << rp_iter->second.read_type
		  << std::endl;
      }

      rp_iter->second.found_pair = true;
      rp_iter->second.read2 = alignment;
      // Only need to do something if class is still unknown
      if (rp_iter->second.read_type == RC_UNKNOWN) {
	if (debug) {
	  std::cerr << "Checking mate " << alignment.Name() << " " << alignment.QueryBases() << std:: endl;
	}

	int32_t data_value;
	ReadType read_type;
	ProcessSingleRead(alignment, chrom_ref_id, locus,
			  &data_value, &read_type);

	if (debug) {
	  std::cerr << "Mate found to be   " << read_type << std:: endl;
	}
	rp_iter->second.read_type = read_type;
	rp_iter->second.data_value = data_value;
      }
      continue; // move on to next read
    }
  
    /* Discard read pair if position is irrelevant */
    if (debug) {
      std::cerr << "Checking for discard" << std::endl;
    }
    if (FindDiscardedRead(alignment, chrom_ref_id, locus)) {
      ReadPair read_pair;
      read_pair.read_type = RC_DISCARD;
      read_pair.read1 = alignment;
      read_pairs->insert(std::pair<std::string, ReadPair>(aln_key, read_pair));
      continue;
    }

    /* Check if read is spanning */
    if (debug) {
      std::cerr << "Checking for spanning" << std::endl;
    }
    int32_t insert_size;
    if (FindSpanningRead(alignment, chrom_ref_id, locus, &insert_size)) {
      ReadPair read_pair;
      read_pair.read_type = RC_SPAN;
      read_pair.read1 = alignment;
      read_pair.data_value = insert_size;
      read_pairs->insert(std::pair<std::string, ReadPair>(aln_key, read_pair));
    }

    /* If read is at all close to the STR,
       perform realignment to determine read type 
       This will find some enclosing, FRRs, and spanning
    */
    if (debug) {
      std::cerr << "Processing single read" << std::endl;
    }
    int32_t data_value;
    ReadType read_type;
    ReadPair read_pair;
    ProcessSingleRead(alignment, chrom_ref_id, locus,
		      &data_value, &read_type);
    read_pair.read_type = read_type;
    read_pair.read1 = alignment;
    read_pair.data_value = data_value;
    read_pairs->insert(std::pair<std::string, ReadPair>(aln_key, read_pair));
  }
  /*  Second pass through reads where only one end processed */
  for (std::map<std::string, ReadPair>::iterator iter = read_pairs->begin();
       iter != read_pairs->end(); iter++) {
    if (iter->second.found_pair || iter->second.read_type != RC_UNKNOWN) {
      continue;
    }
    if (debug) {
      std::cerr << "Attempting to rescue mate " << iter->first << std::endl;
    }
    BamAlignment matepair;
    if (!RescueMate(bamreader, iter->second.read1, &matepair)) {
      continue;
    }
    if (debug) {
      std::cerr << "Found mate for " << iter->first << std::endl;
    }
    iter->second.read2 = matepair;
    int32_t data_value;
    ReadType read_type;
    ProcessSingleRead(matepair, chrom_ref_id, locus,
		      &data_value, &read_type);
    if (debug) {
      std::cerr << "Processed mate, found " << read_type << " " << data_value << std::endl;
    }
    iter->second.read_type = read_type;
    iter->second.data_value = data_value;
  }
  return true;
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
  // Only includes obvious cases, pre/post flank taken care of elsewhere
  bool span1 = alignment.RefID() == chrom_ref_id && alignment.GetEndPosition() <= locus.start &&
    alignment.MateRefID() == chrom_ref_id && alignment.MatePosition() >= locus.end;
  bool span2 = alignment.RefID() == chrom_ref_id && alignment.Position() >= locus.end &&
    alignment.MateRefID() == chrom_ref_id && alignment.MatePosition() <= locus.start-read_length;
  if (span1 || span2) {
    *insert_size = abs(alignment.TemplateLength());
    return true;
  }
  return false;
}

/*
  Process a single read using realignment
  
  Return false if something goes wrong
 */
bool ReadExtractor::ProcessSingleRead(BamAlignment alignment,
				      const int32_t& chrom_ref_id,
				      const Locus& locus,
				      int32_t* data_value,
				      ReadType* read_type) {
  /* If read in vicinity but not close to STR, save for later */
  if (alignment.RefID() == chrom_ref_id &&
      (alignment.Position() > locus.end || alignment.GetEndPosition() < locus.start)) {
    *read_type = RC_UNKNOWN;
    return true;
  }
  int32_t pos, pos_rev;
  int32_t score, score_rev;
  int32_t nCopy, nCopy_rev;
  std::string seq = alignment.QueryBases();
  std::string seq_rev = reverse_complement(seq);
  int32_t read_length = (int32_t)seq.size();
  /* Perform realignment and classification */
  if (!expansion_aware_realign(seq, locus.pre_flank, locus.post_flank, locus.motif,
			       &nCopy, &pos, &score)) {
    return false;
  }
  if (!expansion_aware_realign(seq_rev, locus.pre_flank, locus.post_flank, locus.motif,
			       &nCopy_rev, &pos_rev, &score_rev)) {
    return false;
  }
  if (score_rev > score) {
    nCopy = nCopy_rev;
    pos = pos_rev;
    score = score_rev;
    seq = seq_rev;
  }
  SingleReadType srt;
  if (!classify_realigned_read(seq, locus.motif, pos, nCopy, score,
			       (int32_t)locus.pre_flank.size(), &srt)) {
    return false;
  }
  if (debug) {
    std::cerr << "Processing single read found " << score << " " << pos << " " << srt << std::endl;
  }
  // Process according to guessed read type
  /* Spanning cases */
  // 5.2_filter_spanning_only_core.py#L86 - preflank case
  if ((alignment.MateRefID() == chrom_ref_id) &&
      (alignment.MatePosition() >= locus.end - read_length) &&
      (srt == SR_PREFLANK)) {
    *read_type = RC_SPAN;
    *data_value = abs(locus.start + nCopy*(int32_t)locus.motif.size() - read_length
		      - (alignment.MatePosition() + read_length));
    return true;
  }
  // 5.2_filter_spanning_only_core.py#L92 - postflank case
  if ((alignment.MateRefID() == chrom_ref_id) &&
      (alignment.MatePosition() <= locus.start) &&
      (srt == SR_POSTFLANK)) {
    *read_type = RC_SPAN;
    *data_value = abs(locus.end - nCopy*(int32_t)locus.motif.size()
		      - alignment.MatePosition() + read_length);
    return true;
  }
  /* FRR cases */
  // 5.2_filter_FRR_only_core.py:136
  if (srt == SR_IRR) {
    *read_type = RC_FRR;
    if (alignment.MatePosition() < locus.start) {
      *data_value = locus.start - (alignment.MatePosition()+read_length);
    } else {
      *data_value = alignment.MatePosition() - locus.end;
    }
    return true;
  }
  /* Enclosing cases */
  if (srt == SR_ENCLOSING) {
    *read_type = RC_ENCL;
    *data_value = nCopy;
    return true;
  }
  /* Partially spanning cases - get bound */
  if (srt == SR_PREFLANK || srt == SR_POSTFLANK) {
    *read_type = RC_BOUND;
    *data_value = nCopy;
    return true;
  }
  /* If we get here, we're not sure */
  *read_type = RC_UNKNOWN;
  return true;
}

/*
  Rescue mate if the mate didn't get mapped to the right region
  TODO: this doesn't handle if the mate is unaligned?
 */
bool ReadExtractor::RescueMate(BamCramMultiReader* bamreader,
			       BamAlignment alignment, BamAlignment* matepair) {
  const BamHeader* bam_header = bamreader->bam_header();
  std::string aln_key1 = trim_alignment_name(alignment);
  if (debug) {
    std::cerr << "Looking for mate in " << bam_header->ref_name(alignment.MateRefID()) <<
      " " << alignment.MatePosition() << std::endl;
  }
  bamreader->SetRegion(bam_header->ref_name(alignment.MateRefID()),
		       alignment.MatePosition()-1, alignment.MatePosition()+1);
  BamAlignment aln;
  while (bamreader->GetNextAlignment(aln)) {
    std::string aln_key2 = trim_alignment_name(aln);
    if (debug) {
      std::cerr << "Looking for " << aln_key1 << " found " << aln_key2 << std::endl;
    }
    if (aln_key1 == aln_key2) {
      *matepair = aln;
      return true;
    }
  }
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
