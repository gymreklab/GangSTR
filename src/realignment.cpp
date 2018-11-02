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

#include "src/realignment.h"
#include <algorithm>
#include <sstream>
#include <iostream>

using namespace std;

bool find_longest_stretch(const std::string& seq,
			  const std::string& motif,
			  int32_t* nCopy_stretch,
			  int32_t* nCopy_total){
  int32_t read_len = (int32_t)seq.size();
  int32_t period = (int32_t)motif.size();
  bool motif_found;
  bool on_stretch = false;
  int32_t longest_stretch = 0, current_stretch = 0, total = 0;
  for (int i = 0; i < read_len - period; i++){
    motif_found = true;
    for (int j = 0; j < period; j++){
      if (motif[j] != seq[i + j]){
        on_stretch = false;
        motif_found = false;
        break;
      }
    }
    if (motif_found){
      total++;
      i+=period-1;
      if (on_stretch){
        current_stretch++;
      }
      else{
        on_stretch = true;
        current_stretch = 1;
      }
    }
    if (current_stretch > longest_stretch){
      longest_stretch = current_stretch;
    }
  }
  *nCopy_stretch = longest_stretch;
  *nCopy_total = total;
}


bool expansion_aware_realign(const std::string& seq,
			     const std::string& qual,
			     const std::string& pre_flank,
			     const std::string& post_flank,
			     const std::string& motif,
			     const int32_t& min_match,
			     int32_t* nCopy, 
			     int32_t* start_pos, 
			     int32_t* end_pos, 
			     int32_t* score,
			     FlankMatchState* fm_start,
			     FlankMatchState* fm_end) {

  *fm_start = FM_NOMATCH;
  *fm_end = FM_NOMATCH;
  int32_t read_len = (int32_t)seq.size();
  int32_t period = (int32_t)motif.size();
  int32_t min_nCopy = 0, total_nCopy = 0;
  // Find longest stretch of motif as starting point of our search.
  find_longest_stretch(seq, motif, &min_nCopy, &total_nCopy);
  if (min_nCopy < 2 and total_nCopy < 10){
      *nCopy = 0;
      *score = 0;
      *start_pos = 0;
      *end_pos = 0;
      return true;
  }
  int32_t max_score = 0;
  int32_t second_best_score = 0;
  int32_t max_nCopy = 0;
  int32_t second_best_nCopy = 0;
  int32_t max_start_pos = 0;
  int32_t max_end_pos = 0;
  int32_t current_score = 0;
  int32_t current_start_pos = 0, current_end_pos = 0;
  int32_t current_nCopy, current_num_mismatch;
  int32_t prev_score = 0;
  int32_t prev_prev_score = -1;
  std::string template_sub, sequence_sub;
  MARGIN = 1 * period;
  
  //cerr << min_nCopy << " ";
  for (current_nCopy=min_nCopy; current_nCopy<(int32_t)(read_len/period)+2; current_nCopy++) {
    std::stringstream var_realign_ss;
    var_realign_ss << pre_flank;
    for (int i = 0; i<current_nCopy; i++) {
      var_realign_ss << motif;
    }
    var_realign_ss << post_flank;
    std::string var_realign_string = var_realign_ss.str();
    

    if (!striped_smith_waterman(var_realign_string, seq, qual, &current_start_pos, &current_end_pos, &current_score, &current_num_mismatch)) {
      return false;
    }
    

    // Flank match check
    // Preflank
    if (read_len - current_start_pos - min_match >= 0 &&
	read_len - current_start_pos + min_match <= read_len){ //Full match is possible
      sequence_sub = seq.substr(read_len - current_start_pos - min_match, 2 * min_match);
      template_sub = var_realign_string.substr(read_len - min_match, 2 * min_match);
      if (sequence_sub == template_sub){
	*fm_start = FM_COMPLETE;
      }
      else{
	*fm_start = FM_NOMATCH;
      }
    }
    else{
      *fm_start = FM_NOMATCH;
    }
    // Postflank
    if (read_len - current_start_pos + current_nCopy * period + min_match <= read_len &&
	read_len - current_start_pos + current_nCopy * period - min_match >= 0){ //Full match is possible
      sequence_sub = seq.substr(read_len - current_start_pos + current_nCopy * period - min_match, 2 * min_match);
      template_sub = var_realign_string.substr(read_len + current_nCopy * period - min_match, 2 * min_match);
      
      if (sequence_sub == template_sub){
	*fm_end = FM_COMPLETE;
      }
      else{
	*fm_end = FM_NOMATCH;
      }
    }
    else{
      *fm_end = FM_NOMATCH;
    }
    if (current_score >= max_score) {
      second_best_score = max_score;
      second_best_nCopy = max_nCopy;
      max_score = current_score;
      max_nCopy = current_nCopy;
      max_start_pos = current_start_pos;
      max_end_pos = current_end_pos;
    }
  
    if (*fm_start == FM_COMPLETE && *fm_end == FM_COMPLETE){
      break;
    }
    // Stop if score is relatively high, but lower than max
    if (current_score > 0.7 * SSW_MATCH_SCORE * read_len and 
	current_score <= max_score and
	prev_score == current_score and
	prev_prev_score == prev_score and
	total_nCopy > max_nCopy){
      //max_nCopy--;
      break;
    }
    if (current_score == read_len*SSW_MATCH_SCORE) {
      break;
    }
    prev_prev_score = prev_score;
    prev_score = current_score;
  //cerr << current_nCopy << endl;
  }
  if (max_nCopy < 0.85 * read_len / period and 
      *fm_start == FM_NOMATCH and *fm_end == FM_NOMATCH){
    max_nCopy = 0;
  }
  if (max_nCopy > read_len / period){ // did we overshoot?
    int32_t overshoot = max_nCopy * period - read_len;
    max_end_pos-=overshoot;
    max_nCopy = read_len / period;
  }
  *nCopy = max_nCopy;
  *score = max_score;
  *start_pos = max_start_pos;
  *end_pos = max_end_pos;
    
  return true;
}

bool trace_back(const std::vector<std::vector<int32_t> > score_matrix, 
            const int32_t& start_pos, 
            const int32_t& start_pos_temp,
            const std::string& seq1, 
            const std::string& seq2,
            std::string* seq1_rea){
  sw_move best_move;
  int32_t x = start_pos;
  int32_t y = start_pos_temp;
  if(!next_move(score_matrix, x, y, &best_move)){
    return false;
  }
  while (best_move != SW_END){
    if (best_move == SW_DIAG){

    }
    else if (best_move == SW_UP){

    }
    else{

    }
    if(!next_move(score_matrix, x, y, &best_move)){
      return false;
    }
  }
}

bool next_move(std::vector<std::vector<int32_t> > score_matrix, 
            const int32_t& x, 
            const int32_t& y, 
            sw_move* move){
  int32_t diag = score_matrix.at(x - 1).at(y - 1);
  int32_t up = score_matrix.at(x - 1).at(y);
  int32_t left = score_matrix.at(x).at(y - 1);
  if (diag >= up and diag >= left){
    *move = (diag != 0 ? SW_DIAG : SW_END);
    return true;
  }
  else if (up > diag and up >= left){
    *move = (up != 0 ? SW_UP : SW_END);
    return true;
  }
  else if (left > diag and left > up){
    *move = (left != 0 ? SW_LEFT : SW_END);
    return true;
  }
  else{
    return false;
  }
}

bool striped_smith_waterman(const std::string& ref,
        const std::string& seq,
        const std::string& qual,
        int32_t* pos, int32_t* end, int32_t* score, int32_t* mismatches) {

  // SSW Objects
  StripedSmithWaterman::Aligner* aligner;
  StripedSmithWaterman::Filter* filter;
  StripedSmithWaterman::Alignment* alignment;

  aligner = new StripedSmithWaterman::Aligner(SSW_MATCH_SCORE, 
                              SSW_MISMATCH_SCORE, 
                              SSW_GAP_OPEN, 
                              SSW_GAP_EXTEND);
  filter = new StripedSmithWaterman::Filter;
  alignment = new StripedSmithWaterman::Alignment;
  int32_t maskLen = seq.size() / 2;
  maskLen = maskLen < 15 ? 15 : maskLen;
  maskLen = 15;
  aligner->Align(seq.c_str(), ref.c_str(), (int32_t)ref.size(), *filter, alignment, maskLen);

  *pos = alignment->ref_begin;
  *end = alignment->ref_end;
  *score = alignment->sw_score;
  *mismatches = alignment->mismatches;

  delete aligner;
  delete filter;
  delete alignment;
  return true;
}

static void ssw_PrintAlignment(const StripedSmithWaterman::Alignment& alignment){
  cerr << "===== SSW result =====" << endl;
  cerr << "Best Smith-Waterman score:\t" << alignment.sw_score << endl
       << "Next-best Smith-Waterman score:\t" << alignment.sw_score_next_best << endl
       << "Reference start:\t" << alignment.ref_begin << endl
       << "Reference end:\t" << alignment.ref_end << endl
       << "Query start:\t" << alignment.query_begin << endl
       << "Query end:\t" << alignment.query_end << endl
       << "Next-best reference end:\t" << alignment.ref_end_next_best << endl
       << "Number of mismatches:\t" << alignment.mismatches << endl
       << "Cigar: " << alignment.cigar_string << endl;
  cerr << "======================" << endl;
}

bool smith_waterman(const std::string& seq1,
        const std::string& seq2,
        const std::string& qual,
        int32_t* pos, int32_t* pos_temp, int32_t* score) {
  // The scoring matrix contains an extra row and column for the gap (-), hence
  // the +1 here
  int32_t rows = (int32_t)seq1.size() + 1;
  int32_t cols = (int32_t)seq2.size() + 1;

  // Initialize the scoring matrix
  int32_t current_score;
  int32_t start_pos, start_pos_temp;
  std::vector<std::vector<int32_t> > score_matrix;
  score_matrix.resize(rows, std::vector<int32_t>(cols, 0));
  if (!create_score_matrix(rows, cols, seq1, seq2, qual,
         &score_matrix, &start_pos, &start_pos_temp, &current_score)) {
    return false;
  }
  *pos = start_pos-(int32_t)seq2.size();
  *pos_temp = start_pos_temp - (int32_t)seq1.size();
  *score = current_score;
  return true;
}

/*
  Create a matrix of scores representing trial alignments of the two sequences.
  Sequence alignment can be treated as a graph search problem. This function
  creates a graph (2D matrix) of scores, which are based on trial alignments
  of different base pairs. The path with the highest cummulative score is the
  best alignment.
 */
bool create_score_matrix(const int32_t& rows, const int32_t& cols,
       const std::string& seq1,
       const std::string& seq2,
       const std::string& qual,
       std::vector<std::vector<int32_t> >* score_matrix,
       int32_t* start_pos, int32_t* start_pos_temp, int32_t* current_score) {
  int32_t max_score = 0;
  int32_t max_pos_row = -1; // The row and column of highest score in the matrix
  int32_t max_pos_col = -1;

  for (int32_t i=1; i<rows; i++) {
    for (int32_t j=1; j<cols; j++) {
      if (!calc_score(i, j, seq1, seq2, qual, score_matrix)) {
  return false;
      }
      if (score_matrix->at(i).at(j) > max_score) {
  max_score = score_matrix->at(i).at(j);
  max_pos_row = i;
  max_pos_col = j;
      }
    }
  }
  if (max_pos_row == -1 || max_pos_col == -1) {
    *current_score = 0;
    *start_pos = -1;
    return true;
  }
  *current_score = max_score;
  *start_pos = max_pos_row;
  // *start_pos_temp = max_pos_col;
  return true;
}

/*
  Calculate score for a given x, y position in the scoring matrix.
  The score is based on the up, left, and upper-left neighbors.
 */
bool calc_score(const int32_t& i, const int32_t& j,
    const std::string& seq1, const std::string& seq2,
    const std::string& qual,
    std::vector<std::vector<int32_t> >* score_matrix) {
  int32_t max_score = 0;
  int32_t baseq = int32_t(qual.at(j-1));
  int32_t similarity = (seq1.at(i-1)==seq2.at(j-1)) ? 
    MATCH_SCORE : MISMATCH_SCORE;
  // TODO pass threshold instead of hard code
  // int32_t similarity = (seq1.at(i-1)==seq2.at(j-1)) ? 
  //   MATCH_SCORE : (baseq>45 ? MISMATCH_SCORE : MISMATCH_SCORE / 4);
  int32_t diag_score = score_matrix->at(i-1).at(j-1) + similarity;
  if (diag_score > max_score) {
    max_score = diag_score;
  }
  int32_t up_score = score_matrix->at(i-1).at(j) + GAP_SCORE;
  if (up_score > max_score) {
    max_score = up_score;
  }
  int32_t left_score = score_matrix->at(i).at(j-1) + GAP_SCORE;
  if (left_score > max_score) {
    max_score = left_score;
  }
  score_matrix->at(i).at(j) = max_score;
  return true;
}

bool classify_realigned_read(const std::string& seq,
			     const std::string& motif,
			     const int32_t& start_pos,
			     const int32_t& end_pos,
			     const int32_t& nCopy,
			     const int32_t& score,
			     const int32_t& prefix_length,
			     const int32_t& min_match,
			     const bool& isMapped,
			     const std::string& pre_flank,
			     const std::string& post_flank,
			     const FlankMatchState& fm_start,
			     const FlankMatchState& fm_end,
			     SingleReadType* single_read_class) {

  int32_t i,j, limit;
  bool flank_match, failed_flank_test = false;
  *single_read_class = SR_UNKNOWN;
  // Get coords of the STR
  int32_t start_str = prefix_length;
  int32_t end_str = prefix_length + nCopy*(int32_t)motif.size();

  // Check if read starts in the STR
  bool start_in_str = false;
  bool end_in_str = false;
  if ((start_pos >= start_str-MARGIN) && (start_pos <= end_str+MARGIN)) {
    start_in_str = true;
  }
  if ((end_pos >= start_str-MARGIN) && (end_pos <= end_str+MARGIN)) {
    end_in_str = true;
  }

    if (seq == "gctgctgctgctgctgctgctgctgctgctgctgctgctgctgctgctgctgctgctgctgctgctgctgctgctgctgctgctgctgctgctgctgctgctgctgctgctgctgctgctgctgctgctgctgcagctgctgctgctgcg")
      cerr <<"\n>\n"<< end_pos << " " << end_str <<"\n>\n"<< endl;

  // Check if perfect flanks exist:
  if (fm_start == FM_COMPLETE && fm_end == FM_COMPLETE){
    *single_read_class = SR_ENCLOSING;
    return true;
  }
  else if (fm_start == FM_COMPLETE && end_in_str){
    *single_read_class = SR_PREFLANK;
    return true;
  }
  else if (fm_end == FM_COMPLETE && start_in_str){
    *single_read_class = SR_POSTFLANK;
    return true;
  }

  // TODO:
  /*
    Only check for FRR now, no need for checking for other types of reads.
    Remove extra checks (specially flank checks)
  */

  // Set threshold for match
  int32_t score_threshold = (int32_t)(MATCH_PERC_THRESHOLD*seq.size()*SSW_MATCH_SCORE);

  if (isMapped && (score < score_threshold || nCopy == 0)) {
    *single_read_class = SR_UNKNOWN;
    return true;
  } else if (start_in_str && end_in_str) {
    *single_read_class = SR_UM_POT_IRR;
    return true;
  } 
  /*
  else if (start_in_str && !end_in_str) {
    // Post flank check
    flank_match = true;
    j = 0;
    limit = (end_pos - end_str >= min_match ? end_str - start_pos + min_match - 1 : end_pos - start_pos - 1);
    
    for (i = max(min(end_str - start_pos, (int32_t)seq.size() - 1),0) ; 
       i <= min(limit, (int32_t)seq.size() - 1);
       i++){

     if (seq.at(i)!=post_flank.at(j)){
       flank_match = false;
     }
     j++;   
    }
    if (flank_match){
      *single_read_class = SR_POSTFLANK;
      return true;
    }
    else{
      failed_flank_test = true;
      *single_read_class = SR_UNKNOWN;
    }
  } else if (!start_in_str && end_in_str) {
    // Pre flank check
    flank_match = true;
    j = start_str - start_pos >= min_match ? start_str - min_match : start_pos;
    for (i = min(max(start_str - start_pos - min_match, 0), (int32_t)seq.size())
            ; i <min(start_str - start_pos, (int32_t)seq.size()) ; i++){
      if (seq.at(i)!=pre_flank.at(j)){
        flank_match = false;
      }
      j++;
    }

    if(flank_match){
      *single_read_class = SR_PREFLANK;
      return true;
    }
    else{
      failed_flank_test = true;
      *single_read_class = SR_UNKNOWN;
    }
  } else if (start_pos < start_str && end_pos > end_str && score > score_threshold &&
    (start_str - start_pos <= seq.size() - (end_str - start_str))) {
    *single_read_class = SR_ENCLOSING;

    // Pre flank check
    flank_match = true;
    j = start_str - start_pos >= min_match ? start_str - min_match : start_pos;
    for (i = min(max(start_str - start_pos - min_match, 0), (int32_t)seq.size())
            ; i <min(start_str - start_pos, (int32_t)seq.size()) ; i++){
      // cerr<<seq.at(i);
      if (seq.at(i)!=pre_flank.at(j)){
        flank_match = false;
        // break;
      }
      j++;
    }
    // if (flank_match){
    //   cerr<<" -> PASS!!";
    // }
    // cerr<<endl;
    // j = start_str - start_pos >= min_match ? start_str - min_match : start_pos;
    // for (i = min(max(start_str - start_pos - min_match, 0), (int32_t)seq.size())
    //         ; i <min(start_str - start_pos, (int32_t)seq.size()) ; i++){
    //   cerr<<pre_flank.at(j);
    //   j++;
    // }
    // cerr<<endl;


    // Post flank check
    // flank_match = true;
    if (flank_match){
       j = 0;
       limit = (end_pos - end_str >= min_match ? end_str - start_pos + min_match - 1 : end_pos - start_pos - 1);

       // cerr << (end_pos - end_str >= min_match) << endl;
       // cerr << end_str - start_pos + min_match - 1 << endl;
       // cerr << "Size: "<< seq.size() - 1 << endl;
       for (i = max(min(end_str - start_pos, (int32_t)seq.size() - 1),0); 
          i <= min(limit, (int32_t)seq.size() - 1);
          i++){
	 //cerr<<seq.at(i);
        if (seq.at(i)!=post_flank.at(j)){
          flank_match = false;
          // break;
        }
        j++;   
       }
       
       //if (flank_match){
       //  cerr << " -> PASS!!";
       //}
       //cerr<<endl;
       //j = 0;
       //for (i = min(end_str - start_pos, (int32_t)seq.size() - 1) ; 
       //    i <= min(limit, (int32_t)seq.size() - 1);
       //    i++){
       //  cerr<<post_flank.at(j);
       //  j++;
       // }
    }
    // If either flanks didn't match reference
    if (flank_match){
      *single_read_class = SR_ENCLOSING;
      return true;
    }
    else {
      failed_flank_test = true;
      *single_read_class = SR_UNKNOWN;
    }
  }

  double FRR_slip = 0.1 * prefix_length;
  if (!isMapped or failed_flank_test){ // If isMapped is false, or failed flank test, check if FRR
    if (start_pos < start_str + FRR_slip && start_pos > start_str - FRR_slip &&
      end_pos < end_str + FRR_slip && end_pos > end_str - FRR_slip &&
      nCopy > 0.7 * seq.size() / motif.size() && score > 0.7 * score_threshold){
      // cerr<<"nCopy: " << nCopy<<"\tscore: "<<score<<"/"<<seq.size()*SSW_MATCH_SCORE<<endl;
      *single_read_class = SR_UM_POT_IRR;
      return true;
    }
  }
  */
  // If no other class matches this read pair:
  *single_read_class = SR_UNKNOWN;
  return true; 

}
