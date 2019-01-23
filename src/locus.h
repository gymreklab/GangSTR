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

#ifndef SRC_LOCUS_H__
#define SRC_LOCUS_H__

#include <map>
#include <utility>
#include <string>
#include <vector>
#include <stdint.h>

struct GenomeRegion{
  std::string chrom;
  int start;
  int end;
};

class Locus {
 public:
  Locus();
  virtual ~Locus();
  void Reset();
  std::string chrom;
  int start;
  int end;

  int period;
  std::string motif;
  std::string pre_flank;
  std::string post_flank;

  // Off target loci
  bool offtarget_set;
  std::vector<GenomeRegion> offtarget_regions;
  double offtarget_share;

  // Fill in these fields - separate for each sample
  std::map<std::string,int> allele1;
  std::map<std::string,int> allele2;
  std::map<std::string,int> lob1;
  std::map<std::string,int> hib1;
  std::map<std::string,int> lob2;
  std::map<std::string,int> hib2;
  std::map<std::string,double> a1_se;
  std::map<std::string,double> a2_se;
  std::map<std::string,double> min_neg_lik;
  std::map<std::string,size_t> enclosing_reads;
  std::map<std::string,size_t> spanning_reads;
  std::map<std::string,size_t> frr_reads;
  std::map<std::string,size_t> flanking_reads;
  std::map<std::string,size_t> depth;
  std::map<std::string,bool> called;
  std::map<std::string,std::vector<double> > expansion_probs;
  std::map<std::string, std::map<std::pair<int32_t, int32_t>, double> > grid_likelihoods; // log10 lik
  std::map<std::string, double> q_scores;
  // Likelihood grid
  int32_t grid_min_allele;
  int32_t grid_max_allele;

  // Threshold for calling expansions
  int32_t expansion_threshold;

  // Stutter params
  double stutter_up;
  double stutter_down;
  double stutter_p;
};

#endif  // SRC_LOCUS_H__
