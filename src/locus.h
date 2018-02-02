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

#include <string>

class Locus {
 public:
  Locus();
  virtual ~Locus();

  std::string chrom;
  int start;
  int end;
  int period;
  std::string motif;
  std::string pre_flank;
  std::string post_flank;

  // Fill in these fields
  double insert_size_mean;
  double insert_size_stddev;
  int allele1;
  int allele2;
  int lob1;
  int hib1;
  int lob2;
  int hib2;
  double min_neg_lik;
  size_t enclosing_reads;
  size_t spanning_reads;
  size_t frr_reads;
  size_t flanking_reads;
  size_t depth;

  // Off target loci
  std::string offchrom;
  int offstart;
  int offend;
};

#endif  // SRC_LOCUS_H__
