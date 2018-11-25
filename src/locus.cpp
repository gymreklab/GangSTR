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

#include "src/locus.h"

using namespace std;

Locus::Locus() {
  Reset();
}

void Locus::Reset(){
  chrom = "";
  start = -1;
  end = -1;
  period = -1;

  allele1.clear();
  allele2.clear();
  lob1.clear();
  lob2.clear();
  hib1.clear();
  hib2.clear();
  min_neg_lik.clear();
  enclosing_reads.clear();
  spanning_reads.clear();
  frr_reads.clear();
  flanking_reads.clear();
  depth.clear();
  called.clear();

  offtarget_set = false;
  offtarget_regions.clear();
  offtarget_share = 0.0;
}
Locus::~Locus() {}
