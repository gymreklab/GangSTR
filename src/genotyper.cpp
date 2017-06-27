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

#include <iostream>

#include "src/genotyper.h"

using namespace std;

Genotyper::Genotyper(BamCramMultiReader& _bamreader,
		     RefGenome _refgenome,
		     const Options& _options) {
  refgenome = &_refgenome;
  bamreader = &_bamreader;
  options = &_options;
}

bool Genotyper::ProcessLocus(Locus* locus) {
  // Load all read data
  likelihood_maximizer.Reset();
  if (!read_extractor.ExtractReads(*bamreader, &likelihood_maximizer)) {
    return false;
  }
  // Maximize the likelihood
  int32_t allele1, allele2;
  if (!likelihood_maximizer.OptimizeLikelihood(&allele1, &allele2)) {
    return false;
  }
  // Write VCF Output - TODO
  return true;
}

void Genotyper::Debug() {
  cerr << "testing refgenome" << endl;
  std::string seq;
  refgenome->GetSequence("3", 63898361, 63898392, &seq);
  cerr << seq << endl;
  cerr << "testing bam" << endl;
  bamreader->SetRegion("1", 0, 10000);
  BamAlignment aln;
  if (bamreader->GetNextAlignment(aln)) { // Requires SetRegion was called
    std::string testread = aln.QueryBases();
    cerr << testread << endl;
  } else {
    cerr << "testing bam failed" << endl;
  }
}

Genotyper::~Genotyper() {}
