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

#include <string>
#include <sstream>

#include "src/vcf_writer.h"

using namespace std;

VCFWriter::VCFWriter(const std::string& _vcffile,
		     const std::string& full_command) {
  writer_.open(_vcffile.c_str());
  // Write header
  writer_ << "##fileformat=VCFv4.1" << std::endl;
  writer_ << "##command=" << full_command << std::endl;
  writer_ << "##INFO=<ID=END,Number=1,Type=Integer,Description=\"End position of variant\">" << endl;
  writer_ << "##INFO=<ID=RU,Number=1,Type=String,Description=\"Repeat motif\">" << endl;
  writer_ << "##INFO=<ID=REF,Number=1,Type=Float,Description=\"Reference copy number\">" << endl;
  writer_ << "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">" << endl;
  writer_ << "##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Read Depth\">" << endl;
  writer_ << "##FORMAT=<ID=GB,Number=1,Type=String,Description=\"Genotype given in bp difference from reference\">" << endl;
  writer_ << "##FORMAT=<ID=CI,Number=1,Type=String,Description=\"Confidence intervals\">" << endl;
  writer_ << "##FORMAT=<ID=RC,Number=1,Type=String,Description=\"Number of reads in each class (enclosing, spanning, FRR, bounding\">" << endl;
  writer_ << "##FORMAT=<ID=Q,Number=1,Type=Float,Description=\"Min. negative likelihood\">" << endl;
  writer_ << "##FORMAT=<ID=INS,Number=1,Type=String,Description=\"Insert size mean and stddev\">" << endl;
  writer_ << "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample" << endl;
  writer_.flush();
}

void VCFWriter::WriteRecord(const Locus& locus) {
  int ref_size = (locus.end-locus.start+1)/locus.period;
  stringstream ref_allele;
  stringstream alt_alleles;
  stringstream gt_str;
  for (int i=0; i<ref_size; i++) {
    ref_allele << locus.motif;
  }
  if (locus.allele1 != ref_size) {
    for (int i=0; i<locus.allele1; i++) {
      alt_alleles << locus.motif;
    }
  }
  if (locus.allele2 != ref_size) {
    if (locus.allele1 != ref_size) {
      alt_alleles << ",";
    }
    for (int i=0; i<locus.allele2; i++) {
      alt_alleles << locus.motif;
    }
  }
  if (locus.allele1 == ref_size && locus.allele2 == ref_size) {
    alt_alleles << ".";
    gt_str << "0/0";
  } else if (locus.allele1 == ref_size) {
    gt_str << "0/1";
  } else if (locus.allele2 == ref_size) {
    gt_str << "1/0";
  } else {
    gt_str << "1/2";
  }
  writer_ << locus.chrom << "\t"
	  << locus.start << "\t"
	  << ".\t"
	  << ref_allele.str() << "\t"
	  << alt_alleles.str() << "\t"
	  << "." << "\t"
	  << "." << "\t"
	  << "END=" << locus.end << ";"
	  << "RU=" << locus.motif << ";"
	  << "REF=" << ref_size << "\t"
	  << "GT:DP:GB:CI:RC:Q:INS" << "\t"
	  << gt_str.str() << ":"
	  << locus.depth << ":"
	  << locus.allele1 << "," << locus.allele2 << ":"
	  << locus.lob1 << "-" << locus.hib1 << "," << locus.lob2 << "-" << locus.hib2 << ":"
	  << locus.enclosing_reads << "," << locus.spanning_reads << "," << locus.frr_reads << "," << locus.flanking_reads << ":"
	  << locus.min_neg_lik << ":"
	  << locus.insert_size_mean << "," << locus.insert_size_stddev
	  << endl;
  writer_.flush();
}

VCFWriter::~VCFWriter() {
  writer_.close();
}
