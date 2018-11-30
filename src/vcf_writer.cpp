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

#include <set>
using namespace std;

VCFWriter::VCFWriter(const std::string& _vcffile,
		     const std::string& full_command,
		     SampleInfo& _sample_info) {
  sample_info = &_sample_info;
  writer_.open(_vcffile.c_str());
  // Write header
  writer_ << "##fileformat=VCFv4.1" << std::endl;
  writer_ << "##command=" << full_command << std::endl;
  writer_ << "##INFO=<ID=END,Number=1,Type=Integer,Description=\"End position of variant\">" << endl;
  writer_ << "##INFO=<ID=RU,Number=1,Type=String,Description=\"Repeat motif\">" << endl;
  writer_ << "##INFO=<ID=REF,Number=1,Type=Float,Description=\"Reference copy number\">" << endl;
  writer_ << "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">" << endl;
  writer_ << "##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Read Depth\">" << endl;
  writer_ << "##FORMAT=<ID=REPCN,Number=1,Type=String,Description=\"Genotype given in number of copies of the repeat motif\">" << endl;
  writer_ << "##FORMAT=<ID=REPCI,Number=1,Type=String,Description=\"Confidence intervals\">" << endl;
  writer_ << "##FORMAT=<ID=RC,Number=1,Type=String,Description=\"Number of reads in each class (enclosing, spanning, FRR, bounding)\">" << endl;
  writer_ << "##FORMAT=<ID=Q,Number=1,Type=Float,Description=\"Min. negative likelihood\">" << endl;
  writer_ << "##FORMAT=<ID=INS,Number=1,Type=String,Description=\"Insert size mean and stddev\">" << endl;
  writer_ << "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT";
  const std::set<std::string> rg_samples = sample_info->GetSamples();
  sample_names.clear();
  for (std::set<std::string>::iterator it=rg_samples.begin();
       it != rg_samples.end(); it++) {
    writer_ << "\t" << *it;
    sample_names.push_back(*it);
  }
  writer_ << endl;
  writer_.flush();
}

void VCFWriter::WriteRecord(Locus& locus) {
  std::string locus_motif = locus.motif;
  // Set REF
  int refsize = (locus.end-locus.start+1)/locus.period;
  stringstream ref_allele;
  for (int i=0; i<refsize; i++) {
    ref_allele << locus.motif;
  }
  // Set ALT alleles
  stringstream alt_alleles;
  std::set<int> alt_allele_lengths;
  for (size_t i=0; i<sample_names.size(); i++) {
    std::string samp = sample_names[i];
    if (locus.allele1[samp] != refsize) {
      alt_allele_lengths.insert(locus.allele1[samp]);
    }
    if (locus.allele2[samp] != refsize) {
      alt_allele_lengths.insert(locus.allele2[samp]);
    }
  } 
  std::vector<int> alt_allele_lengths_v(alt_allele_lengths.begin(), alt_allele_lengths.end());

  std::map<int, int>alt_length_to_ind;
  if (alt_allele_lengths.size() == 0) {
    alt_alleles << ".";
  } else {
    for (int i=0; i<alt_allele_lengths_v[0]; i++) {
      alt_alleles << locus.motif;
    }
    alt_length_to_ind[alt_allele_lengths_v[0]] = 1;
    for (size_t i=1; i<alt_allele_lengths_v.size(); i++) {
      alt_length_to_ind[alt_allele_lengths_v[i]] = i+1;
      for (int j=0; j<alt_allele_lengths_v[i]; j++) {
	alt_alleles << locus.motif;
      }
    }
  }

  // Write locus info
  writer_ << locus.chrom << "\t"
	  << locus.start << "\t"
	  << ".\t"
	  << ref_allele.str() << "\t"
	  << alt_alleles.str() << "\t"
	  << "." << "\t"
	  << "." << "\t"
	  << "END=" << locus.end << ";"
	  << "RU=" << locus.motif << ";"
	  << "REF=" << refsize << "\t"
	  << "GT:DP:REPCN:REPCI:RC:Q:INS";
  
  // Write info for each sample
  stringstream gt_str;
  int period = static_cast<int>(locus_motif.size());
  for (size_t i=0; i<sample_names.size(); i++) {
    std::string samp = sample_names[i];
    if (!locus.called[samp]) {
      writer_ << "\t.";
      continue;
    }
    gt_str.str("");
    gt_str << alt_length_to_ind[locus.allele1[samp]] << "/" << alt_length_to_ind[locus.allele2[samp]];
    writer_ << "\t"
	    << gt_str.str() << ":"
	    << locus.depth[samp] << ":"
      //	    << (locus.allele1[samp]-refsize)*period << "," <<  (locus.allele2[samp]-refsize)*period << ":"
      //	    << (locus.lob1[samp]-refsize)*period << "-" 
      //	    << (locus.hib1[samp]-refsize)*period << "," 
      //	    << (locus.lob2[samp]-refsize)*period << "-" 
      //	    << (locus.hib2[samp]-refsize)*period << ":"
	    << locus.allele1[samp] << "," <<  locus.allele2[samp] << ":"
      	    << locus.lob1[samp] << "-" 
	    << locus.hib1[samp] << "," 
	    << locus.lob2[samp] << "-" 
	    << locus.hib2[samp] << ":"
      	    << locus.enclosing_reads[samp] << "," << locus.spanning_reads[samp] << "," << locus.frr_reads[samp] << "," << locus.flanking_reads[samp] << ":"
	    << locus.min_neg_lik[samp] << ":"
	    << sample_info->GetInsertMean(samp) << "<" << sample_info->GetInsertSdev(samp);
  }
  writer_ << endl;
  writer_.flush();
}

VCFWriter::~VCFWriter() {
  writer_.close();
}
