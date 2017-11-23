// Modified from https://github.com/tfwillems/HipSTR/blob/3a7c6f043d822a74df3ad42be00afc35b91f12fa/src/stringops.h
#ifndef STRING_OPS_H_
#define STRING_OPS_H_

#include <string>
#include <vector>

void split_by_delim(const std::string &s, char delim, 
		    std::vector<std::string>& substrings);

std::string uppercase(std::string str);
std::string lowercase(std::string str);

bool string_starts_with(std::string& s, std::string prefix);

bool string_ends_with(std::string& s, std::string suffix); 

bool orderByLengthAndSequence(const std::string& s1, const std::string s2);

int length_suffix_match(std::string& s1, std::string& s2);

std::string reverse_complement(std::string nucs);
char complement(const char nucleotide);

#endif
