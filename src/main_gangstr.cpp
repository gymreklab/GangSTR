#include <getopt.h>
#include <stdlib.h>

#include <iostream>
#include <sstream>

#include "src/common.h"

using namespace std;

void show_help() {
  std::stringstream help_msg;
  help_msg << "\nGangSTR [OPTIONS] "
	   << "--bam <file1[,file2,...]> "
	   << "--ref <reference.fa> "
	   << "--loci <strref.tab> "
	   << "--out <outprefix> "
	   << "\n\nOptions:\n"
	   << "-h,--help      display this help screen\n"
	   << "-v,--verbose   print out useful progress messages\n"
	   << "This program takes in aligned reads in BAM format\n"
	   << "and genotypes a reference set of STRs\n\n";
  cerr << help_msg.str();
  exit(1);
}

void parse_commandline_options(int argc, char* argv[]) {
  enum LONG_OPTIONS {
    OPT_BAMFILES,
    OPT_REFFA,
    OPT_LOCI,
    OPT_OUT,
    OPT_HELP,
    OPT_VERBOSE,
    OPT_VERSION
  };
  static struct option long_options[] = {
    {"--bam", 1, 0, OPT_BAMFILES},
    {"--ref", 1, 0, OPT_REFFA},
    {"--loci", 1, 0, OPT_LOCI},
    {"--out", 1, 0, OPT_OUT},
    {"--help", 0, 0, OPT_HELP},
    {"--verbose", 0, 0, OPT_VERBOSE},
    {"--version", 0, 0, OPT_VERSION},
    {NULL, no_argument, NULL, 0},
  };
  int ch;
  int option_index = 0;
  ch = getopt_long(argc, argv, "hv?",
                   long_options, &option_index);
  while (ch != -1) {
    switch (ch) {
    case OPT_BAMFILES:
      //TODO
      break;
    case OPT_REFFA:
      //TODO
      break;
    case OPT_LOCI:
      //TODO
      break;
    case OPT_OUT:
      //TODO
      break;
    case OPT_HELP:
    case 'h':
      show_help();
    case OPT_VERBOSE:
    case 'v':
      //TODO
      break;
    case OPT_VERSION:
      cerr << _GIT_VERSION << endl;
      exit(0);
    case '?':
      show_help();
    default:
      show_help();
    };
  };
  // Leftover arguments are errors
  if (optind < argc) {
    PrintMessageDieOnError("Unnecessary leftover arguments", ERROR);
  }
  // Perform other error checking TODO
}

int main(int argc, char* argv[]) {
  parse_commandline_options(argc, argv);
}
