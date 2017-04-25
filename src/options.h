#ifndef SRC_OPTIONS_H__
#define SRC_OPTIONS_H__

#include <vector>
#include <string>

class Options {
 public:
  Options();
  virtual ~Options();

  // User defined options
  std::vector<std::string> bamfiles;
  std::string reffa;
  std::string locifile;
  std::string outprefix;
  bool verbose;
};

#endif  // SRC_OPTIONS_H__
