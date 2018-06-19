
#ifndef UTIL_REPOZ_HH
#define UTIL_REPOZ_HH
#include <fstream>

inline void repoz( std::fstream& ntape ){
  ntape.clear();
  ntape.seekg(0);
}

#endif
