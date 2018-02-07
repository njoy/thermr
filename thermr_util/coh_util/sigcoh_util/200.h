#include <iostream>
#include <vector>


auto readBraggParameters( std::vector<double>& fl, int l, double econ ){
  // bragg parameters already read from endf6
  int nl = fl[5];
  double blast = 0, scon = 1, lenext = fl[9];
  //enext=sigfig(enext,7,-1);
  for ( int i = 1; i <= nl; ++i ){
    l = 1 + 2 * (i-1);
    fl[l-1] = fl[l+8-1] * econ;
    fl[l+1-1] = fl[l+9-1] - blast;
    blast = fl[l+9-1];
  } 
}
