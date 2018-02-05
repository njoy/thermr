#include <iostream>
#include <vector>


auto readBraggParameters( std::vector<double> fl, int l, double econ ){
  int k;
  

 // bragg parameters already read from endf6
   k=fl[6];
   int nl=k;
   double blast=0;
   double scon=1;
   double lenext=fl[9];
   //enext=sigfig(enext,7,-1);
   for ( int i = 1; i < nl; ++i ){
      l=1+2*(i-1);
      fl[l]=fl[l+8]*econ;
      fl[l+1]=fl[l+9]-blast;
      blast=fl[l+9];
    } 
}
