#include <iostream>
#include "general_util/sigfig.h"



//template <typename Range, typename Float>
template <typename Float>
auto e_ep_mu( Float T, Float& teff, Float& teff2 ){
  // should only be here if iform = 0
  std::vector<double> egrid { 1.e-5, 1.78e-5, 2.5e-5, 3.5e-5, 5.0e-5, 7.0e-5, 
    1.e-4, 1.26e-4, 1.6e-4, 2.0e-4, 2.53e-4, 2.97e-4, 3.5e-4, 4.2e-4, 5.06e-4, 
    6.15e-4, 7.5e-4, 8.7e-4, 1.012e-3, 1.23e-3, 1.5e-3, 1.8e-3, 2.03e-3, 
    2.277e-3, 2.6e-3, 3e-3, 3.5e-3, 4.048e-3, 4.5e-3, 5e-3, 5.6e-3, 6.325e-3, 
    7.2e-3, 8.1e-3, 9.108e-3, 0.01, 0.01063, 0.0115, 0.012397, 0.0133, 0.01417, 
    0.015, 0.016192, 0.0182, 0.0199, 0.020493, 0.0215, 0.0228, 0.0253, 0.028, 
    0.030613, 0.0338, 0.0365, 0.0395, 0.042757, 0.0465, 0.05, 0.056925, 0.0625, 
    0.069, 0.075, 0.081972, 0.09, 0.096, 0.1035, 0.111573, 0.12, 0.128, 0.1355, 
    0.145728, 0.16, 0.172, 0.184437, 0.2, 0.2277, 0.2510392, 0.2705304, 
    0.2907501, 0.3011332, 0.3206421, 0.3576813, 0.39, 0.4170351, 0.45, 
    0.5032575, 0.56, 0.625, 0.7, 0.78, 0.86, 0.95, 1.05, 1.16, 1.28, 1.42, 1.55, 
    1.7, 1.855, 2.02, 2.18, 2.36, 2.59, 2.855, 3.12, 3.42, 3.75, 4.07, 4.46, 
    4.9, 5.35, 5.85, 6.4, 7.0, 7.65, 8.4, 9.15, 9.85, 10.0 };

  
  Float kb = 8.6173303E-5;
  Float tev = kb*T;
  teff  *= kb;
  teff2 *= kb;
  // write ( in the 6222 section of thermr output )
  // CONTIO
  // za  awr   0  5  1  0    (this is funky bc the 5 shows up as 1 in the output)
  // TAB1IO
  // 1    1   -1  1  1  2  
  // 2    2  
  // 1e-5 1 emax 1
  // TAB2IO
  // T    0    3  1  1  nne
  // nne  2
  
  int ie;
  Float enow, ep;

  int imax = 20;
  std::vector<double> esi(nne+1), xsi(nne+1), 
  ubar(egrid.size()), p2(egrid.size()), p3(egrid.size()), x(imax);

  // loop over given incident energy grid
  while (true){
    std::cout << " --- 305 --- " << std::endl;
    ie = 0;

    while (true){
      std::cout << " --- 310 --- " << std::endl;
      ie = ie+1;
      enow = egrid[ie-1];
      if ( T > 3000.0 ){
        std::cout << "do temp approx" << std::endl;
        Float tone = 0.0253/kb;
        Float elo = egrid[0];
        enow = elo*exp(log(enow/elo)*log((T/tone)*egrid[egrid.size()-1]/elo)/log(egrid[egrid.size()-1]/elo));
      }
      enow = sigfig(enow,8,0);
      esi[ie-1] = enow;
      xsi[ie-1] = 0.0;
      ubar[ie-1] = 0.0;
      p2[ie-1] = 0.0;
      p3[ie-1] = 0.0;
      ep = 0.0;
      x[0] = ep;
      return;

    }


  } 

}
