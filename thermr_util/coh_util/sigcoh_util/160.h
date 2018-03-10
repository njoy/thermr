#include <iostream>
#include <vector>
#include "form.h"
#include "175.h"
#include "165.h"

auto tausq( int m1, int m2, int m3, double c1, double c2 ){
  double twopis = 4 * M_PI * M_PI;
  double tausqVal = (c1*(m1*m1+m2*m2+m1*m2)+(m3*m3*c2))*twopis;
  return tausqVal;
}

auto do160( int lat, double& w1, double& w2, double& w3, int l1, int l2, int l3,
  int& k, double& c1, double& c2, double& t2, 
  std::vector<double>& wrk, double& wint, int nw, double& ulim, double& recon,
  double& tsqx, double& eps ){
   double tsq, tau, w, f;
   std::cout << "160" << std::endl;
   tsq = tausq( l1, -l2, l3, c1, c2 );
   // This should return and we should go to beginning of i3 loop
   if (tsq <= 0 or tsq > ulim ){ return 1; }
   //if (tsq <= 0 or tsq > ulim ){ end_175_180_185( wek, k, recon, ulim ); }
   tau = sqrt(tsq);
   w = exp(-tsq*t2*wint)*w1*w2*w3/tau;
   f = w * form(lat,l1,-l2,l3);
   //if (k > 0 and tsq > tsqx) go to 165
   if (k > 0 and tsq > tsqx) { 
     bool continueLoop= do165( k, wrk, recon, ulim, f, tsq, eps, nw );
     if ( continueLoop ){ return 1; }
   }
   k = k + 1;
   if ((2*k) > nw) { std::cout << "oh no, sigcoh! Storage exceeded" << std::endl; } 
   wrk[2*k-2] = tsq;
   wrk[2*k-1] = f;
   return 1;

}

