#include <iostream>
#include <vector>

auto tausq( int m1, int m2, int m3, double c1, double c2 ){
  double twopis = 4 * M_PI * M_PI;
  double tausqVal = (c1*(m1*m1+m2*m2+m1*m2)+(m3*m3*c2))*twopis;
  return tausqVal;
}

auto do160( int lat, double w1, double w2, double w3, int l1, int l2, int l3,
  double tsq, int k, double f ){
   tsq = tausq( l1, -l2, l3, c1, c2 );
   // if (tsq <= 0 or tsq > ulim ) go to 175
   tau = sqrt(tsq);
   w = exp(-tsq*t2*wint)*w1*w2*w3/tau;
   f = w * form(lat,l1,-l2,l3);
   // if (k > 0 and tsq > tsqx) go to 165
   k = k + 1;
   if ((2*k).gt.nw) std::cout << "oh no, sigcoh! Storage exceeded" << std::endl; } 
   wrk[2*k-1] = tsq;
   wrk[2*k] = f;

}

