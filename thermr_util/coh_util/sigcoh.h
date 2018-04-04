#include <iostream>
#include <vector>
#include "sigcoh_util/form.h"
#include "sigcoh_util/terp.h"


auto tausq( int m1, int m2, int m3, double c1, double c2 ){
  return (c1*(m1*m1+m2*m2+m1*m2)+(m3*m3*c2))*4*M_PI*M_PI;
}

void swapVals( double& a, double& b ){
  double c = a; a = b; b = c;
}
  
auto computeCrossSections( double e, std::vector<double>& fl, 
    std::vector<double>& s, double emax, double scon, double recon, int nl, 
    std::vector<double> p, int k ){
  // compute cross sections at this energy
  // 210 continue
   double re = 1.0/e;
   double elim;
   for ( int il = 0; il < nl; ++il ){
      s[il]=0;
   }
   int last=0;
   for ( int i = 1; i <= k; ++i ){
      double tsq=fl[2*i-1-1];
      elim=tsq*recon;
      // if (elim >= e) exit
      double f=fl[2*i-1];
      if (e > emax) f=0;
      double u=1-2*elim*re;
      int lmax=nl-1;
      // call legndr(u,p,lmax)
      for ( int il = 0; il < nl; ++il ){
         s[il]=s[il]+f*p[il];
      }
      if (i == k) last=1;
   }
   for ( int il = 0; il < nl; ++il ){
      s[il]=s[il]*scon*re;
   }
   if (last == 1) elim=emax;
   if (elim > emax) elim=emax;

   return;

}




auto sigcoh( double e, double enext, std::vector<double> s, int nl, int lat, 
  double temp, double emax, int natom ){
 /*-------------------------------------------------------------------
  * Compute the first nl Legendre components of the coherent scatter-
  * ing at energy e from lattice type lat.  Here enext is the next
  * Bragg edge.  Initialize if e=0.  A list of reciprocal lattice
  * shells and weights is precomputed and stored for use at all e.
  * Long, closely-spaced shells are grouped together to speed up the
  * calculation.
  *       lat=1  graphite
  *       lat=2  be
  *       lat=3  beo
  *       lat=10 read from endf6
  * nl returns no. of Bragg edges on initialization call.
  *-------------------------------------------------------------------
  */
  int nw, k, i1m, i1, l1, i2m, i2, l2, i3m, i3, l3, l, i, j,
    il, lmax, last;
  double amne, econ, tsqx, a, c, amsc, scoh, wal2, wint, x, w1, w2, w3, tsq, 
    tau, w, f, st, sf, blast, re, t2, ulim, phi, elim, u, twopis, c1, c2, scon;
  std::vector<double> p(6,0);
  int nd = 10;
  std::vector<double> dwf1 { 2.1997, 2.7448, 3.2912, 3.8510, 4.4210, 4.9969, 
    6.1624, 7.3387, 9.6287, 11.992 },
    dwf2 { 3.16663, 3.88842, 4.62944, 5.40517, 6.19880, 7.0042, 8.63665, 
      10.2865, 0, 0 },
    dwf3 {  2.153, 2.6374, 3.1348, 3.6513, 4.1798, 4.7164, 5.8052, 6.9068, 
      0, 0 },
    tmp { 296, 400, 500, 600, 700, 800, 1000, 1200, 1600, 2000 };

  
  double gr1 = 2.4573e-8, gr2 = 6.700e-8, gr3 = 12.011e0, gr4 = 5.50e0, 
    be1 = 2.2856e-8, be2 = 3.5832e-8, be3 = 9.01e0, be4 = 7.53e0, 
    beo1 = 2.695e-8, beo2 = 4.39e-8, beo3 = 12.5e0, beo4 = 1.e0, 
    sqrt3 = 1.732050808e0, cw = 0.658173e-15, eps = 0.05e0, 
    zero = 0, hbar = 1.05457266e-27, amu = 1.6605402e-24, amassn = 1.008664904,
    ev = 1.60217733e-12;

  // save k,recon,scon

  // initialize.
  // go to 210 
  if (e > 0) {
//  computeCrossSections( e, fl, s, emax, scon, recon, nl, p, k );
  }
 
  twopis = (4*M_PI*M_PI);
  amne = amassn * amu;
  econ = ev * 8 * ( amne / hbar ) / hbar;
  tsqx = econ / 20;

  if (lat == 1){
    // graphite constants
    a = gr1; c = gr2; amsc = gr3; scoh = gr4 / natom; 
    wal2 = terp(tmp,dwf1,nd,temp,2);
  }
  else if (lat == 2) {
    // beryllium constants
    a = be1; c = be2; amsc = be3; scoh = be4/natom;
    wal2 = terp(tmp,dwf2,nd,temp,2);
  }
  else if (lat == 3){
    // beryllium oxide constants
    a = beo1; c = beo2; amsc = beo3; scoh = beo4/natom;
    wal2 = terp(tmp,dwf3,nd,temp,2);
  } 
  else {
    std::cout << "OH NO! Error over here. Illegal lat value" << std::endl;
    // call error('sigcoh','illegal lat.',' ')
  }
  c1 = 4 / ( 3 * a * a );
  c2 = 1 / ( c * c );
  scon = scoh * ( 16 * M_PI*M_PI )/( 2 * a * a * c * sqrt3 * econ );
  wint = cw * amsc * wal2;
  t2 = hbar / ( 2.0 * amu * amsc );
  ulim = econ * emax;
  nw = 10000;
  std::vector<double> wrk(nw,0.0);


  // compute and sort lattice factors.
  phi = ulim / twopis;
  i1m = a * std::pow(phi,0.5) + 1;
  k = 0;

  for ( int i1 = 1; i1 <= i1m; ++i1 ){
    l1 = i1 - 1;
    i2m = 0.5 * ( l1 + sqrt( 3 * ( a * a * phi - l1 * l1 ) ) ) + 1;;
    for ( int i2 = i1; i2 <= i2m; ++i2 ){
      l2 = i2 - 1;
      x = phi - c1 * ( l1 * l1 + l2 * l2 - l1 * l2 );
      i3m = c * std::pow(x,0.5) + 1;
      for ( int i3 = 1; i3 <= i3m; ++i3 ){
        l3 = i3 - 1;
        w1 = (l1 == l2) ? 1 : 2;

        if      (l1 == 0 and l2 == 0) { w2 = 0.5; }
        else if (l1 == 0 or  l2 == 0) { w2 = 1;   }
        else                          { w2 = 2;   } 
        
        w3 = (l3 == 0) ? 1 : 2;
         
        for ( double&& l2: { l2, -l2 } ){
          tsq = tausq(l1,l2,l3,c1,c2);
          if (tsq > 0 and tsq <= ulim ){
            tau = std::pow(tsq,0.5);
            w = exp(-tsq*t2*wint)*w1*w2*w3/tau;
            f = w * form( lat, l1, l2, l3 );

            if (k > 0 and tsq > tsqx) {
              for ( int i = 1; i <= k; ++i ){
                if (tsq < wrk[2*i-1-1] or tsq >= (1+eps)*wrk[2*i-1-1]){
                  if ( i == k ){ 
                    k = k + 1;
                    if ((2*k) > nw) std::cout << "call error('sigcoh','storage exceeded.',' ')" << std::endl;
                    if ((2*k) > nw) { return wrk; } 
                    wrk[2*k-2] = tsq;
                    wrk[2*k-1] = f;
                    break;
                  }
                  continue; 
                }
                wrk[2*i-1]=wrk[2*i-1]+f;
                break;
              }
            }
            else {
              k = k + 1;
              if ((2*k) > nw) { std::cout << "oh no, sigcoh! Storage exceeded" << std::endl; } 
              if ((2*k) > nw) { return wrk; } 
              wrk[2*k-2] = tsq;
              wrk[2*k-1] = f;
            }

          } 
        }
      }
    } 
  }

  for ( int i = 1; i <= k - 1; ++i ){
    for ( int j = i + 1; j <= k; ++j ){
      if (wrk[2*j-1-1] < wrk[2*i-1-1]) {
        swapVals( wrk[2*j-2], wrk[2*i-2] );
        swapVals( wrk[2*j-1], wrk[2*i-1] );
      } 
    }
  }
  k += 1;
  wrk[2*k-2]=ulim;
  wrk[2*k-1]=wrk[2*k-2-1];
  nw=2*k;
  enext=wrk[1-1]/econ;
  nl=k;
  return wrk;
  // return wrk as fl

}

















