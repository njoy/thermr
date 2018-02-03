#include <iostream>
#include <vector>

auto computeCrossSections( double e, std::vector<double> fl, 
  std::vector<double> s, double emax, double scon, double recon, int nl,
  std::vector<double> p, int k ){
   // compute cross sections at this energy
   // 210 continue
   double re = 1 / e;
   for ( size_t il = 0; il < nl; ++il ){ 
     s[il] = 0;
   } 

   double elim;
   int last = 0;
   for ( size_t i = 1; i <= k; ++i ){
      double tsq = fl[2*i-1-1];
      elim = tsq * recon;
      if (elim >= e) return;
      double f = fl[2*i-1];
      if (e > emax) f = 0;
      double u = 1 - 2*elim*re;
      int lmax = nl - 1;
      //call legndr(u,p,lmax)
      for ( size_t il = 0; il < nl; ++il ){
         s[il]=s[il]+f*p[il];
      }
      if (i == k) last = 1;
   }
   for ( size_t il = 0; il < nl; ++il ){
      s[il] = s[il]*scon*re;
   }
   if (last == 1) elim=emax;
   if (elim > emax) elim=emax;
   //enext = sigfig(elim,7,-1);
   //if (e.gt.sigfig(enext,7,-1)) enext=sigfig(elim,7,+1)
   return;
}



