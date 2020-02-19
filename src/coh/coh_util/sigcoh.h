#include "coh/coh_util/sigcoh_util/form.h"
#include "coh/coh_util/sigcoh_util/terp.h"
#include "coh/coh_util/sigcoh_util/legndr.h"
#include "general_util/sigfig.h"
#include "generalTools/constants.h"
#include <iostream>
#include <range/v3/all.hpp>

template <typename Range, typename Float>
bool finish( int& k, const Float& f, const Float& tau_sq, Range& vec1, Range& vec2 ){
  k++;
  if ((unsigned) k > 10000){ return true; } //std::cout << "storage exceeded" << std::endl;
  if ((unsigned) k > vec1.size()){
    vec1.resize(vec1.size()*2);
    vec2.resize(vec2.size()*2);
  }
  vec1[k-1] = tau_sq; vec2[k-1] = f;
  return false;
}


template <typename Float>
auto tausq( int m1, int m2, int m3, Float c1, Float c2 ){
  /* This function computes the value tau^2, which is defined in the General
   * Atomics HEXSCAT documentation. On pg. 62 of the HEXSCAT pdf document
   * (which is in HEXSCAT Appendix, Section 1: Formulation), we have that 
   *     tau^2/4pi^2  = (4/3a^2) * ( l1^2 + l2^2 + l1*l2 ) + l3 / c^2
   */
  return (c1*(m1*m1+m2*m2+m1*m2)+c2*(m3*m3))*4*M_PI*M_PI;
}

template <typename Float>
void swapVals( Float& a, Float& b ){
  Float c = a; a = b; b = c;
}
  

template <typename Range, typename Float>
auto prepareBraggEdges( int lat, Float temp, Float emax, int natom, Range& vec1, Range& vec2 ){
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
  int l1;
  Float w, f;
  
  // These are lattice factors (borrowed from HEXSCAT)
  Float 
    gr1 = 2.4573e-8,//http://www.phy.ohiou.edu/~asmith/NewATOMS/HOPG.pdf
    gr2 = 6.700e-8, //http://www.phy.ohiou.edu/~asmith/NewATOMS/HOPG.pdf 
    be1 = 2.2856e-8,//http://periodictable.com/Properties/A/LatticeConstants.html
    be2 = 3.5832e-8,//http://periodictable.com/Properties/A/LatticeConstants.html 
    beo1 = 2.695e-8,//https://link.springer.com/chapter/10.1007%2F10681719_737
    beo2 = 4.39e-8; //https://link.springer.com/chapter/10.1007%2F10681719_737
                    // II-VI and I-VII Compounds; Semimagnetic Compounds
  // These are masses
  Float gr3  = 12.011e0, 
        be3  = 9.01e0, 
        beo3 = 12.5e0;  // Mass of BeO is actually 25, but apparently we 
                        // divide by 2 because I suppose avg mass per atom
  // These are the characteristic coherent cross sections for hte material.
  // These first appear in Eq. 222 on pg. 166.
  Float gr4  = 5.50, // pg. 18 Neutron Physics Karl-Heinrich Beckurts, Karl Wirtz
        be4  = 7.53, // pg. 18 Neutron Physics Karl-Heinrich Beckurts, Karl Wirtz
        beo4 = 1.0;

  Float cw = 0.658173e-15;

  // Temperatures interpolated over when trying to get correct Debye-Waller
  // Coefficient.
  Range temps {296, 400, 500, 600, 700, 800, 1000, 1200, 1600, 2000};

  // Debye Waller to interpolate over
  // This appears in Eq. 220 of the manual, as 'W' in the exponential term
  // This also appears in the General Atomics HEXSCAT report as the integral
  // shown below (as wal2 description)
  Range dwf ( temps.size() );
  Float a, c, amsc, scoh;
  if (lat == 1){       // GRAPHITE
    a = gr1; c = gr2; amsc = gr3; scoh = gr4/natom; 
    dwf = {2.1997,2.7448,3.2912,3.8510,4.4210,4.9969,6.1624,7.3387,9.6287,11.992};      
  }
  else if (lat == 2) { // BERYLLIUM
    a = be1; c = be2; amsc = be3; scoh = be4/natom;
    dwf = {3.16663,3.88842,4.62944,5.40517,6.19880,7.0042,8.63665,10.2865,0,0};  
  }
  else if (lat == 3){  // BERYLLIUM OXIDE
    a = beo1; c = beo2; amsc = beo3; scoh = beo4/natom;
    dwf = {2.153,2.6374,3.1348,3.6513,4.1798,4.7164,5.8052,6.9068,0,0};
  } 
  else {
    std::cout << "OH NO! Error over here. Illegal lat value" << std::endl;
    throw std::exception();
  }


 /* If this is the first entry (E=0) for an ENDF-III type material, the 
  * appropriate lattice constants are selected and the Debye-Waller coefficient 
  * is obtained for the desired temperature by interpolation. Then the 
  * reciprocal lattice wave vectors and structure factors are computed, 
  * sorted into shells, and stored for later use.
  */
  Float econ = ev * 8 * (massNeutron)/(hbar*hbar*1e4);
  Float tsqx = econ / 20;
  Float scon = scoh * ( 16.0 * M_PI*M_PI )/( 2.0 * a * a * c * sqrt(3) * econ );

  // wal2 is (supposed to be) equal to
  //
  //       1   |' f(w)        
  //      ---  |  ----  coth( w / kb2 T ) dw
  //       M  ,|   w   
  //
  // this seems to be the case from pg. 5 of the General Atomics HEXSCAT code,
  // in the table of descriptions for input values into the original HEXSCAT.
  //
  // In this case, the main equation we're trying to compute ( Eq. 1 in the
  // GA HEXSCAT documentation, on pg. 1 ) can be computed using an exponential 
  // term 
  //                exp[ ( -hbar2^2 tau^2 / 2 ) * wal2 ]

  Float wal2 = terp(temps,dwf,temp,2);
  Float wint = cw*amsc*wal2,
          c1 = 4.0/(3.0*a*a),
          c2 = 1.0/(c*c),
          t2 = hbar*1e4/(2.0*amu3*amsc);
  // This is the tau^2 value that corresponds to the maximum considered energy.
  // Calculated according to Eq. 223.
  Float ulim = emax * ev3 * 8.0 * massNeutron3/ (hbar3*hbar3*1e4);
  Float phi = ulim / ( 4 * M_PI * M_PI );  // phi = ( tau_max/2pi )^2
  int k = 0;

  // l1 --> 0 : a * tau_max / 2pi + 1
  // l1max = alpha * sqrt(phi), on pg 63 of the HEXSCAT document pdf. 
  int i1m = a * sqrt(phi) + 1;

  for ( int l1 = 0; l1 < i1m; ++l1 ){
    int i2m = 0.5 * ( l1 + sqrt( 3 * ( a * a * phi - l1 * l1 ) ) ) + 1;

    for ( int l2 = l1; l2 < i2m; ++l2 ){
      int i3m = c * sqrt(phi - c1*( l1*l1 + l2*l2 - l1*l2 )) + 1;

      for ( int l3 = 0; l3 < i3m; ++l3 ){

        Float w1 = (l1 == l2) ? 1 : 2; // M1 on pg.3 of HEXSCAT appendix
        Float w2 = (l1 == 0 and l2 == 0) ? 0.5 
                 : (l1 == 0 or l2 == 0 ) ? 1.0 : 2.0; // M2 in HEXSCAT appendix
        Float w3 = (l3 == 0) ? 1 : 2;  // M3 on pg.3 of HEXSCAT appendix
         
        // Consider +/- l2 because of Eq. 4 on pg. 4 of the HEXSCAT appendix.
        for ( Float&& l2: { l2, -l2 } ){
          Float tau_sq = tausq(l1,l2,l3,c1,c2);

          if (tau_sq > 0 and tau_sq <= ulim ){
            Float tau = sqrt(tau_sq);

            // w1*w2*w3 --> weighting factor M, as dfined on pg. 3 of the 
            // General Atomics HEXSCAT appendix
            w = exp(-tau_sq*t2*wint)*w1*w2*w3/tau;
            f = w * form( lat, l1, l2, l3 );

            if (k > 0 and tau_sq > tsqx) {
              for ( int i = 1; i <= k; ++i ){

                // As tau_i gets large, the values of tau_i get more and more 
                // closely spaced together, so a range of tau values can be 
                // lumped together to give a single effective tau_i and f_i.
                // This uses a 5% (eps) grouping factor.
                if (tau_sq < vec1[i] or tau_sq >= 1.05*vec1[i]){
                  if ( i == k ){ 
                    if (finish(k,f,tau_sq,vec1,vec2)){ 
                      vec1.resize(k);
                      vec2.resize(k);
                      return scon;
                    }
                    break;
                  }
                }
                else {                        
                  if ( i > int(vec2.size()) ){ vec1.resize(vec1.size()*2); 
                                               vec2.resize(vec2.size()*2); 
                  }
                  vec2[i] += f;
                  break;
                }
              }
            }
            else {
              if (finish(k,f,tau_sq,vec1,vec2)){ 
                vec1.resize(k);
                vec2.resize(k);
                return scon;
              }
            }
          } 
        }
      }
    } 
  }

  for ( int i = 0; i < k - 1; ++i ){
    for ( int j = i + 1; j < k; ++j ){
      if (vec1[j] < vec1[i]) {
        swapVals(vec1[j], vec1[i]);
        swapVals(vec2[j], vec2[i]);
      } 

    }
  }

  ++k;
  vec1[k-1] = ulim;
  vec2[k-1] = vec2[k-2];
  vec1.resize(k);
  vec2.resize(k);
  return scon; 
  // Scon contains material specific constants

}



template <typename Range, typename Float>
auto computeCrossSections( Float e, Range& vec1, Range& vec2, Float emax, 
  Float scon, Float recon, Range& s ){
  // compute cross sections at this energy
   Float elim;
   for ( Float& sVal : s ){ sVal = 0.0; }
   int nl = s.size();
   Range p(nl,0.0);
   int last = 0;

   for ( size_t i = 0; i < vec1.size(); ++i ){
      Float tau_sq=vec1[i];
      elim = tau_sq*recon;

      if (elim >= e) { break; }
      Float f = ( e > emax ) ? 0.0 : vec2[i];
      Float u = 1.0-2.0*elim/e;
      // u here is equal to fl for l = 1 (P1 component).
      // This is defined in the General Atomics HEXSCAT paper, in Part 1 
      // Formulation. If l == 0, fl = 1. But if l == 1, then
      //            fl = 1 - tau^2 lambda^2 / 8 pi^2
      //    which simplifies to 
      //                 1 - tau^2 hbar2^2 / 4 m_n E
      legndr(u,p,nl-1);
      for ( int il = 0; il < nl; ++il ){
         s[il] += f*p[il];
      }
      if (i == vec1.size()-1) { last = 1; }
   }
   for ( int il = 0; il < nl; ++il ){
      s[il] *= scon/e;
   }
   if (last == 1 or elim > emax ) { elim=emax; }

   Float enext = sigfig(elim,7,-1);
   if (e > sigfig(enext,7,-1)){
     enext = sigfig(elim,7,+1);
   }
   return enext;
}



