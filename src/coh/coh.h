#include "coh/coh_util/upstk.h"
#include "coh/coh_util/sigcoh.h"
#include <vector>
#include <range/v3/all.hpp>
#include "general_util/sigfig.h"

template <typename Float, typename Range>
auto findLocation( const Range& Egrid, const Float& bragg, int j = 0){
  for (size_t i = j; i < Egrid.size()-1; ++i){
    if (Egrid[i+1] > bragg){
      return i;
    }
  }
  return Egrid.size();
}

template <typename Float, typename Range>
Float addPoint( const Float& x1, const Float& x2, Range& finalE, Range& finalXS,
  Range& vec1, Range& vec2, const Float& emax, const Float& scon, const Float& recon, 
  int nbragg, int& counter, const Float& tol ){ 

  using std::max;

  Range s1 (6,0.0), s2(6,0.0);
  auto enext1 = computeCrossSections( x1, vec1, vec2, emax, scon, recon, s1, nbragg );
  auto enext2 = computeCrossSections( x2, vec1, vec2, emax, scon, recon, s2, nbragg );

  Float xm = (x1+x2)*0.5;

  if ( abs(x1-x2) < 3.0e-5*xm){   // The energies are sufficiently close to each
      finalE[++counter] = x2;     // to each other so that we don't need to 
      finalXS[ counter] = s2[0];  // worry about putting in an intermeidate val
      return enext2; 
  }

  Range sm(6,0.0);
  auto enextM = computeCrossSections( xm, vec1, vec2, emax, scon, recon, sm, nbragg );

  Float ym = (s1[0]+s2[0])*0.5;
  Float test = max(tol*sm[0],1e-6);
  if (abs(sm[0]-ym) <= test ){
    // Can reasonably interpolate here, just add in the right boundary withou
    // any more heartache
    finalE[++counter] = x2;
    finalXS[ counter] = s2[0];
    return enext2;
  }

  else {
    // Add in any intermediate points that we need between left bound and middle
    addPoint( x1, xm, finalE, finalXS, vec1,  vec2, emax, scon, recon, nbragg, counter, tol );
    // Add in any intermediate points that we need between middle and right bound
    return addPoint( xm, x2, finalE, finalXS, vec1,  vec2, emax, scon, recon, nbragg, counter, tol );
  }
}



template <typename Float, typename Range>
auto coh( const Float& temp, int lat, const Float& emax, int numAtoms, const Range& Egrid, const Float tol){
  std::vector<Float> vec1 (5000,0.0), vec2 (5000,0.0);
  auto out = prepareBraggEdges(lat,temp,emax,numAtoms,vec1,vec2);
  int nbragg = std::get<0>(out);
  Float scon = std::get<1>(out);

  Float recon = 5.1803120897E-20;
  auto enext = vec1[0]*recon;

  Range finalE  (100,0.0),
        finalXS (100,0.0);

  finalE[0]  = Egrid[0];
  finalXS[0] = 0.0;

  int counter = 0, i_old = 0;
 
  while (finalE[counter] <= emax) {
    if (counter > 0.8*finalE.size()){ 
        finalE.resize(finalE.size()*2);
        finalXS.resize(finalXS.size()*2);
    }
    int i = findLocation(Egrid,enext);
    if ( i > i_old ){
      for ( int k = i_old+1; k < i+1; ++k ){
        enext = addPoint( finalE[counter], Egrid[k], finalE, finalXS,
                  vec1, vec2, emax, scon, recon, nbragg, counter, tol );
      }
    }
    enext = addPoint( finalE[counter], enext, finalE, finalXS, vec1, 
                      vec2, emax, scon, recon, nbragg, counter, tol );
    i_old = i;
  }

  Float a = sigfig(finalE[counter-1],7,0);
  Float c = sigfig(emax,7,0);

  if (abs(c-a)/a > 1e-10 ){
      ++counter;
      Range s(6,0.0);
      computeCrossSections( emax, vec1, vec2, emax, scon, recon, s, nbragg );
      finalE[counter]    = finalE[counter-1];
      finalXS[counter]   = finalXS[counter-1];
      finalE[counter-1]  = emax;
      finalXS[counter-1] = s[0];
  }
  
  finalE[counter+1] = Egrid[Egrid.size()-1];
  finalE.resize(counter+2);
  finalXS[counter+1] = 0.0;
  finalXS.resize(counter+2);

  return std::make_tuple(finalE,finalXS);

}



