#include "coh/coh_util/upstk.h"
#include "coh/coh_util/sigcoh.h"
#include <vector>
#include <range/v3/all.hpp>
#include "general_util/sigfig.h"
#include "generalTools/constants.h"

template <typename Float, typename Range>
auto findLocation( const Range& Egrid, const Float& bragg, int j = 0){
  for (size_t i = j; i < Egrid.size()-1; ++i){
    if ( bragg <= Egrid[i+1] ){ return i; }
  }
  return Egrid.size()-1;
}

template <typename Range>
auto findWhichPointsWeAdded( const Range& finalE, const Range& Egrid ){
  Range addedE (finalE.size(),0.0);
  int addedE_counter = 0;
  for (size_t i = 0; i < finalE.size(); ++i){
    bool newValue = true;
    for ( size_t j = 0; j < Egrid.size(); ++j ){
      if ( abs(finalE[i]-Egrid[j])/finalE[i] < 2e-6 ){
        newValue = false;
        break;
      }
    }
    if (newValue) { addedE[addedE_counter++] = finalE[i]; }
  }
  addedE.resize(addedE_counter);
  return addedE;
}




template <typename Float, typename Range>
Float addPoint( const Float& x1, const Float& x2, Range& finalE, Range& finalXS,
  Range& vec1, Range& vec2, const Float& emax, const Float& scon, const Float& recon, 
  int& counter, const Float& tol ){ 
  using std::max;
  Range s1 (6,0.0), s2(6,0.0), sMid(6,0.0);

                                 computeCrossSections( x1, vec1, vec2, emax, scon, recon, s1 );
  auto eNextAfterRightBoundary = computeCrossSections( x2, vec1, vec2, emax, scon, recon, s2 );

  if (counter+5 > int(finalE.size())){ finalE.resize( finalE.size() *2);
                                       finalXS.resize(finalXS.size()*2); }

  Float xm = (x1+x2)*0.5;
  if ( abs(x1-x2) < 3.0e-5*xm){   // The energies are sufficiently close to each
      finalE[++counter] = x2;     // to each other so that we don't need to 
      finalXS[ counter] = s2[0];  // worry about putting in an intermeidate val
      return eNextAfterRightBoundary;
  }

  computeCrossSections( xm, vec1, vec2, emax, scon, recon, sMid );

  Float ym = (s1[0]+s2[0])*0.5;
  if (abs(sMid[0]-ym) <= max(tol*sMid[0],1e-6) ){
    // Can reasonably interpolate here, just add in the right boundary without
    // any more heartache
    finalE[++counter] = x2;
    finalXS[ counter] = s2[0];
    return eNextAfterRightBoundary;
  }

  else {
    // Add in any intermediate points that we need between left bound and middle
    addPoint( x1, xm, finalE, finalXS, vec1,  vec2, emax, scon, recon, counter, tol );
    // Add in any intermediate points that we need between middle and right bound
    return addPoint( xm, x2, finalE, finalXS, vec1,  vec2, emax, scon, recon, counter, tol );
  }
}

template <typename Float, typename Range>
auto coh( const Float& temp, int lat, const Float& emax, int numAtoms, const Range& Egrid, const Float tol){
  using std::abs;
  std::vector<Float> vec1 (50,0.0), vec2 (50,0.0);
  Float scon  = prepareBraggEdges(lat,temp,emax,numAtoms,vec1,vec2);
  Float recon = (hbar*hbar)/(ev*8*(massNeutron))*1e4; 
  auto enext = vec1[0]*recon;

  Range finalE  (100,0.0), 
        finalXS (100,0.0); 
  finalE[0]  = Egrid[0];

  int counter = 0, i_old = 0;// count2 = 0;

  while (finalE[counter] <= emax) {
    // Where is out next bragg edge? 
    int i = findLocation(Egrid,enext); // Egrid[i] <= Enext <= Egrid[i+1]

    // i_old is the index of where our old bragg peak would fit in Egrid. If 
    // i_old == i, then our bragg peaks are right next to each other and we 
    // can just add them in without putting in an intermediate Egrid point. 
    // Otherwise we need to loop through all intermediate Egrid points and add
    // them in
    for ( int k = i_old+1; k < i+1; ++k ){
      enext = addPoint( finalE[counter], Egrid[k], finalE, finalXS, vec1, vec2, 
                         emax, scon, recon, counter, tol );
    }
    // Add in the next bragg peak to our finalE
    enext = addPoint( finalE[counter], enext, finalE, finalXS, vec1, vec2, emax, 
                      scon, recon, counter, tol );
    i_old = i;
  }

  Float a = sigfig(finalE[counter-1],7,0),
        c = sigfig(emax,7,0);

  if (abs(c-a)/a > 1e-10 ){
    ++counter;
    Range s(6,0.0);
    computeCrossSections( emax, vec1, vec2, emax, scon, recon, s );
    finalE[ counter] = finalE[ counter-1]; finalE[ counter-1] = emax;
    finalXS[counter] = finalXS[counter-1]; finalXS[counter-1] = s[0];
  }
  
  finalE[counter+1]  = Egrid[Egrid.size()-1]; finalE.resize( counter+2);
  finalXS[counter+1] = 0.0;                   finalXS.resize(counter+2);


  Range addedE = findWhichPointsWeAdded( finalE, Egrid );

  finalE.resize( counter+2);
  finalXS.resize(counter+2);

  return std::make_tuple(finalE,finalXS,addedE);

}



