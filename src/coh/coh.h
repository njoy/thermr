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
  int nbragg, int& finalE_counter, const Float& tol, bool check ){
  Range s1 (6,0.0), s2(6,0.0);

  auto enext1 = computeCrossSections( x1, vec1, vec2, emax, scon, recon, s1, nbragg );
  auto enext2 = computeCrossSections( x2, vec1, vec2, emax, scon, recon, s2, nbragg );

  Float xm = (x1+x2)*0.5;
  if ( abs(x1-x2) < 3.0e-5*xm){
      finalE[++finalE_counter] = x2;
      finalXS[finalE_counter] = s2[0];
      //std::cout << x2 << "    " << s2[0] << std::endl;
  //std::cout << std::endl;
      return enext2; 
  }


  Range sm(6,0.0);
  //std::cout << std::endl;
  //std::cout << std::endl;
  auto enextM = computeCrossSections( xm, vec1, vec2, emax, scon, recon, sm, nbragg );

  auto test1 = tol*sm[0];
  if (test1 < 1e-6) { test1 = 1e-6; }
  auto ym = (s1[0]+s2[0])*0.5;

  if ( check == false ){
    finalE[++finalE_counter] = x2;
    finalXS[finalE_counter] = s2[0];
    //std::cout << x2 << "    " << s2[0] << std::endl;
  //std::cout << std::endl;
    return enext2;
  }

  if (abs(sm[0]-ym) <= test1 ){
    finalE[++finalE_counter] = x2;
    finalXS[finalE_counter] = s2[0];
    //std::cout << x2 << "    " << s2[0] << std::endl;
  //std::cout << std::endl;
    return enext2;
  }
  else {
    addPoint( x1, xm, finalE, finalXS, vec1,  vec2, emax, scon, recon, nbragg, finalE_counter, tol, check );
  //std::cout << std::endl;
    return addPoint( xm, x2, finalE, finalXS, vec1,  vec2, emax, scon, recon, nbragg, finalE_counter, tol , check);
  }
}



template <typename Float, typename Range>
auto coh( const Float& temp, int lat, const Float& emax, int numAtoms, const Range& Egrid, const Float tol){
  std::vector<Float> vec1 (5000,0.0), vec2 (5000,0.0);
  auto out = prepareBraggEdges(lat,temp,emax,numAtoms,vec1,vec2);
  int nbragg = std::get<0>(out);
  Float scon = std::get<1>(out);

  //std::cout << vec2[0] << "   " << vec2[1] << "   " << vec2[2] << std::endl;
  Float recon = 5.1803120897E-20;




  Range braggs(int(5.5*nbragg),0.0);
  auto enext = vec1[0]*recon;

  Range finalE (1000,0.0);
  Range finalXS (1000,0.0);
  finalE[0] = Egrid[0];
  finalXS[0] = 0.0;
  int finalE_counter = 0;
  int i_old = 0;
 
  for ( int ibragg = 0; ibragg < 5.5*nbragg; ++ibragg ){
    if (finalE[finalE_counter] > emax){ break; }
    int i = findLocation(Egrid,enext);
    if ( i > i_old ){
      for ( int k = i_old+1; k < i+1; ++k ){
        enext = addPoint( finalE[finalE_counter], Egrid[k], finalE, finalXS,
                  vec1, vec2, emax, 
                  scon, recon, nbragg, finalE_counter, tol, false );
      }
    }
    enext = addPoint( finalE[finalE_counter], enext, finalE, finalXS, vec1, vec2, emax, 
              scon, recon, nbragg, finalE_counter, tol, true );
    i_old = i;
  }
  if ( sigfig(finalE[finalE_counter],7,0) - sigfig(emax,7,0)  > 1e-6 ){
      finalE[finalE_counter+1] = finalE[finalE_counter];
      finalXS[finalE_counter+1] = finalXS[finalE_counter];
      Range s(6,0.0);
      computeCrossSections( emax, vec1, vec2, emax, scon, recon, s, nbragg );
      finalE[finalE_counter] = emax;
      finalXS[finalE_counter] = s[0];
      ++finalE_counter;
  }



  
  finalE[finalE_counter+1] = Egrid[Egrid.size()-1];
  finalE.resize(finalE_counter+2);
  finalXS[finalE_counter+1] = 0.0;
  finalXS.resize(finalE_counter+2);
  std::cout.precision(15);
  for ( int j = 0; j < finalE_counter+2; ++j ){
      //finalXS[j] = sigfig(finalXS[j],7,0);
      std::cout << finalE[j] << "    " << finalXS[j] << std::endl;
  }

  return std::make_tuple(finalE,finalXS);

}




/*
template <typename Float, typename Range>
auto coh( const Float& temp, int lat, const Float& emax, int numAtoms, 
  const Range& Egrid, const Float& tol){
  std::vector<Float> vec1 (5000,0.0), vec2 (5000,0.0);
  auto out = prepareBraggEdges(lat,temp,emax,numAtoms,vec1,vec2);
  //std::cout << vec1[0] << "   " << vec1[1] << "   " << vec1[2] << std::endl;
  //std::cout << vec1[3] << "   " << vec1[4] << "   " << vec1[5] << std::endl;
  //std::cout << vec2[0] << "   " << vec2[1] << "   " << vec2[2] << std::endl;
  //std::cout << vec2[3] << "   " << vec2[4] << "   " << vec2[5] << std::endl;
  int nbragg = std::get<0>(out);
  Float scon = std::get<0>(out);
  Range x(5,0.0);//, y(5,0.0), z(5,0.0)

  int ix  = 1;
  int j   = 0;
  int iex = 0;
  int nlt = 5;
  if ( int(Egrid.size()) < nlt ){ nlt = int(Egrid.size()); }

  Range finalE ( 1000, 0.0 );
  int final_i = 0;

  Range stk ( 1000, 0.0 );
  Float recon = 5.1803120897E-20;
  Float enext = vec1[0]*recon;
  enext = sigfig(enext,7,-1);
  size_t i1 = 1;

  while (true){
    //std::cout << " --- 100 --- " << std::endl;
    for ( i1 = 1; i1 <= Egrid.size(); ++i1 ){
      ++iex;
      //std::cout << iex << "    " << Egrid[iex-1] << std::endl;
      //std::cout << iex << std::endl;
      if ( Egrid[iex-1] > (1.0+1e-10)*enext ){ break; }
      //finalE[iex-1] = Egrid[iex-1];
      ++final_i;
      finalE[final_i-1] = Egrid[iex-1];
      x[0] = Egrid[iex-1];
      ++j;
    }

    //std::cout << " --- 105 --- " << std::endl;
    ++ix;
    x[ix-1] = Egrid[iex-1];
    //std::cout << Egrid[iex-1] << std::endl;
    if (ix >= nlt){ break; }
  }
  //std::cout << (x|ranges::view::all) << std::endl;

  //std::cout << iex << "    " << ix << std::endl;
  Range ej(20,0.0);
  Range s(6,0.0);
  Float e;
  //std::cout << enext << std::endl;
  //std::cout << iex << std::endl;
  finalE[final_i-1] = Egrid[final_i-1];
  finalE[final_i-0] = enext;
  //std::cout << (finalE|ranges::view::all) << std::endl;
  e = enext;
  enext = computeCrossSections( enext, vec1, vec2, emax, scon, recon, s, nbragg );
  //std::cout << enext << "   " << s[0] << std::endl;
  stk[0] = e;
  stk[1] = s[0];
  //std::cout << enext << std::endl;
  //std::cout << (finalE|ranges::view::all) << std::endl;
  //std::cout << (s|ranges::view::all) << std::endl;
    
  int i = 1;
  while ( true ){
    std::cout << " --- 120 --- " << std::endl;
    e = enext;
    enext = computeCrossSections( enext, vec1, vec2, emax, scon, recon, s, nbragg );
    //std::cout << enext << "  " << s[0] << std::endl;
    upstk(e,s,stk,1,i,2);
    //std::cout << (stk|ranges::view::all) << std::endl;
    if ( e > 5e-4 ){ return; }
    for ( int ix = 0; ix < 5; ++ix ){
      std::cout << x[ix] << std::endl;
      if (x[ix] >= stk[0+(i-1)*2]*(1.0+1e-10)){
          std::cout << " --- 135 --- " << std::endl;

          if (x[ix]>= stk[0+(i-2)*2]*(1.0-1e-10)){
            std::cout << " --- 140 --- " << std::endl;
            if ( i == 20 ){ std::cout << "go to 160" << std::endl; }
            auto xm = (stk[0+(i-2)*2] + stk[0+(i-1)*2])*0.5;
            xm = sigfig(xm,7,0);
            if ( (stk[0+(i-2)*2] - stk[0+(i-1)*2]) < 3e-5*xm ){ 
              std::cout << " --- 160 --- " << std::endl; 
              ++j;
              ej[0] = stk[0+(i-1)*2];
  
              std::cout << " --- 170 --- " << std::endl; 
              if (ej[0] <= x[2]*(1.0+1e-10) or iex == int(Egrid.size())){
                std::cout << " --- 190 --- " << std::endl;
  
                if (stk[0+(i-1)*2] > emax*(1.0+1e-10)){
                    std::cout << " go to 230 (we're basically finished)" << std::endl;
                    return;
                }
                std::cout << " want to add " << ej[0] << " xs in other reactions" << std::endl;
                --i;
                if ( i > 1 ){
                    std::cout << " go to 125 " << std::endl;
                }
                break; // go to 120
            }
            return;
          }
          std::cout << "continue with 140 shit" << std::endl;

          //std::cout << stk[0+(i-2)*2] << "   " << stk[0+(i-1)*2] << "   " << i << std::endl;
          return;
        }
        e = x[ix];
        enext = computeCrossSections( e, vec1, vec2, emax, scon, recon, s, nbragg );
        upstk(e,s,stk,1,i,2);


        break;
      }
    }
  }
  return;
  std::cout << tol << std::endl;
}
*/
