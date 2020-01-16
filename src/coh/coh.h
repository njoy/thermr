#include "coh/coh_util/upstk.h"
#include "coh/coh_util/sigcoh.h"




template <typename Float, typename Range>
auto coh( const Float& temp, int lat, const Float& emax, int numAtoms, const Range& Egrid){
  std::vector<Float> vec1 (5000,0.0), vec2 (5000,0.0);
  auto out = prepareBraggEdges(lat,temp,emax,numAtoms,vec1,vec2);
  //std::cout << vec2[0] << "   " << vec2[1] << "   " << vec2[2] << std::endl;
  //std::cout << vec2[3] << "   " << vec2[4] << "   " << vec2[5] << std::endl;
  int nbragg = std::get<0>(out);
  Float scon = std::get<0>(out);

  int ix  = 1;
  int j   = 0;
  int iex = 0;
  int nlt = 5;
  if ( int(Egrid.size()) < nlt ){ nlt = int(Egrid.size()); }


  Float recon = 5.1803120897E-20;
  Float enext = vec1[0]*recon;
  size_t i = 1;

  while (true){
    std::cout << " --- 100 --- " << std::endl;
    for ( i = 1; i <= Egrid.size(); ++i ){
      ++iex;
      if ( Egrid[i-1] > (1.0+1e-10)*enext ){ break; }
      ++j;
    }

    std::cout << " --- 105 --- " << std::endl;
    ++ix;
    if (ix >= nlt){ break; }
  }

  std::cout << iex << "    " << ix << std::endl;
  Range s(6,0.0);
  enext = computeCrossSections( enext, vec1, vec2, emax, scon, recon, s, nbragg );
  std::cout << (s|ranges::view::all) << std::endl;
    

  
}
