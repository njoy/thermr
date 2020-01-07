#include <iostream>
//#include "coh/coh_util/sigcoh_util/terp.h"
//#include "iel/iel_util/terpa.h"
//#include <cmath>



template <typename Float>
auto getIncohElasticEquiprobableAngles( const Float& E, const Float& DebyeWaller,
  const Float& sigma_b, int numAtoms, int numAngles){
  
  Float sigma = sigma_b/(numAtoms*4.0*DebyeWaller*E)*(1.0-exp(-4.0*DebyeWaller*E));
  Float u = -1.0, unow, x2;
  for ( int n = 1; n <= numAngles; ++n ){
    auto x2 = exp(-2.0*E*DebyeWaller*(1.0-u));
    unow = 1.0+1.0/(2.0*E*DebyeWaller)*log((1.0-exp(-4.0*E*DebyeWaller))/n + x2);
    std::cout << x2 << std::endl;
  }

}




