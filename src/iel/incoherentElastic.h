#include <iostream>

template <typename Float, typename Range>
auto getIncohElasticEquiprobableAngles( const Float& E, const Float& dwf,
  const Float& sigma_b, int nAtoms, Range& equiprobCosines){
  
  int nAngles = equiprobCosines.size();
  Float sigma = sigma_b/(nAtoms*4.0*dwf*E)*(1.0-exp(-4.0*dwf*E)),
        lastCosine = -1.0, thisCosine, x2, mu_i, c2;

  c2 = 2.0*E*dwf;
  Float exp_4EW = exp(-4.0*E*dwf);

  for ( int i = 0; i < nAngles; ++i ){
    x2 = exp(-c2*(1.0-lastCosine));
    thisCosine = 1.0+(1.0/c2)*log((1-exp_4EW)/nAngles + x2);
    mu_i = nAngles/c2*(exp(-c2*(1.0-thisCosine))*(c2*thisCosine-1.0)-x2*(c2*lastCosine-1.0))/(1.0-exp_4EW);
    equiprobCosines[i] = mu_i;
    lastCosine = thisCosine;
  }
  // Return the scattering cross section sigma(E)
  return sigma_b/(2.0*nAtoms) * (1.0-exp_4EW)/c2;

}




