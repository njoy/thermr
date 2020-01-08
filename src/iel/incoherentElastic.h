#include <iostream>

template <typename Float, typename Range>
auto getIncohElasticDataSingleEnergy( const Float& E, const Float& dwf,
  Range& equiprobCosines){
  
  int nAngles = equiprobCosines.size();
  Float lastCosine = -1.0, thisCosine, exp_2EW_cos, mu_i,
        exp_4EW = exp(-4.0*E*dwf),
        inv_2EW = 1.0/(2.0*E*dwf);

  for ( int i = 0; i < nAngles; ++i ){
    exp_2EW_cos = exp( -2.0*E*dwf*(1.0-lastCosine) );
    thisCosine  = 1.0+inv_2EW*log((1-exp_4EW)/nAngles + exp_2EW_cos );

    mu_i = nAngles * inv_2EW * 
           ( exp(-2.0*E*dwf*(1.0-thisCosine))*(2.0*E*dwf*thisCosine-1.0) - 
             exp_2EW_cos*(2.0*E*dwf*lastCosine-1.0) ) / (1.0-exp_4EW);
    equiprobCosines[i] = mu_i;
    lastCosine = thisCosine;
  }
  // Return the scattering cross section sigma(E)
}


template <typename Float>
Float getIncohElasticXS( const Float& sigma_b, const Float& dwf, 
  const Float& E, int nAtoms ){
  return sigma_b/(2.0*nAtoms) * (1.0-exp(-4.0*E*dwf)) / ( 2.0*E*dwf ); 
}



/*

template <typename Float, typename Range>
auto getIncohElasticDataAllEnergies( const Range& Egrid, const Float& Emax, 
  const Float& sigma_b, int nAtoms, int nBins ){
  
  for 

}

*/






