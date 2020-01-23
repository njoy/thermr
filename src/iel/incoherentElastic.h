#include <iostream>
#include "coh/coh_util/sigcoh_util/terp.h"

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
}


template <typename Float>
Float getIncohElasticXS( const Float& sigma_b, const Float& dwf, 
  const Float& E, int nAtoms ){
  return sigma_b/(2.0*nAtoms) * (1.0-exp(-4.0*E*dwf)) / ( 2.0*E*dwf ); 
}


template <typename Range>
auto getOnRightGrid( const Range& currentX, const Range& currentY, const Range& desiredX ){
  Range desiredY(desiredX.size(),0.0);
  for (size_t i = 0; i < desiredX.size()-1; ++i){
    desiredY[i] = terp(currentX,currentY,desiredX[i],3);
  }
  desiredY[desiredX.size()-1] = 0.0;
  return desiredY;
}


template <typename Range, typename Float>
auto getIncohElasticXSGrid_withInterpolation( const Float& sigma_b, const Float& debyeWaller, 
  int numAtoms, const Range& currentX, const Range& desiredX ){
  Range currentY(currentX.size(),0.0);
  for ( size_t i = 0; i < currentX.size(); ++i ){
    currentY[i] = getIncohElasticXS( sigma_b, debyeWaller, currentX[i], numAtoms );
  }
  return getOnRightGrid(currentX,currentY,desiredX);
}  


/*
 // This isn't really going to be used right now because although lovely and 
 // beautiful, the answers it returns differ significantly from those provided
 // by legacy THERMR
template <typename Range, typename Float>
auto getIncohElasticXSGrid_directly( const Float& sigma_b, const Float& debyeWaller, 
  int numAtoms, const Range& desiredX ){
  Range currentY(desiredX.size(),0.0);
  for ( size_t i = 0; i < desiredX.size(); ++i ){
    currentY[i] = getIncohElasticXS( sigma_b, debyeWaller, desiredX[i], numAtoms );
  }
  currentY[desiredX.size()-1] = 0.0;
  return currentY;
}  
*/







