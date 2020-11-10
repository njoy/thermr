#include "ENDFtk.hpp"
using namespace njoy::ENDFtk;
using ContinuumEnergyAngle  = section::Type<6>::ContinuumEnergyAngle;
using IncoherentElastic = section::Type<7,2>::IncoherentElastic;
using ThermalScatteringData = section::Type<6>::ContinuumEnergyAngle::ThermalScatteringData;
using Variant = section::Type< 6 >::ContinuumEnergyAngle::Variant;
using ReactionProduct = section::Type< 6 >::ReactionProduct;


template <typename Range>
auto search( Range xRange, double x, int i, int left, int right ){
  if ( xRange[i] <= x and x <= xRange[i+1] ){ return i; }
  if ( x > xRange[i] ){ left  = i; }
  else                { right = i; }
  i = (left + right)*0.5;
  return search(xRange,x,i,left,right);
}



double interpolate( std::vector<double> xVec, std::vector<double> yVec, 
                    const double& x ){
  int len = xVec.size();
  if (xVec.size() == 1 ){ return yVec[0];     }
  if ( x <= xVec[0]     ){ return yVec[0];     } 
  if ( x >= xVec[len-1] ){ return yVec[len-1]; }
  int index = search(xVec, x, int(len*0.5), 0, len);
  double b = yVec[index];
  double m = (yVec[index+1]-yVec[index])/(xVec[index+1]-xVec[index]);
  return m*(x-xVec[index])+b;
}









auto getIncohElasticDataSingleEnergy( const double& E, const double& dwf,
  std::vector<double>& equiprobCosines ){
  
  int nAngles = equiprobCosines.size()-2;
  equiprobCosines[0] = E;
  double lastCosine = -1.0, thisCosine, exp_2EW_cos, mu_i,
         exp_4EW = exp(-4.0*E*dwf),
         inv_2EW = 1.0/(2.0*E*dwf);

  for ( int i = 0; i < nAngles; ++i ){
    exp_2EW_cos = exp( -2.0*E*dwf*(1.0-lastCosine) );
    thisCosine  = 1.0+inv_2EW*log((1-exp_4EW)/nAngles + exp_2EW_cos );

    mu_i = nAngles * inv_2EW * 
           ( exp(-2.0*E*dwf*(1.0-thisCosine))*(2.0*E*dwf*thisCosine-1.0) - 
             exp_2EW_cos*(2.0*E*dwf*lastCosine-1.0) ) / (1.0-exp_4EW);
    equiprobCosines[i+2] = mu_i;
    lastCosine = thisCosine;
  }
}




auto incoherentElastic( ContinuumEnergyAngle law, int nbin, double debyeWallerFactor ){
      std::vector<double> esi;
      for (const auto& subsection : law.subsections()){
        esi.push_back(std::get<ThermalScatteringData>(subsection).energy());
      }

      std::vector<std::vector<double>> equiProbCosinesVec;
      std::vector<double> equiprobCosines(nbin+2,1.0);

      for ( const double& E : esi ){
        getIncohElasticDataSingleEnergy( E, debyeWallerFactor, equiprobCosines );
        equiProbCosinesVec.push_back(equiprobCosines);
      }

      int n2 = nbin+2;
      auto firstSCR = equiProbCosinesVec[0];
      ThermalScatteringData chunky( esi[0], n2, std::move(firstSCR) );
      std::vector<Variant> chunks {chunky};
      for ( size_t j = 1; j < esi.size(); ++j){
        auto scratch = equiProbCosinesVec[j];
        ThermalScatteringData chunky( esi[j], n2, std::move(scratch) );
        chunks.push_back(chunky);
      }


      return chunks;
}






/*
auto incoherentElastic( section::Type<7,2>::IncoherentElastic incoh_law,
  const double& temp, const double& emax, int nbin, section::Type<6> chunk ){
      double dwf = 0.0;

      std::vector<double> temperatures = incoh_law.temperatures();
      std::vector<double> debyeWaller  = incoh_law.debyeWallerValues();
        
      if ( temp < 0.9*temperatures[0] or 
           temp > 1.1*temperatures[temperatures.size()-1] ){
        std::cout << "oh no put an error here" << std::endl;
      }
      else {
        dwf = interpolate(temperatures, debyeWaller, temp);
      }

      //auto chunk = section6Vec[0];

      auto products = chunk.products();
      //auto products = section6Vec[0].products();
      auto law = std::get<ContinuumEnergyAngle>(products[0].distribution());
      auto subsections = law.subsections();

      std::vector<double> esi;
      for (const auto& subsection : subsections){
        esi.push_back(std::get<ThermalScatteringData>(subsection).energy());
      }
      std::vector<std::vector<double>> equiProbCosinesVec;
      std::vector<double> equiprobCosines(nbin);

      for ( const double& E : esi ){
        equiprobCosines.resize(nbin);
        getIncohElasticDataSingleEnergy( E, dwf, equiprobCosines );
        equiprobCosines.insert(equiprobCosines.begin(), 1);
        equiprobCosines.insert(equiprobCosines.begin(), E);
        equiProbCosinesVec.push_back(equiprobCosines);
      }

      int n2 = nbin+2;
      auto firstSCR = equiProbCosinesVec[0];
      ThermalScatteringData chunky( esi[0], n2, std::move(firstSCR) );
      std::vector<Variant> chunks {chunky};
      for ( size_t j = 1; j < esi.size(); ++j){
        auto scratch = equiProbCosinesVec[j];
        ThermalScatteringData chunky( esi[j], n2, std::move(scratch) );
        chunks.push_back(chunky);
      }

      long lep = 1;
      ContinuumEnergyAngle continuumChunk( lep, 
        {(long) esi.size()}, {2}, std::move(chunks) );

      //int jp = 0, lct = 1;
      std::vector<ReactionProduct> iel_products = 
        {ReactionProduct({ 1., 1, -1, 1, {2}, {2}, { 1.e-5, emax }, 
                         { 1., 1. }}, continuumChunk )};

      return iel_products;
}

*/

