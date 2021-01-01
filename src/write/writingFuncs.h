


section::Type<6> prepareMF6_inelastic( const std::vector<double> incidentEnergies,
  std::vector<std::vector<double>> totalSCR, int nbin, double emax, int mtref, 
  double za, double pendf_awr){

  std::vector<Variant> chunks;
  for ( size_t j = 0; j < incidentEnergies.size(); ++j){
    std::vector<double> scratch = totalSCR[j];
    ThermalScatteringData chunk(incidentEnergies[j], nbin+2, std::move(scratch));
    chunks.push_back(chunk);
  }

  ContinuumEnergyAngle continuumChunk(1, {(long) incidentEnergies.size()}, 
                                     {2}, std::move( chunks ) );
 
  std::vector<ReactionProduct> products = 
    {ReactionProduct({ 1., 1, -1, 1, {2}, {2}, { 1.e-5, emax }, 
                     { 1., 1. }}, continuumChunk )};
  section::Type<6> inelastic(mtref, za, pendf_awr, 0, 1, std::move(products));
  return inelastic;
}




void prepareMF3_inelastic( std::vector<double>& MF3_energies, 
  std::vector<double>& MF3_XS, section::Type<3> MF3_2, double emax,
  std::vector<double> initialEnergies, std::vector<double> xsi, int iinc ) {

  std::vector<double> desiredEnergies; 
  std::vector<double> finalXS;

  if (iinc == 1){
    initialEnergies = MF3_2.energies();
    xsi             = MF3_2.crossSections();
  }

  for (auto energy : MF3_2.energies()){
    if (energy >= emax or std::fabs((energy-emax)/emax) < 1e-10){ 
        if (energy > emax*(1+1e-10)){ 
            energy = emax; 
        }
        desiredEnergies.emplace_back(energy); 
        if (std::fabs(energy-emax) < 1e-10*emax){ 
          desiredEnergies.emplace_back(1.00001*energy); 
        }
        break; 
    }
    desiredEnergies.emplace_back(energy);
  }
  for ( const auto& energy : desiredEnergies ){
    finalXS.emplace_back(terp(initialEnergies,xsi,energy,5));
  }
  desiredEnergies.emplace_back(2e7);
  finalXS.emplace_back(0.0);

  if (desiredEnergies[desiredEnergies.size()-2] > emax){ 
      finalXS[desiredEnergies.size()-2] = 0.0; }

  MF3_energies = desiredEnergies;
  MF3_XS       = finalXS;
}




section::Type<3> prepareMF3_cohElastic( std::vector<double> MF3CohElasticEnergies,
  std::vector<double> MF3CohElasticCrossSections, int mtref, int za, 
  double pendf_awr, double temp, std::vector<double>& MF3_energies, 
  std::vector<double>& MF3_XS ){

  std::vector<double> MF3_energies_New = MF3CohElasticEnergies;
  std::vector<double> MF3_XS_New;

  // Write out bragg peaks
  section::Type<3> MF3_cohEl( mtref+1, za, pendf_awr, temp, 0.0, 0,
                              std::vector<long>(1,long(MF3CohElasticEnergies.size())),
                              std::vector<long>(1,2),
                              std::move( MF3CohElasticEnergies), 
                              std::move( MF3CohElasticCrossSections) );


  // Write ensure that other cross sections are on same grid as bragg
  int numSec3Energies = MF3_energies.size();

  for (const auto& E : MF3_energies_New){
    int begin = 0, end = 0;
    for ( int j = 0; j < numSec3Energies; ++j){
      if (MF3_energies[j] >= E){ 
        if (j-3 < 0){
          begin = 0;
          end   = 5;
        }
        else if (j+3 >= numSec3Energies ){
          end   = numSec3Energies - 1;
          begin = end - 5;
        }
        else {
          begin = j-2;
          end   = j+3;
        }
        break; 
      }
    }
    int interpOrder = (begin == 0 or end == numSec3Energies - 1) ? 3 : 4,
        nl          = (begin == 0 or end == numSec3Energies - 1) ? 4 : 5;

    std::vector<double> 
      temporary1 (MF3_energies.begin()+begin,MF3_energies.begin()+end),
      temporary2 (MF3_XS.begin()      +begin,MF3_XS.begin()      +end);

    MF3_XS_New.push_back(terp(temporary1,temporary2,E,interpOrder,nl));
  }
   MF3_XS_New[MF3_XS_New.size()-2] = 0.0;
   MF3_XS_New[MF3_XS_New.size()-1] = 0.0;

   MF3_energies = MF3_energies_New;
   MF3_XS       = MF3_XS_New;

  return MF3_cohEl;
}






