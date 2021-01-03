#include "ENDFtk.hpp"

using ReactionProduct = section::Type< 6 >::ReactionProduct;
using CoherentElastic = section::Type<7,2>::CoherentElastic;
using Unknown = section::Type<6>::Unknown;


auto prepareMF6_cohElastic( int nendf, int matde, double emax, int mtref, 
  double za, double pendf_awr, std::vector<DirectoryRecord>& index, 
  std::vector<section::Type<6>>& section6Vec ){

  tree::Tape<std::string> leaprTape ( utility::slurpFileToMemory("tape" +
                                      std::to_string(nendf)));

  file::Type<7> MF7 = leaprTape.material(matde).front().file(7).parse<7>();

  auto cohElasticInfo = std::get<CoherentElastic>((MF7.MT(2_c)).scatteringLaw());

  // Write the MF6 coherent inelastic section
  std::vector<ReactionProduct> products {ReactionProduct(
    { 1., 1, -cohElasticInfo.numberBraggEdges(), 0, {2}, {2}, {1.e-5, emax}, 
    {1., 1.}}, Unknown() )};
 
  section::Type<6> cohElastic(mtref+1, za, pendf_awr, 0, 1, std::move(products));
  index.emplace_back( 6, mtref+1, cohElastic.NC(), 0 );
  section6Vec.push_back( std::move(cohElastic) );
  return cohElasticInfo;

}





double prepareMF6_incohElastic(int nendf, int matde, double temp,
  double emax, int nbin, int mtref, double za, double pendf_awr,
  std::vector<section::Type<6>>& section6Vec, std::ostream& error,
  std::vector<DirectoryRecord>& index ){

  tree::Tape<std::string> leaprTape ( utility::slurpFileToMemory("tape" +
                                      std::to_string(nendf)));

  file::Type<7> MF7 = leaprTape.material(matde).front().file(7).parse<7>();

  njoy::ENDFtk::section::Type<7,2> leapr_MT2 = MF7.MT(2_c);
  auto incoh_law = std::get<IncoherentElastic>(leapr_MT2.scatteringLaw());
  auto chunk = section6Vec[0];
  
  std::vector<double> temperatures = incoh_law.temperatures(),
                      debyeWaller  = incoh_law.debyeWallerValues();
          
  if ( temp < 0.9*temperatures[0] or 
       temp > 1.1*temperatures[temperatures.size()-1] ){
    error << "Cannot obtain Debye-Waller factor " <<
             "- Temperature out of interpolation range" << std::endl;
  }

  double debyeWallerFactor = interpolate(temperatures, debyeWaller, temp);

  auto law = std::get<ContinuumEnergyAngle>( 
                               chunk.reactionProducts()[0].distribution());

  auto chunks = incoherentElastic( law, nbin, debyeWallerFactor );

  ContinuumEnergyAngle continuumChunk( long(1), 
    {(long) chunks.size()}, {2}, std::move(chunks) );

  std::vector<ReactionProduct> iel_products = 
                   { ReactionProduct({ 1., 1, -1, 1, {2}, {2}, { 1.e-5, emax }, 
                                     { 1., 1. }}, continuumChunk ) };

  section::Type<6> incElastic(mtref+1, za, pendf_awr, 0, 1, 
                              std::move(iel_products));

  index.emplace_back( 6, mtref+1, incElastic.NC(), 0 );
  section6Vec.push_back( std::move(incElastic) );

  return debyeWallerFactor;

}


template <typename Tuple>
void prepareMF6_E_Ep_Mu( const Tuple& E_Ep_Mu_output, int nbin, double emax, 
  int mtref, double za, double pendf_awr, std::vector<DirectoryRecord>& index, 
  std::vector<section::Type<6>>& section6Vec ){

  std::vector<double>              incidentEnergies = std::get<0>(E_Ep_Mu_output);
  std::vector<std::vector<double>> totalSCR         = std::get<1>(E_Ep_Mu_output);


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

  index.emplace_back( 6, mtref, inelastic.NC(), 0 );

  section6Vec.push_back( std::move(inelastic) );
}



void prepareMF6_E_Mu_Ep( double emax, LaboratoryAngleEnergy labAngleEnergy,
  int mtref, double za, double pendf_awr, std::vector<DirectoryRecord>& index,
  std::vector<section::Type<6>>& section6Vec ){

  std::vector<ReactionProduct> products = { ReactionProduct( 
    { 1., 1, -1, 7, {2}, {2}, { 1.e-5, emax }, { 1., 1. }}, 
    std::move(labAngleEnergy))};

  section::Type<6> inelastic(mtref, za, pendf_awr, 0, 1, std::move(products));
  index.emplace_back( 6, mtref, inelastic.NC(), 0 );
  section6Vec.push_back( std::move(inelastic) );

}



void prepareMF3_inelastic( std::vector<double>& MF3_energies, 
  std::vector<double>& MF3_XS, section::Type<3> MF3_2, double emax,
  std::vector<double> initialEnergies, std::vector<double> xsi, int iinc,
  int mtref, double za, double pendf_awr, double temp ) {

  std::vector<double> desiredEnergies, finalXS;

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

  //return MF3;

}



section::Type<3> prepareMF3_incohElastic( std::vector<double> inelasticEnergies,
  std::vector<double> MF3_energies, double boundXS, double debyeWallerFactor,
  int mtref, int natom, double za, double pendf_awr, double temp ){

  double c1 = boundXS/(2*natom);
  std::vector<double> MF3_XS_oldGrid;
  for (auto& energy : inelasticEnergies){
    double c2 = debyeWallerFactor*energy*2;
    MF3_XS_oldGrid.push_back(c1/c2*(1.0-exp(-2*c2)));
  }

  std::vector<double> MF3_XS;
  for (auto& energy : MF3_energies){
    MF3_XS.push_back(terp(inelasticEnergies,MF3_XS_oldGrid,energy,3,
                          inelasticEnergies.size()));
  } 





  MF3_XS[MF3_XS.size()-2] = 0.0;
  MF3_XS[MF3_XS.size()-1] = 0.0;
  section::Type<3> MF3_incohEl( mtref+1, za, pendf_awr, temp, 0.0, 0,
                              std::vector<long>(1,long(MF3_energies.size())),
                              std::vector<long>(1,2),
                              std::move( MF3_energies), 
                              std::move( MF3_XS) );

  return MF3_incohEl;

}



section::Type<3> prepareMF3_cohElastic( std::tuple<std::vector<double>,
  std::vector<double>,std::vector<double>> cohElasticOutput, int mtref, int za, 
  double pendf_awr, double temp, std::vector<double>& MF3_energies, 
  std::vector<double>& MF3_XS ){

  std::vector<double> MF3CohElasticEnergies = std::get<0>(cohElasticOutput),
                      MF3CohElasticXS       = std::get<1>(cohElasticOutput);
  std::vector<double> MF3_energies_New = MF3CohElasticEnergies,
                      MF3_XS_New;

  // Write out bragg peaks
  section::Type<3> MF3_cohEl( mtref+1, za, pendf_awr, temp, 0.0, 0,
                              std::vector<long>(1,long(MF3CohElasticEnergies.size())),
                              std::vector<long>(1,2),
                              std::move( MF3CohElasticEnergies), 
                              std::move( MF3CohElasticXS ) );

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


void writeTape( std::vector<Material> materials, tree::Tape<std::string> pendf,
                int nout ){

  Tape tape( TapeIdentification(std::string(pendf.TPID().text()),1), 
             std::move(materials) );
  std::string buffer;
  auto materialOutput = std::back_inserter( buffer );
  tape.print( materialOutput );
  std::string name = "tape"+std::to_string(nout);
  std::ofstream out(name);
  out << buffer;
  out.close();
  }
  



