#define CATCH_CONFIG_MAIN
#include "catch.hpp"
#include "thermr.cpp"
#include <range/v3/all.hpp>
#include "generalTools/testing.h"
#include <typeinfo>

using namespace njoy::ENDFtk;
using Tabulated = section::Type<7,4>::Tabulated;
using ContinuumEnergyAngle = section::Type< 6 >::ContinuumEnergyAngle;
using LaboratoryAngleEnergy = section::Type< 6 >::LaboratoryAngleEnergy;
using ThermalScatteringData = section::Type< 6 >::ContinuumEnergyAngle::ThermalScatteringData;


void checkCohElastic( section::Type<6> my_elastic_section, 
                      section::Type<6> good_elastic_section ){
  REQUIRE( my_elastic_section.ZA() == good_elastic_section.ZA() );
  REQUIRE( my_elastic_section.AWR() == good_elastic_section.AWR() );
  REQUIRE( my_elastic_section.MT() == good_elastic_section.MT() );

  auto good_products = good_elastic_section.products();
  auto good_law      = std::get<Unknown>(good_products[0].distribution());
  auto my_products = my_elastic_section.products();
  auto my_law = std::get<Unknown>(my_products[0].distribution() );
  REQUIRE( my_law.LAW() == good_law.LAW() );
  auto good_multiplicity = good_products[0].multiplicity();
  auto   my_multiplicity =   my_products[0].multiplicity();

  checkVec(good_multiplicity.energies(),my_multiplicity.energies());
  checkVec(good_multiplicity.multiplicities(), my_multiplicity.multiplicities() );
  checkVec(good_multiplicity.boundaries(), my_multiplicity.boundaries() );
  checkVec(good_multiplicity.interpolants(), my_multiplicity.interpolants() );
}




void check_E_mu_Ep( LaboratoryAngleEnergy good_law, LaboratoryAngleEnergy my_law ){
    auto good_energies = good_law.angularDistributions();
    auto   my_energies =   my_law.angularDistributions();

    REQUIRE( good_energies[0].energy() == Approx(my_energies[0].energy()) );
    auto   my_cosines = my_energies[0].energyDistributions();
    auto good_cosines = good_energies[0].energyDistributions();

    REQUIRE( good_energies.size() == my_energies.size() );

    for (size_t ienergy = 0; ienergy < good_energies.size(); ++ienergy) {

      checkVec( good_energies[0].interpolants(), my_energies[0].interpolants() );
      checkVec( good_energies[0].boundaries(), my_energies[0].boundaries() );

      REQUIRE( good_cosines.size() == my_cosines.size() );
      for (size_t icos = 0; icos < good_cosines.size(); ++icos){
        REQUIRE( good_cosines[icos].cosine() == Approx(my_cosines[icos].cosine()) );
        checkVec( good_cosines[icos].energies(),       my_cosines[icos].energies() );
        checkVec( good_cosines[icos].probabilities(),  my_cosines[icos].probabilities() );
        checkVec( good_cosines[icos].interpolants(),   my_cosines[icos].interpolants() );
        checkVec( good_cosines[icos].boundaries(),     my_cosines[icos].boundaries() );
      }
    }
}







template <typename Section>
void checkContinuumEnergyAngle(Section good_inelastic, Section my_section){
  auto good_products = good_inelastic.products();
  auto good_law      = std::get<ContinuumEnergyAngle>(good_products[0].distribution());

  REQUIRE( my_section.ZA() == good_inelastic.ZA() );
  REQUIRE( my_section.AWR() == good_inelastic.AWR() );
  REQUIRE( my_section.MT() == good_inelastic.MT() );

  auto my_products = my_section.products();
  auto my_law = std::get< ContinuumEnergyAngle >(my_products[0].distribution() );
  REQUIRE( my_law.LAW() == good_law.LAW() );
  REQUIRE( my_law.LEP() == good_law.LEP() );
  REQUIRE( my_law.NE() == good_law.NE() ); // number inc. E points on egrid
  REQUIRE( my_law.NR() == good_law.NR() );

  auto good_multiplicity = good_products[0].multiplicity();
  auto   my_multiplicity =   my_products[0].multiplicity();

  checkVec(good_multiplicity.energies(), my_multiplicity.energies() );
  checkVec(good_multiplicity.multiplicities(), my_multiplicity.multiplicities() );
  checkVec(good_multiplicity.boundaries(), my_multiplicity.boundaries() );
  checkVec(good_multiplicity.interpolants(), my_multiplicity.interpolants() );

  auto good_energies = good_law.subsections();
  auto   my_energies =   my_law.subsections();
  REQUIRE( good_energies.size() == my_energies.size() );
  
  for ( size_t i = 0; i < my_energies.size(); ++i ){ // Incoming energies
    auto my_entry   =   my_energies[i];
    auto good_entry = good_energies[i];
    auto good_subsection = std::get<ThermalScatteringData>(good_entry);
    auto   my_subsection = std::get<ThermalScatteringData>(  my_entry);

    REQUIRE( good_subsection.LANG()   == my_subsection.LANG() );
    REQUIRE( good_subsection.LTT()    == my_subsection.LTT() );
    REQUIRE( good_subsection.energy() == Approx(my_subsection.energy()).epsilon(1e-6) );
    REQUIRE( good_subsection.NW()     == my_subsection.NW() );
    REQUIRE( good_subsection.N2()     == my_subsection.N2() );

    auto good_data = good_subsection.data();
    auto   my_data =   my_subsection.data();

    checkVec(good_subsection.energies(),my_subsection.energies());
    checkVec(good_subsection.data(),    my_subsection.data());
    checkVec(good_subsection.PP(),my_subsection.PP());

    REQUIRE( good_subsection.cosines().size() == my_subsection.cosines().size() );
    for (size_t i = 0; i < good_subsection.cosines().size(); ++i){
      checkVec(good_subsection.cosines()[i],my_subsection.cosines()[i]);
    }
  }
}



TEST_CASE( "thermr" ){
    /*
  GIVEN( "NJOY Test 9 - H in H2O Example" ){
    
    njoy::njoy21::lipservice::iRecordStream<char> iss( std::istringstream(
      "24 23 25/\n"
      "101 125 4 1 2 0 0 1 222 2/\n"
      "296./\n"
      ".5 .625/" ) );
    njoy::njoy21::lipservice::THERMR thermr(iss);
    nlohmann::json jsonTHERMR(thermr);

    njoy::ENDFtk::syntaxTree::Tape<std::string> 
      pendfTape(njoy::utility::slurpFileToMemory("h2o_tape23")),
      leaprTape(njoy::utility::slurpFileToMemory("h2o_tape24")),
      rightTape(njoy::utility::slurpFileToMemory("h2o_tape25"));


    auto out = finalTHERMR( jsonTHERMR, leaprTape, pendfTape);
   
    auto my_section = out[0].MT(int(jsonTHERMR["mtref"])); 

    njoy::ENDFtk::file::Type<6> MF6 = 
        rightTape.materialNumber(125).front().fileNumber(6).parse<6>();
    auto good_inelastic  = MF6.MT(222);

    checkContinuumEnergyAngle(good_inelastic, my_section);


  } // GIVEN 
  GIVEN( "ENDF-B/VIII.0 Be" ){
    
    njoy::njoy21::lipservice::iRecordStream<char> iss( std::istringstream(
     " 24 23 25/\n"
     " 26 425 4 1 2 1 0 1 222 2/\n"
     " 296./\n"
     " .5 .2/" ));
    njoy::njoy21::lipservice::THERMR thermr(iss);
    nlohmann::json jsonTHERMR(thermr);

    njoy::ENDFtk::syntaxTree::Tape<std::string> 
      pendfTape(njoy::utility::slurpFileToMemory("be_tape23")),
      leaprTape(njoy::utility::slurpFileToMemory("be_tape24")),
      rightTape(njoy::utility::slurpFileToMemory("be_tape25"));


    auto out = finalTHERMR( jsonTHERMR, leaprTape, pendfTape);
   
    auto my_section = out[0].MT(int(jsonTHERMR["mtref"])); 

    njoy::ENDFtk::file::Type<6> MF6 = 
        rightTape.materialNumber(425).front().fileNumber(6).parse<6>();
    auto good_inelastic  = MF6.MT(222);

    checkContinuumEnergyAngle(good_inelastic, my_section);


  } // GIVEN 


  GIVEN( "NJOY Test 9 - H in H2O Example (Multiple Temps)" ){
    
    njoy::njoy21::lipservice::iRecordStream<char> iss( std::istringstream(
     "24 23 25/\n"
     "101 125 4 2 2 0 0 1 222 2/\n"
     "296. 500.0/\n"
     "0.5 .05/" ));
    njoy::njoy21::lipservice::THERMR thermr(iss);
    nlohmann::json jsonTHERMR(thermr);

    njoy::ENDFtk::syntaxTree::Tape<std::string> 
      pendfTape(njoy::utility::slurpFileToMemory("h2oMultT_tape23")),
      leaprTape(njoy::utility::slurpFileToMemory("h2oMultT_tape24")),
      rightTape(njoy::utility::slurpFileToMemory("h2oMultT_tape25"));

    auto out = finalTHERMR( jsonTHERMR, leaprTape, pendfTape);
   
    auto my_section = out[0].MT(int(jsonTHERMR["mtref"])); 

    std::vector<njoy::ENDFtk::file::Type<6>> MF6_variousTemps;
    for (auto& material : rightTape.materialNumber(125)){
      MF6_variousTemps.push_back(material.fileNumber(6).parse<6>());
    }


    for (size_t t = 0; t < jsonTHERMR["tempr"].size(); ++t){
      auto   my_inelastic_section = out[t].MT(int(jsonTHERMR["mtref"])); 
      auto good_inelastic_section = MF6_variousTemps[t].MT(int(jsonTHERMR["mtref"]));
      checkContinuumEnergyAngle(good_inelastic_section, my_inelastic_section);
    }


  } // GIVEN 


  GIVEN( "NJOY Test 9 - H in H2O Example (E mu E')" ){
    njoy::njoy21::lipservice::iRecordStream<char> iss( std::istringstream(
     "24 23 25/\n"
     "101 125 4 1 2 0 1 1 222 2/\n"
     "296./\n"
     "0.5 .625/" ));
    njoy::njoy21::lipservice::THERMR thermr(iss);
    nlohmann::json jsonTHERMR(thermr);

    njoy::ENDFtk::syntaxTree::Tape<std::string> 
      pendfTape(njoy::utility::slurpFileToMemory("h2oEmuEp_tape23")),
      leaprTape(njoy::utility::slurpFileToMemory("h2oEmuEp_tape24")),
      rightTape(njoy::utility::slurpFileToMemory("h2oEmuEp_tape25"));

    auto out = finalTHERMR( jsonTHERMR, leaprTape, pendfTape);
   
    auto my_section = out[0].MT(int(jsonTHERMR["mtref"])); 

    njoy::ENDFtk::file::Type<6> MF6 = 
        rightTape.materialNumber(125).front().fileNumber(6).parse<6>();
    auto good_inelastic  = MF6.MT(222);



    auto good_products = good_inelastic.products();
    auto good_law      = std::get<LaboratoryAngleEnergy>(good_products[0].distribution());

    auto my_products = my_section.products();
    auto my_law = std::get<LaboratoryAngleEnergy>(my_products[0].distribution() );

    REQUIRE( my_section.AWR() == good_inelastic.AWR() );
    REQUIRE( my_section.MT() == good_inelastic.MT() );

    check_E_mu_Ep( good_law, my_law );


  } // GIVEN 
  */


  GIVEN( "H in ZrH" ){

    WHEN( "1 temperature is provided" ){
    
      njoy::njoy21::lipservice::iRecordStream<char> iss( std::istringstream(
        "24 23 25/\n"
        "7 125 4 1 2 12 0 1 222 2/\n"
        "296./\n"
        ".5 .625/" ));
      njoy::njoy21::lipservice::THERMR thermr(iss);
      nlohmann::json jsonTHERMR(thermr);

      njoy::ENDFtk::syntaxTree::Tape<std::string> 
        pendfTape(njoy::utility::slurpFileToMemory("zrh_tape23")),
        leaprTape(njoy::utility::slurpFileToMemory("zrh_tape24")),
        rightTape(njoy::utility::slurpFileToMemory("zrh_tape25"));


      auto out = finalTHERMR( jsonTHERMR, leaprTape, pendfTape);
  
      auto   my_elastic_section = out[0].MT(int(jsonTHERMR["mtref"])+1); 
      auto   my_inelastic_section = out[0].MT(int(jsonTHERMR["mtref"])); 
  
      njoy::ENDFtk::file::Type<6> MF6 = 
          rightTape.materialNumber(int(jsonTHERMR["matdp"])).front().fileNumber(6).parse<6>();
      auto good_elastic_section = MF6.MT(int(jsonTHERMR["mtref"])+1);
      auto good_inelastic_section = MF6.MT(int(jsonTHERMR["mtref"]));
  
  
      checkContinuumEnergyAngle(  good_elastic_section,   my_elastic_section);
      checkContinuumEnergyAngle(good_inelastic_section, my_inelastic_section);
  

    } // WHEN

    WHEN( "multiple temperatures are provided" ){
    
  //    njoy::njoy21::lipservice::iRecordStream<char> iss( std::istringstream(
  //      "24 23 25/\n"
  //      "7 125 4 1 2 12 0 1 222 2/\n"
  //      "296./\n"
  //      ".5 .625/" ));
  //    njoy::njoy21::lipservice::THERMR thermr(iss);
  //    nlohmann::json jsonTHERMR(thermr);
//
//      njoy::ENDFtk::syntaxTree::Tape<std::string> 
//        pendfTape(njoy::utility::slurpFileToMemory("zrh_tape23")),
 //       leaprTape(njoy::utility::slurpFileToMemory("zrh_tape24")),
  //      rightTape(njoy::utility::slurpFileToMemory("zrh_tape25"));
//
//
 //     auto out = finalTHERMR( jsonTHERMR, leaprTape, pendfTape);
 // 
  //    auto   my_elastic_section = out[0].MT(int(jsonTHERMR["mtref"])+1); 
   //   auto   my_inelastic_section = out[0].MT(int(jsonTHERMR["mtref"])); 
 // 
  //    njoy::ENDFtk::file::Type<6> MF6 = 
//          rightTape.materialNumber(int(jsonTHERMR["matdp"])).front().fileNumber(6).parse<6>();
 //     auto good_elastic_section = MF6.MT(int(jsonTHERMR["mtref"])+1);
  //    auto good_inelastic_section = MF6.MT(int(jsonTHERMR["mtref"]));
 // 
  
  //    checkContinuumEnergyAngle(  good_elastic_section,   my_elastic_section);
   //   checkContinuumEnergyAngle(good_inelastic_section, my_inelastic_section);
  

    } // WHEN



  } // GIVEN 


} // TEST CASE

