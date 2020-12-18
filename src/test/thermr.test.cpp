#define CATCH_CONFIG_MAIN
#include "catch.hpp"
#include "thermr.cpp"
#include <range/v3/all.hpp>
#include "generalTools/testing.h"
#include <typeinfo>
#include <stdio.h>

using namespace njoy::ENDFtk;
using Tabulated = section::Type<7,4>::TabulatedFunctions;
using ContinuumEnergyAngle = section::Type< 6 >::ContinuumEnergyAngle;
using LaboratoryAngleEnergy = section::Type< 6 >::LaboratoryAngleEnergy;
using ThermalScatteringData = section::Type< 6 >::ContinuumEnergyAngle::ThermalScatteringData;



void checkFile3( section::Type<3> good_MF3, section::Type<3> my_MF3 ){

  REQUIRE( good_MF3.MT() == my_MF3.MT() );
  REQUIRE( good_MF3.ZA() == my_MF3.ZA() );
  REQUIRE( good_MF3.AWR() == Approx(my_MF3.AWR()) );
  REQUIRE( good_MF3.LR() == my_MF3.LR() );
  REQUIRE( good_MF3.QM() == Approx(my_MF3.QM()) );
  REQUIRE( good_MF3.QI() == Approx(my_MF3.QI()) );
  REQUIRE( good_MF3.NP() == my_MF3.NP() );
  REQUIRE( good_MF3.NR() == my_MF3.NR() );
  checkVec( good_MF3.interpolants(), my_MF3.interpolants() );
  checkVec( good_MF3.boundaries(), my_MF3.boundaries() );
  checkVec( good_MF3.energies(), my_MF3.energies() );
  checkVec( good_MF3.crossSections(), my_MF3.crossSections() );
  REQUIRE( good_MF3.NC() == my_MF3.NC() );

}




void checkCohElastic( section::Type<6> my_elastic_section, 
                      section::Type<6> good_elastic_section ){
  REQUIRE( my_elastic_section.ZA() == good_elastic_section.ZA() );
  REQUIRE( my_elastic_section.AWR() == good_elastic_section.AWR() );
  REQUIRE( my_elastic_section.MT() == good_elastic_section.MT() );

  auto good_products = good_elastic_section.reactionProducts();
  auto good_law      = std::get<Unknown>(good_products[0].distribution());
  auto my_products = my_elastic_section.reactionProducts();
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

    REQUIRE( good_energies[0].incidentEnergy() == Approx(my_energies[0].incidentEnergy()) );
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
  auto good_products = good_inelastic.reactionProducts();
  auto good_law      = std::get<ContinuumEnergyAngle>(good_products[0].distribution());

  REQUIRE( my_section.ZA() == good_inelastic.ZA() );
  REQUIRE( my_section.AWR() == good_inelastic.AWR() );
  REQUIRE( my_section.MT() == good_inelastic.MT() );

  auto my_products = my_section.reactionProducts();
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

  auto good_energies = good_law.distributions();
  auto   my_energies =   my_law.distributions();
  REQUIRE( good_energies.size() == my_energies.size() );
  
  for ( size_t i = 0; i < my_energies.size(); ++i ){ // Incoming energies
    auto my_entry   =   my_energies[i];
    auto good_entry = good_energies[i];
    auto good_subsection = std::get<ThermalScatteringData>(good_entry);
    auto   my_subsection = std::get<ThermalScatteringData>(  my_entry);

    REQUIRE( good_subsection.LANG()   == my_subsection.LANG() );
    REQUIRE( good_subsection.LTT()    == my_subsection.LTT() );
    REQUIRE( good_subsection.incidentEnergy() == Approx(my_subsection.incidentEnergy()).epsilon(1e-6) );
    REQUIRE( good_subsection.NW()     == my_subsection.NW() );
    REQUIRE( good_subsection.N2()     == my_subsection.N2() );

    auto good_data = good_subsection.data();
    auto   my_data =   my_subsection.data();

    checkVec(good_subsection.energies(),my_subsection.energies());
    checkVec(good_subsection.data(),    my_subsection.data());
    checkVec(good_subsection.PP(),my_subsection.PP());

    REQUIRE( good_subsection.cosines().size() == my_subsection.cosines().size() );
    for (size_t j = 0; j < good_subsection.cosines().size(); ++j){
      checkVec(good_subsection.cosines()[j],my_subsection.cosines()[j]);
    }
  }
}





TEST_CASE( "thermr" ){
    /*
  */
  GIVEN( "NJOY Test 9 - H in H2O Example" ){
    
    njoy::njoy21::lipservice::iRecordStream<char> iss( std::istringstream(
      "24 23 25/\n"
      "101 125 4 1 2 0 0 1 222 2/\n"
      "296./\n"
      ".5 .625/" ) );
    njoy::njoy21::lipservice::THERMR thermr(iss);
    nlohmann::json jsonTHERMR(thermr);

    rename( "h2o_tape23" , "tape23");
    rename( "h2o_tape24" , "tape24");

    finalTHERMR( jsonTHERMR );//, std::cout, std::cerr );


    tree::Tape<std::string> goodTape(utility::slurpFileToMemory("h2o_tape25"));
    tree::Tape<std::string> myTape  (utility::slurpFileToMemory("tape25"    ));

    file::Type<6> good_MF6 = goodTape.material(125).front().file(6).parse<6>();
    file::Type<6> my_MF6   =   myTape.material(125).front().file(6).parse<6>();

    section::Type<6> good_inelastic  = good_MF6.MT(int(jsonTHERMR["mtref"]));
    section::Type<6> my_inelastic    =   my_MF6.MT(int(jsonTHERMR["mtref"]));

    checkContinuumEnergyAngle(good_inelastic, my_inelastic);

    file::Type<3> good_file3 = goodTape.material(125).front().file(3).parse<3>();
    file::Type<3> my_file3   =   myTape.material(125).front().file(3).parse<3>();
    section::Type<3> good_MF3 = good_file3.MT(int(jsonTHERMR["mtref"]));
    section::Type<3>   my_MF3 =   my_file3.MT(int(jsonTHERMR["mtref"]));

    checkFile3( good_MF3, my_MF3 );
    rename( "tape23" , "h2o_tape23");
    rename( "tape24" , "h2o_tape24");



  } // GIVEN 

  GIVEN( "ENDF-B/VIII.0 Be" ){
    
    njoy::njoy21::lipservice::iRecordStream<char> iss( std::istringstream(
     " 24 23 25/\n"
     " 26 425 4 1 2 2 0 1 222 2/\n"
     " 296./\n"
     " .5 .2/" ));
    njoy::njoy21::lipservice::THERMR thermr(iss);
    nlohmann::json jsonTHERMR(thermr);

    rename( "be_tape23" , "tape23");
    rename( "be_tape24" , "tape24");

    finalTHERMR( jsonTHERMR );//, std::cout, std::cerr );

    tree::Tape<std::string> goodTape(utility::slurpFileToMemory("be_tape25"));
    tree::Tape<std::string> myTape  (utility::slurpFileToMemory("tape25"    ));

    file::Type<6> good_MF6 = goodTape.material(425).front().file(6).parse<6>();
    file::Type<6> my_MF6   =   myTape.material(425).front().file(6).parse<6>();

    section::Type<6> good_inelastic  = good_MF6.MT(int(jsonTHERMR["mtref"]));
    section::Type<6> my_inelastic    =   my_MF6.MT(int(jsonTHERMR["mtref"]));

    checkContinuumEnergyAngle(good_inelastic, my_inelastic);

    file::Type<3> good_file3 = goodTape.material(425).front().file(3).parse<3>();
    file::Type<3> my_file3   =   myTape.material(425).front().file(3).parse<3>();
    section::Type<3> good_MF3 = good_file3.MT(int(jsonTHERMR["mtref"]));
    section::Type<3>   my_MF3 =   my_file3.MT(int(jsonTHERMR["mtref"]));

    section::Type<3> good_MF3Bragg = good_file3.MT(int(jsonTHERMR["mtref"])+1);
    section::Type<3>   my_MF3Bragg =   my_file3.MT(int(jsonTHERMR["mtref"])+1);

    checkFile3( good_MF3, my_MF3 );
    checkFile3( good_MF3Bragg, my_MF3Bragg );

    rename( "tape23" , "be_tape23");
    rename( "tape24" , "be_tape24");

  } // GIVEN 


  GIVEN( "NJOY Test 9 - H in H2O Example (Multiple Temps)" ){
    
    njoy::njoy21::lipservice::iRecordStream<char> iss( std::istringstream(
     "24 23 25/\n"
     "101 125 4 2 2 0 0 1 222 2/\n"
     "296. 500.0/\n"
     "0.5 .05/" ));
    njoy::njoy21::lipservice::THERMR thermr(iss);
    nlohmann::json jsonTHERMR(thermr);

    rename( "h2oMultT_tape23" , "tape23");
    rename( "h2oMultT_tape24" , "tape24");

    finalTHERMR( jsonTHERMR );//, std::cout, std::cerr );

    tree::Tape<std::string> goodTape(utility::slurpFileToMemory("h2oMultT_tape25"));
    tree::Tape<std::string> myTape  (utility::slurpFileToMemory("tape25"    ));

    file::Type<6> good_MF6 = goodTape.material(125).front().file(6).parse<6>();
    file::Type<6> my_MF6   =   myTape.material(125).front().file(6).parse<6>();

    section::Type<6> good_inelastic  = good_MF6.MT(int(jsonTHERMR["mtref"]));
    section::Type<6> my_inelastic    =   my_MF6.MT(int(jsonTHERMR["mtref"]));

    checkContinuumEnergyAngle(good_inelastic, my_inelastic);

    file::Type<3> good_file3 = goodTape.material(125).front().file(3).parse<3>();
    file::Type<3> my_file3   =   myTape.material(125).front().file(3).parse<3>();
    section::Type<3> good_MF3 = good_file3.MT(int(jsonTHERMR["mtref"]));
    section::Type<3>   my_MF3 =   my_file3.MT(int(jsonTHERMR["mtref"]));

    checkFile3( good_MF3, my_MF3 );

    rename( "tape23" , "h2oMultT_tape23");
    rename( "tape24" , "h2oMultT_tape24");



  } // GIVEN 


  GIVEN( "NJOY Test 9 - H in H2O Example (E mu E')" ){
    njoy::njoy21::lipservice::iRecordStream<char> iss( std::istringstream(
     "24 23 25/\n"
     "101 125 4 1 2 0 1 1 222 2/\n"
     "296./\n"
     "0.5 .625/" ));
    njoy::njoy21::lipservice::THERMR thermr(iss);
    nlohmann::json jsonTHERMR(thermr);


    rename( "h2oEmuEp_tape23" , "tape23");
    rename( "h2oEmuEp_tape24" , "tape24");

    finalTHERMR( jsonTHERMR );//, std::cout, std::cerr );
    std::cout << "finished " << std::endl;


    tree::Tape<std::string> goodTape(utility::slurpFileToMemory("h2oEmuEp_tape25"));
    tree::Tape<std::string> myTape  (utility::slurpFileToMemory("tape25"    ));

    file::Type<6> good_MF6 = goodTape.material(125).front().file(6).parse<6>();
    file::Type<6> my_MF6   =   myTape.material(125).front().file(6).parse<6>();


    section::Type<6> good_inelastic  = good_MF6.MT(int(jsonTHERMR["mtref"]));
    section::Type<6> my_inelastic    =   my_MF6.MT(int(jsonTHERMR["mtref"]));

    auto good_products = good_inelastic.reactionProducts();
    auto good_law      = std::get<LaboratoryAngleEnergy>(good_products[0].distribution());

    auto my_products = my_inelastic.reactionProducts();
    auto my_law = std::get<LaboratoryAngleEnergy>(my_products[0].distribution() );

    REQUIRE( my_inelastic.AWR() == good_inelastic.AWR() );
    REQUIRE( my_inelastic.MT() == good_inelastic.MT() );

    check_E_mu_Ep( good_law, my_law );

    file::Type<3> good_file3 = goodTape.material(125).front().file(3).parse<3>();
    file::Type<3> my_file3   =   myTape.material(125).front().file(3).parse<3>();
    section::Type<3> good_MF3 = good_file3.MT(int(jsonTHERMR["mtref"]));
    section::Type<3>   my_MF3 =   my_file3.MT(int(jsonTHERMR["mtref"]));

    checkFile3( good_MF3, my_MF3 );


    rename( "tape23" , "h2oEmuEp_tape23");
    rename( "tape24" , "h2oEmuEp_tape24");


  } // GIVEN 


  GIVEN( "H in ZrH" ){

    WHEN( "1 temperature is provided" ){
    
      njoy::njoy21::lipservice::iRecordStream<char> iss( std::istringstream(
        "24 23 25/\n"
        "7 125 4 1 2 12 0 1 222 2/\n"
        "296./\n"
        ".5 .625/" ));
      njoy::njoy21::lipservice::THERMR thermr(iss);
      nlohmann::json jsonTHERMR(thermr);

      rename( "zrh_tape23" , "tape23");
      rename( "zrh_tape24" , "tape24");

      finalTHERMR( jsonTHERMR );//, std::cout, std::cerr );

      tree::Tape<std::string> goodTape(utility::slurpFileToMemory("zrh_tape25"));
      tree::Tape<std::string> myTape  (utility::slurpFileToMemory("tape25"    ));

      file::Type<6> good_MF6 = goodTape.material(125).front().file(6).parse<6>();
      file::Type<6> my_MF6   =   myTape.material(125).front().file(6).parse<6>();

      section::Type<6> good_inelastic  = good_MF6.MT(int(jsonTHERMR["mtref"]));
      section::Type<6> my_inelastic    =   my_MF6.MT(int(jsonTHERMR["mtref"]));

      checkContinuumEnergyAngle(good_inelastic, my_inelastic);

      file::Type<3> good_file3 = goodTape.material(125).front().file(3).parse<3>();
      file::Type<3> my_file3   =   myTape.material(125).front().file(3).parse<3>();
      section::Type<3> good_MF3 = good_file3.MT(int(jsonTHERMR["mtref"]));
      section::Type<3>   my_MF3 =   my_file3.MT(int(jsonTHERMR["mtref"]));

      checkFile3( good_MF3, my_MF3 );

      rename( "tape23" , "zrh_tape23");
      rename( "tape24" , "zrh_tape24");

    } // WHEN

    WHEN( "multiple temperatures are provided" ){
    
      njoy::njoy21::lipservice::iRecordStream<char> iss( std::istringstream(
        "24 23 25/\n"
        "7 125 4 3 2 12 0 1 222 2/\n"
        "300. 900. 1200./\n"
        "0.5 .25/" ));
      njoy::njoy21::lipservice::THERMR thermr(iss);
      nlohmann::json jsonTHERMR(thermr);


      rename( "zrhMultT_tape23" , "tape23");
      rename( "zrhMultT_tape24" , "tape24");

      finalTHERMR( jsonTHERMR );//, std::cout, std::cerr );

      tree::Tape<std::string> goodTape(utility::slurpFileToMemory("zrhMultT_tape25"));
      tree::Tape<std::string> myTape  (utility::slurpFileToMemory("tape25"    ));

      file::Type<6> good_MF6 = goodTape.material(125).front().file(6).parse<6>();
      file::Type<6> my_MF6   =   myTape.material(125).front().file(6).parse<6>();

      section::Type<6> good_inelastic  = good_MF6.MT(int(jsonTHERMR["mtref"]));
      section::Type<6> my_inelastic    =   my_MF6.MT(int(jsonTHERMR["mtref"]));

      checkContinuumEnergyAngle(good_inelastic, my_inelastic);

      file::Type<3> good_file3 = goodTape.material(125).front().file(3).parse<3>();
      file::Type<3> my_file3   =   myTape.material(125).front().file(3).parse<3>();
      section::Type<3> good_MF3 = good_file3.MT(int(jsonTHERMR["mtref"]));
      section::Type<3>   my_MF3 =   my_file3.MT(int(jsonTHERMR["mtref"]));

      checkFile3( good_MF3, my_MF3 );

      rename( "tape23" , "zrhMultT_tape23");
      rename( "tape24" , "zrhMultT_tape24");

        
    } // WHEN

  } // GIVEN 


} // TEST CASE

