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
  //REQUIRE( my_elastic_section.ZA() == good_elastic_section.ZA() );
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
void checkInelastic(Section good_inelastic, Section my_section){
  auto good_products = good_inelastic.products();
  auto good_law      = std::get<ContinuumEnergyAngle>(good_products[0].distribution());

  //REQUIRE( my_section.ZA() == good_inelastic.ZA() );
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
  GIVEN( "NJOY Test 9 - H in H2O Example" ){
    std::string endfFile = njoy::utility::slurpFileToMemory("h2o_tape25");
    njoy::ENDFtk::syntaxTree::Tape<std::string > tape(endfFile);


    njoy::ENDFtk::file::Type<6> MF6 = tape.materialNumber(1301).front().fileNumber(6).parse<6>();
    auto good_inelastic  = MF6.MT(222);
    std::string leaprOut = njoy::utility::slurpFileToMemory("h2o_tape24");
    njoy::ENDFtk::syntaxTree::Tape<std::string > leaprTape(leaprOut);

    int matde = 101, matdp = 1301, nbin  = 4, iinc  = 2, icoh  = 0, 
        iform = 0,   natom = 1,    mtref = 222;
 
    std::vector<double> temps {296.0};
    double  tol = 0.5, emax = 0.625; 

    auto out = thermr(matde, matdp, nbin, iinc, icoh, iform, natom, mtref, temps,
           tol, emax, leaprTape);
   
    auto my_section = out[0].MT(mtref); 
    checkInelastic(good_inelastic, my_section);

  } // GIVEN

  GIVEN( "ENDF-B/VIII.0 Be" ){
    std::string endfFile = njoy::utility::slurpFileToMemory("be_tape25");
    njoy::ENDFtk::syntaxTree::Tape<std::string > tape(endfFile);


    njoy::ENDFtk::file::Type<6> MF6 = tape.materialNumber(425).front().fileNumber(6).parse<6>();
    auto good_inelastic_section = MF6.MT(222);
    auto good_elastic_section   = MF6.MT(223);
    std::string leaprOut = njoy::utility::slurpFileToMemory("be_tape24");
    njoy::ENDFtk::syntaxTree::Tape<std::string > leaprTape(leaprOut);

    int matde = 26, matdp = 425, nbin  = 4, iinc  = 2, icoh  = 2, 
        iform = 0,   natom = 1,  mtref = 222;
 
    std::vector<double> temps {296.0};
    double  tol = 0.5, emax = 0.2; 

    auto out = thermr(matde, matdp, nbin, iinc, icoh, iform, natom, mtref, temps,
           tol, emax, leaprTape);
    auto my_inelastic_section = out[0].MT(mtref); 
    checkInelastic(good_inelastic_section, my_inelastic_section);

    auto my_elastic_section   = out[1].MT(mtref+1); 

    checkCohElastic( my_elastic_section, good_elastic_section );


  } // GIVEN


  GIVEN( "NJOY Test 9 - H in H2O Example (500 K)" ){

    std::string leaprOut = njoy::utility::slurpFileToMemory("h2o500K_tape24");
    njoy::ENDFtk::syntaxTree::Tape<std::string > leaprTape(leaprOut);

    std::string endfFile = njoy::utility::slurpFileToMemory("h2o500K_tape25");
    njoy::ENDFtk::syntaxTree::Tape<std::string > tape(endfFile);

    std::vector<njoy::ENDFtk::file::Type<6>> MF6_variousTemps;
    for (auto& material : tape.materialNumber(1301)){
      MF6_variousTemps.push_back(material.fileNumber(6).parse<6>());
    }

    int matde = 101, matdp = 1301, nbin  = 4, iinc  = 2, icoh  = 0, 
        iform = 0,   natom = 1,    mtref = 222;
 
    std::vector<double> temps {500.0};
    double  tol = 0.5, emax = 0.05; 

    auto out = thermr(matde, matdp, nbin, iinc, icoh, iform, natom, mtref, temps,
           tol, emax, leaprTape);
   
    auto   my_inelastic_section = out[0].MT(mtref); 
    auto good_inelastic_section = MF6_variousTemps[0].MT(mtref);
    checkInelastic(good_inelastic_section, my_inelastic_section);

  } // GIVEN





  GIVEN( "NJOY Test 9 - H in H2O Example (Multiple Temps)" ){

    std::string leaprOut = njoy::utility::slurpFileToMemory("h2oMultT_tape24");
    njoy::ENDFtk::syntaxTree::Tape<std::string > leaprTape(leaprOut);

    std::string endfFile = njoy::utility::slurpFileToMemory("h2oMultT_tape25");
    njoy::ENDFtk::syntaxTree::Tape<std::string > tape(endfFile);

    std::vector<njoy::ENDFtk::file::Type<6>> MF6_variousTemps;
    for (auto& material : tape.materialNumber(1301)){
      MF6_variousTemps.push_back(material.fileNumber(6).parse<6>());
    }

    int matde = 101, matdp = 1301, nbin  = 4, iinc  = 2, icoh  = 0, 
        iform = 0,   natom = 1,    mtref = 222;
 
    std::vector<double> temps {296.0,500.0};
    double  tol = 0.5, emax = 0.05; 

    auto out = thermr(matde, matdp, nbin, iinc, icoh, iform, natom, mtref, temps,
           tol, emax, leaprTape);
   
    for (size_t t = 0; t < temps.size(); ++t){
      auto   my_inelastic_section = out[t].MT(mtref); 
      auto good_inelastic_section = MF6_variousTemps[t].MT(mtref);
      checkInelastic(good_inelastic_section, my_inelastic_section);
    }
  } // GIVEN

  GIVEN( "NJOY Test 9 - H in H2O Example (E mu E')" ){
    std::string endfFile = njoy::utility::slurpFileToMemory("h2oEmuEp_tape25");
    njoy::ENDFtk::syntaxTree::Tape<std::string > tape(endfFile);


    njoy::ENDFtk::file::Type<6> MF6 = tape.materialNumber(1301).front().fileNumber(6).parse<6>();
    auto good_inelastic  = MF6.MT(222);
    std::string leaprOut = njoy::utility::slurpFileToMemory("h2oEmuEp_tape24");
    njoy::ENDFtk::syntaxTree::Tape<std::string > leaprTape(leaprOut);
    int matde = 101, matdp = 1301, nbin  = 4, iinc  = 2, icoh  = 0, 
        iform = 1,   natom = 1,    mtref = 222;
 
    std::vector<double> temps {296.0};
    double  tol = 0.5, emax = 0.625; 

    auto out = thermr(matde, matdp, nbin, iinc, icoh, iform, natom, mtref, temps,
           tol, emax, leaprTape);
   
    auto my_section = out[0].MT(mtref); 

    auto good_products = good_inelastic.products();
    auto good_law      = std::get<LaboratoryAngleEnergy>(good_products[0].distribution());

    auto my_products = my_section.products();
    auto my_law = std::get<LaboratoryAngleEnergy>(my_products[0].distribution() );

    REQUIRE( my_section.AWR() == good_inelastic.AWR() );
    REQUIRE( my_section.MT() == good_inelastic.MT() );


    check_E_mu_Ep( good_law, my_law );

  } // GIVEN


} // TEST CASE 



