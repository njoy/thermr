#define CATCH_CONFIG_MAIN
#include "catch.hpp"
#include "thermr.cpp"
#include <range/v3/all.hpp>
#include "correct_h2o_info/tape24.h"
#include "generalTools/testing.h"
#include "correct_h2o_info/correct_6222_altered_egrid.h"
#include "correct_be_info/correct_6222.h"
#include "correct_be_info/tape24.h"
#include <typeinfo>

using namespace njoy::ENDFtk;
using Tabulated = section::Type<7,4>::Tabulated;
using ContinuumEnergyAngle = section::Type< 6 >::ContinuumEnergyAngle;
using ThermalScatteringData = section::Type< 6 >::ContinuumEnergyAngle::ThermalScatteringData;

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
    //REQUIRE( good_subsection.NW()     == my_subsection.NW() );
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
    njoy::ENDFtk::file::Type<7> MF7 = 
      leaprTape.materialNumber(101).front().fileNumber(7).parse<7>();
    int matde = 101, matdp = 1301, nbin  = 4, iinc  = 2, icoh  = 0, 
        iform = 0,   natom = 1,    mtref = 222;
 
    std::vector<double> temps {296.0};
    double  tol = 0.5, emax = 0.625; 

    //auto out = thermr(matde, matdp, nbin, iinc, icoh, iform, natom, mtref, temps,
    //       tol, emax, MF7);
   
    //auto my_section = out[0].MT(mtref); // Because I'm currently returning an optional from thermr
    //checkInelastic(good_inelastic, my_section);



  } // GIVEN
    /*
  */

  GIVEN( "ENDF-B/VIII.0 Be" ){
    std::string endfFile = njoy::utility::slurpFileToMemory("be_tape25");
    njoy::ENDFtk::syntaxTree::Tape<std::string > tape(endfFile);


    njoy::ENDFtk::file::Type<6> MF6 = tape.materialNumber(425).front().fileNumber(6).parse<6>();
    auto good_inelastic_section = MF6.MT(222);
    auto good_elastic_section   = MF6.MT(223);
    std::string leaprOut = njoy::utility::slurpFileToMemory("be_tape24");
    njoy::ENDFtk::syntaxTree::Tape<std::string > leaprTape(leaprOut);
    njoy::ENDFtk::file::Type<7> MF7 = 
      leaprTape.materialNumber(26).front().fileNumber(7).parse<7>();


    /*
    std::cout << good_elastic.ZA() << std::endl;
    std::cout << good_elastic.AWR() << std::endl;
    std::cout << good_elastic.JP() << std::endl;
    std::cout << good_elastic.LCT() << std::endl;
    std::cout << good_elastic.NK() << std::endl;
    std::cout << good_elastic.MT() << std::endl;
    std::cout << std::endl;
    */
    //auto products = good_elastic_section.products();
    /*
    //std::cout << products.size() << std::endl;
    std::cout << products[0].ZAP() << std::endl;
    std::cout << products[0].AWP() << std::endl;
    std::cout << products[0].LIP() << std::endl;
    std::cout << products[0].LAW() << std::endl;
    */
    //auto multiplicity = products[0].multiplicity();
    /*
    std::cout << std::endl;
    std::cout << (multiplicity.energies()|ranges::view::all) << std::endl;
    //std::cout << (multiplicity.boundaries()|ranges::view::all) << std::endl;
    std::cout << (multiplicity.multiplicities()|ranges::view::all) << std::endl;
    //std::cout << (multiplicity.interpolants()|ranges::view::all) << std::endl;

    */

    int matde = 26, matdp = 425, nbin  = 4, iinc  = 2, icoh  = 2, 
        iform = 0,   natom = 1,  mtref = 222;
 
    std::vector<double> temps {296.0};
    double  tol = 0.5, emax = 0.2; 

    //auto out = thermr(matde, matdp, nbin, iinc, icoh, iform, natom, mtref, temps,
    //       tol, emax, MF7);
    //auto my_inelastic_section = out[0].MT(mtref); 
    //std::cout << "going into test" << std::endl;
    //checkInelastic(good_inelastic_section, my_inelastic_section);
    //std::cout << "out of test" << std::endl;
    //auto my_elastic_section = out[0].MT(mtref+1); 
    //auto good_elastic_section = out[1].MT(mtref+1); 
    //std::cout << my_elastic_section.ZA() << std::endl;

  //REQUIRE( my_section.ZA() == good_inelastic.ZA() );
  //REQUIRE( my_elastic_section.AWR() == good_elastic_section.AWR() );
  //REQUIRE( my_elastic_section.MT() == good_elastic_section.MT() );

  //auto good_products = good_elastic_section.products();
  //auto good_law      = std::get<ContinuumEnergyAngle>(good_products[0].distribution());
  //auto my_products = my_elastic_section.products();
  //auto my_law = std::get< ContinuumEnergyAngle >(my_products[0].distribution() );
  //REQUIRE( my_law.LAW() == good_law.LAW() );
  //REQUIRE( my_law.LEP() == good_law.LEP() );
  //REQUIRE( my_law.NE() == good_law.NE() ); // number inc. E points on egrid
  //REQUIRE( my_law.NR() == good_law.NR() );

    /*

  */


  } // GIVEN

  GIVEN( "NJOY Test 9 - H in H2O Example (Multiple Temps)" ){
    std::string endfFile = njoy::utility::slurpFileToMemory("h2oMultT_tape25");
    njoy::ENDFtk::syntaxTree::Tape<std::string > tape(endfFile);

    /*

    std::vector<njoy::ENDFtk::file::Type<6>> MF6_variousTemps;
    for (auto& material : tape.materialNumber(1301)){
      //njoy::ENDFtk::file::Type<6> thisMF6 = material.fileNumber(6).parse<6>();
      MF6_variousTemps.push_back(material.fileNumber(6).parse<6>());
    }

    auto good_inelastic  = MF6.MT(222);
    std::string leaprOut = njoy::utility::slurpFileToMemory("h2oMultT_tape24");
    njoy::ENDFtk::syntaxTree::Tape<std::string > leaprTape(leaprOut);
    njoy::ENDFtk::file::Type<7> MF7 = 
      leaprTape.materialNumber(101).front().fileNumber(7).parse<7>();
      */

    /*
    auto good_products = good_inelastic.products();
    auto good_law      = std::get<ContinuumEnergyAngle>(good_products[0].distribution());
    auto good_multiplicity = good_products[0].multiplicity();
    std::cout << (good_multiplicity.energies()|ranges::view::all)  << std::endl;
    std::cout << (good_multiplicity.multiplicities()|ranges::view::all)  << std::endl;
     
    auto good_energies = good_law.subsections();
    auto good_subsection = std::get<ThermalScatteringData>(good_energies[0]);

    std::cout << good_subsection.energy() << std::endl;
    std::cout << good_law.LAW() << std::endl;
    std::cout << good_law.LEP() << std::endl;
    std::cout << good_law.NE() << std::endl;
    std::cout << good_law.NR() << std::endl;
    */




    //REQUIRE( my_section.ZA() == good_inelastic.ZA() );
    //REQUIRE( my_section.AWR() == good_inelastic.AWR() );
    //REQUIRE( my_section.MT() == good_inelastic.MT() );

    //REQUIRE( my_law.LAW() == good_law.LAW() );
    //REQUIRE( my_law.LEP() == good_law.LEP() );
    //REQUIRE( my_law.NE() == good_law.NE() ); // number inc. E points on egrid
    //REQUIRE( my_law.NR() == good_law.NR() );


    int matde = 101, matdp = 1301, nbin  = 4, iinc  = 2, icoh  = 0, 
        iform = 0,   natom = 1,    mtref = 222;
 
    std::vector<double> temps {296.0,500.0};
    double  tol = 0.5, emax = 0.625; 

    //auto out = thermr(matde, matdp, nbin, iinc, icoh, iform, natom, mtref, temps,
    //       tol, emax, MF7);
   
    //auto my_section = out[0].MT(mtref); 
    //checkInelastic(good_inelastic, my_section);

    /*
    */


  } // GIVEN

  GIVEN( "NJOY Test 9 - H in H2O Example (E mu E')" ){
    std::string endfFile = njoy::utility::slurpFileToMemory("h2oEmuEp_tape25");
    njoy::ENDFtk::syntaxTree::Tape<std::string > tape(endfFile);


    njoy::ENDFtk::file::Type<6> MF6 = tape.materialNumber(1301).front().fileNumber(6).parse<6>();
    auto good_inelastic  = MF6.MT(222);
    std::string leaprOut = njoy::utility::slurpFileToMemory("h2oEmuEp_tape24");
    njoy::ENDFtk::syntaxTree::Tape<std::string > leaprTape(leaprOut);
    njoy::ENDFtk::file::Type<7> MF7 = 
      leaprTape.materialNumber(101).front().fileNumber(7).parse<7>();
    int matde = 101, matdp = 1301, nbin  = 4, iinc  = 2, icoh  = 0, 
        iform = 1,   natom = 1,    mtref = 222;
 
    std::vector<double> temps {296.0};
    double  tol = 0.5, emax = 0.625; 

    //auto out = thermr(matde, matdp, nbin, iinc, icoh, iform, natom, mtref, temps,
    //       tol, emax, MF7);
   
    //auto my_section = out[0].MT(mtref); // Because I'm currently returning an optional from thermr
    //checkInelastic(good_inelastic, my_section);



  } // GIVEN


} // TEST CASE 



