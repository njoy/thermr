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
void checkLAW(Section good_section, Section my_section){
  auto good_products = good_section.products();
  auto good_law      = std::get<ContinuumEnergyAngle>(good_products[0].distribution());

  //REQUIRE( my_section.ZA() == good_section.ZA() );
  REQUIRE( my_section.AWR() == good_section.AWR() );
  REQUIRE( my_section.MT() == good_section.MT() );

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
    auto good_section  = MF6.MT(222);
    std::string leaprOut = njoy::utility::slurpFileToMemory("h2o_tape24");
    njoy::ENDFtk::syntaxTree::Tape<std::string > leaprTape(leaprOut);
    njoy::ENDFtk::file::Type<7> MF7 = 
      leaprTape.materialNumber(101).front().fileNumber(7).parse<7>();
    int matde = 101, matdp = 1301, nbin  = 4, iinc  = 2, icoh  = 0, 
        iform = 0,   natom = 1,    mtref = 222;
 
    std::vector<double> temps {296.0};
    double  tol = 0.5, emax = 0.625; 

    auto out = thermr(matde, matdp, nbin, iinc, icoh, iform, natom, mtref, temps,
           tol, emax, MF7);
   
    auto my_section = *out; // Because I'm currently returning an optional from thermr
    checkLAW(good_section, my_section);



  } // GIVEN

  GIVEN( "ENDF-B/VIII.0 Be" ){
    std::string endfFile = njoy::utility::slurpFileToMemory("be_tape25");
    njoy::ENDFtk::syntaxTree::Tape<std::string > tape(endfFile);


    njoy::ENDFtk::file::Type<6> MF6 = tape.materialNumber(425).front().fileNumber(6).parse<6>();
    auto good_section  = MF6.MT(222);
    std::string leaprOut = njoy::utility::slurpFileToMemory("be_tape24");
    njoy::ENDFtk::syntaxTree::Tape<std::string > leaprTape(leaprOut);
    njoy::ENDFtk::file::Type<7> MF7 = 
      leaprTape.materialNumber(26).front().fileNumber(7).parse<7>();



    int matde = 26, matdp = 425, nbin  = 4, iinc  = 2, icoh  = 2, 
        iform = 0,   natom = 1,  mtref = 222;
 
    std::vector<double> temps {296.0};
    double  tol = 0.5, emax = 0.2; 

    //std::cout << good_section.ZA() << std::endl;
    auto out = thermr(matde, matdp, nbin, iinc, icoh, iform, natom, mtref, temps,
           tol, emax, MF7);
   

    auto my_section = *out; // Because I'm currently returning an optional from thermr
    checkLAW(good_section, my_section);



  } // GIVEN


} // TEST CASE 




/*
TEST_CASE( "thermr trial " ){
    std::string coolTape25 = njoy::utility::slurpFileToMemory("tape25");
    std::cout << coolTape25 << std::endl;
} // TEST CARE 
*/


/*
TEST_CASE( "thermr" ){
  std::cout.precision(15);
  GIVEN( "H in H2O Example" ){
    WHEN( "" ){

      auto begin = correct_6222.begin(), end = correct_6222.end();
      long lineNumber = 1;
      StructureDivision division1(begin,end,lineNumber);
      njoy::ENDFtk::file::Type<6> good_MF6(division1,begin,end,lineNumber);
      auto good_section   = good_MF6.MT(222);
      auto good_products = good_section.products();
      auto good_law = std::get< ContinuumEnergyAngle >(good_products[0].distribution() );

      begin = h2o_leapr_out.begin(), end = h2o_leapr_out.end();
      lineNumber = 1;
      StructureDivision division(begin,end,lineNumber);
      njoy::ENDFtk::file::Type<7> h2o_MF7(division,begin,end,lineNumber);
      
      int matde = 101, matdp = 1301, nbin  = 4, iinc  = 2, icoh  = 0, 
          iform = 0,   natom = 1,    mtref = 222;
 
      std::vector<double> temps {296.0};
      double  tol = 0.5, emax = 0.2; 

      auto out = thermr(matde, matdp, nbin, iinc, icoh, iform, natom, mtref, temps,
             tol, emax, h2o_MF7);
   
      auto my_section = *out; // Because I'm currently returning an optional from thermr

      //std::cout << typeid(my_section).name() << std::endl;
      //std::cout << typeid(good_section).name() << std::endl;

      REQUIRE( my_section.ZA() == good_section.ZA() );
      REQUIRE( my_section.AWR() == good_section.AWR() );
      REQUIRE( my_section.MT() == good_section.MT() );

      auto my_products = my_section.products();
      auto my_law = std::get< ContinuumEnergyAngle >(my_products[0].distribution() );
      REQUIRE( my_law.LAW() == good_law.LAW() );
      REQUIRE( my_law.LEP() == good_law.LEP() );
      REQUIRE( my_law.NE() == good_law.NE() );
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

        checkVec(good_subsection.data(),    my_subsection.data());
        checkVec(good_subsection.energies(),my_subsection.energies());
        checkVec(good_subsection.PP(),my_subsection.PP());

        REQUIRE( good_subsection.cosines().size() == my_subsection.cosines().size() );
        for (size_t i = 0; i < good_subsection.cosines().size(); ++i){
          checkVec(good_subsection.cosines()[i],my_subsection.cosines()[i]);
        }
      }



 
    } // WHEN
  } // GIVEN
  */

  //GIVEN( "ENDF-B/VIII.0 Be Input" ){
      /*
      auto begin = correct_be.begin(), end = correct_be.end();
      long lineNumber = 1;
      StructureDivision division1(begin,end,lineNumber);
      njoy::ENDFtk::file::Type<6> good_MF6(division1,begin,end,lineNumber);
      auto good_section   = good_MF6.MT(222);
      auto good_products = good_section.products();
      auto good_law = std::get< ContinuumEnergyAngle >(good_products[0].distribution() );

      begin = be_leapr_out.begin(), end = be_leapr_out.end();
      lineNumber = 1;
      StructureDivision division(begin,end,lineNumber);
      njoy::ENDFtk::file::Type<7> be_MF7(division,begin,end,lineNumber);
      
      int matde = 26,  matdp = 425,  nbin  = 4, iinc  = 2, icoh  = 2, 
          iform = 0,   natom = 1,    mtref = 222;
 
      matdp = 4009; // really i need to read this in from tape23
      std::vector<double> temps {296.0};
      double  tol = 0.5, emax = 0.2; 

      //auto out = thermr(matde, matdp, nbin, iinc, icoh, iform, natom, mtref, temps,
      //       tol, emax, be_MF7);
   
      //auto section = *out;
      //std::cout << typeid(out).name() << std::endl;
      //std::cout << out->ZA() << std::endl;
      //auto my_section = *out; // Because I'm currently returning an optional from thermr


      //REQUIRE( my_section.ZA() == good_section.ZA() );
      //REQUIRE( my_section.AWR() == good_section.AWR() );
      //REQUIRE( my_section.MT() == good_section.MT() );

      //auto my_products = my_section.products();
      */
      /*
      auto my_law = std::get< ContinuumEnergyAngle >(my_products[0].distribution() );
      REQUIRE( my_law.LAW() == good_law.LAW() );
      REQUIRE( my_law.LEP() == good_law.LEP() );
      REQUIRE( my_law.NE() == good_law.NE() );
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

        checkVec(good_subsection.data(),    my_subsection.data());
        checkVec(good_subsection.energies(),my_subsection.energies());
        checkVec(good_subsection.PP(),my_subsection.PP());

        REQUIRE( good_subsection.cosines().size() == my_subsection.cosines().size() );
        for (size_t i = 0; i < good_subsection.cosines().size(); ++i){
          checkVec(good_subsection.cosines()[i],my_subsection.cosines()[i]);
        }
      }

      */



//  } // GIVEN
//} // TEST CASE









      /*
      std::string buffer;
      auto output = std::back_inserter(buffer);
      MF3.print(output,27);
      std::cout << buffer << std::endl;
      */
 
