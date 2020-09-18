#define CATCH_CONFIG_MAIN
#include "catch.hpp"
#include "thermr.cpp"
#include <range/v3/all.hpp>
//#include "correct_h2o_info/tape23.h"
#include "correct_h2o_info/tape24.h"
#include "generalTools/testing.h"
//#include "correct_h2o_info/correct_6222.h"
#include "correct_h2o_info/correct_6222_altered_egrid.h"
#include <typeinfo>

using namespace njoy::ENDFtk;
using Tabulated = section::Type<7,4>::Tabulated;
using ContinuumEnergyAngle = section::Type< 6 >::ContinuumEnergyAngle;
using ThermalScatteringData = section::Type< 6 >::ContinuumEnergyAngle::ThermalScatteringData;

TEST_CASE( "thermr" ){
  GIVEN( "" ){
    WHEN( "" ){

      auto begin = correct_6222.begin(), end = correct_6222.end();
      long lineNumber = 1;
      StructureDivision division1(begin,end,lineNumber);
      njoy::ENDFtk::file::Type<6> good_MF6(division1,begin,end,lineNumber);
      auto good_section = good_MF6.MT(222);
      auto good_products = good_section.products();
      auto good_law = std::get< ContinuumEnergyAngle >(good_products[0].distribution() );

        /*
        std::cout << "INCIDENT ENERGY    " << incidentEnergy << std::endl;
        std::cout << "PP                 " << pp << std::endl;
        */
        //std::cout << (outgoingEnergies|ranges::view::all) << std::endl;
        //std::cout << (cosines|ranges::view::all) << std::endl;

        /*
        std::cout << outgoingEnergies[0] << "  " << (cosines[0]|ranges::view::all) << std::endl;
        std::cout << outgoingEnergies[1] << "  " << (cosines[1]|ranges::view::all) << std::endl;
        std::cout << outgoingEnergies[2] << "  " << (cosines[2]|ranges::view::all) << std::endl;
        std::cout << outgoingEnergies[3] << "  " << (cosines[3]|ranges::view::all) << std::endl;
        break;
        */

      //}
        /*
      */

       
    
      begin = h2o_leapr_out.begin(), end = h2o_leapr_out.end();
      lineNumber = 1;
      StructureDivision division(begin,end,lineNumber);
      njoy::ENDFtk::file::Type<7> h2o_MF7(division,begin,end,lineNumber);
      
      int matde = 101;
      int matdp = 1301;
      int nbin  = 4;
      //int ntemp = 1;
      int iinc  = 2;
      int icoh  = 0; 
      int iform = 0; 
      int natom = 1;
      int mtref = 222;
      //int iprint = 2;
 
      std::vector<double> temps {296.0};
      double  tol = 0.5;
      double emax = 0.2; 

      auto out = thermr(matde, matdp, nbin, iinc, icoh, iform, natom, mtref, temps,
             tol, emax, h2o_MF7);
   
      auto my_section = *out;

      //std::cout << typeid(my_section).name() << std::endl;
      //std::cout << typeid(good_section).name() << std::endl;

      REQUIRE( my_section.ZA() == good_section.ZA() );
      REQUIRE( my_section.AWR() == good_section.AWR() );
      REQUIRE( my_section.MT() == good_section.MT() );

      //auto my_section = my_MF6;//.MT(222);
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
      
      std::cout.precision(15);
      for ( size_t i = 0; i < my_energies.size(); ++i ){ // Incoming energies
          auto my_entry   =   my_energies[i];
          auto good_entry = good_energies[i];
          std::cout << "i    " << i << std::endl;
          auto good_subsection = std::get<ThermalScatteringData>(good_entry);
          auto   my_subsection = std::get<ThermalScatteringData>(  my_entry);

          std::cout << "   A " << std::endl;
          REQUIRE( good_subsection.LANG() == my_subsection.LANG() );
          REQUIRE( good_subsection.LTT() == my_subsection.LTT() );
          REQUIRE( good_subsection.energy() == Approx(my_subsection.energy()).epsilon(1e-6) );
          std::cout << "   B " << std::endl;
          REQUIRE( good_subsection.NW() == my_subsection.NW() );
          REQUIRE( good_subsection.N2() == my_subsection.N2() );

          std::cout << "   C " << std::endl;
          auto good_data = good_subsection.data();
          auto   my_data =   my_subsection.data();

          std::cout << "   D " << std::endl;
          checkVec(good_subsection.data(),    my_subsection.data());
          checkVec(good_subsection.energies(),my_subsection.energies());
          checkVec(good_subsection.PP(),my_subsection.PP());
          REQUIRE( good_subsection.cosines().size() == my_subsection.cosines().size() );
          for (size_t i = 0; i < good_subsection.cosines().size(); ++i){
            checkVec(good_subsection.cosines()[i],my_subsection.cosines()[i]);
          }
      }


      //std::cout << products[0].multiplicity().energies().size() << std::endl;
      //std::cout << (products[0].multiplicity().multiplicities() | ranges::view::all )  << std::endl;
      //std::cout << (good_products[0].multiplicity().boundaries() | ranges::view::all )  << std::endl;
      //std::cout << (good_products[0].multiplicity().interpolants() | ranges::view::all )  << std::endl;






      //auto my_products = mf6.products();


      /*
      REQUIRE( 1001. == Approx( products[0].ZAP() ) );
      REQUIRE( 1001. == Approx( products[0].productIdentifier() ) );
      REQUIRE( 0.9986234 == Approx( products[0].AWP() ) );
      REQUIRE( 0.9986234 == Approx( products[0].productMassRatio() ) );
      REQUIRE( 0 == products[0].LIP() );
      REQUIRE( 0 == products[0].productModifierFlag() );
      REQUIRE( 1 == products[0].LAW() );
      */


      /*
      begin = correct_6222_IncE_1.begin(), end = correct_6222_IncE_1.end();
      lineNumber = 1;
      StructureDivision division1(begin,end,lineNumber);
      ThermalScatteringData  chunk(division1,begin,end,lineNumber);
      */


      /*
      std::string buffer;
      auto output = std::back_inserter(buffer);
      MF3.print(output,27);
      std::cout << buffer << std::endl;
      */
 
 
    } // WHEN
  } // GIVEN
} // TEST CASE
