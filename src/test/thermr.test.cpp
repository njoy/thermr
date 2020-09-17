#define CATCH_CONFIG_MAIN
#include "catch.hpp"
#include "thermr.cpp"
#include <range/v3/all.hpp>
//#include "correct_h2o_info/tape23.h"
#include "correct_h2o_info/tape24.h"
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
      njoy::ENDFtk::file::Type<6> MF6(division1,begin,end,lineNumber);
      auto section = MF6.MT(222);
      auto products = section.products();
      auto law = std::get< ContinuumEnergyAngle >(products[0].distribution() );

      //std::cout << products[0].multiplicity().energies().size() << std::endl;
      //std::cout << (products[0].multiplicity().multiplicities() | ranges::view::all )  << std::endl;
      std::cout << (products[0].multiplicity().boundaries() | ranges::view::all )  << std::endl;
      std::cout << (products[0].multiplicity().interpolants() | ranges::view::all )  << std::endl;


      auto energies = law.subsections();
      
      for ( const auto& entry : energies ){ // Incoming energies
          std::cout.precision(15);
        //auto lang = std::visit( [] (const auto& variant) {return variant.LANG();}, entry );
        //std::cout << lang << std::endl;
        auto subsection = std::get<ThermalScatteringData>(entry);
        //std::cout << subsection.LANG() << std::endl;
        //std::cout << subsection.LANG() << std::endl;
        //std::cout << subsection.LTT() << std::endl;
        //std::cout << subsection.energy() << std::endl;
        //std::cout << subsection.NW() << std::endl;
        //std::cout << subsection.N2() << std::endl;
        /*
        std::cout << law.LEP() << std::endl;
        std::cout << (law.boundaries()|ranges::view::all) << std::endl;
        std::cout << (law.interpolants()|ranges::view::all) << std::endl;
        std::cout << std::endl;
        std::cout << std::endl;
        */
        auto incidentEnergy   = subsection.energy();
        auto outgoingEnergies = subsection.energies();
        auto cosines = subsection.cosines();
        auto pp = subsection.PP();

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

      }
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

      thermr(matde, matdp, nbin, iinc, icoh, iform, natom, mtref, temps,
             tol, emax, h2o_MF7);
   

      /*
      begin = correct_6222_IncE_1.begin(), end = correct_6222_IncE_1.end();
      lineNumber = 1;
      StructureDivision division1(begin,end,lineNumber);
      ThermalScatteringData  chunk(division1,begin,end,lineNumber);
      */




      /*
     
             */
      
      /*
      std::string buffer;
      auto output = std::back_inserter(buffer);
      MF3.print(output,27);
      std::cout << buffer << std::endl;
      */
 
 
        //std::string fileName = "/Users/ameliajo/thermr/src/test/tape24";
        //std::string fileName = "/Users/ameliajo/thermr/src/test/tape24_smaller_alpha_beta_grid";
        //int mat = 101;
        //int iform = 0;
        //int iinc = 2;
        //int nbin = 8;
        //int natom = 1;
        //std::vector<double> temperatures {296.0};
        //double tol = 0.05, emax = 0.625;
        //thermr(mat,fileName,iform,iinc,nbin,temperatures,tol,emax,natom);
        //REQUIRE( 0 == 0 );
    } // WHEN
  } // GIVEN
} // TEST CASE
