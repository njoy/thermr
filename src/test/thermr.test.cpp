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
        /*
      auto begin = h2o_nin.begin(), end = h2o_nin.end();
      long lineNumber = 1;
      StructureDivision division(begin,end,lineNumber);
      njoy::ENDFtk::file::Type<3> MF3(division,begin,end,lineNumber);
      //njoy::ENDFtk::section::Type<3,1> MT1 = MF3.MT(1_c);
      auto MT1 = MF3.MT(1_c);
      auto MT2 = MF3.MT(2_c);
      std::cout << MT1.ZA() << std::endl; 
      std::cout << MT1.AWR() << std::endl; 
      */
      auto begin = correct_6222.begin(), end = correct_6222.end();
      long lineNumber = 1;
      StructureDivision division1(begin,end,lineNumber);
      njoy::ENDFtk::file::Type<6> MF6(division1,begin,end,lineNumber);
      auto section = MF6.MT(222);
      auto products = section.products();
      //std::cout << products.size() << std::endl;

      auto law = std::get< ContinuumEnergyAngle >(products[0].distribution() );

      //std::cout << "LAW    " << law.LAW() << std::endl;
      //std::cout << "NW     " << law.NW()  << std::endl;
      //std::cout << "N2     " << law.N2()  << std::endl;
      //std::cout << "N2     " << law.N2()  << std::endl;
      //auto data = law.

        /*
      */



      /*
  CHECK( 5 == chunk.LTT() );
  CHECK( 1e-5 == Approx( chunk.energy() ) );

  CHECK( 18 == chunk.NW() );
  CHECK( 6 == chunk.N2() );

  auto data = chunk.data();
  CHECK( 18 == data.size() );
  CHECK( 0. == Approx( data[0] ) );
  CHECK( 0. == Approx( data[1] ) );
  CHECK( 0. == Approx( data[2] ) );
  CHECK( 0. == Approx( data[3] ) );
  CHECK( 0. == Approx( data[4] ) );
  CHECK( 0. == Approx( data[5] ) );
  CHECK( 9.999999e-6 == Approx( data[6] ) );
  CHECK( 9.477167e+1 == Approx( data[7] ) );
  CHECK( -5.379121e-1 == Approx( data[8] ) );
  CHECK( 0.21062848 == Approx( data[9] ) );
  CHECK( 0.70490082 == Approx( data[10] ) );
  CHECK( 9.552579e-1 == Approx( data[11] ) );
  CHECK( 1.265100e-1 == Approx( data[12] ) );
  CHECK( 0. == Approx( data[13] ) );
  CHECK( 0. == Approx( data[14] ) );
  CHECK( 0. == Approx( data[15] ) );
  CHECK( 0. == Approx( data[16] ) );
  CHECK( 0. == Approx( data[17] ) );

  CHECK( 3 == chunk.NEP() );
  CHECK( 3 == chunk.numberEnergies() );

  auto energies = chunk.energies();
  CHECK( 3 == energies.size() );
  CHECK( 0. == Approx( energies[0] ) );
  CHECK( 9.999999e-6 == Approx( energies[1] ) );
  CHECK( 1.265100e-1 == Approx( energies[2] ) );

  auto pp = chunk.PP();
  CHECK( 3 == pp.size() );
  CHECK( 0. == Approx( pp[0] ) );
  CHECK( 9.477167e+1 == Approx( pp[1] ) );
  CHECK( 0. == Approx( pp[2] ) );

  auto cosines = chunk.cosines();
  CHECK( 3 == cosines.size() );
  CHECK( 4 == cosines[0].size() );
  CHECK( 0. == Approx( cosines[0][0] ) );
  CHECK( 0. == Approx( cosines[0][1] ) );
  CHECK( 0. == Approx( cosines[0][2] ) );
  CHECK( 0. == Approx( cosines[0][3] ) );
  CHECK( 4 == cosines[1].size() );
  CHECK( -5.379121e-1 == Approx( cosines[1][0] ) );
  CHECK( 0.21062848 == Approx( cosines[1][1] ) );
  CHECK( 0.70490082 == Approx( cosines[1][2] ) );
  CHECK( 9.552579e-1 == Approx( cosines[1][3] ) );
  CHECK( 4 == cosines[2].size() );
  CHECK( 0. == Approx( cosines[2][0] ) );
  CHECK( 0. == Approx( cosines[2][1] ) );
  CHECK( 0. == Approx( cosines[2][2] ) );
  CHECK( 0. == Approx( cosines[2][3] ) );

  CHECK( 4 == chunk.NC() );
  */





        /*
      auto energies = law.subsections();
      
      for ( const auto& entry : energies ){ // Incoming energies
          std::cout.precision(15);
        //auto lang = std::visit( [] (const auto& variant) {return variant.LANG();}, entry );
        //std::cout << lang << std::endl;
        auto subsection = std::get<ThermalScatteringData>(entry);
        //std::cout << subsection.LANG() << std::endl;
        std::cout << subsection.LANG() << std::endl;
        std::cout << subsection.LTT() << std::endl;
        std::cout << subsection.energy() << std::endl;
        std::cout << subsection.NW() << std::endl;
        std::cout << subsection.N2() << std::endl;
        std::cout << std::endl;
        std::cout << std::endl;
        auto incidentEnergy   = subsection.energy();
        auto outgoingEnergies = subsection.energies();
        auto cosines = subsection.cosines();
        auto pp = subsection.PP();

        //std::cout << pp << std::endl;
        //std::cout << "INCIDENT ENERGY    " << incidentEnergy << std::endl;
        //std::cout << (outgoingEnergies|ranges::view::all) << std::endl;
        //std::cout << (cosines|ranges::view::all) << std::endl;

        //std::cout << outgoingEnergies[0] << "  " << (cosines[0]|ranges::view::all) << std::endl;
        //std::cout << outgoingEnergies[1] << "  " << (cosines[1]|ranges::view::all) << std::endl;
        //std::cout << outgoingEnergies[2] << "  " << (cosines[2]|ranges::view::all) << std::endl;
        //std::cout << outgoingEnergies[3] << "  " << (cosines[3]|ranges::view::all) << std::endl;
        break;

      }
      */

       



      //njoy::ENDFtk::section::Type<3,1> MT1 = MF3.MT(1_c);
      //auto MT1 = MF3.MT(1_c);
      //auto MT2 = MF3.MT(2_c);
      //std::cout << MT1.ZA() << std::endl; 
      //std::cout << MT1.AWR() << std::endl; 
 
    
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
