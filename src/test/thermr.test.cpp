#define CATCH_CONFIG_MAIN
#include "catch.hpp"
#include "thermr.cpp"
#include <range/v3/all.hpp>
#include "correct_h2o_info/tape23.h"
#include "correct_h2o_info/tape24.h"
#include <typeinfo>

using Tabulated = section::Type<7,4>::Tabulated;

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
    
      auto begin = h2o_leapr_out.begin(), end = h2o_leapr_out.end();
      long lineNumber = 1;
      StructureDivision division(begin,end,lineNumber);
      njoy::ENDFtk::file::Type<7> h2o_MF7(division,begin,end,lineNumber);
      
      int matde = 101;
      int matdp = 1301;
      int nbin  = 8;
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
