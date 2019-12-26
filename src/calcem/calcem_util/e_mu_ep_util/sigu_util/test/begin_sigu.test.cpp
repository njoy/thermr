#define CATCH_CONFIG_MAIN
#include "catch.hpp"
#include "calcem/calcem_util/e_mu_ep_util/sigu_util/begin_sigu.h"
#include <range/v3/all.hpp>
#include "generalTools/testing.h" 


TEST_CASE( "begin sigu (113,116)" ){
  GIVEN( "inputs" ){
    int jbeta = -7, lat = 1, iinc = 2, 
        lasym = 0;

    std::cout << std::setprecision(10) ;
    double tev = 1.5e-5, az = 11.9, tevz = 2.53e-1, 
      az2 = 0.0, teff2 = 0.0, sb = 5.53, sb2 = 0.0,
      teff = 6.14e-2, tolin = 5e-2, u,
      e = 1.0000000474974513E-003;

    std::vector<double> alpha { 1.1, 2.2, 3.3, 4.5, 5.8 },
      beta { 0.1, 0.2, 1.3, 1.4, 2.5, 2.6, 3.7 }, x(20,0.0), y(20,0.0);
    std::vector<double> sab(alpha.size()*beta.size());
    for ( size_t i = 0; i < alpha.size(); ++i ){
      for ( size_t j = 0; j < beta.size(); ++j ){
        sab[i*beta.size()+j] = 0.01*((j+1) + 0.1*(i+1));
      } 
    } 

    std::vector<double> cosines { -1.0, -0.8, -0.5, 0.1, 0.9 },
      correct_xVals(5), correct_yVals(5), correct_x(20,0.0), correct_y(20,0.0);

    WHEN( "lat = 1" ){
      correct_xVals = { 7.1395952E-4, 7.3838244E-4, 7.7664539E-4, 
                        2.6299999E-2, 2.6299999E-2 },
      correct_yVals = { 67960510168.6, 32016311225.0, 9925375996.58, 0.0, 0.0 };
      THEN("x, y, and jbeta are correctly returned for each scattering cosine"){
        for (size_t i = 0; i < cosines.size(); ++i){
          correct_x[0] = correct_xVals[i];
          correct_y[0] = correct_yVals[i];
          do_113_116( jbeta, lat, x, y, e, tev, tevz, cosines[i], alpha, beta, sab, 
                      az, lasym, teff, sb, sb2, iinc);
          REQUIRE( ranges::equal(x, correct_x, equal) );
          REQUIRE( ranges::equal(y, correct_y, equal) );
          REQUIRE( jbeta == 1 );
        }
      } // THEN
    } // WHEN

    WHEN( "lat = 0" ){
      lat = 0;
      correct_xVals = { 7.1395952E-4, 7.3838244E-4, 7.7664539E-4, 9.4450005E-4, 
                        9.4450005E-4 };
      correct_yVals = { 157.26389831, 167.22962548, 185.61792063, 252.19972257, 
                        1223152.9616};
      THEN("x, y, and jbeta are correctly returned for each scattering cosine"){
        for (size_t i = 0; i < cosines.size(); ++i){
          correct_x[0] = correct_xVals[i];
          correct_y[0] = correct_yVals[i];
          do_113_116( jbeta, lat, x, y, e, tev, tevz, cosines[i], alpha, beta, sab, 
                      az, lasym, teff, sb, sb2, iinc);
          REQUIRE( ranges::equal(x, correct_x, equal) );
          REQUIRE( ranges::equal(y, correct_y, equal) );
          REQUIRE( jbeta == -7 );
        }
      } // THEN

    } // WHEN





  } // GIVEN

  GIVEN( "inputs 2" ){
    int jbeta = -7, lat, iinc = 2, lasym = 0;

    std::cout << std::setprecision(10) ;
    double e = 1.0e-5, tev = 2.5507297688E-2, az = 0.99917,
      tevz = 2.53E-2, az2 = 0.0, teff2 = 0.0, sb = 163.72792237360667, sb2 = 0.0,
      teff = 0.12044192657731301, tolin = 5e-2, u = -1.0;

    std::vector<double> alpha { 1.1, 2.2, 3.3, 4.5, 5.8 },
      beta { 0.1, 0.2, 1.3, 1.4, 2.5, 2.6, 3.7 }, x(20,0.0), y(20,0.0);

    std::vector<double> sab(alpha.size()*beta.size());
    for ( size_t i = 0; i < alpha.size(); ++i ){
      for ( size_t j = 0; j < beta.size(); ++j ){
        sab[i*beta.size()+j] = 0.01*((j+1) + 0.1*(i+1));
      } 
    } 

    std::vector<double> cosines { -1.0, -0.8, -0.5, 0.1, 0.9 },
      correct_xVals(5), correct_yVals(5), correct_x(20,0.0), correct_y(20,0.0);


    //std::vector<double> 
    //correct_xVals  { 1.7236804E-12, 2.6945097E-12, 6.9119573E-12, 2.54E-3, 2.54E-3 },
    //correct_yVals  { 71.0333342974, 88.81240724978, 142.244019088, 
    //                              168116.849423, 177887.8702032};
                     
    //std::vector<double> correct_x (20,0.0); 
    //std::vector<double> correct_y (20,0.0); 


    WHEN( "lat = 1" ){
      lat = 1;
      correct_xVals = { 1.7236804E-12, 2.6945097E-12, 6.9119573E-12, 2.54E-3, 2.54E-3 },
      correct_yVals = { 71.0333342974, 88.81240724978, 142.244019088, 
                                    168116.849423, 177887.8702032};
 
      THEN("x, y, and jbeta are correctly returned for each scattering cosine"){
        for (size_t i = 0; i < cosines.size(); ++i){
          correct_x[0] = correct_xVals[i];
          correct_y[0] = correct_yVals[i];
          do_113_116( jbeta, lat, x, y, e, tev, tevz, cosines[i], alpha, beta, sab, 
                      az, lasym, teff, sb, sb2, iinc);
          REQUIRE( ranges::equal(x, correct_x, equal) );
          REQUIRE( ranges::equal(y, correct_y, equal) );
          REQUIRE( jbeta == 1 );
        }
      } // THEN
    } // WHEN


    WHEN( "lat = 1" ){
      lat = 0;
      cosines = { -1.0, -0.8, -0.5, 0.1, 0.9 },
      correct_xVals = { 1.7236804E-12, 2.6945097E-12, 6.9119573E-12, 2.5607298E-3, 
                        2.5607298E-3 };
      correct_yVals = { 71.323686538, 89.175432324, 142.8254485, 168734.045036, 
                        178497.453965};

      THEN("x, y, and jbeta are correctly returned for each scattering cosine"){
        for (size_t i = 0; i < cosines.size(); ++i){
          correct_x[0] = correct_xVals[i];
          correct_y[0] = correct_yVals[i];
          do_113_116( jbeta, lat, x, y, e, tev, tevz, cosines[i], alpha, beta, sab, 
                      az, lasym, teff, sb, sb2, iinc);
          REQUIRE( ranges::equal(x, correct_x, equal) );
          REQUIRE( ranges::equal(y, correct_y, equal) );
          REQUIRE( jbeta == 1 );
        }
      } // THEN
    } // WHEN
  } // GIVEN
} // TEST CASE
