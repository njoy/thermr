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

    u = -1.0;
    do_113_116( jbeta, lat, x, y, e, tev, tevz, u, alpha, beta, sab, az, lasym, 
                teff, sb, sb2, iinc);
    std::vector<double> correct_x (20,0.0); correct_x[0] = 7.1395952E-4;
    std::vector<double> correct_y (20,0.0); correct_y[0] = 63376459056.429459;  
    correct_x[0] = 7.1395952000007134E-004;
    correct_y[0] = 67960510168.646820;
    REQUIRE( ranges::equal(x, correct_x, equal) );
    REQUIRE( ranges::equal(y, correct_y, equal) );
    REQUIRE( jbeta == 1 );

    u = -0.8;
    do_113_116( jbeta, lat, x, y, e, tev, tevz, u, alpha, beta, sab, az, lasym, 
                teff, sb, sb2, iinc);
    correct_x[0] = 7.3838244000007369E-004;
    correct_y[0] = 32016311225.079830;
    REQUIRE( ranges::equal(x, correct_x, equal) );
    REQUIRE( ranges::equal(y, correct_y, equal) );
    REQUIRE( jbeta == 1 );

 

    u = -0.5;
    do_113_116( jbeta, lat, x, y, e, tev, tevz, u, alpha, beta, sab, az, lasym, 
                teff, sb, sb2, iinc);
    correct_x[0] = 7.7664539000007755E-004;
    correct_y[0] = 9925375996.5867958;
    REQUIRE( ranges::equal(x, correct_x, equal) );
    REQUIRE( ranges::equal(y, correct_y, equal) );
    REQUIRE( jbeta == 1 );

   
    u = 0.1;
    do_113_116( jbeta, lat, x, y, e, tev, tevz, u, alpha, beta, sab, az, lasym, 
                teff, sb, sb2, iinc);
    correct_x[0] = 2.6299999000002627E-002;
    correct_y[0] = 0.0;
    REQUIRE( ranges::equal(x, correct_x, equal) );
    REQUIRE( ranges::equal(y, correct_y, equal) );
    REQUIRE( jbeta == 1 );

    u = 0.9;
    do_113_116( jbeta, lat, x, y, e, tev, tevz, u, alpha, beta, sab, az, lasym, 
                teff, sb, sb2, iinc);
    correct_x[0] = 2.6299999000002627E-002;
    correct_y[0] = 0.0;
    REQUIRE( ranges::equal(x, correct_x, equal) );
    REQUIRE( ranges::equal(y, correct_y, equal) );
    REQUIRE( jbeta == 1 );






  } // GIVEN
    /*

  */

  GIVEN( "inputs 2" ){
    int jbeta = -7, lat = 1, iinc = 2, lasym = 0;

    std::cout << std::setprecision(10) ;
    double e = 1.0e-5, tev = 2.5507297687999999E-002, az = 0.99917000000000000,
      tevz = 2.5300000000000000E-002, az2 = 0.0, teff2 = 0.0, sb = 163.72792237360667, sb2 = 0.0,
      teff = 0.12044192657731301, tolin = 5e-2, u = -1.0;

  std::vector<double> alpha { 1.1, 2.2, 3.3, 4.5, 5.8 },
    beta { 0.1, 0.2, 1.3, 1.4, 2.5, 2.6, 3.7 }, x(20,0.0), y(20,0.0);


  std::vector<double> sab(alpha.size()*beta.size());
  for ( size_t i = 0; i < alpha.size(); ++i ){
    for ( size_t j = 0; j < beta.size(); ++j ){
      sab[i*beta.size()+j] = 0.01*((j+1) + 0.1*(i+1));
    } 
  } 


    u = -1.0;
    do_113_116( jbeta, lat, x, y, e, tev, tevz, u, alpha, beta, 
        sab, az, lasym, teff, sb, sb2, iinc);
    std::vector<double> correct_x (20,0.0); correct_x[0] = 1.7236804E-12;
    std::vector<double> correct_y (20,0.0); correct_y[0] = 71.033334297426663;  
    REQUIRE( ranges::equal(x, correct_x, equal) );
    REQUIRE( ranges::equal(y, correct_y, equal) );
    REQUIRE( 1 == jbeta );

    u = -0.8;
    do_113_116( jbeta, lat, x, y, e, tev, tevz, u, alpha, beta, 
        sab, az, lasym, teff, sb, sb2, iinc);
    correct_x[0] = 2.6945097000002692E-012;
    correct_y[0] = 88.812407249785352;  
    REQUIRE( ranges::equal(x, correct_x, equal) );
    REQUIRE( ranges::equal(y, correct_y, equal) );
    REQUIRE( 1 == jbeta );

    u = -0.5;
    do_113_116( jbeta, lat, x, y, e, tev, tevz, u, alpha, beta, 
        sab, az, lasym, teff, sb, sb2, iinc);
    correct_x[0] = 6.9119573000006907E-012;
    correct_y[0] = 142.24401908809887;  
    REQUIRE( ranges::equal(x, correct_x, equal) );
    REQUIRE( ranges::equal(y, correct_y, equal) );
    REQUIRE( 1 == jbeta );

    u = 0.1; 
    do_113_116( jbeta, lat, x, y, e, tev, tevz, u, alpha, beta, 
        sab, az, lasym, teff, sb, sb2, iinc);
    correct_x[0] = 2.5400000000002539E-3;
    correct_y[0] = 168116.84942347769;
    REQUIRE( ranges::equal(x, correct_x, equal) );
    REQUIRE( ranges::equal(y, correct_y, equal) );
    REQUIRE( 1 == jbeta );

    u = 0.9; 
    do_113_116( jbeta, lat, x, y, e, tev, tevz, u, alpha, beta, 
        sab, az, lasym, teff, sb, sb2, iinc);
    correct_x[0] = 2.5400000000002539E-3;
    correct_y[0] = 177887.87020327282;
    REQUIRE( ranges::equal(x, correct_x, equal) );
    REQUIRE( ranges::equal(y, correct_y, equal) );
    REQUIRE( 1 == jbeta );






  } // GIVEN

  /*


  GIVEN( "inputs 3" ){
    int jbeta = -80, lat = 1, iinc = 2, 
        lasym = 0;

    std::cout << std::setprecision(10) ;
    double e = 1.0e-6, tev = 1.5e-4, az = 11.9,
      tevz = 2.2e-4, az2 = 0.0, teff2 = 0.0, sb = 5.53, sb2 = 0.0,
      teff = 6.14e-2, tolin = 5e-2, u = 0.1, root1 = 0.00092700235;

    std::vector<double> alpha(40),beta(80),x(20,0.0), y(20,0.0);
    for ( int i = 0; i < 40;     ++i ){ alpha[i] = 0.1*i + i%6 + 0.001; }
    for ( int i = 0; i < 80;     ++i ){ beta[i]  = 0.2*i + i%4 + 0.025; }

  std::vector<double> sab(alpha.size()*beta.size());
  for ( size_t i = 0; i < alpha.size(); ++i ){
    for ( size_t j = 0; j < beta.size(); ++j ){
      sab[i*beta.size()+j] = 0.01*((j+1) + 0.1*(i+1));
    } 
  } 


    do_113_116( jbeta, lat, x, y, e, tev, tevz, u, alpha, beta, 
        sab, az, lasym, teff, sb, sb2, iinc);

   // REQUIRE( 6.5e-6 == Approx(x[0]).epsilon(1e-6) );
   // REQUIRE( 46226.360425 == Approx(y[0]).epsilon(1e-6) ); 

    for ( size_t i = 1; i < x.size(); ++i ){ 
    //  REQUIRE( 0 == Approx(x[i]).epsilon(1e-6) ); 
    //  REQUIRE( 0 == Approx(y[i]).epsilon(1e-6) ); 
    }
 
    //REQUIRE( 1 == jbeta );

  } // GIVEN



*/
} // TEST CASE


