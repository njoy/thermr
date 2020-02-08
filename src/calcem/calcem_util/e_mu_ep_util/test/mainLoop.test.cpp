#include "catch.hpp"
#include "calcem/calcem_util/e_mu_ep_util/mainLoop.h"
#include "generalTools/testing.h"


TEST_CASE( " mu-E' ordering " ){
  int imax = 20, lat = 0, iinc = 2, lasym = 0;
  double tev = 2.5507297688e-2, tol = 0.05;
  std::vector<double> 
  alphas { 1.1, 2.2, 3.3, 4.5, 5.8 },
  betas { 0.1, 0.2, 1.3, 1.4, 2.5, 2.6, 3.7 },
  sab {-0.18259619, -0.30201347, -3.93654779, -3.98809174, -4.33545607, 
  -4.39515402, -5.88934921, -0.76225291, -0.81658341, -3.14161459, -3.30566186, 
  -3.90554652, -3.96233362, -5.23696660, -1.19182884, -1.23155471, -2.79610565, 
  -2.95633099, -3.74989225, -3.80837585, -4.93373911, -1.58342860, -1.61310713, 
  -2.71233943, -2.84291608, -3.69699590, -3.75199349, -4.77433858, -1.96121202, 
  -1.98720663, -2.78454600, -2.88531460, -3.71288120, -3.77142141, -4.71158392 };

  double az = 0.99917, sigma_b = 163.72792237, sigma_b2 = 0.0, teff = 0.120441926;
  std::vector<double> eVec, correctEnergies;
  double enow = 1e-5;

  std::vector<double> correct_uj;




  enow = 1e-6; 
  auto out = do_530_etc(enow,tev,tol,lat,iinc,lasym,alphas,betas,sab,az,sigma_b,sigma_b2,teff);
  REQUIRE( 6325.1941714771147 == Approx(std::get<0>(out)).epsilon(1e-6) );
  REQUIRE( 2.8371734942764979E-003 == Approx(std::get<1>(out)).epsilon(1e-6) );
  correct_uj = {-1, -0.75, -0.625, -0.5, -0.375, -0.25, -0.125, 0, 0.125, 0.25, 0.375, 0.5, 0.625, 0.75, 1 };
  for (size_t i = 0; i < correct_uj.size(); ++i){
    REQUIRE( correct_uj[i] == Approx(std::get<2>(out)[i]).epsilon(1e-6) );
  }




  enow = 1e-5;
  out = do_530_etc(enow,tev,tol,lat,iinc,lasym,alphas,betas,sab,az,sigma_b,sigma_b2,teff);
  REQUIRE( 2004.20546598 == Approx(std::get<0>(out)).epsilon(1e-6) );
  REQUIRE( 0.00903406640 == Approx(std::get<1>(out)).epsilon(1e-6) );
  correct_uj = {-1, -0.75, -0.625, -0.5, -0.375, -0.25, -0.125, 0, 0.125, 0.25, 0.375, 0.5, 0.625, 0.75, 1 };
  for (size_t i = 0; i < correct_uj.size(); ++i){
    REQUIRE( correct_uj[i] == Approx(std::get<2>(out)[i]).epsilon(1e-6) );
  }



  enow = 1e-4;
  out = do_530_etc(enow,tev,tol,lat,iinc,lasym,alphas,betas,sab,az,sigma_b,sigma_b2,teff);
  REQUIRE( 648.0310824 == Approx(std::get<0>(out)).epsilon(1e-6) );
  REQUIRE( 0.031418177 == Approx(std::get<1>(out)).epsilon(1e-6) );
  correct_uj = {-1, -0.75, -0.625, -0.5, -0.375, -0.25, -0.125, 0, 0.125, 0.25, 0.375, 0.5, 0.625, 0.75, 1 };
  for (size_t i = 0; i < correct_uj.size(); ++i){
    REQUIRE( correct_uj[i] == Approx(std::get<2>(out)[i]).epsilon(1e-6) );
  }


  enow = 1e-2;
  out = do_530_etc(enow,tev,tol,lat,iinc,lasym,alphas,betas,sab,az,sigma_b,sigma_b2,teff);
  REQUIRE( 123.58909310679307 == Approx(std::get<0>(out)).epsilon(1e-6) );
  REQUIRE( 9.0734245819064419E-002 == Approx(std::get<1>(out)).epsilon(1e-6) );
  correct_uj = {-1, -0.75, -0.625, -0.5, -0.375, -0.25, -0.125, 0, 0.125, 0.25, 0.375, 0.5, 0.625, 0.75, 1 };
  for (size_t i = 0; i < correct_uj.size(); ++i){
    REQUIRE( correct_uj[i] == Approx(std::get<2>(out)[i]).epsilon(1e-6) );
  }


  enow = 1.0;
  out = do_530_etc(enow,tev,tol,lat,iinc,lasym,alphas,betas,sab,az,sigma_b,sigma_b2,teff);
  REQUIRE( 40.227979467296890 == Approx(std::get<0>(out)).epsilon(1e-6) );
  REQUIRE( 0.54810870138245560  == Approx(std::get<1>(out)).epsilon(1e-6) );
  correct_uj = {-1, -0.75, -0.625, -0.5, -0.375, -0.25, -0.125, 0, 0.125, 0.25, 0.375, 0.5, 0.625, 0.75, 1 };
  for (size_t i = 0; i < correct_uj.size(); ++i){
    REQUIRE( correct_uj[i] == Approx(std::get<2>(out)[i]).epsilon(1e-6) );
  }





  enow = 1e-2;


  tol = 1.0;
  out = do_530_etc(enow,tev,tol,lat,iinc,lasym,alphas,betas,sab,az,sigma_b,sigma_b2,teff);
  REQUIRE( 122.34746244407418 == Approx(std::get<0>(out)).epsilon(1e-6) );
  REQUIRE( 2.4551863335822984E-002 == Approx(std::get<1>(out)).epsilon(1e-6) );
  correct_uj = {-1, -0.75, -0.625, -0.5, -0.375, -0.25, -0.125, 0, 0.125, 0.25, 0.375, 0.5, 0.625, 0.75, 1 };
  for (size_t i = 0; i < correct_uj.size(); ++i){
    REQUIRE( correct_uj[i] == Approx(std::get<2>(out)[i]).epsilon(1e-6) );
  }

  tol = 5.0;
  out = do_530_etc(enow,tev,tol,lat,iinc,lasym,alphas,betas,sab,az,sigma_b,sigma_b2,teff);
  REQUIRE( 143.75464230251364 == Approx(std::get<0>(out)).epsilon(1e-6) );
  REQUIRE( 3.8463983799896351E-002 == Approx(std::get<1>(out)).epsilon(1e-6) );
  correct_uj = {-1, -0.75, -0.625, -0.5, -0.375, -0.25, -0.125, 0, 0.125, 0.25, 0.375, 0.5, 0.625, 0.75, 1 };
  for (size_t i = 0; i < correct_uj.size(); ++i){
    REQUIRE( correct_uj[i] == Approx(std::get<2>(out)[i]).epsilon(1e-6) );
  }

  tol = 10.0;
  out = do_530_etc(enow,tev,tol,lat,iinc,lasym,alphas,betas,sab,az,sigma_b,sigma_b2,teff);
  REQUIRE( 145.79479928728830 == Approx(std::get<0>(out)).epsilon(1e-6) );
  REQUIRE( 5.0169921218852623E-002 == Approx(std::get<1>(out)).epsilon(1e-6) );
  correct_uj = {-1, -0.75, -0.625, -0.5, -0.375, -0.25, -0.125, 0, 0.125, 0.25, 0.375, 0.5, 0.625, 0.75, 1 };
  for (size_t i = 0; i < correct_uj.size(); ++i){
    REQUIRE( correct_uj[i] == Approx(std::get<2>(out)[i]).epsilon(1e-6) );
  }



} // TEST CASE







TEST_CASE( "E-mu-E' ordering " ){
  int imax = 20, lat = 0, iinc = 2, lasym = 0;
  double tev = 2.5507297688e-2, tol = 0.05;
  std::vector<double> 
  alphas { 1.1, 2.2, 3.3, 4.5, 5.8 },
  betas { 0.1, 0.2, 1.3, 1.4, 2.5, 2.6, 3.7 },
  sab {-0.18259619, -0.30201347, -3.93654779, -3.98809174, -4.33545607, 
  -4.39515402, -5.88934921, -0.76225291, -0.81658341, -3.14161459, -3.30566186, 
  -3.90554652, -3.96233362, -5.23696660, -1.19182884, -1.23155471, -2.79610565, 
  -2.95633099, -3.74989225, -3.80837585, -4.93373911, -1.58342860, -1.61310713, 
  -2.71233943, -2.84291608, -3.69699590, -3.75199349, -4.77433858, -1.96121202, 
  -1.98720663, -2.78454600, -2.88531460, -3.71288120, -3.77142141, -4.71158392 };

  double az = 0.99917, sigma_b = 163.72792237, sigma_b2 = 0.0, teff = 0.120441926;
  std::vector<double> eVec, correctEnergies;

  double enow = 1e-2;

  mu_ep(enow,tev,tol,lat,iinc,lasym,alphas,betas,sab,az,sigma_b,sigma_b2,teff);

}



