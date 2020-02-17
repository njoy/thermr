#include "catch.hpp"
#include "calcem/calcem_util/e_mu_ep.h"
#include <range/v3/all.hpp>
#include "generalTools/testing.h"
#include <range/v3/all.hpp>


TEST_CASE( "main E mu E' function" ){
  int imax = 20, lat = 0, iinc = 2, lasym = 0;
  double tev = 2.5507297688e-2, tol = 2.0;
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
  std::vector<double> eVec;

  GIVEN( "temperature is reasonable (doesn't need scaling)" ){
    eVec = {1E-5, 1.78E-5, 2.5E-5, 3.5E-5, 5E-5, 7E-5, 1E-4, 1.26E-4, 1.6E-4, 2E-4};
    e_mu_ep(eVec,tev,tol,lat,iinc,lasym,alphas,betas,sab,az,sigma_b,sigma_b2,teff);
  } // GIVEN 
} // TEST CASE 






