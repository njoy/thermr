#include "catch.hpp"
#include "calcem/calcem_util/e_ep_mu_util/e_ep_mu.h"
#include <range/v3/all.hpp>
#include "generalTools/testing.h"
#include <range/v3/all.hpp>




TEST_CASE( "do we need a midpoint" ){ 
  std::vector<double> x(20,0.0), y(20*65,0.0); 
  
  std::vector<double> beginningOfY { 195158.09939262934, -0.87238772679496734, -0.61817875518802279, -0.36545660682035069, -0.11418564421716945, 0.13568062577950979, 0.38420180765640111, 0.63141488906433418, 0.87738226704862166 };

  x[0] = 2.54E-3;
  
  for (size_t i = 0; i < beginningOfY.size(); ++i){ 
    y[i*x.size()+0] = beginningOfY[i]; 
  }



  double xm = 1.27e-3;
  std::vector<double> s { 235586.48147822541, -0.86806816809249565, -0.60677324533223020, -0.34926974510882292, -9.5577908119815064E-002, 0.15436521665672856, 0.40054354740949549, 0.64302993119542795, 0.88181200416435501, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };

  int i = 2;
  int nl = 9;
  double tol = 1.0;
  REQUIRE( needMidpoint(x,y,xm,i,nl,s,tol) == false );
  tol = 1e-2;
  REQUIRE( needMidpoint(x,y,xm,i,nl,s,tol) == true  );
  tol = 5e-1;
  REQUIRE( needMidpoint(x,y,xm,i,nl,s,tol) == true  );
  tol = 0.56;
  REQUIRE( needMidpoint(x,y,xm,i,nl,s,tol) == true  );
  tol = 0.59;
  REQUIRE( needMidpoint(x,y,xm,i,nl,s,tol) == false );


} // TEST CASE



TEST_CASE( "do 330" ){ 
  std::vector<double> x(20,0.0), y(20*65,0.0); 
  
  std::vector<double> beginningOfY { 194040.33954861510, -0.74537422649273466, -0.24002340213077786, 0.25973645751508195, 0.75430496529005719 };

  x[0] = 2.5607298E-3;
  double enow = 1e-5;
  
  for (size_t i = 0; i < beginningOfY.size(); ++i){ 
    y[i*x.size()+0] = beginningOfY[i]; 
  }

  int i = 2, lat = 0, iinc = 2, lasym = 0;
  double tev = 2.5507297687999999E-2, tol = 5.0000000000000003E-002;
  std::vector<double> alphas { 1.1, 2.2, 3.3, 4.5, 5.8 },
                       betas { 0.1, 0.2, 1.3, 1.4, 2.5, 2.6, 3.7 },
                       sab {  
-0.18259619450343811, 
 -0.30201347821403990, 
  -3.9365477946016552, 
  -3.9880917447839432, 
  -4.3354560748618214, 
  -4.3951540230337844, 
  -5.8893492113551460, 
 -0.76225291555860053, 
 -0.81658341547833779, 
  -3.1416145973110434, 
  -3.3056618636120363, 
  -3.9055465231884647, 
  -3.9623336249041428, 
  -5.2369666084323878, 
  -1.1918288416393834, 
  -1.2315547119120109, 
  -2.7961056502182311, 
  -2.9563309997796829, 
  -3.7498922512606674, 
  -3.8083758567455739, 
  -4.9337391161748112, 
  -1.5834286060899199, 
  -1.6131071358084847, 
  -2.7123394367355611, 
  -2.8429160837626020, 
  -3.6969959035560183, 
  -3.7519934969692779, 
  -4.7743385860598142, 
  -1.9612120279318532, 
  -1.9872066311269845, 
  -2.7845460064772798, 
  -2.8853146068676256, 
  -3.7128812054039422, 
  -3.7714214156555439, 
  -4.7115839227552696 };
  double az =   0.99917000000000000, sigma_b = 163.72792237360667, sigma_b2 =  0.0000000000000000, teff =   0.12044192657731301;
  int nnl = -5, nl = 5;
  std::vector<double> p2(118,0.0), p3(118,0.0), xsi(99,0.0);
  double xlast = 0.0, ylast = 0.0;

  /*
  std::cout << y[0*x.size()+0] << std::endl;
  std::cout << y[0*x.size()+1] << std::endl;
  std::cout << y[0*x.size()+2] << std::endl;
  std::cout << std::endl;
  std::cout << y[1*x.size()+0] << std::endl;
  std::cout << y[1*x.size()+1] << std::endl;
  std::cout << y[1*x.size()+2] << std::endl;
  */

  do_330(enow,x,y,i,tev,tol,lat,iinc,lasym,alphas,betas,sab,az,sigma_b,sigma_b2,teff,nnl,nl);

} // TEST CASE




/*

TEST_CASE( "313" ) {
  int lat = 1, jbeta = -7, iskip = 0;
  double enow = 1e-5, ep = 0.0, tev = 2.5507297688E-2;
  std::vector<double> betas { 0.1, 0.2, 1.3, 1.4, 2.5, 2.6, 3.7 }, x(20,0.0);


  GIVEN( "Various jbeta values, fixed enow" ){ 
    std::vector<int> jbetaVec(15), final_jbeta(15);
    std::vector<double> final_ep(15);

    jbetaVec    = ranges::view::iota(-7,8);
    final_jbeta = {1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 3, 4, 5, 6, 7};
    final_ep    = {2.54e-3, 2.54e-3, 2.54e-3, 2.54e-3, 2.54e-3, 2.54e-3, 2.54e-3, 
                   2.54e-3, 2.54e-3, 5.07e-3, 3.29e-2, 3.543e-2, 6.326e-2, 
                   6.579e-2, 9.362e-2};
    enow = 1e-5;
    for ( size_t i = 0; i < jbetaVec.size(); ++i ){
      jbeta = jbetaVec[i];
      ep = do_313( lat, jbeta, enow, betas, x, tev );
      REQUIRE( final_jbeta[i] == jbeta );
      REQUIRE( final_ep[i] == Approx(ep).epsilon(1e-6) );
    }


    final_jbeta = {1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 3, 4, 5, 6, 7};
    final_ep    = {3.53e-3, 3.53e-3, 3.53e-3, 3.53e-3, 3.53e-3, 3.53e-3, 3.53e-3, 
                   3.53e-3, 3.53e-3, 6.06e-3, 3.389e-2, 3.642e-2, 6.425e-2, 
                   6.678e-2, 9.461e-2};
    enow = 1e-3;
    for ( size_t i = 0; i < jbetaVec.size(); ++i ){
      jbeta = jbetaVec[i];
      ep = do_313( lat, jbeta, enow, betas, x, tev );
      REQUIRE( final_jbeta[i] == jbeta );
      REQUIRE( final_ep[i] == Approx(ep).epsilon(1e-6) );
    } 

    final_jbeta = {-7, -6, -5, -4, -3, -2, -1, 1, 1, 2, 3, 4, 5, 6, 7 };
    final_ep    = {6.39000153e-3, 3.4220001e-2, 3.6750001e-2, 6.4580001e-2, 
                   6.7110001e-2, 9.4940001e-2, 9.7470001e-2, 0.10253,  0.10253, 
                   0.10506, 0.13289, 0.13542, 0.16325, 0.16578, 0.19361 };
    enow = 1e-1;
    for ( size_t i = 0; i < jbetaVec.size(); ++i ){
      jbeta = jbetaVec[i];
      //std::cout << jbeta << std::endl;
      ep = do_313( lat, jbeta, enow, betas, x, tev );
      REQUIRE( final_jbeta[i] == jbeta );
      REQUIRE( final_ep[i] == Approx(ep).epsilon(1e-6) );
    } 

  } // GIVEN



} // TEST CASE
*/




/*

TEST_CASE( "E-E'-mu" ){
  REQUIRE( true );

  std::vector<double> alphas { 1.1, 2.2, 3.3, 4.5, 5.8 },
                       betas { 0.1, 0.2, 1.3, 1.4, 2.5, 2.6, 3.7 };
  std::vector<double> sab(alphas.size()*betas.size());
  for ( size_t i = 0; i < alphas.size(); ++i ){
    for ( size_t j = 0; j < betas.size(); ++j ){
      sab[i*betas.size()+j] = 0.01*((j+1) + 0.1*(i+1));
    } 
  } 
  double az = 0.99917, tol = 5e-2;
  int lasym = 0, lat = 1, iinc = 2, jmax = 55550;
  int nne = 88, nnl = -9, nl = 9;
  double T = 296.0, teff = 1397.671, teff2 = 0.0, sigma_b = 163.72792237360667, sigma_b2 = 0.0;

  std::cout.precision(15);

  e_ep_mu( T, teff, teff2, jmax, nne, nnl, nl, tol, sigma_b, sigma_b2, az, lasym, lat, iinc, alphas, betas, sab );


}
*/
