#include "catch.hpp"
#include "calcem/calcem_util/e_ep_mu_util/e_ep_mu.h"
#include <range/v3/all.hpp>
#include "generalTools/testing.h"
#include <range/v3/all.hpp>

template <typename RangeInt, typename Range, typename Float>
void checkFirstEprime( const RangeInt& jbetaVec, int lat, const Float& enow, const Range& betas, const Range& x, const Float& tev, const RangeInt& finalJbeta, const Range& finalEp ){
  int jbeta;
  Float ep;
  for ( size_t i = 0; i < jbetaVec.size(); ++i ){
    jbeta = jbetaVec[i];
    ep = findFirstEprime( lat, jbeta, enow, betas, x, tev );
    REQUIRE( finalJbeta[i] == jbeta );
    REQUIRE( finalEp[i] == Approx(ep).epsilon(1e-6) );
  }
}




/*

TEST_CASE( "do we need a midpoint" ){ 
  std::vector<double> x(20,0.0), y(20*65,0.0); 
  
  std::vector<double> initialY { 195158.09939262934, -0.87238772679496734, -0.61817875518802279, -0.36545660682035069, -0.11418564421716945, 0.13568062577950979, 0.38420180765640111, 0.63141488906433418, 0.87738226704862166 };

  x[0] = 2.54E-3;
  
  for (size_t i = 0; i < initialY.size(); ++i){ 
    y[i*x.size()+0] = initialY[i]; 
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


*/

/*
TEST_CASE( "do 330" ){ 

  int imax = 20;
  std::vector<double> x(imax,0.0), y(imax*65,0.0), initialY { 194040.3395, 
                     -0.745374226, -0.240023402, 0.259736457, 0.754304965 };
  x[0] = 2.5607298E-3;
  for (size_t i = 0; i < initialY.size(); ++i){ y[i*x.size()+0] = initialY[i]; }

  int i = 2, lat = 0, iinc = 2, lasym = 0;
  double tev = 2.5507297688e-2, tol = 5.0E-2;
  double enow;
  std::vector<double> alphas { 1.1, 2.2, 3.3, 4.5, 5.8 },
                       betas { 0.1, 0.2, 1.3, 1.4, 2.5, 2.6, 3.7 },
                       sab { -0.18259619, -0.30201347, -3.93654779, -3.98809174, 
  -4.33545607, -4.39515402, -5.88934921, -0.76225291, -0.81658341, -3.14161459, 
  -3.30566186, -3.90554652, -3.96233362, -5.23696660, -1.19182884, -1.23155471, 
  -2.79610565, -2.95633099, -3.74989225, -3.80837585, -4.93373911, -1.58342860, 
  -1.61310713, -2.71233943, -2.84291608, -3.69699590, -3.75199349, -4.77433858, 
  -1.96121202, -1.98720663, -2.78454600, -2.88531460, -3.71288120, -3.77142141, 
  -4.71158392 };

  double az = 0.99917, sigma_b = 163.72792237, sigma_b2 = 0.0, teff = 0.120441926;
  int nnl = -5, nl = 5, jbeta = 1,j = 0;


  double ulast, u2last, u3last;
  std::vector<double> correctX(imax), pdf(imax), mu1(imax), mu2(imax), mu3(imax), mu4(imax);


  GIVEN( "Small initial energy (E=1e-5 eV)" ){
    enow = 1e-5;

    THEN( "Returned x values, y vector, and moment values are correct" ){
      auto out = do_330(enow,x,y,i,j,tev,tol,lat,iinc,lasym,alphas,betas,sab,az,sigma_b,sigma_b2,teff,nnl,nl,jbeta);

      ulast  = std::get<0>(out); REQUIRE( 2609.673 == Approx(ulast ).epsilon(1e-6));
      u2last = std::get<1>(out); REQUIRE(-6650.153 == Approx(u2last).epsilon(1e-6));
      u3last = std::get<2>(out); REQUIRE( 358.7471 == Approx(u3last).epsilon(1e-6));

      correctX = { 2.560729E-3, 1.920547E-3, 1.280364E-3, 3.200912E-4, 1.600456E-4, 
      8.002281E-5, 4.001140E-5, 2.000570E-5, 1.000285E-5, 5.001426E-6, 2.500713E-6, 
      1.250356E-6, 6.251782E-7, 3.125891E-7, 1.562945E-7, 7.814728E-8, 3.907364E-8, 
      1.953682E-8, 9.768411E-9, 0.0 };
      pdf = { 194040.3, 213357.4, 234616.8, 270904.3, 278088.3, 281762.2, 
           284435.9, 284849.7, 281681.0, 201641.8, 142771.2, 100595.0, 71050.27, 
           50112.44, 35395.39, 25015.29, 17684.03, 12502.95, 8840.380, 0.0};
      mu1 = {     -0.7453742, -0.7422468, -0.7375231, -0.7158896, -0.6992124, 
      -0.6768329, -0.6449507, -0.6021970, -0.5463158, -0.6019580, -0.6442864, 
      -0.6756062, -0.6972283, -0.7129893, -0.7239285, -0.7315982, -0.7369996, 
      -0.7408188, -0.7435070, 0.00000000 };
      mu2 = {      -0.24002340, -0.23307917, -0.22265027, -0.17496084, -0.13820280, 
      -0.08941394, -0.02005056,  0.07379239,  0.19662745,  0.07431822, -0.01860502, 
      -0.08672420, -0.13381656, -0.16855131, -0.19263496, -0.20951361, -0.22139824, 
      -0.22977689, -0.23570553,  0.0 };
      mu3 = {    0.2597364, 0.2667551, 0.2772268, 0.3250109, 0.3618243, 0.4093513, 
      0.4769398, 0.5740460, 0.6965196, 0.5745724, 0.4783592, 0.4120090, 0.3662356, 
      0.3314537, 0.3073638, 0.2904852, 0.2786010, 0.2702225, 0.2643070, 0.0 };
      mu4 = {    0.7543049, 0.7574966, 0.7623170, 0.7840282, 0.8005880, 0.8224674, 
      0.8530495, 0.8969429, 0.9528127, 0.8971790, 0.8537006, 0.8236852, 0.8025615, 
      0.7869531, 0.7760527, 0.7683948, 0.7629976, 0.7591794, 0.7564790, 0.0 };

      REQUIRE( 1  == i );
      REQUIRE( 20 == j );
      REQUIRE( 2  == jbeta );

      REQUIRE( ranges::equal(x,correctX,equal) );
      std::vector<std::vector<double>> correct { pdf, mu1, mu2, mu3, mu4 };
      for ( size_t i = 0; i < correct.size(); ++i ){
        for ( size_t j = 0; j < pdf.size(); ++j ){
          REQUIRE( correct[i][j] == Approx(y[i*x.size()+j]).epsilon(1e-6) );
        }
      }
      for ( size_t i = correct.size()*x.size(); i < y.size(); ++i){
        REQUIRE( 0.0 == Approx(y[i]).epsilon(1e-6) );
      }

    } // THEN
  } // GIVEN



  GIVEN( "Medium initial energy (E=1e-3 eV)" ){
    enow = 1e-3;

    THEN( "Returned x values, y vector, and moment values are correct" ){

      auto out = do_330(enow,x,y,i,j,tev,tol,lat,iinc,lasym,alphas,betas,sab,az,sigma_b,sigma_b2,teff,nnl,nl,jbeta);

      ulast  = std::get<0>(out); REQUIRE( 2993.879 == Approx(ulast ).epsilon(1e-6));
      u2last = std::get<1>(out); REQUIRE(-240.7431 == Approx(u2last).epsilon(1e-6));
      u3last = std::get<2>(out); REQUIRE( 427.3755 == Approx(u3last).epsilon(1e-6));
  

      correctX = { 2.5607298E-3, 2.5607297E-3, 2.5607295E-3, 6.4018245E-4, 
      1.6004562E-4, 8.0022810E-5, 4.0011405E-5, 2.0005703E-5, 1.0002852E-5, 
      5.0014260E-6, 2.5007130E-6, 1.2503565E-6, 6.2517825E-7, 3.1258913E-7, 
      1.5629457E-7, 7.8147285E-8, 3.9073643E-8, 1.9536822E-8, 9.7684110E-9, 0.0 };
      pdf = { 194040.3, 22928.53, 22928.53, 22004.47, 10421.57, 7271.812, 
      5109.760, 3602.049, 2543.152, 1796.915, 1270.129, 897.9471, 634.8844, 
      448.9098, 317.4197, 224.4469, 158.7070, 112.2224, 79.35316, 0.0 };
      mu1 = {      -0.7453742, -0.6639123, -0.6639123, -0.5966787, -0.6825695, 
      -0.7041852,  -0.7182603, -0.7278000, -0.7343954, -0.7390024, -0.7422431, 
      -0.7445198,  -0.7461270, -0.7472623, -0.7480646, -0.7486316, -0.7490325, 
      -0.7493159,  -0.7495163, 0.0 };
      mu2 = {      -0.24002340, -0.06436183, -0.06436181, 0.08567835, -0.10222948, 
      -0.14962484, -0.18042466, -0.20130111, -0.21574559, -0.22584518, -0.23293962, 
      -0.23794773, -0.24148212, -0.24397868, -0.24574300, -0.24699014, -0.24787184, 
      -0.24849522, -0.24893599, 0.0 };
      mu3 = {    0.2597364, 0.4272364, 0.4272364, 0.5853635, 0.3965669, 0.3496719, 
      0.3192034, 0.2985083, 0.2841580, 0.2741063, 0.2670555, 0.2620485, 0.2585157, 
      0.2560201, 0.2542564, 0.2530096, 0.2521280, 0.2515047, 0.2510640, 0.0 };
      mu4 = {    0.7543049, 0.8233353, 0.8233353, 0.9011784, 0.8153225, 0.7947936, 
      0.7812307, 0.7719456, 0.7654775, 0.7609340, 0.7577044, 0.7554555, 0.7538610, 
      0.7527317, 0.7519323, 0.7513667, 0.7509666, 0.7506835, 0.7504834, 0.0 };

      REQUIRE( 1  == i );
      REQUIRE( 43 == j );
      REQUIRE( 2  == jbeta );


      REQUIRE( ranges::equal(x,correctX,equal) );
      std::vector<std::vector<double>> correct { pdf, mu1, mu2, mu3, mu4 };
      for ( size_t i = 0; i < correct.size(); ++i ){
        for ( size_t j = 0; j < pdf.size(); ++j ){
          REQUIRE( correct[i][j] == Approx(y[i*x.size()+j]).epsilon(1e-6) );
        }
      }
      for ( size_t i = correct.size()*x.size(); i < y.size(); ++i){
        REQUIRE( 0.0 == Approx(y[i]).epsilon(1e-6) );
      }
  
    } // THEN 
  } // GIVEN

} // TEST CASE

*/

/*

TEST_CASE( "313" ) {
  int lat = 1, jbeta = -7;
  double enow, ep, tev = 2.5507297688E-2;
  std::vector<double> betas { 0.1, 0.2, 1.3, 1.4, 2.5, 2.6, 3.7 }, x(20,0.0);


  GIVEN( "Various jbeta values, fixed enow" ){ 
    std::vector<int> jbetaVec(15), finalJbeta(15);
    std::vector<double> finalEp(15);

    jbetaVec    = ranges::view::iota(-7,8);

    enow = 1e-5;
    finalJbeta = {1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 3, 4, 5, 6, 7};
    finalEp    = {2.54e-3, 2.54e-3, 2.54e-3, 2.54e-3, 2.54e-3, 2.54e-3, 2.54e-3, 
                   2.54e-3, 2.54e-3, 5.07e-3, 3.29e-2, 3.543e-2, 6.326e-2, 
                   6.579e-2, 9.362e-2};
    checkFirstEprime( jbetaVec, lat, enow, betas, x, tev, finalJbeta, finalEp );


    enow = 1e-3;
    finalJbeta = {1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 3, 4, 5, 6, 7};
    finalEp    = {3.53e-3, 3.53e-3, 3.53e-3, 3.53e-3, 3.53e-3, 3.53e-3, 3.53e-3, 
                  3.53e-3, 3.53e-3, 6.06e-3, 3.389e-2, 3.642e-2, 6.425e-2, 
                  6.678e-2, 9.461e-2};
    checkFirstEprime( jbetaVec, lat, enow, betas, x, tev, finalJbeta, finalEp );


    enow = 1e-1;
    finalJbeta = {-7, -6, -5, -4, -3, -2, -1, 1, 1, 2, 3, 4, 5, 6, 7 };
    finalEp    = {6.39000153e-3, 3.4220001e-2, 3.6750001e-2, 6.4580001e-2, 
                  6.7110001e-2, 9.4940001e-2, 9.7470001e-2, 0.10253,  0.10253, 
                  0.10506, 0.13289, 0.13542, 0.16325, 0.16578, 0.19361 };
    checkFirstEprime( jbetaVec, lat, enow, betas, x, tev, finalJbeta, finalEp );

  } // GIVEN



} // TEST CASE

*/









TEST_CASE( "do 330 (and some things around it)" ){ 

  int imax = 20;
  std::vector<double> x(imax,0.0), y(imax*65,0.0), initialY { 194040.3395, 
                     -0.745374226, -0.240023402, 0.259736457, 0.754304965 };
  x[0] = 2.5607298E-3;
  for (size_t i = 0; i < initialY.size(); ++i){ y[i*x.size()+0] = initialY[i]; }

  int lat = 0, iinc = 2, lasym = 0;
  double tev = 2.5507297688e-2, tol = 5.0E-2;
  double enow = 1e-5;
  std::vector<double> alphas { 1.1, 2.2, 3.3, 4.5, 5.8 },
                       betas { 0.1, 0.2, 1.3, 1.4, 2.5, 2.6, 3.7 },
                       sab { -0.18259619, -0.30201347, -3.93654779, -3.98809174, 
  -4.33545607, -4.39515402, -5.88934921, -0.76225291, -0.81658341, -3.14161459, 
  -3.30566186, -3.90554652, -3.96233362, -5.23696660, -1.19182884, -1.23155471, 
  -2.79610565, -2.95633099, -3.74989225, -3.80837585, -4.93373911, -1.58342860, 
  -1.61310713, -2.71233943, -2.84291608, -3.69699590, -3.75199349, -4.77433858, 
  -1.96121202, -1.98720663, -2.78454600, -2.88531460, -3.71288120, -3.77142141, 
  -4.71158392 };

  double az = 0.99917, sigma_b = 163.72792237, sigma_b2 = 0.0, teff = 0.120441926;
  int nnl = -5, nl = 5, jbeta = 1,j = 0;


  do_330_extra(enow,x,y,j,tev,tol,lat,iinc,lasym,alphas,betas,sab,az,sigma_b,sigma_b2,teff,nnl,nl,jbeta);




} // TEST CASE






















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

