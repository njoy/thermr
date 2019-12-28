#define CATCH_CONFIG_MAIN
#include "catch.hpp"
#include "calcem/calcem_util/e_ep_mu_util/sigl_util/beginningLoop.h"
#include <range/v3/all.hpp>
#include "generalTools/testing.h"


//auto equal = [](auto x, auto y, double tol = 1e-6){return x == Approx(y).epsilon(tol);};
auto equal2 = [](auto x, auto y, double tol = 1e-5){return x == Approx(y).epsilon(tol);};

/*
TEST_CASE( "110" ){
  std::vector<double> x(20,0.0), y(20,0.0);
  x[0] =  1.0; x[2] = -1.0;
  y[0] =  2.5; y[2] =  4.0;
  
  int i = 3;
  double e = 1e-5, ep = 1e-4, tev = 2.55e-2;
  
  std::vector<double> alpha { 1.1, 2.2, 3.3, 4.5, 5.8 },
                       beta { 0.1, 0.2, 1.3, 1.4, 2.5, 2.6, 3.7 };
  std::vector<double> sab(alpha.size()*beta.size());
  for ( size_t i = 0; i < alpha.size(); ++i ){
    for ( size_t j = 0; j < beta.size(); ++j ){
      sab[i*beta.size()+j] = 0.01*((j+1) + 0.1*(i+1));
    } 
  } 
  double az = 0.99917, tevz = 2.53e-2, sigma_b = 163.72792237360667, 
         sigma_b2 = 0.0, teff = 0.120441926577313, tol = 2.5e-2, ymax = 1e-3;
  int lasym = 0, lat = 1, iinc = 2; 

  do_110(i,x,y,e,ep,tev,alpha,beta,sab,az,tevz,lasym,lat,sigma_b,sigma_b2,teff,iinc,tol,ymax);
  std::vector<double> 
    correct_x { 1.0, 0.0, -0.5, -0.75, -0.875, -0.9375, -0.96875, -0.984375, 
    -0.9921875, -0.99609375, -0.99804687, -0.99902343, -0.99951171, -0.99975585, 
    -0.99987792, -0.99993896, -0.99996948, -0.99998474, -0.99999237, -1.0},
    correct_y { 2.5, 0.0, 143662.33773, 136248.865859, 132948.05965, 
    131385.1453676, 130624.0785461, 130248.4714971, 130061.878966, 129968.882923, 
    129922.4597632, 129899.2668299, 129887.6750200, 129881.880278, 129878.983198, 
    129877.5346130, 129876.8103382, 129876.4482054, 129876.267140, 4.0};

  REQUIRE(ranges::equal(correct_x, x, equal));
  REQUIRE(ranges::equal(correct_y, y, equal));


}
*/



/*
TEST_CASE( "170-180" ){
  double fract = 312.0, sum = 0.0, muLeft = -1.0, xsLeft = 1252.9281013622765, xil = 2.0100502512564842; 
  int i = 4, j = 0;

  std::vector<double> 
    muVec  { 1.0, 0.99, -5.0E-3, -0.5025, -1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 },
    xsVec { 1250.5627281224217, 1250.5869704560703, 1252.1124859871622, 1252.5658596359865, 1252.9281013622765, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};

  WHEN( "170 -> leave " ){
    THEN( "nothing changed" ){ 
      auto xn = do_170_175_180( fract, sum, xsVec, muVec, xsLeft, muLeft, i, j, xil );
      REQUIRE( fract   == Approx(312.0).epsilon(1e-6) );
      REQUIRE( sum     == Approx(  0.0).epsilon(1e-6) );
      REQUIRE( muLeft  == Approx( -1.0).epsilon(1e-6) );
      REQUIRE( xsLeft  == Approx(1252.92810136).epsilon(1e-6) );
      REQUIRE( i == 4 );
      REQUIRE( j == 1 ); 
      REQUIRE(xn == Approx(-0.7509833168).epsilon(1e-6));
    } // THEN
  } // WHEN 


  WHEN( "170 -> 175 -> 180 -> leave " ){
    xsLeft *= 0.001;
    THEN( "nothing changed" ){ 
      auto xn = do_170_175_180( fract, sum, xsVec, muVec, xsLeft, muLeft, i, j, xil );
      REQUIRE( fract   == Approx(312.0).epsilon(1e-6) );
      REQUIRE( sum     == Approx(  0.0).epsilon(1e-6) );
      REQUIRE( muLeft  == Approx( -1.0).epsilon(1e-6) );
      REQUIRE( xsLeft  == Approx(1.25292810136).epsilon(1e-6) );
      REQUIRE( i == 4 );
      REQUIRE( j == 1 ); 
      REQUIRE(xn == Approx( -0.5025).epsilon(1e-6));
    } // THEN
  } // WHEN 
 
} // TEST CASE

*/









    







TEST_CASE( "dothings" ){
  double e = 1e-2, ep = 1e-3, tev = 0.025, tolin = 5e-2, az = 0.99917, 
         sigma_b = 4.0, sigma_b2 = 0.0, teff = 0.12;
  int nL = -9, lat = 1, iinc = 2, lasym = 0;
  std::vector<double> alphas { 1.1, 2.2, 3.3, 4.5, 5.8 },
                       betas { 0.1, 0.2, 1.3, 1.4, 2.5, 2.6, 3.7 };
  std::vector<double> sab(alphas.size()*betas.size());
  for ( size_t i = 0; i < alphas.size(); ++i ){
    for ( size_t j = 0; j < betas.size(); ++j ){
      sab[i*betas.size()+j] = 0.01*((j+1) + 0.1*(i+1));
    } 
  } 

  std::vector<double> correct_s, s;

  /*
  */
  ep = 1e-1;
  s = sigl(ep,e,tev,tolin,nL,lat,iinc,alphas,betas,sab,az,lasym,sigma_b,sigma_b2,teff);
  correct_s = {71.844290065748709, -0.48647654336295154, -0.28577865744470676, -8.8254959343803965E-2, 0.10934278668830079, 0.30704325144135147, 0.50477104242836623, 0.70270321344353548, 0.90176859041415269, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
  REQUIRE(ranges::equal(correct_s, s, equal));


  ep = 1e-2;


  //e = 9.9999997764825821E-3;
  //ep = 9.9999997764825821E-3;
  //sigma_b = 4.0;
  //teff = 0.11999999731779099;
  //tev = 2.5000000372529030E-002;


  s = sigl(ep,e,tev,tolin,nL,lat,iinc,alphas,betas,sab,az,lasym,sigma_b,sigma_b2,teff);
  correct_s = { 220.16164823349109, -0.14479545028098159, 0.16598934410207883, 0.40051344934325717, 0.59855092490738082, 0.75607791889432496, 0.87463820097482681, 0.95401954612035467, 0.99332099672089436, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
  REQUIRE(ranges::equal(correct_s, s, equal));

  ep = 1e-3;
  correct_s = {61.173606284237110, -0.87503424644581240, -0.62526004804699387, -0.37535265378315336, -0.12546954443858765, 0.12447684514433568, 0.37451189318707867, 0.62447532003983108, 0.87522679160684391, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
  s = sigl(ep,e,tev,tolin,nL,lat,iinc,alphas,betas,sab,az,lasym,sigma_b,sigma_b2,teff);
  REQUIRE(ranges::equal(correct_s, s, equal));


  ep = 1e-4;
  s = sigl(ep,e,tev,tolin,nL,lat,iinc,alphas,betas,sab,az,lasym,sigma_b,sigma_b2,teff);
  correct_s = {19.661973513357161, -0.87500800464553019, -0.62509139062387742, -0.37511937787300242, -0.12515611411904232,  0.12483210150811240,  0.37484817883993737, 0.62485204314697129, 0.87507351378686860, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
  REQUIRE(ranges::equal(correct_s, s, equal));

  ep = 1e-5;
  s = sigl(ep,e,tev,tolin,nL,lat,iinc,alphas,betas,sab,az,lasym,sigma_b,sigma_b2,teff);
  correct_s = {6.2278467641892856, -0.87500211825460317, -0.62502992284537084, -0.37503857945425939, -0.12505006715664307, 0.12494679820216802, 0.37495231013203850, 0.62495468630262829, 0.87502358904232658, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
  REQUIRE(ranges::equal(correct_s, s, equal));

  /*

  */


  //std::cout << std::setprecision(15) << ep << std::endl;
  //std::cout << (s|ranges::view::all) << std::endl;

  


}
/*
*/















/*
TEST_CASE( "110 120 130" ){
  std::vector<double> x(20,0.0);
  std::vector<double> y(20,0.0);
  x[0] = 1.0; x[1] = 0.99; x[2] = -1.0;
  y[0] = 1.35700e5; y[1] = 1.35701e5; y[2] = 1.35809e5;

  std::vector<double> alpha { 1.1, 2.2, 3.3, 4.5, 5.8 };
  std::vector<double> beta { 0.1, 0.2, 1.3, 1.4, 2.5, 2.6, 3.7 };
    

  std::vector<double> sab(alpha.size()*beta.size());
  for ( size_t i = 0; i < alpha.size(); ++i ){
    for ( size_t j = 0; j < beta.size(); ++j ){
      sab[i*beta.size()+j] = 0.01*((j+1) + 0.1*(i+1));
    } 
  } 

  int lasym = 0, lat = 1, iinc = 2, nlmax = 65, nl = 10, i = 3, 
      j = 0, nbin = 8;

  double e = 1.0e-6, ep = 1.2e-4, tev = 1.5e-4, bbm = 0.0, az = 11.9,
    tevz = 2.2e-4, sb = 5.53, sb2 = 0.0,
    teff = 6.14e-2, tolin = 5e-2, sigmin = 1.0e-32, 
    s1bb = 1.1369180380, tol = 2.5e-2, xtol = 1.0e-5, 
    seep = 91200.0, yl = 13500, ymax = 13500, fract = 0.0, xl = -1.0, 
    eps = -1.0e-3;


  GIVEN( "inputs 1" ){

    THEN( "110-->110, 110-->120, 120-->110, 120-->130, not many iterations" ){


      auto out = do_110_120_130_for_sigl( i, x, y, e, ep, tev, tevz, alpha, 
        beta, sab, az, lasym, teff, lat, sb, sb2, iinc, nl, 
        sigmin, s, nbin, fract, xl, j, ymax, yl, tol, xtol );



      //ymax = adaptiveLinearization( x, y, e, ep, tev, tevz, alpha, beta, sab, 
      //  bbm, az, lasym, teff, lat, sb, sb2, iinc, eps, seep, s1bb );

      double gral = std::get<0>(out);
      double sum = std::get<1>(out);
      REQUIRE( 0 == Approx(gral).epsilon(1e-6) ); 
      REQUIRE( 0 == Approx(sum).epsilon(1e-6) ); 

      REQUIRE( 3 == Approx(i).epsilon(1e-6) ); 
      REQUIRE( 0 == Approx(j).epsilon(1e-6) ); 
      REQUIRE( 9 == Approx(nbin).epsilon(1e-6) ); 
      REQUIRE(-1 == Approx(xl).epsilon(1e-6) );
      REQUIRE( 30174.6306224 == Approx(fract).epsilon(1e-6) );
      REQUIRE( 135700.0 == Approx(yl).epsilon(1e-6) );
      //REQUIRE( 135829.6496457 == Approx(ymax).epsilon(1e-6) );

      std::vector<double> correctX = { 1.0, 0.99, -1.0, -0.005, -1.0 },
        correctY = { 135757.913, 135758.3455, 135829.6496, 135797.21918, 135809.0 };

      for ( size_t i = 0; i < x.size(); ++i ){ 
        if ( i < 5 ){ REQUIRE( correctX[i] == Approx(x[i]).epsilon(1e-6) ); }
        else        { REQUIRE( 0.0         == Approx(x[i]).epsilon(1e-6) ); }
      }


      for ( size_t i = 0; i < y.size(); ++i ){ 
        if ( i < 5 ){ REQUIRE( correctY[i] == Approx(y[i]).epsilon(1e-6) ); }
        else        { REQUIRE( 0.0         == Approx(y[i]).epsilon(1e-6) ); }
      }

      REQUIRE( 271571.6756021== Approx(s[0]).epsilon(1e-6) );
      for ( size_t i = 1; i < s.size(); ++i ){ 
        REQUIRE( 0.0 == Approx(s[i]).epsilon(1e-6) ); 
      }


    } // THEN
  } // GIVEN

*/

  /*

  GIVEN( "inputs 2" ){
    y[0] = 1.00000e5; y[1] = 1.00001e5; y[2] = 1.00009e5;

    e = 1.0e-3; ep = 1.2e-2; tev = 1.5e-1;
      tevz = 2.2e-1; 
      teff = 6.14e-0; 


    THEN( "110-->110, 110-->120, 120-->110, 120-->120, 120-->130, many iterations" ){


      auto out = do_110_120_130_for_sigl( i, x, y, e, ep, tev, tevz, alpha, 
        beta, sab, az, lasym, teff, lat, sb, sb2, iinc, nl, 
        sigmin, s, nbin, fract, xl, j, ymax, yl, tol, xtol );



      //ymax = adaptiveLinearization( x, y, e, ep, tev, tevz, alpha, beta, sab, 
      //  bbm, az, lasym, teff, lat, sb, sb2, iinc, eps, seep, s1bb );

      double gral = std::get<0>(out);
      double sum = std::get<1>(out);
      REQUIRE( 0 == Approx(gral).epsilon(1e-6) ); 
      REQUIRE( 0 == Approx(sum).epsilon(1e-6) ); 

      REQUIRE( 3 == Approx(i).epsilon(1e-6) ); 
      REQUIRE( 0 == Approx(j).epsilon(1e-6) ); 
      REQUIRE( 9 == Approx(nbin).epsilon(1e-6) ); 
      REQUIRE(-1 == Approx(xl).epsilon(1e-6) );
      REQUIRE( 120.57844468516407 == Approx(fract).epsilon(1e-6) );
      REQUIRE( 100000.0 == Approx(yl).epsilon(1e-6) );
      //REQUIRE( 538.71588696219601 == Approx(ymax).epsilon(1e-6) );

      std::vector<double> correctX = { 1.0, 0.99, -1.0, 0.99125, 0.990625, 
        0.9903125, 0.99015625, 0.99007813, 0.99003907, 0.99001954, 0.99000977, 
        0.99, -0.99902832, -0.99951416, -0.99975708, -0.99987854, -0.99993927, 
        -0.99996963, -0.99998481, -1.0 },
      correctY = { 461.24225912, 464.2445488, 538.7158869, 463.8735389, 
        464.0591955, 464.151910, 464.198238, 464.2213947, 464.23297, 464.238758, 
        464.2416537, 100001.0, 538.74710638, 538.7314976, 538.72369255, 
        538.7197898, 538.71783, 538.7168628, 538.71637506, 100009.0 };

      for ( size_t i = 0; i < x.size(); ++i ){ 
        REQUIRE( correctX[i] == Approx(x[i]).epsilon(1e-6) );
      }


      for ( size_t i = 0; i < y.size(); ++i ){ 
        REQUIRE( correctY[i] == Approx(y[i]).epsilon(1e-6) );
      }

      REQUIRE( 1085.2060021664768 == Approx(s[0]).epsilon(1e-6) );
      for ( size_t i = 1; i < s.size(); ++i ){ 
        REQUIRE( 0.0 == Approx(s[i]).epsilon(1e-6) ); 
      }


    } // THEN
  } // GIVEN

  GIVEN( "inputs 2" ){
    x[0] = 1.00; x[1] = 0.99; x[2] = -1.00;
    y[0] = 0.00; y[1] = 0.00; y[2] = 0.00;

    e = 1.0e-3; ep = 1.2e-2; tev = 1.5e-5;
      tevz = 2.2e-1; 
      teff = 6.14e-3; 
      xl =-1.0;
      yl = 0.0;
      i = 3;
      ymax = 1e-3;


    THEN( "110-->110, 110-->120, 120-->110, 120-->120, 120-->130, many iterations" ){


      //auto out = do_110_120_130_for_sigl(i, x, y, e, ep, tev, tevz, alpha, beta, sab, bbm, az, 
      //lasym, teff, lat, sb, sb2, iinc, nl, sigmin, s, 
      //  nbin, fract, xl, j, ymax, eps, seep, yl, s1bb, tol, xtol);

      auto out = do_110_120_130_for_sigl( i, x, y, e, ep, tev, tevz, alpha, 
        beta, sab, az, lasym, teff, lat, sb, sb2, iinc, nl, 
        sigmin, s, nbin, fract, xl, j, ymax, yl, tol, xtol );



      double gral = std::get<0>(out);
      double sum = std::get<1>(out);
      REQUIRE( 0 == Approx(gral).epsilon(1e-6) ); 
      REQUIRE( 0 == Approx(sum).epsilon(1e-6) ); 

      REQUIRE( 0 == Approx(i).epsilon(1e-6) ); 
      REQUIRE( 0 == Approx(j).epsilon(1e-6) ); 
      REQUIRE( 8 == Approx(nbin).epsilon(1e-6) ); 
      REQUIRE( 1.0 == Approx(xl).epsilon(1e-6) );
      REQUIRE( 0.0 == Approx(fract).epsilon(1e-6) );
      REQUIRE( 0 == Approx(yl).epsilon(1e-6) );
      REQUIRE( 1e-3 == Approx(ymax).epsilon(1e-6) );
      

      std::vector<double> correctX = { 1.0, 0.99, 0.4925, -5.0E-3, -1.0 };

      for ( size_t i = 0; i < x.size(); ++i ){ 
        if ( i < correctX.size() ){
          REQUIRE( correctX[i] == Approx(x[i]).epsilon(1e-6) );
        }
        else {
          REQUIRE( 0.0 == Approx(x[i]).epsilon(1e-6) );
        }
      }


      for ( size_t i = 0; i < y.size(); ++i ){ 
        REQUIRE( 0.0 == Approx(y[i]).epsilon(1e-6) );
      }

      for ( size_t i = 0; i < s.size(); ++i ){ 
        REQUIRE( 0.0 == Approx(s[i]).epsilon(1e-6) ); 
      }


    } // THEN
  } // GIVEN

} // TEST CASE

  */

