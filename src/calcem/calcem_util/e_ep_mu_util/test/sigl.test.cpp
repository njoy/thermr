#define CATCH_CONFIG_MAIN
#include "catch.hpp"
#include "calcem/calcem_util/e_ep_mu_util/sigl.h"
#include <range/v3/all.hpp>
#include "generalTools/testing.h"






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
         sigma_b2 = 0.0, teff = 0.120441926577313, tol = 2.5e-2;
  int lasym = 0, lat = 1, iinc = 2; 

  double xsMax = maxOf4Vals( y[0], y[1], y[2], 0.001 );
  addMidpointsRight(i,x,y,e,ep,tev,alpha,beta,sab,az,lasym,lat,sigma_b,sigma_b2,teff,iinc,tol,xsMax);
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




TEST_CASE( "Get pdf value" ){
  double e, tev = 0.025, tolin = 5e-2*0.5, az = 0.99917, sigma_b = 4.0, sigma_b2 = 0.0, teff = 0.12;
  int lat = 1, iinc = 2, lasym = 0;
  bool equiprobableBins = true;
  std::vector<double> alphas { 1.1, 2.2, 3.3, 4.5, 5.8 },
                       betas { 0.1, 0.2, 1.3, 1.4, 2.5, 2.6, 3.7 };
  std::vector<double> sab(alphas.size()*betas.size());
  for ( size_t i = 0; i < alphas.size(); ++i ){
    for ( size_t j = 0; j < betas.size(); ++j ){
      sab[i*betas.size()+j] = 0.01*((j+1) + 0.1*(i+1));
    } 
  } 

  std::vector<double> finalEnergies, correctPDFVals;

  e = 1.0;
  finalEnergies  = { 1e-2, 0.1, 0.8, 1.0, 2.0 },
  correctPDFVals = { 0.316342, 0.8001197, 1.0023345, 3.0989839, 0.0 };
  for ( size_t i = 0; i < finalEnergies.size(); ++i ){
    REQUIRE( getPDF(finalEnergies[i],e,tev,tolin,lat,iinc,alphas,betas,sab,az,
     lasym,sigma_b,sigma_b2,teff) == Approx(correctPDFVals[i]).epsilon(1e-6) );
  }

  e = 1e-2;
  correctPDFVals = { 71.84429, 220.1616, 61.17360, 19.66197, 6.227846 };
  finalEnergies  = { 1e-1, 1e-2, 1e-3, 1e-4, 1e-5 };
  for ( size_t i = 0; i < finalEnergies.size(); ++i ){
    REQUIRE( getPDF(finalEnergies[i],e,tev,tolin,lat,iinc,alphas,betas,sab,az,
     lasym,sigma_b,sigma_b2,teff) == Approx(correctPDFVals[i]).epsilon(1e-6) );
  }

  e  = 1e-1; 
  correctPDFVals = { 262.937411, 1.02706201, 0.32551817, 0.10296110 };
  finalEnergies = { 1e-2, 1e-3, 1e-4, 1e-5 };
  for ( size_t i = 0; i < finalEnergies.size(); ++i ){
    REQUIRE( getPDF(finalEnergies[i],e,tev,tolin,lat,iinc,alphas,betas,sab,az,
     lasym,sigma_b,sigma_b2,teff) == Approx(correctPDFVals[i]).epsilon(1e-6) );
  }




} // TEST CASE


TEST_CASE( "sigl" ){

  GIVEN ( "equiprobable angle bins are requested" ){
  double e = 1e-2, ep, tev = 0.025, tolin = 5e-2, az = 0.99917, 
         sigma_b = 4.0, sigma_b2 = 0.0, teff = 0.12;
  int lat = 1, iinc = 2, lasym = 0, nbin = 8;
  bool equiprobableBins = true;
  std::vector<double> alphas { 1.1, 2.2, 3.3, 4.5, 5.8 },
                       betas { 0.1, 0.2, 1.3, 1.4, 2.5, 2.6, 3.7 };
  std::vector<double> sab(alphas.size()*betas.size());
  for ( size_t i = 0; i < alphas.size(); ++i ){
    for ( size_t j = 0; j < betas.size(); ++j ){
      sab[i*betas.size()+j] = 0.01*((j+1) + 0.1*(i+1));
    } 
  } 

  std::vector<double> correct_s, s;

  WHEN ( "2 bins requested" ){
      std::vector<std::vector<double>> correctMuVecs;
      std::vector<double> finalEnergies;

      nbin = 2;

      e = 1.0;
      finalEnergies = { 1e-2, 1e-1, 5e-1, 1.0, 2.0 };
      correctMuVecs = { 
        {-0.289542, 0.581023 },
        {-0.458477, 0.527666 },
        { 0.300582, 0.783071 },
        { 0.914921, 0.999029 },
        { 0.000000, 0.000000 } };
      for ( size_t i = 0; i < correctMuVecs.size(); ++i ){
        s = sigl(finalEnergies[i],e,tev,tolin,lat,iinc,alphas,betas,sab,az,
                 lasym,sigma_b,sigma_b2,teff,nbin,equiprobableBins);
      }

      e = 1e-2;
      correctMuVecs = {
        {  0.0000000, 0.0000000 },
        { -0.1877993, 0.6040790 },
        {  0.2550645, 0.8945141 },
        { -0.5002791, 0.4996727 },
        { -0.5000937, 0.4999014 },
        { -0.5000301, 0.4999693 } };
      finalEnergies = { 1e0, 1e-1, 1e-2, 1e-3, 1e-4, 1e-5 };
      for ( size_t i = 0; i < correctMuVecs.size(); ++i ){
        ep = finalEnergies[i];
        s = sigl(ep,e,tev,tolin,lat,iinc,alphas,betas,sab,az,lasym,sigma_b,
                 sigma_b2,teff,nbin,equiprobableBins);
        REQUIRE(ranges::equal(correctMuVecs[i], s, equal));
      }





  } // WHEN



  WHEN( "8 bins are requested" ){ 
      
    nbin = 8;
    std::vector<std::vector<double>> correctMuVecs;
    std::vector<double> finalEnergies;

    e = 1e-2;

    correctMuVecs = {
    {-0.486476,-0.285778,-0.088254, 0.109342,0.307043,0.504771,0.702703,0.901768},
    {-0.144795, 0.165989, 0.400513, 0.598550,0.756077,0.874638,0.954019,0.993320},
    {-0.875034,-0.625260,-0.375352,-0.125469,0.124476,0.374511,0.624475,0.875226},
    {-0.875008,-0.625091,-0.375119,-0.125156,0.124832,0.374848,0.624852,0.875073},
    {-0.875002,-0.625029,-0.375038,-0.125050,0.124946,0.374952,0.624954,0.875023} };
    
    finalEnergies = { 1e-1, 1e-2, 1e-3, 1e-4, 1e-5 };
    for ( size_t i = 0; i < correctMuVecs.size(); ++i ){
      s = sigl(finalEnergies[i],e,tev,tolin,lat,iinc,alphas,betas,sab,az,lasym,
               sigma_b,sigma_b2,teff,nbin,equiprobableBins);
      REQUIRE(ranges::equal(correctMuVecs[i], s, equal));
    }


    e = 1e-1;
    correctMuVecs = {
    { 0.164498, 0.843460, 0.891257, 0.928531,0.957995,0.978053,0.991724,0.998927},
    {-0.486476,-0.285778,-0.088255, 0.109342,0.307043,0.504771,0.702703,0.901768},
    {-0.862936,-0.593610,-0.331200,-0.075487,0.173851,0.417021,0.654333,0.885978},
    {-0.871318,-0.615340,-0.361387,-0.109464,0.140475,0.388454,0.634476,0.878524},
    {-0.873847,-0.621967,-0.370717,-0.120097,0.129896,0.379267,0.628014,0.876137} };
    finalEnergies = { 1e-1,1e-2,1e-3,1e-4,1e-5 };
    for ( size_t i = 0; i < correctMuVecs.size(); ++i ){
      s = sigl(finalEnergies[i],e,tev,tolin,lat,iinc,alphas,betas,sab,az,lasym,
               sigma_b,sigma_b2,teff,nbin,equiprobableBins);
      REQUIRE(ranges::equal(correctMuVecs[i], s, equal));
    }




    } // WHEN






  } // GIVEN

  GIVEN( "legendre expansion is requested (in fortran, nlin (nL) is a positive number)" ){
    WHEN( "8 legendre order requested" ){
      double e = 1e-2, ep = 1e-3, tev = 0.025, tolin = 5e-2, az = 0.99917, 
             sigma_b = 4.0, sigma_b2 = 0.0, teff = 0.12;
      int nL = 9, lat = 1, iinc = 2, lasym = 0, nbin = 8;
      bool equiprobableBins = false;
      std::vector<double> alphas { 1.1, 2.2, 3.3, 4.5, 5.8 },
                           betas { 0.1, 0.2, 1.3, 1.4, 2.5, 2.6, 3.7 };
      std::vector<double> sab(alphas.size()*betas.size());
      for ( size_t i = 0; i < alphas.size(); ++i ){
        for ( size_t j = 0; j < betas.size(); ++j ){
          sab[i*betas.size()+j] = 0.01*((j+1) + 0.1*(i+1));
        } 
      } 
  
      std::vector<double> correct_s, s, finalEnergies;
      std::vector<std::vector<double>> correctMuVecs;

      e = 1e-2;
      correctMuVecs = {
      { 0.208139,  -0.126102,  0.031542, 9.98544E-3,-5.09514E-2, 3.52565E-3,
       -2.33062E-2,-4.73643E-2 },
      { 0.574789,   0.212627, 0.147169, 0.16253311, 0.10777361, 7.40515E-2, 
        8.82367E-2, 7.03711E-2 },
      {-3.03203E-4,-7.81019E-3, 1.80460E-4, -2.53228E-2, 2.75738E-4, -4.86267E-2, 
        1.89076E-5,-6.30910E-2 },
      {-9.61312E-5,-7.80410E-3, 6.30817E-5, -2.54786E-2, 8.61904E-5, -4.87435E-2, 
        5.19336E-6,-6.29823E-2 },
      {-3.04130E-5,-7.80909E-3, 2.05935E-5, -2.55233E-2, 2.74940E-5, -4.87740E-2, 
        1.59765E-6,-6.29532E-2 } };

      finalEnergies = { 1e-1,1e-2,1e-3,1e-4,1e-5 };
      for ( size_t i = 0; i < correctMuVecs.size(); ++i ){
        s = sigl(finalEnergies[i],e,tev,tolin,lat,iinc,alphas,betas,sab,az,lasym,
                 sigma_b,sigma_b2,teff,nbin,equiprobableBins);
        REQUIRE(ranges::equal(correctMuVecs[i], s, equal));
      }



    } // WHEN 

    WHEN( "20 legendre values requested" ){
      double e = 1e-1, ep = 1e-3, tev = 0.025, tolin = 5e-2, az = 0.99917, 
             sigma_b = 4.0, sigma_b2 = 0.0, teff = 0.12;
      int nL = 20, lat = 1, iinc = 2, lasym = 0, nbin = 19;
      bool equiprobableBins = false;
      std::vector<double> alphas { 1.1, 2.2, 3.3, 4.5, 5.8 },
                           betas { 0.1, 0.2, 1.3, 1.4, 2.5, 2.6, 3.7 };
      std::vector<double> sab(alphas.size()*betas.size());
      for ( size_t i = 0; i < alphas.size(); ++i ){
        for ( size_t j = 0; j < betas.size(); ++j ){
          sab[i*betas.size()+j] = 0.01*((j+1) + 0.1*(i+1));
        } 
      } 

      std::cout.precision(15);
      std::vector<double> correct_s, s;

      ep = 1e-2;
      s = sigl(ep,e,tev,tolin,lat,iinc,alphas,betas,sab,az,lasym,sigma_b,sigma_b2,teff,nbin,equiprobableBins);
      correct_s = {0.20813984053303056, -0.12206101754333731, 3.5420004472892841E-2, 1.7480367993305570E-2, -3.6371260166596767E-2, 1.7879061157533079E-2, -2.4264323114832629E-3, -2.0463296371717440E-2, 7.3565359220208737E-3, -1.1596846287971821E-2, -1.6803829040721665E-2, -1.4323126894478649E-3, -1.7576947245690058E-2, -1.5796110004199532E-2, -6.8284006823239249E-3, -1.8536020538439134E-2, -1.1281889099425624E-2, -4.9976566135597825E-3, -1.1095096314347460E-2 };
      REQUIRE(ranges::equal(correct_s, s, equal));

      ep = 1e-3;
      correct_s = {3.3493796107324496E-2, -4.8953137143630671E-4, 2.3096900805850973E-4, -4.6380853622502959E-3, 6.2623868662719379E-4, -9.5663377863375778E-3, 1.1012852030670870E-3, -1.5646270367429341E-2, 1.4123251432655790E-3, -2.1631932255922436E-2, 1.1017003870489143E-3, -2.5372240677067729E-2, -4.2123297915736685E-4, -2.3813679486400695E-2, -3.5273490108372208E-3, -1.4395083152475549E-2, -7.6107554563814789E-3, 3.0046226065482571E-3, -1.0329665857216333E-2 };
      s = sigl(ep,e,tev,tolin,lat,iinc,alphas,betas,sab,az,lasym,sigma_b,sigma_b2,teff,nbin,equiprobableBins);
      REQUIRE(ranges::equal(correct_s, s, equal));


      ep = 1e-4;
      s = sigl(ep,e,tev,tolin,lat,iinc,alphas,betas,sab,az,lasym,sigma_b,sigma_b2,teff,nbin,equiprobableBins);
      correct_s = {1.0552476053064960E-2, -1.2967602067123080E-3, 7.2655858771655812E-5, -4.6044971736832925E-3, 1.9784923415354719E-4, -9.5134807727040738E-3, 3.4780849750069495E-4, -1.5602844354506646E-2, 4.4720947349806764E-4, -2.1670240451877584E-2, 3.5347924496352198E-4, -2.5569650753882744E-2, -1.2139318096246812E-4, -2.4335016577193738E-2, -1.0996235955583397E-3, -1.5193920685933497E-2, -2.4054450617538416E-3, 2.3617840955264813E-3, -3.3148292912540955E-3 };      
      REQUIRE(ranges::equal(correct_s, s, equal));


      ep = 1e-5;
      s = sigl(ep,e,tev,tolin,lat,iinc,alphas,betas,sab,az,lasym,sigma_b,sigma_b2,teff,nbin,equiprobableBins);
      correct_s = {3.3357650219262638E-3, -1.4108638799838619E-3, -1.0344411850295859E-4, -4.8382963236342710E-3, -2.5974813987816522E-4, -9.8781732397069336E-3, -2.5265773306037899E-4, -1.5914769193947603E-2, -1.2274888117982594E-4, -2.1863378409880618E-2, 3.7182953701513782E-5, -2.5512244416752403E-2, 2.3139913076728955E-4, -2.3944331921534098E-2, 1.7567663511084200E-4, -1.4740814554808142E-2, -2.6219538052183952E-4, 2.7226745391690536E-3, -7.0088884419827573E-4 };
      REQUIRE(ranges::equal(correct_s, s, equal));


    } // WHEN 

  } // GIVEN 

}

/*
*/

/*

*/
