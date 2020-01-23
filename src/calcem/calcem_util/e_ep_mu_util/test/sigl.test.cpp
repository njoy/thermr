#define CATCH_CONFIG_MAIN
#include "catch.hpp"
#include "calcem/calcem_util/e_ep_mu_util/sigl.h"
#include <range/v3/all.hpp>
#include "generalTools/testing.h"


auto populate_SAB( int a_size, int b_size ){
  std::vector<double> sab(a_size*b_size);
  for ( int i = 0; i < a_size; ++i ){
    for ( int j = 0; j < b_size; ++j ){
      sab[i*b_size+j] = 0.01*((j+1) + 0.1*(i+1));
    } 
  } 
  return sab;
}







TEST_CASE( "adaptively refine between two points in muVec grid" ){
  GIVEN( "physical constants, S(a,b) grid, and energies (E->E')" ){
    int lasym = 0, lat = 1, iinc = 2; 
    double az = 0.99917,    tevz = 2.53e-2, sigma_b = 163.7279223, sigma_b2 = 0,
         teff = 0.120441926, tol = 2.5e-2,      tev = 2.55e-2;
    std::vector<double> alphas { 1.1, 2.2, 3.3, 4.5, 5.8 },
                         betas { 0.1, 0.2, 1.3, 1.4, 2.5, 2.6, 3.7 };
    auto sab = populate_SAB( alphas.size(), betas.size() );

    double e = 1e-5, ep = 1e-4;

    WHEN( "initial mu vector is populated only with extreme cosines" ){
      std::vector<double> muVec(20,0.0), xsVec(20,0.0), correct_muVec(20), correct_xsVec(20);
      muVec[0] =  1.0; muVec[2] = -1.0;
      xsVec[0] =  2.5; xsVec[2] =  4.0;
      int i = 3;

      THEN( "the cosine vector (and corresponding xs vector) are adaptively refined" ){
        double xsMax = maxOf4Vals( xsVec[0], xsVec[1], xsVec[2], 0.001 );
        addMidpointsRight(i,muVec,xsVec,e,ep,tev,alphas,betas,sab,az,lasym,lat,sigma_b,sigma_b2,teff,iinc,tol,xsMax);
        correct_muVec = { 1, 0, -0.5, -0.75, -0.875, -0.9375, -0.96875, -0.984375, 
        -0.9921875, -0.9960937, -0.9980468, -0.9990234, -0.9995117, -0.9997558, 
        -0.9998779, -0.9999389, -0.9999694, -0.9999847, -0.9999923, -1.0};
        correct_xsVec = { 2.5, 0.0, 143662.3, 136248.8, 132948.0, 131385.1, 130624.0, 
        130248.4, 130061.8, 129968.8, 129922.4, 129899.2, 129887.6, 129881.8, 
        129878.9, 129877.5, 129876.8, 129876.4, 129876.2, 4.0};
        REQUIRE(ranges::equal(correct_muVec, muVec, equal));
        REQUIRE(ranges::equal(correct_xsVec, xsVec, equal));
      } // THEN 
    } // WHEN 
  } // GIVEN
} // TEST CASE



   

TEST_CASE( "Get pdf value" ){
  GIVEN( "material constants and energies (E->E')" ){
    double e, tev = 0.025, tolin = 5e-2*0.5, az = 0.99917, sigma_b = 4.0, sigma_b2 = 0.0, teff = 0.12;
    int lat = 1, iinc = 2, lasym = 0;
    bool equiprobableBins = true;
    std::vector<double> alphas { 1.1, 2.2, 3.3, 4.5, 5.8 },
                         betas { 0.1, 0.2, 1.3, 1.4, 2.5, 2.6, 3.7 };
  
    auto sab = populate_SAB( alphas.size(), betas.size() );
  
    WHEN( "provide various E and E' values" ){
      std::vector<double> eVec {1e-2,1e-1,1e0};
      std::vector<std::vector<double>> 
        epVecs  { { 1e-1, 1e-2, 1e-3, 1e-4, 1e-5 },
                  { 1e-2, 1e-3, 1e-4, 1e-5       },
                  { 1e-2, 0.1, 0.8, 1.0, 2.0     } },
        pdfVecs { { 71.84429, 220.1616, 61.17360, 19.66197, 6.227846 },
                  { 262.937411, 1.02706201, 0.32551817, 0.10296110   },
                  { 0.316342, 0.8001197, 1.0023345, 3.0989839, 0.0   } };
  
      THEN( "integrate the incoherent inelastic scattering xs over all mu" ){
        for ( size_t i = 0; i < eVec.size(); ++i ){
          for ( size_t j = 0; j < epVecs[i].size(); ++j ){
            REQUIRE( getPDF(epVecs[i][j],eVec[i],tev,tolin,lat,iinc,alphas,
                            betas,sab,az,lasym,sigma_b,sigma_b2,teff) == 
                     Approx(pdfVecs[i][j]).epsilon(1e-6) );
          }
        }
      } // THEN
    } // WHEN
  } // GIVEN
} // TEST CASE





TEST_CASE( "sigl" ){

  double e = 1e-2, tev = 0.025, tolin = 5e-2, az = 0.99917, sigma_b = 4.0, 
  sigma_b2 = 0.0, teff = 0.12;
  int lat = 1, iinc = 2, lasym = 0, nbin = 8;

  std::vector<double> alphas { 1.1, 2.2, 3.3, 4.5, 5.8 },
                       betas { 0.1, 0.2, 1.3, 1.4, 2.5, 2.6, 3.7 };

  auto sab = populate_SAB( alphas.size(), betas.size() );



  GIVEN ( "equiprobable angle bins are requested" ){

    bool equiprobableBins = true;
    std::vector<double> correct_s, s;

  WHEN ( "2 bins requested" ){
      nbin = 2;

      std::vector<double> epVec;
      std::vector<std::vector<double>> correctMu;

      e = 1e-2;
      epVec = { 1e0, 1e-1, 1e-2, 1e-3, 1e-4, 1e-5 };
      correctMu = { { 0.0000000, 0.0000000 }, {-0.1877993, 0.6040790 },
                    { 0.2550645, 0.8945141 }, {-0.5002791, 0.4996727 },
                    {-0.5000937, 0.4999014 }, {-0.5000301, 0.4999693 } };
      for ( size_t i = 0; i < correctMu.size(); ++i ){
        s = sigl(epVec[i],e,tev,tolin,lat,iinc,alphas,betas,sab,az,lasym,sigma_b,
                 sigma_b2,teff,nbin,equiprobableBins);
        REQUIRE(ranges::equal(correctMu[i], s, equal));
      }


      e = 1.0;
      epVec = { 1e-2, 1e-1, 5e-1, 1.0, 2.0 };
      correctMu = { {-0.458477, 0.527666}, {-0.289542, 0.581023}, 
                    { 0.300582, 0.783071}, { 0.914921, 0.999029}, { 0, 0 } };
      for ( size_t i = 0; i < correctMu.size(); ++i ){
        s = sigl(epVec[i],e,tev,tolin,lat,iinc,alphas,betas,sab,az,lasym,sigma_b,
                 sigma_b2,teff,nbin,equiprobableBins);
        REQUIRE(ranges::equal(correctMu[i], s, equal));
      }

  } // WHEN



  WHEN( "8 bins are requested" ){ 
      
    nbin = 8;
    std::vector<std::vector<double>> correctMu;
    std::vector<double> epVec;

    e = 1e-2;

    correctMu= {
    {-0.486476,-0.285778,-0.088254, 0.109342,0.307043,0.504771,0.702703,0.901768},
    {-0.144795, 0.165989, 0.400513, 0.598550,0.756077,0.874638,0.954019,0.993320},
    {-0.875034,-0.625260,-0.375352,-0.125469,0.124476,0.374511,0.624475,0.875226},
    {-0.875008,-0.625091,-0.375119,-0.125156,0.124832,0.374848,0.624852,0.875073},
    {-0.875002,-0.625029,-0.375038,-0.125050,0.124946,0.374952,0.624954,0.875023} };
    
    epVec = { 1e-1, 1e-2, 1e-3, 1e-4, 1e-5 };
    for ( size_t i = 0; i < correctMu.size(); ++i ){
      s = sigl(epVec[i],e,tev,tolin,lat,iinc,alphas,betas,sab,az,lasym,
               sigma_b,sigma_b2,teff,nbin,equiprobableBins);
      REQUIRE(ranges::equal(correctMu[i], s, equal));
    }


    e = 1e-1;
    correctMu = {
    { 0.164498, 0.843460, 0.891257, 0.928531,0.957995,0.978053,0.991724,0.998927},
    {-0.486476,-0.285778,-0.088255, 0.109342,0.307043,0.504771,0.702703,0.901768},
    {-0.862936,-0.593610,-0.331200,-0.075487,0.173851,0.417021,0.654333,0.885978},
    {-0.871318,-0.615340,-0.361387,-0.109464,0.140475,0.388454,0.634476,0.878524},
    {-0.873847,-0.621967,-0.370717,-0.120097,0.129896,0.379267,0.628014,0.876137} };
    epVec = { 1e-1,1e-2,1e-3,1e-4,1e-5 };
    for ( size_t i = 0; i < correctMu.size(); ++i ){
      s = sigl(epVec[i],e,tev,tolin,lat,iinc,alphas,betas,sab,az,lasym,
               sigma_b,sigma_b2,teff,nbin,equiprobableBins);
      REQUIRE(ranges::equal(correctMu[i], s, equal));
    }




    } // WHEN
  } // GIVEN

  GIVEN( "legendre expansion is requested (in fortran, nlin (nL) is a positive number)" ){
    bool equiprobableBins = false;
    WHEN( "8 legendre order requested" ){

      nbin = 8;
      std::vector<double> correct_s, s, epVec;
      std::vector<std::vector<double>> correctMu;

      e = 1e-2;
      correctMu = {
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

      epVec = { 1e-1,1e-2,1e-3,1e-4,1e-5 };
      for ( size_t i = 0; i < correctMu.size(); ++i ){
        s = sigl(epVec[i],e,tev,tolin,lat,iinc,alphas,betas,sab,az,lasym,
                 sigma_b,sigma_b2,teff,nbin,equiprobableBins);
        REQUIRE(ranges::equal(correctMu[i], s, equal));
      }
    } // WHEN 

    WHEN( "19 legendre values requested" ){
      nbin = 19;

      std::vector<double> correct_s, s;

      e = 1e-1;

      std::vector<double> epVec { 1e-5, 1e-4, 1e-3, 1e-2 };

      double ep = 1e-2;
      s = sigl(ep,e,tev,tolin,lat,iinc,alphas,betas,sab,az,lasym,sigma_b,sigma_b2,teff,nbin,equiprobableBins);
      correct_s = { 0.2081398, -0.1220610,  0.035420004, 0.01748036,-0.03637126, 
      0.01787906, -2.426432E-3,-0.02046329, 7.356535E-3,-0.01159684,-0.01680382, 
     -1.432312E-3,-0.01757694, -0.01579611,-6.828400E-3,-0.01853602,-0.01128188, 
     -4.997656E-3,-.01109509 };
      REQUIRE(ranges::equal(correct_s, s, equal));

      ep = 1e-3;
      correct_s = {  3.349379E-2, -4.895313E-4,  2.309690E-4, -4.638085E-3, 
       6.262386E-4, -9.566337E-3,  1.101285E-3, -1.564627E-2,  1.412325E-3,
      -2.163193E-2,  1.101700E-3, -2.537224E-2, -4.212329E-4, -2.381367E-2,
      -3.527349E-3, -1.439508E-2, -7.610755E-3,  3.004622E-3, -1.032966E-2 };
      s = sigl(ep,e,tev,tolin,lat,iinc,alphas,betas,sab,az,lasym,sigma_b,sigma_b2,teff,nbin,equiprobableBins);
      REQUIRE(ranges::equal(correct_s, s, equal));

      ep = 1e-4;
      s = sigl(ep,e,tev,tolin,lat,iinc,alphas,betas,sab,az,lasym,sigma_b,sigma_b2,teff,nbin,equiprobableBins);
      correct_s = {  1.055247E-2, -1.296760E-3,  7.265585E-5, -4.604497E-3, 
       1.978492E-4, -9.513480E-3,  3.478084E-4, -1.560284E-2,  4.472094E-4,
      -2.167024E-2,  3.534792E-4, -2.556965E-2, -1.213931E-4, -2.433501E-2,
      -1.099623E-3, -1.519392E-2, -2.405445E-3,  2.361784E-3, -3.314829E-3 };      
      REQUIRE(ranges::equal(correct_s, s, equal));


      ep = 1e-5;
      s = sigl(ep,e,tev,tolin,lat,iinc,alphas,betas,sab,az,lasym,sigma_b,sigma_b2,teff,nbin,equiprobableBins);
      correct_s = {  3.335765E-3, -1.410863E-3, -1.034441E-4, -4.838296E-3, 
      -2.597481E-4, -9.878173E-3, -2.526577E-4, -1.591476E-2, -1.227488E-4, 
      -2.186337E-2,  3.718295E-5, -2.551224E-2,  2.313991E-4, -2.394433E-2, 
      1.7567663E-4, -1.474081E-2, -2.621953E-4,  2.722674E-3, -7.008888E-4 };
      REQUIRE(ranges::equal(correct_s, s, equal));


    } // WHEN 

  } // GIVEN 

}

