#define CATCH_CONFIG_MAIN
#include "catch.hpp"
#include "inelastic/e_ep_mu_util/sigl.h"
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

    std::vector<double> boundXsVec {sigma_b, sigma_b2};
    double e = 1e-5, ep = 1e-4;

    WHEN( "initial mu vector is populated only with extreme cosines" ){
      std::vector<double> muVec(20,0.0), xsVec(20,0.0), correct_muVec(20), 
                                                        correct_xsVec(20);
      muVec[0] =  1.0; muVec[2] = -1.0;
      xsVec[0] =  2.5; xsVec[2] =  4.0;
      int i = 3;

      THEN( "the cosine vector (and corresponding xs vector) are adaptively refined" ){
        double xsMax = maxOf4Vals( xsVec[0], xsVec[1], xsVec[2], 0.001 );
        addMidpointsRight(i,muVec,xsVec,e,ep,tev,alphas,betas,sab,az,lasym,lat,
                          boundXsVec,teff,iinc,tol,xsMax);
        correct_muVec = { 1.0, 0.0, -0.5, -0.75, -0.875, -0.93750, -0.96875, 
          -0.98437500, -0.99218750, -0.99609375, -0.99804687, -0.99902343, 
          -0.99951171, -0.99975585, -0.99987792, -0.99993896, -0.99996948, 
          -0.99998474, -0.99999237,  -1.0 };
        correct_xsVec = { 2.5, 0.0, 9690.6394250, 9692.1192405, 9692.8043549, 
           9693.1347393, 9693.2970539, 9693.3775111, 9693.4175671, 9693.4375521, 
           9693.4475340, 9693.4525222, 9693.4550157, 9693.4562623, 9693.4568855, 
           9693.4571971, 9693.4573529, 9693.4574309, 9693.4574698, 4.0 };
        REQUIRE(ranges::equal(correct_muVec, muVec, equal));
        REQUIRE(ranges::equal(correct_xsVec, xsVec, equal));
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

  std::vector<double> boundXsVec {sigma_b, sigma_b2};
  std::vector<double> epVec, s, correctPDFs;
  std::vector<std::vector<double>> correctMu;
  double pdfVal;

  GIVEN ( "equiprobable angle bins are requested" ){

    bool equiprobableBins = true;

    WHEN ( "2 angular bins requested" ){
      nbin = 2;
      std::vector<double> s(nbin,0.0);

      THEN( "mu bins are correctly generated for each E, E' combination" ){
        e = 1e-2;
        epVec = { 1e0, 1e-1, 1e-2, 1e-3, 1e-4, 1e-5 };
        correctMu = { { 0.0000000, 0.0000000 }, {-0.1877993, 0.6040790 },
                      { 0.0000000, 0.0000000 }, {-0.5002791, 0.4996727 },
                      {-0.5000937, 0.4999014 }, {-0.5000301, 0.4999693 } };
        correctPDFs = { 0.0, 71.84429, 0.0, 61.173606,  19.661973, 6.2278467 };

        for ( size_t i = 0; i < correctPDFs.size(); ++i ){
          pdfVal = sigl(epVec[i],e,tev,tolin,lat,iinc,alphas,betas,sab,az,lasym,
                   boundXsVec,teff,s,equiprobableBins);
          REQUIRE(ranges::equal(correctMu[i], s, equal));
          REQUIRE( pdfVal == Approx(correctPDFs[i]) );
        }

        e = 1.0;
        epVec = { 1e-2, 1e-1, 5e-1, 1.0, 2.0 };
        correctMu = { {-0.458477, 0.527666}, {-0.289542, 0.581023}, 
                      { 0.300582, 0.783071}, { 0.626336, 0.879973}, { 0, 0 } };
        correctPDFs = { 0.3163424, 0.8001197, 0.9932267, 0.437514, 0.0 };
        for ( size_t i = 0; i < correctMu.size(); ++i ){
          pdfVal =  sigl(epVec[i],e,tev,tolin,lat,iinc,alphas,betas,sab,az,lasym,
                         boundXsVec,teff,s,equiprobableBins);
          REQUIRE(ranges::equal(correctMu[i], s, equal));
          REQUIRE( pdfVal == Approx(correctPDFs[i]) );
        }
      } // THEN
    } // WHEN




    WHEN ( "4 angular bins requested" ){
      nbin = 4;
      std::vector<double> s(nbin,0.0);

      AND_WHEN( "we have room temperature water inputs" ){
        sigma_b  = 163.727922;  sigma_b2 = 0.0;      tev   = 2.5507297688e-2;
        teff     = 0.1204419;   az       = 0.99917;  
        lasym    = 0;           lat      = 0;        iinc  = 2;
        boundXsVec = {sigma_b, sigma_b2};

        sab = { -0.1825961, -0.3020134, -3.9365477, -3.9880917, -4.3354560, 
        -4.3951540, -5.8893492, -0.7622529, -0.8165834, -3.1416145, -3.3056618, 
        -3.9055465, -3.9623336, -5.2369666, -1.1918288, -1.2315547, -2.7961056, 
        -2.9563309, -3.7498922, -3.8083758, -4.9337391, -1.5834286, -1.6131071, 
        -2.7123394, -2.8429160, -3.6969959, -3.7519934, -4.7743385, -1.9612120, 
        -1.9872066, -2.7845460, -2.8853146, -3.7128812, -3.7714214, -4.7115839 };

        
        THEN( "mu bins are correctly generated for each E, E' combination" ){
          std::cout.precision(15);
          e = 1e-2;
          epVec = { 1e-5, 1e-4, 1e-3, 1e-2 };
          correctMu = { 
                        { -0.7525806, -0.25574937, 0.24416911, 0.74731460 },
                        { -0.7576714, -0.26752370, 0.23168669, 0.74130696 },
                        { -0.7673533, -0.29204855, 0.20221922, 0.72482378 },
                        { -0.5204081,  0.22722905, 0.71478052, 0.95924769 }
          };

          for ( size_t i = 0; i < correctMu.size(); ++i ){
            pdfVal = sigl(epVec[i],e,tev,tolin,lat,iinc,alphas,betas,sab,az,lasym,
                     boundXsVec,teff,s,equiprobableBins);
            REQUIRE(ranges::equal(correctMu[i], s, equal));
          }



          e = 1e-2;
          tolin = 1e-2;
          epVec = { 1e-5, 1e-4, 1e-3, 1e-2 };
          correctMu = {
              { -0.7525806, -0.25574937, 0.24416911, 0.74731460 },
              { -0.7576714, -0.26752370, 0.23168669, 0.74130696 },
              { -0.7673533, -0.29204855, 0.20221922, 0.72482378 },
              { -0.5226632,  0.22242174, 0.70996733, 0.95789197 } 
          };

          for ( size_t i = 0; i < correctMu.size(); ++i ){
            pdfVal = sigl(epVec[i],e,tev,tolin,lat,iinc,alphas,betas,sab,az,lasym,
                          boundXsVec,teff,s,equiprobableBins);
            REQUIRE(ranges::equal(correctMu[i], s, equal));
          }


        } // THEN

      } // AND WHEN

    } // WHEN

    WHEN( "8 bins are requested" ){ 
      
      nbin = 8;
      std::vector<double> s(nbin,0.0);

      THEN( "mu bins are correctly generated for each E, E' combination" ){
        epVec = { 1e-1, 1e-2, 1e-3, 1e-4, 1e-5 };

        e = 1e-2;
        correctMu= {
        { -0.48647654, -0.28577865, -8.82549593E-2, 0.10934278, 0.30704325, 
           0.50477104,  0.70270321,  0.90176859},
        {  0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0},
        { -0.875034246, -0.62526004, -0.37535265, -0.125469544, 0.124476845, 
           0.374511893,  0.62447532,  0.87522679},
        { -0.875008004, -0.62509139, -0.37511937, -0.125156114, 0.124832101, 
           0.374848178,  0.62485204,  0.87507351 },
        { -0.875002118, -0.62502992, -0.37503857, -0.125050067, 0.124946798, 
           0.374952310,  0.62495468,  0.87502358 } };
        correctPDFs = { 71.844290, 0.0, 61.1736062, 19.6619735, 6.22784676 };
    
        for ( size_t i = 0; i < correctMu.size(); ++i ){
          pdfVal = sigl(epVec[i],e,tev,tolin,lat,iinc,alphas,betas,sab,az,lasym,
                        boundXsVec,teff,s,equiprobableBins);
          REQUIRE(ranges::equal(correctMu[i], s, equal));
          REQUIRE( pdfVal == Approx(correctPDFs[i]) );
        }

        e = 1e-1;
        correctMu = {
          { -0.874774735, -0.651125119, -0.456915386, -0.287144697, 
            -0.137871601, -5.6328329E-3, 0.112095602,  0.217506519 },
          { -0.486476543, -0.285778657, -8.825495934E-2, 0.109342786, 
             0.307043251,  0.504771042,  0.702703213,  0.901768590 } ,
          { -0.862935950, -0.593610548, -0.331200729, -7.54874931E-2, 
             0.173851635,  0.417021713,  0.654333678,  0.885978062 },
          { -0.87131822,  -0.61534055,  -0.36138752,  -0.10946478, 
             0.14047594,   0.38845396,   0.63447693,   0.878524053 },
          { -0.87384791, -0.62196710, -0.37071696, -0.12009711, 0.12989681, 
             0.37926704, 0.62801413, 0.876137228 } };
        correctPDFs = { 2.39519541, 262.9374112,  1.0270620, 0.32551817, 0.1029611 };

        for ( size_t i = 0; i < correctMu.size(); ++i ){
          pdfVal = sigl(epVec[i],e,tev,tolin,lat,iinc,alphas,betas,sab,az,lasym,
                        boundXsVec,teff,s,equiprobableBins);
          REQUIRE(ranges::equal(correctMu[i], s, equal));
          REQUIRE( pdfVal == Approx(correctPDFs[i]) );
        }

      } // THEN
    } // WHEN
  } // GIVEN

  /*
  GIVEN( "legendre expansion is requested" ){
    bool equiprobableBins = false;

    WHEN( "8 legendre order requested" ){
      nbin = 8;
      std::vector<double> s(nbin,0.0);

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
        pdfVal = sigl(epVec[i],e,tev,tolin,lat,iinc,alphas,betas,sab,az,lasym,
                      boundXsVec,teff,s,equiprobableBins);
        //REQUIRE(ranges::equal(correctMu[i], s, equal));
      }
    } // WHEN 

    WHEN( "19 legendre values requested" ){
      nbin = 19;
      std::vector<double> s(nbin,0.0);

      THEN( "mu bins are correctly generated for each E, E' combination" ){
        e = 1e-1;

        epVec = { 1e-5, 1e-4, 1e-3, 1e-2 };
        correctMu = { 
          { 3.335765E-3, -1.410863E-3, -1.034441E-4, -4.838296E-3, -2.597481E-4, 
           -9.878173E-3, -2.526577E-4, -1.591476E-2, -1.227488E-4, -2.186337E-2,  
            3.718295E-5, -2.551224E-2,  2.313991E-4, -2.394433E-2,  1.756766E-4, 
           -1.474081E-2, -2.621953E-4,  2.722674E-3, -7.008888E-4 },
          { 1.055247E-2, -1.296760E-3,  7.265585E-5, -4.604497E-3,  1.978492E-4, 
           -9.513480E-3,  3.478084E-4, -1.560284E-2,  4.472094E-4, -2.167024E-2,  
            3.534792E-4, -2.556965E-2, -1.213931E-4, -2.433501E-2, -1.099623E-3, 
           -1.519392E-2, -2.405445E-3,  2.361784E-3, -3.314829E-3 },      
          { 3.349379E-2, -4.895313E-4,  2.309690E-4, -4.638085E-3,  6.262386E-4, 
           -9.566337E-3,  1.101285E-3, -1.564627E-2,  1.412325E-3, -2.163193E-2,  
            1.101700E-3, -2.537224E-2, -4.212329E-4, -2.381367E-2, -3.527349E-3, 
           -1.439508E-2, -7.610755E-3,  3.004622E-3, -1.032966E-2 },
          { 0.2081398,  -0.1220610,    0.035420004,  0.01748036,  -0.03637126, 
            0.01787906, -2.426432E-3, -0.02046329,   7.356535E-3, -0.01159684, 
           -0.01680382, -1.432312E-3, -0.01757694,  -0.01579611,  -6.828400E-3, 
           -0.01853602, -0.01128188,  -4.997656E-3, -0.01109509 } };

        for ( size_t i = 0; i < correctMu.size(); ++i ){
          pdfVal = sigl(epVec[i],e,tev,tolin,lat,iinc,alphas,betas,sab,az,lasym,
                        boundXsVec,teff,s,equiprobableBins);
          REQUIRE(ranges::equal(correctMu[i], s, equal));
        }
      } // THEN 
    } // WHEN 
  } // GIVEN 
*/
} // TEST CASE

