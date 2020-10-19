#define CATCH_CONFIG_MAIN
#include "catch.hpp"
#include "calcem/calcem_util/e_mu_ep_util/sigu.h"
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
    std::vector<double> boundXsVec {sb,sb2};

    std::vector<double> alpha { 1.1, 2.2, 3.3, 4.5, 5.8 },
      beta { 0.1, 0.2, 1.3, 1.4, 2.5, 2.6, 3.7 }, x(20,0.0), y(20,0.0);
    std::vector<double> sab(alpha.size()*beta.size());
    for ( size_t i = 0; i < alpha.size(); ++i ){
      for ( size_t j = 0; j < beta.size(); ++j ){
        sab[i*beta.size()+j] = 0.01*((j+1) + 0.1*(i+1));
      } 
    } 

    std::vector<double> cosines {-1.0,-0.8,-0.5,0.1,0.9}, correct_xVals(5), 
                        correct_yVals(5), correct_x(20,0.0), correct_y(20,0.0);

    WHEN( "lat = 1" ){
      correct_xVals = { 7.1395952E-4, 7.3838244E-4, 7.7664539E-4, 
                        2.6299999E-2, 2.6299999E-2 },
      correct_yVals = { 67960510168.6, 32016311225.0, 9925375996.58, 0.0, 0.0 };
      THEN("x, y, and jbeta are correctly returned for each scattering cosine"){
        for (size_t i = 0; i < cosines.size(); ++i){
          correct_x[0] = correct_xVals[i];
          correct_y[0] = correct_yVals[i];
          initializeEpXS( jbeta, lat, x, y, e, tev, tevz, cosines[i], alpha, beta, sab, 
                      az, lasym, teff, boundXsVec, iinc);
          break;
          //checkVec(x,correct_x);
          //checkVec(y,correct_y);
          //REQUIRE( ranges::equal(x, correct_x, equal) );
          //REQUIRE( ranges::equal(y, correct_y, equal) );
          REQUIRE( jbeta == 1 );
        }
      } // THEN
    } // WHEN

    /*
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
          initializeEpXS( jbeta, lat, x, y, e, tev, tevz, cosines[i], alpha, beta, sab, 
                      az, lasym, teff, boundXsVec, iinc);
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
    std::vector<double> boundXsVec {sb,sb2};

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
          initializeEpXS( jbeta, lat, x, y, e, tev, tevz, cosines[i], alpha, beta, sab, 
                      az, lasym, teff, boundXsVec, iinc);
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
          initializeEpXS( jbeta, lat, x, y, e, tev, tevz, cosines[i], alpha, beta, sab, 
                      az, lasym, teff, boundXsVec, iinc);
          REQUIRE( ranges::equal(x, correct_x, equal) );
          REQUIRE( ranges::equal(y, correct_y, equal) );
          REQUIRE( jbeta == 1 );
        }
      } // THEN
    } // WHEN
  } // GIVEN
} // TEST CASE




TEST_CASE( "sigu" ){

  int lasym = 0, lat = 1, iinc = 2, nemax = 60;
  double e = 1e-5, tev = 2.55e-2, az = 0.99917, sb = 4.0, sb2 = 0.0, teff = 0.12, tolin = 5e-2, u;
  std::vector<double> boundXsVec {sb,sb2};

  std::vector<double> alphas { 1.1, 2.2, 3.3, 4.5, 5.8 },
                       betas { 0.1, 0.2, 1.3, 1.4, 2.5, 2.6, 3.7 };
  std::vector<double> sab(alphas.size()*betas.size());
  for ( size_t i = 0; i < alphas.size(); ++i ){
    for ( size_t j = 0; j < betas.size(); ++j ){
      sab[i*betas.size()+j] = 0.01*((j+1) + 0.1*(i+1));
    } 
  } 
  std::cout.precision(15);

  GIVEN( "various scattering cosines" ){
    WHEN( "u = -1.0" ){
      u = -1.0;
      std::vector<double> s1(nemax,0.0),s2(nemax,0.0),s3(2*nemax,0.0),
      correct_s1 { 393.910081, 0, 4.309201E-13, 8.618402E-13, 1.7236804E-12, 
      1.9342239E-8, 3.8682755E-8, 7.7363787E-8, 1.5472585E-7, 3.0944997E-7, 
      6.1889821E-7, 1.2377947E-6, 2.4755877E-6, 4.9511736E-6, 9.9023454E-6, 
      1.9804689E-5, 3.9609376E-5, 7.9218751E-5, 1.5843750E-4, 3.1687500E-4, 
      6.3375000E-4, 1.2675000E-3, 2.5350000E-3, 3.8025000E-3, 4.4362500E-3, 
      4.7531250E-3, 4.9115625E-3, 4.9907813E-3, 5.0303907E-3, 5.0501954E-3, 
      5.0600977E-3, 5.0650489E-3, 5.0675245E-3, 5.0687623E-3, 5.0693812E-3, 
      5.0696906E-3, 5.0698453E-3, 5.0699227E-3, 5.0699614E-3, 5.0699807E-3, 
      5.0699904E-3, 5.0699952E-3, 5.0699976E-3, 5.0699988E-3, 5.0699994E-3, 
      5.0699997E-3, 5.0699999E-3, 5.0700000E-3, 1.2027500E-2, 1.8985000E-2, 
      3.2900000E-2, 3.5430000E-2, 6.3260000E-2, 6.5790000E-2, 7.9705000E-2, 
      8.6662500E-2, 9.0141250E-2, 9.1880625E-2, 9.2750313E-2, 9.3185157E-2, 
      9.3402579E-2, 9.3511290E-2, 9.3565645E-2, 9.3592823E-2, 9.3606412E-2, 
      9.3613206E-2, 9.3616603E-2, 9.3618302E-2, 9.3619151E-2, 9.3619576E-2, 
      9.3619788E-2, 9.3619894E-2, 9.3619947E-2, 9.3619974E-2, 9.3619987E-2, 
      9.3619994E-2, 9.3619997E-2, 9.3619999E-2, 9.362E-2, 0, 0, 0, 0, 0, 0, 0, 0, 
      0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
      correct_s2 { 78.0, 0, 0.86812825908, 1.2276132005, 1.7358961666, 
      176.21088440, 244.91981576, 338.16164285, 462.73065256, 625.71756699, 
      833.24881808, 1088.5346227, 1389.5569827, 1727.2900621, 2085.6787026, 
      2444.1229634, 2781.8564282, 3082.3091937, 3335.5023797, 3537.7410097, 
      3689.1586656, 3790.0173521, 3836.0912359, 3832.4067140, 3823.7633604, 
      3818.4357145, 3815.5700435, 3814.0916618, 3813.3416363, 3812.9639802, 
      3812.7745004, 3812.6795972, 3812.6321052, 3812.6083491, 3812.5964686, 
      3812.5905286, 3812.5875585, 3812.5860724, 3812.5853294, 3812.5849588, 
      3812.5847726, 3812.5846804, 3812.5846343, 3812.5846113, 3812.5845998, 
      3812.5845940, 3812.5845902, 1629.3457424, 2160.1978761, 2361.7853083, 
      2435.4047613, 2429.3914507, 1901.8474047, 1864.3607137, 1607.0867876, 
      1457.6959400, 1380.6181996, 1341.7985483, 1322.3552138, 1312.6295853, 
      1307.7663214, 1305.3346438, 1304.1188132, 1303.5108869, 1303.2069240, 
      1303.0549538, 1302.9789687, 1302.9409649, 1302.9219743, 1302.9124677, 
      1302.9077257, 1302.9053546, 1302.9041691, 1302.9035652, 1302.9032744, 
      1302.9031178, 1302.9030507, 1302.9030060, 12.967700530, 0, 0, 0, 0, 0, 0, 
      0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
  
      sigu(e,u,tev,alphas,betas,sab,tolin,az,iinc,lat,lasym,boundXsVec,teff,s3);
    
      for ( size_t i = 0; i < correct_s1.size(); ++i ){
        REQUIRE( correct_s1[i] == Approx(s3[2*i]).epsilon(1e-6) );
        REQUIRE( correct_s2[i] == Approx(s3[2*i+1]).epsilon(1e-6) );

      }

    } // WHEN 
  
    WHEN( "u = -0.2" ){
      u = -0.2;
      std::vector<double> s1(nemax,0.0),s2(nemax,0.0),s3(2*nemax,0.0),
      correct_s1 { 396.2647985, 0, 1.0993286E-11, 2.1986572E-11, 4.3973143E-11, 
      1.9384488E-8, 3.8725003E-8, 7.7406032E-8, 1.5476809E-7, 3.0949221E-7, 
      6.1894044E-7, 1.2378369E-6, 2.4756298E-6, 4.9512157E-6, 9.9023875E-6, 
      1.9804731E-5, 3.9609418E-5, 7.9218792E-5, 1.5843754E-4, 3.1687503E-4, 
      6.3375002E-4, 1.2675000E-3, 2.5350000E-3, 3.8025000E-3, 4.4362500E-3, 
      4.7531250E-3, 4.9115625E-3, 4.9907813E-3, 5.0303907E-3, 5.0501954E-3, 
      5.0600977E-3, 5.0650489E-3, 5.0675245E-3, 5.0687623E-3, 5.0693812E-3, 
      5.0696906E-3, 5.0698453E-3, 5.0699227E-3, 5.0699614E-3, 5.0699807E-3, 
      5.0699904E-3, 5.0699952E-3, 5.0699976E-3, 5.0699988E-3, 5.0699994E-3, 
      5.0699997E-3, 5.0699999E-3, 5.0700000E-3, 1.2027500E-2, 1.8985000E-2, 
      3.2900000E-2, 3.5430000E-2, 6.3260000E-2, 6.5790000E-2, 7.9705000E-2, 
      8.6662500E-2, 9.0141250E-2, 9.1880625E-2, 9.2750313E-2, 9.3185157E-2, 
      9.3402579E-2, 9.3511290E-2, 9.3565645E-2, 9.3592823E-2, 9.3606412E-2, 
      9.3613206E-2, 9.3616603E-2, 9.3618302E-2, 9.3619151E-2, 9.3619576E-2, 
      9.3619788E-2, 9.3619894E-2, 9.3619947E-2, 9.3619974E-2, 9.3619987E-2, 
      9.3619994E-2, 9.3619997E-2, 9.3619999E-2, 9.3620000E-2, 0, 0, 0, 0, 0, 0, 
      0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
      correct_s2 { 78.0, 0, 4.38478329, 6.20047815, 8.76771402, 182.388489, 
      256.632392, 360.355873, 504.181921, 701.187934, 965.586062, 1308.62820, 
      1730.20135, 2207.78372, 2692.59709, 3125.70252, 3467.07443, 3709.87555, 
      3870.25560, 3969.98766, 4025.91244, 4046.18048, 4028.53443, 3995.00495, 
      3976.33544, 3966.74328, 3961.89795, 3959.46446, 3958.24517, 3957.63491, 
      3957.32963, 3957.17695, 3957.10060, 3957.06242, 3957.04333, 3957.03379, 
      3957.02902, 3957.02663, 3957.02544, 3957.02484, 3957.02454, 3957.02439, 
      3957.02432, 3957.02428, 3957.02426, 3957.02425, 3957.02425, 1629.18740, 
      2160.05899, 2361.66359, 2435.32438, 2429.30824, 1901.75904, 1864.27301, 
      1607.00689, 1457.62248, 1380.54783, 1341.72979, 1322.28727, 1312.56206, 
      1307.69900, 1305.26743, 1304.05165, 1303.44375, 1303.13980, 1302.98784, 
      1302.91185, 1302.87385, 1302.85486, 1302.84536, 1302.84061, 1302.83824, 
      1302.83706, 1302.83645, 1302.83616, 1302.83601, 1302.83594, 1302.83590, 
      13.0751407, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
  
      sigu(e,u,tev,alphas,betas,sab,tolin,az,iinc,lat,lasym,boundXsVec,teff,s3);
  
      for ( size_t i = 0; i < correct_s1.size(); ++i ){
        REQUIRE( correct_s1[i] == Approx(s3[2*i]).epsilon(1e-6) );
        REQUIRE( correct_s2[i] == Approx(s3[2*i+1]).epsilon(1e-6) );

      }
    } // WHEN

    WHEN( "u = 0.0" ){
      u = 0.0;
      std::vector<double> s1(nemax,0.0),s2(nemax,0.0),s3(2*nemax,0.0),
      correct_s1 { 396.95156, 0, 9.6893315E-9, 1.9378663E-8, 3.8757325E-8, 
      7.7514650E-8, 1.5502930E-7, 3.1005860E-7, 6.2011720E-7, 1.2402344E-6, 
      2.4804688E-6, 4.9609375E-6, 9.9218750E-6, 1.9843750E-5, 3.9687500E-5, 
      7.9375000E-5, 1.5875000E-4, 3.1750000E-4, 6.3500000E-4, 1.2700000E-3, 
      2.5400000E-3, 3.8050000E-3, 4.4375000E-3, 4.7537500E-3, 4.9118750E-3, 
      4.9909375E-3, 5.0304688E-3, 5.0502344E-3, 5.0601172E-3, 5.0650586E-3, 
      5.0675293E-3, 5.0687647E-3, 5.0693824E-3, 5.0696912E-3, 5.0698456E-3, 
      5.0699228E-3, 5.0699614E-3, 5.0699807E-3, 5.0699904E-3, 5.0699952E-3, 
      5.0699976E-3, 5.0699988E-3, 5.0699994E-3, 5.0699997E-3, 5.0699999E-3, 
      5.0700000E-3, 1.2027500E-2, 1.8985000E-2, 3.2900000E-2, 3.5430000E-2, 
      6.3260000E-2, 6.5790000E-2, 7.9705000E-2, 8.6662500E-2, 9.0141250E-2, 
      9.1880625E-2, 9.2750313E-2, 9.3185157E-2, 9.3402579E-2, 9.3511290E-2, 
      9.3565645E-2, 9.3592823E-2, 9.3606412E-2, 9.3613206E-2, 9.3616603E-2, 
      9.3618302E-2, 9.3619151E-2, 9.3619576E-2, 9.3619788E-2, 9.3619894E-2, 
      9.3619947E-2, 9.3619974E-2, 9.3619987E-2, 9.3619994E-2, 9.3619997E-2, 
      9.3619999E-2, 9.3620000E-2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
      0, 0, 0, 0, 0, 0, 0, 0 },
      correct_s2 { 76.0, 0, 130.140379, 183.957183, 259.903258, 366.850496, 
      516.819051, 725.370941, 1010.72857, 1389.35860, 1864.57645, 2408.22326, 
      2951.04475, 3409.25922, 3735.74749, 3937.68007, 4049.76753, 4105.37501, 
      4125.57645, 4118.70564, 4081.17580, 4038.89659, 4017.36231, 4006.56371, 
      4001.16087, 3998.45898, 3997.10796, 3996.43244, 3996.09468, 3995.92580, 
      3995.84136, 3995.79913, 3995.77802, 3995.76747, 3995.76219, 3995.75955, 
      3995.75823, 3995.75758, 3995.75724, 3995.75708, 3995.75700, 3995.75696, 
      3995.75694, 3995.75693, 3995.75692, 1629.14609, 2160.02328, 2361.63247, 
      2435.30429, 2429.28744, 1901.73691, 1864.25104, 1606.98687, 1457.60408, 
      1380.53020, 1341.71256, 1322.27025, 1312.54514, 1307.68214, 1305.25059, 
      1304.03482, 1303.42693, 1303.12298, 1302.97102, 1302.89504, 1302.85704, 
      1302.83805, 1302.82854, 1302.82380, 1302.82143, 1302.82024, 1302.81964, 
      1302.81935, 1302.81919, 1302.81912, 1302.81908, 13.1022021, 0, 0, 0, 
      0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
  
      sigu(e,u,tev,alphas,betas,sab,tolin,az,iinc,lat,lasym,boundXsVec,teff,s3);
  
      for ( size_t i = 0; i < correct_s1.size(); ++i ){
        REQUIRE( correct_s1[i] == Approx(s3[2*i]).epsilon(1e-6) );
        REQUIRE( correct_s2[i] == Approx(s3[2*i+1]).epsilon(1e-6) );
      }
    } // WHEN

    WHEN( "u = 0.5" ){
      u = 0.5;
      std::vector<double> s1(nemax,0.0),s2(nemax,0.0),s3(2*nemax,0.0),
      correct_s1 { 398.9380160, 0, 9.6893315E-9, 1.9378663E-8, 3.8757325E-8, 
      7.7514650E-8, 1.5502930E-7, 3.1005860E-7, 6.2011720E-7, 1.2402344E-6, 
      2.4804688E-6, 4.9609375E-6, 9.9218750E-6, 1.9843750E-5, 3.9687500E-5, 
      7.9375000E-5, 1.5875000E-4, 3.1750000E-4, 6.3500000E-4, 1.2700000E-3, 
      2.5400000E-3, 3.8050000E-3, 4.4375000E-3, 4.7537500E-3, 4.9118750E-3, 
      4.9909375E-3, 5.0304688E-3, 5.0502344E-3, 5.0601172E-3, 5.0650586E-3, 
      5.0675293E-3, 5.0687647E-3, 5.0693824E-3, 5.0696912E-3, 5.0698456E-3, 
      5.0699228E-3, 5.0699614E-3, 5.0699807E-3, 5.0699904E-3, 5.0699952E-3, 
      5.0699976E-3, 5.0699988E-3, 5.0699994E-3, 5.0699997E-3, 5.0699999E-3, 
      5.0700000E-3, 1.2027500E-2, 1.8985000E-2, 3.2900000E-2, 3.5430000E-2, 
      6.3260000E-2, 6.5790000E-2, 7.9705000E-2, 8.6662500E-2, 9.0141250E-2, 
      9.1880625E-2, 9.2750313E-2, 9.3185157E-2, 9.3402579E-2, 9.3511290E-2, 
      9.3565645E-2, 9.3592823E-2, 9.3606412E-2, 9.3613206E-2, 9.3616603E-2, 
      9.3618302E-2, 9.3619151E-2, 9.3619576E-2, 9.3619788E-2, 9.3619894E-2, 
      9.3619947E-2, 9.3619974E-2, 9.3619987E-2, 9.3619994E-2, 9.3619997E-2, 
      9.3619999E-2, 9.3620000E-2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
      0, 0, 0, 0, 0, 0, 0, 0 },
      correct_s2 {76.0, 0, 132.2128157, 188.1375303, 268.3594339, 384.012189, 
      551.757763, 796.593140, 1155.23014, 1676.67620, 2405.34898, 3310.45735, 
      4173.39148, 4692.06610, 4827.23521, 4759.75492, 4635.59662, 4514.75702, 
      4410.85339, 4318.55264, 4222.68663, 4155.45952, 4125.82928, 4111.62895, 
      4104.65534, 4101.19750, 4099.47551, 4098.61622, 4098.18699, 4097.97249, 
      4097.86526, 4097.81165, 4097.78485, 4097.77145, 4097.76475, 4097.76140, 
      4097.75972, 4097.75889, 4097.75846, 4097.75826, 4097.75815, 4097.75810, 
      4097.75807, 4097.75806, 4097.75805, 1629.03954, 2159.93218, 2361.55340, 
      2435.25405, 2429.23544, 1901.68147, 1864.19602, 1606.93672, 1457.55800, 
      1380.48605, 1341.66943, 1322.22763, 1312.50278, 1307.63990, 1305.20842, 
      1303.99268, 1303.38481, 1303.08087, 1302.92891, 1302.85293, 1302.81493, 
      1302.79594, 1302.78644, 1302.78169, 1302.77932, 1302.77814, 1302.77753, 
      1302.77724, 1302.77709, 1302.77702, 1302.77697, 13.1702091, 0, 0, 0, 
      0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}; 
  
      sigu(e,u,tev,alphas,betas,sab,tolin,az,iinc,lat,lasym,boundXsVec,teff,s3);
    
      for ( size_t i = 0; i < correct_s1.size(); ++i ){
        REQUIRE( correct_s1[i] == Approx(s3[2*i]).epsilon(1e-6) );
        REQUIRE( correct_s2[i] == Approx(s3[2*i+1]).epsilon(1e-6) );

      }
    } // WHEN
  
  
    WHEN( "u = 0.9" ){
      u = 0.9;
      std::vector<double> s1(nemax,0.0),s2(nemax,0.0),s3(2*nemax,0.0),
      correct_s1 { 401.1051298, 0, 9.6893315E-9, 1.9378663E-8, 3.8757325E-8, 
      7.7514650E-8, 1.5502930E-7, 3.1005860E-7, 6.2011720E-7, 1.2402344E-6, 
      2.4804688E-6, 4.9609375E-6, 7.4414063E-6, 9.9218750E-6, 1.9843750E-5, 
      3.9687500E-5, 7.9375000E-5, 1.5875000E-4, 3.1750000E-4, 6.3500000E-4, 
      1.2700000E-3, 2.5400000E-3, 3.8050000E-3, 4.4375000E-3, 4.7537500E-3, 
      4.9118750E-3, 4.9909375E-3, 5.0304688E-3, 5.0502344E-3, 5.0601172E-3, 
      5.0650586E-3, 5.0675293E-3, 5.0687647E-3, 5.0693824E-3, 5.0696912E-3, 
      5.0698456E-3, 5.0699228E-3, 5.0699614E-3, 5.0699807E-3, 5.0699904E-3, 
      5.0699952E-3, 5.0699976E-3, 5.0699988E-3, 5.0699994E-3, 5.0699997E-3, 
      5.0699999E-3, 5.0700000E-3, 1.2027500E-2, 1.8985000E-2, 3.2900000E-2, 
      3.5430000E-2, 6.3260000E-2, 6.5790000E-2, 7.9705000E-2, 8.6662500E-2, 
      9.0141250E-2, 9.1880625E-2, 9.2750313E-2, 9.3185157E-2, 9.3402579E-2, 
      9.3511290E-2, 9.3565645E-2, 9.3592823E-2, 9.3606412E-2, 9.3613206E-2, 
      9.3616603E-2, 9.3618302E-2, 9.3619151E-2, 9.3619576E-2, 9.3619788E-2, 
      9.3619894E-2, 9.3619947E-2, 9.3619974E-2, 9.3619987E-2, 9.3619994E-2, 
      9.3619997E-2, 9.3619999E-2, 9.3620000E-2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
      0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, 
      correct_s2 {77.0, 0, 133.9439956, 191.6950879, 275.7527925, 399.6227679, 
      585.461500, 871.654745, 1329.60870, 2104.23041, 3513.53749, 6165.71807, 
      8245.97898, 9331.69917, 8793.64976, 7084.87261, 5991.50897, 5346.28423, 
      4948.15187, 4687.25293, 4501.21013, 4347.12061, 4256.30937, 4219.18772, 
      4201.86297, 4193.45126, 4189.30235, 4187.24148, 4186.21438, 4185.70164, 
      4185.44548, 4185.31745, 4185.25345, 4185.22145, 4185.20545, 4185.19745, 
      4185.19345, 4185.19145, 4185.19045, 4185.18995, 4185.18970, 4185.18958, 
      4185.18952, 4185.18949, 4185.18947, 4185.18946, 1628.95068, 2159.85734, 
      2361.48880, 2435.21387, 2429.19384, 1901.63704, 1864.15191, 1606.89650, 
      1457.52106, 1380.45067, 1341.63485, 1322.19346, 1312.46881, 1307.60604, 
      1305.17461, 1303.95890, 1303.35104, 1303.04711, 1302.89515, 1302.81918, 
      1302.78118, 1302.76219, 1302.75268, 1302.74794, 1302.74557, 1302.74438, 
      1302.74378, 1302.74349, 1302.74333, 1302.74327, 1302.74322, 13.2249793, 
      0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
  
      sigu(e,u,tev,alphas,betas,sab,tolin,az,iinc,lat,lasym,boundXsVec,teff,s3);
  
      for ( size_t i = 0; i < correct_s1.size(); ++i ){
        REQUIRE( correct_s1[i] == Approx(s3[2*i]).epsilon(1e-6) );
        REQUIRE( correct_s2[i] == Approx(s3[2*i+1]).epsilon(1e-6) );

      }

    } // WHEN
  } // GIVEN

  GIVEN( "varying initial energies" ){
    u = -0.2;
    WHEN( "e = 1e-6" ){
      e = 1e-6;
      std::vector<double> s1(nemax,0.0),s2(nemax,0.0),s3(2*nemax,0.0),
      correct_s1 { 1254.9218718, 0, 1.099328E-12, 2.198657E-12, 4.397314E-12, 
      1.9310579E-8, 3.8616761E-8, 7.7229124E-8, 1.5445385E-7, 3.0890330E-7, 
      6.1780220E-7, 1.2356000E-6, 2.4711957E-6, 4.9423870E-6, 9.8847697E-6, 
      1.9769535E-5, 3.9539066E-5, 7.9078127E-5, 1.5815625E-4, 3.1631250E-4, 
      6.3262500E-4, 1.2652500E-3, 2.5305000E-3, 3.7957500E-3, 4.4283750E-3, 
      4.7446875E-3, 4.9028438E-3, 4.9819219E-3, 5.0214610E-3, 5.0412305E-3, 
      5.0511153E-3, 5.0560577E-3, 5.0585289E-3, 5.0597645E-3, 5.0603823E-3, 
      5.0606912E-3, 5.0608456E-3, 5.0609228E-3, 5.0609614E-3, 5.0609807E-3, 
      5.0609904E-3, 5.0609952E-3, 5.0609976E-3, 5.0609988E-3, 5.0609994E-3, 
      5.0609997E-3, 5.0609999E-3, 5.0610000E-3, 1.2018500E-2, 1.8976000E-2, 
      3.2891000E-2, 3.5421000E-2, 6.3251000E-2, 6.5781000E-2, 7.9696000E-2, 
      8.6653500E-2, 9.0132250E-2, 9.1871625E-2, 9.2741313E-2, 9.3176157E-2, 
      9.3393579E-2, 9.3502290E-2, 9.3556645E-2, 9.3583823E-2, 9.3597412E-2, 
      9.3604206E-2, 9.3607603E-2, 9.3609302E-2, 9.3610151E-2, 9.3610576E-2, 
      9.3610788E-2, 9.3610894E-2, 9.3610947E-2, 9.3610974E-2, 9.3610987E-2, 
      9.3610994E-2, 9.3610997E-2, 9.3610999E-2, 9.3611000E-2, 0, 0, 0, 0, 0, 0, 
      0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
      correct_s2 { 78.0, 0, 13.8621011, 19.6022587, 27.7183467, 1772.43549, 
      2458.52006, 3371.05410, 4537.76992, 5939.40927, 7477.30734, 8978.30497, 
      10265.8546, 11246.7284, 11929.7097, 12379.4378, 12667.2685, 12848.6060, 
      12959.5179, 13019.7818, 13035.3969, 12997.7179, 12877.3240, 12745.0214, 
      12677.6859, 12643.9036, 12626.9964, 12618.5400, 12614.3112, 12612.1968, 
      12611.1395, 12610.6108, 12610.3465, 12610.2143, 12610.1483, 12610.1152, 
      12610.0987, 12610.0905, 12610.0863, 12610.0843, 12610.0832, 12610.0827, 
      12610.0825, 12610.0823, 12610.0823, 12610.0822, 12610.0822, 5147.25293, 
      6828.05843, 7466.38825, 7700.07005, 7681.12148, 6013.41052, 5894.89431, 
      5081.46862, 4609.12507, 4365.41687, 4242.67455, 4181.19733, 4150.44618, 
      4135.06917, 4127.38052, 4123.53622, 4121.61404, 4120.65295, 4120.17244, 
      4119.93218, 4119.81202, 4119.75198, 4119.72192, 4119.70692, 4119.69943, 
      4119.69568, 4119.69377, 4119.69285, 4119.69235, 4119.69214, 4119.69200, 
      41.4077054, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
  
      sigu(e,u,tev,alphas,betas,sab,tolin,az,iinc,lat,lasym,boundXsVec,teff,s3);
    
      for ( size_t i = 0; i < correct_s1.size(); ++i ){
        REQUIRE( correct_s1[i] == Approx(s3[2*i]).epsilon(1e-6) );
        REQUIRE( correct_s2[i] == Approx(s3[2*i+1]).epsilon(1e-6) );

      }

    } // WHEN 
 
    WHEN( "e = 1e-4" ){
      e = 1e-4;
      std::vector<double> s1(nemax,0.0),s2(nemax,0.0),s3(2*nemax,0.0),
      correct_s1 { 124.88344935, 0, 1.099328E-10, 2.198657E-10, 4.397314E-10, 
      3.9807404E-8, 7.9175076E-8, 1.5791042E-7, 3.1538110E-7, 6.3032247E-7, 
      1.2602052E-6, 2.5199706E-6, 5.0395014E-6, 1.0078563E-5, 2.0156687E-5, 
      4.0312935E-5, 8.0625430E-5, 1.6125042E-4, 3.2250040E-4, 6.4500037E-4, 
      1.2900003E-3, 2.5800002E-3, 3.8700001E-3, 4.5150001E-3, 4.8375001E-3, 
      4.9987501E-3, 5.0793751E-3, 5.1196876E-3, 5.1398438E-3, 5.1499219E-3, 
      5.1549610E-3, 5.1574805E-3, 5.1587403E-3, 5.1593702E-3, 5.1596851E-3, 
      5.1598426E-3, 5.1599213E-3, 5.1599607E-3, 5.1599804E-3, 5.1599902E-3, 
      5.1599951E-3, 5.1599976E-3, 5.1599988E-3, 5.1599994E-3, 5.1599997E-3, 
      5.1599999E-3, 5.1600000E-3, 1.2117500E-2, 1.9075000E-2, 3.2990000E-2, 
      3.5520000E-2, 6.3350000E-2, 6.5880000E-2, 7.9795000E-2, 8.6752500E-2, 
      9.0231250E-2, 9.1970625E-2, 9.2840313E-2, 9.3275157E-2, 9.3492579E-2, 
      9.3601290E-2, 9.3655645E-2, 9.3682823E-2, 9.3696412E-2, 9.3703206E-2, 
      9.3706603E-2, 9.3708302E-2, 9.3709151E-2, 9.3709576E-2, 9.3709788E-2, 
      9.3709894E-2, 9.3709947E-2, 9.3709974E-2, 9.3709987E-2, 9.3709994E-2, 
      9.3709997E-2, 9.3709999E-2, 9.3710000E-2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
      0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
      correct_s2 { 77.0, 0, 1.39039714, 1.96614629, 2.78020563, 26.3531200, 
      37.0982545, 52.2515719, 73.5485962, 103.349484, 144.776589, 201.757319, 
      278.759802, 379.804463, 506.209016, 653.132044, 807.200426, 949.556628, 
      1064.57455, 1146.30766, 1196.72498, 1219.44080, 1220.08136, 1217.65266, 
      1216.04915, 1215.17061, 1214.71414, 1214.48183, 1214.36468, 1214.30586, 
      1214.27639, 1214.26164, 1214.25426, 1214.25057, 1214.24872, 1214.24780, 
      1214.24734, 1214.24711, 1214.24699, 1214.24693, 1214.24690, 1214.24689, 
      1214.24688, 1214.24688, 1214.24688, 1214.24688, 519.799730, 685.659904, 
      748.622707, 771.188609, 769.209074, 601.835764, 589.956876, 508.483908, 
      461.195382, 436.800140, 424.514302, 418.360927, 415.283021, 413.743934, 
      412.974378, 412.589604, 412.397213, 412.301018, 412.252924, 412.228877, 
      412.216850, 412.210840, 412.207831, 412.206330, 412.205580, 412.205205, 
      412.205014, 412.204922, 412.204872, 412.204851, 412.204837, 4.11432357, 
      0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
  
      sigu(e,u,tev,alphas,betas,sab,tolin,az,iinc,lat,lasym,boundXsVec,teff,s3);
    
      for ( size_t i = 0; i < correct_s1.size(); ++i ){
        REQUIRE( correct_s1[i] == Approx(s3[2*i]).epsilon(1e-6) );
        REQUIRE( correct_s2[i] == Approx(s3[2*i+1]).epsilon(1e-6) );

      }

    } // WHEN 

    WHEN( "e = 1e-2" ){
      e = 1e-2;
      std::vector<double> s1(nemax,0.0),s2(nemax,0.0),s3(2*nemax,0.0),
      correct_s1 { 15.036463452, 0, 1.0993286E-8, 2.1986571E-8, 4.3973142E-8, 
      9.5583302E-7, 1.8676929E-6, 3.6914126E-6, 7.3388521E-6, 1.4633731E-5, 
      2.9223488E-5, 5.8403002E-5, 1.1676203E-4, 2.3348009E-4, 4.6691621E-4, 
      9.3378844E-4, 1.8675329E-3, 3.7350219E-3, 4.6687664E-3, 4.9022026E-3, 
      4.9313822E-3, 4.9386771E-3, 4.9395890E-3, 4.9398170E-3, 4.9399310E-3, 
      4.9399880E-3, 4.9399952E-3, 4.9400023E-3, 4.9400165E-3, 4.9400450E-3, 
      4.9405009E-3, 4.9423246E-3, 4.9459720E-3, 4.9605617E-3, 5.0189207E-3, 
      5.1356387E-3, 5.6025109E-3, 7.4699998E-3, 1.2530000E-2, 1.5060000E-2, 
      4.2890000E-2, 4.5420000E-2, 7.3250000E-2, 7.5780000E-2, 8.9695000E-2, 
      9.6652500E-2, 0.1001312500, 0.1018706300, 0.1027403200, 0.1031751600, 
      0.1033925800, 0.1035012900, 0.1035556500, 0.1035828300, 0.1035964200, 
      0.1036032100, 0.1036066100, 0.1036083100, 0.1036091600, 0.1036095800, 
      0.1036097900, 0.1036099000, 0.1036099500, 0.1036099800, 0.1036099900, 
      0.1036100000, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
      0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, 
      correct_s2 { 65.0, 0, 0.1008434, 0.1426141, 0.2016868, 0.9403072, 
      1.3143937, 1.8478061, 2.6052421, 3.6783940, 5.1968298, 7.3429621, 
      10.372033, 14.637030, 20.614582, 28.916063, 40.236303, 55.119614, 
      60.670009, 61.927750, 62.081753, 62.120144, 62.124941, 62.126140, 
      62.126739, 62.127039, 62.127077, 78.249647, 78.249689, 78.249774, 
      78.251133, 78.256567, 78.267426, 78.310721, 78.481719, 78.813489, 
      80.014906, 83.298183, 85.931060, 89.015868, 88.019219, 87.069434, 
      64.782413, 63.338579, 53.965359, 48.729252, 46.060350, 44.722878, 
      44.054517, 43.720569, 43.553669, 43.470240, 43.428527, 43.407671, 
      43.397244, 43.392034, 43.389426, 43.388121, 43.387469, 43.387147, 
      43.386986, 43.386902, 43.386863, 43.386840, 43.386832, 0.3683613, 0, 0, 
      0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 
      0, 0, 0, 0, 0, 0, 0, 0 };
  
      sigu(e,u,tev,alphas,betas,sab,tolin,az,iinc,lat,lasym,boundXsVec,teff,s3);
    
      for ( size_t i = 0; i < correct_s1.size(); ++i ){
        REQUIRE( correct_s1[i] == Approx(s3[2*i]).epsilon(1e-6) );
        REQUIRE( correct_s2[i] == Approx(s3[2*i+1]).epsilon(1e-6) );

      }

    } // WHEN 
 
    WHEN( "e = 1e-0" ){
      e = 1e-0;
      std::vector<double> s1(nemax,0.0),s2(nemax,0.0),s3(2*nemax,0.0),
      correct_s1 { 0.3267346590, 0, 1.0993286E-6, 2.1986572E-6, 4.3973143E-6, 
      6.1417312E-5, 1.1843731E-4, 2.3247731E-4, 4.6055731E-4, 9.1671731E-4, 
      1.8290373E-3, 3.6536773E-3, 7.3029572E-3, 1.4601517E-2, 2.9198636E-2, 
      5.8392874E-2, 0.11678135, 0.23355830, 0.46711220, 0.58388915, 
      0.70066610, 0.75905458, 0.81744305, 0.87583153, 0.93422000, 
      0.93675000, 0.96458000, 0.96711000, 0.99494000, 0.99747000, 
      1.00253000, 1.00506000, 1.01897500, 1.03289000, 1.03542000, 
      1.04933500, 1.06325000, 1.06578000, 1.07969500, 1.09361000 },
      correct_s2 { 39.0, 0, 1.7065742E-3, 2.4132472E-3, 3.4124188E-3, 
      1.2737714E-2, 1.7676777E-2, 2.4741426E-2, 3.4773219E-2, 4.8951511E-2, 
      6.8909261E-2, 9.6858627E-2, 0.1356552801, 0.1885250673, 0.2574542065, 
      0.3369416467, 0.3952344752, 0.3477950125, 0.1468085916, 8.2364473E-2, 
      4.3525429E-2, 3.1085044E-2, 2.1980034E-2, 1.5404491E-2, 1.0710479E-2, 
      1.0541368E-2, 8.8407164E-3, 8.6997621E-3, 7.2839529E-3, 7.1667508E-3, 
      6.4162973E-3, 5.8383688E-3, 3.4732670E-3, 2.0655301E-3, 1.8792206E-3, 
      1.1171059E-3, 6.6384237E-4, 6.0388366E-4, 3.5871988E-4, 2.1301905E-4 };//, 
  
      sigu(e,u,tev,alphas,betas,sab,tolin,az,iinc,lat,lasym,boundXsVec,teff,s3);
    
      for ( size_t i = 0; i < correct_s1.size(); ++i ){
        REQUIRE( correct_s1[i] == Approx(s3[2*i]).epsilon(1e-6) );
        REQUIRE( correct_s2[i] == Approx(s3[2*i+1]).epsilon(1e-6) );

      }

    } // WHEN 
 
  } // GIVEN



  GIVEN( "Realistic-ish example" ){
      
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
    int imax = 20, lat = 0, iinc = 2, lasym = 0;
    int nbin = 2;
    std::vector<double> boundXsVec {sigma_b,sigma_b2};


    u = -1.0;
    e = 1e-5;
    int nemax = 1000;
    std::vector<double> s1(nemax,0.0),s2(nemax,0.0),s3(2*nemax,0.0);

    WHEN( " " ){
      sigu(e,u,tev,alphas,betas,sab,tol,az,iinc,lat,lasym,boundXsVec,teff,s3);
      std::vector<double> 
      correct_s1 { 1961.39420, 0.0, 4.3092008E-13, 8.6184015E-13, 1.7236803E-12, 
      1.9500396E-8, 3.8999068E-8, 7.7996412E-8, 1.5599110E-7, 3.1198047E-7, 
      6.2395921E-7, 1.2479167E-6, 2.4958317E-6, 4.9916616E-6, 9.9833214E-6, 
      1.9966641E-5, 3.9933281E-5, 7.9866561E-5, 1.5973312E-4, 3.1946623E-4, 
      6.3893245E-4, 1.2778649E-3, 2.5557298E-3, 5.1114595E-3, 5.1115666E-3, 
      5.1116736E-3, 5.1118877E-3, 5.1123158E-3, 5.1131721E-3, 5.1148846E-3, 
      5.1183097E-3, 5.1251598E-3, 5.1388600E-3, 5.1662604E-3, 5.2210612E-3, 
      5.3306629E-3, 5.5498662E-3, 5.9882729E-3, 6.8650862E-3, 8.6187128E-3, 
      0.0121259660, 0.0156332200, 0.0191404730, 0.0226477270, 0.0261549800, 
      0.0296622340, 0.0331694870, 0.0357202170, 0.0497492310, 0.0637782440, 
      0.0663289740, 0.0803579880, 0.0873724950, 0.0943870010, 0.0 }, 
      correct_s2 { 53.0, 0.0, 29.350193, 41.503882, 58.688237, 5981.4375, 8313.5919, 
      11478.328, 15706.137, 21237.514, 28280.161, 36942.720, 47156.023, 58611.460, 
      70756.269, 82867.650, 94178.612, 103979.62, 111609.69, 116279.17, 116709.48, 
      110675.18, 94945.727, 67280.399, 48540.369, 48539.760, 48538.540, 48536.102, 
      48531.225, 48521.472, 48501.971, 48462.986, 48385.084, 48229.555, 47919.602, 
      47304.192, 46091.829, 43743.764, 39364.949, 31851.393, 20992.188, 14047.844, 
      9558.0062, 6607.6166, 4636.2229, 3248.6850, 2341.1063, 2272.6275, 2184.8212, 
      1574.4957, 1463.6744, 832.26113, 576.51928, 380.64264, 0.0 };
      for ( size_t i = 0; i < correct_s1.size(); ++i ){
        REQUIRE( correct_s1[i] == Approx(s3[2*i]).epsilon(1e-6) );
        REQUIRE( correct_s2[i] == Approx(s3[2*i+1]).epsilon(1e-6) );

      }


    } // WHEN
  } // GIVEN



  GIVEN( "high tolerance, energy, and cosine (debugging)" ){
    int imax = 20, lat = 0, iinc = 2, lasym = 0;
    double tev = 2.5507297688e-2, tol = 5.0;
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
    std::vector<double> correctEnergies;
    std::vector<double> correctSCR;
    std::vector<double> boundXsVec {sigma_b,sigma_b2};

    double enow = 8.6e-1;
    double u = 1.0;

    int nemax = 5000;
    std::vector<double> s1(nemax,0.0),s2(nemax,0.0), s3(nemax*2,0.0);
  
    auto out = sigu(enow,u,tev,alphas,betas,sab,tol,az,iinc,lat,lasym,boundXsVec,teff,s3);
  
    std::vector<double> correctYU_first_110 { 44.6334456177, 54.0, 0, 0, 
    0.765623, 5.979297, 0.79368103, 20.80756, 0.79623176, 20.28317, 0.82428978, 
    1.692838, 0.82684051, 0.7111114, 0.84086953, 22.83815, 0.84788404, 187.8350, 
    0.85489854, 2614.586, 0.85489855, 0, 0.85489856, 0, 0.85489858, 0, 0.85489862, 
    0, 0.85489870, 0, 0.85489886, 0, 0.85489917, 0, 0.85489979, 0, 0.85490104, 
    0, 0.85490353, 0, 0.85490851, 0, 0.85491847, 0, 0.85493840, 0, 0.85497826, 
    0, 0.85505797, 0, 0.85521739, 0, 0.85553623, 0, 0.85617391, 0, 0.85744927, 
    0, 0.86255073, 0, 0.86382610, 0, 0.86446378, 0, 0.86478262, 0, 0.86494204, 
    0, 0.86502175, 0, 0.86506161, 0, 0.86508154, 0, 0.86509150, 0, 0.86509648, 
    0, 0.86509897, 0, 0.86510022, 0, 0.86510084, 0, 0.86510115, 0, 0.86510131, 
    0, 0.86510139, 0, 0.86510143, 0, 0.86510145, 0, 0.86510146, 2153.377, 
    0.87211597, 117.9912, 0.87913048, 10.89082, 0.89315949, 0.1927072, 0.89571022,
    0.4177208, 0.92376824, 1.712638, 0.92631897, 1.591015, 0.954377, 0.1488137 };
    for (size_t i = 0; i < correctYU_first_110.size(); ++i){
      REQUIRE( correctYU_first_110[i] == Approx(out[i]).epsilon(1e-6) );
    }
    for (size_t i = correctYU_first_110.size(); i < out.size(); ++i){
      REQUIRE( 0.0 == Approx(out[i]).epsilon(1e-6) );
    }



  } // GIVEN

} // TEST CASE
*/
}}
