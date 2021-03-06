#define CATCH_CONFIG_MAIN
#include "catch.hpp"
#include "incoherentElastic/incoherentElastic.h"
#include <range/v3/all.hpp>
#include "generalTools/testing.h"

TEST_CASE( "get equiprobable scattering angles" ){
  GIVEN( "an initial energy (i.e. final energy), constants, and cosine" ){
    double debyeWaller = 10.0, sigma_b = 2.5, xs;
    int numAtoms = 1, numAngles = 8;
    std::vector<double> cosines  (numAngles+2, 1.0),
                        energies {1e-5, 1e-4, 1e-3, 1e-2, 1e-1};
    std::vector<std::vector<double>> correctCosinesVec {
    {1e-5,1.0,-0.874977,-0.624939,-0.374914,-0.1249021,0.125097,0.375085,0.625060,0.875022},
    {1e-4,1.0,-0.874770,-0.624395,-0.374145,-0.1240207,0.125979,0.375853,0.625603,0.875228},
    {1e-3,1.0,-0.872682,-0.618908,-0.366416,-0.1151928,0.134774,0.383499,0.630992,0.877266},
    {1e-2,1.0,-0.849361,-0.559611,-0.285741,-0.0260998,0.220719,0.455923,0.680557,0.895529},
    {1e-1,1.0,-0.327071, 0.193319, 0.435175, 0.5967253,0.718423,0.816144,0.897817,0.967984},
    };
    THEN( "the returned cross sections are correct for each E and cosine" ){
      for ( size_t i = 0; i < energies.size(); ++i ){
        getIncohElasticDataSingleEnergy(energies[i],debyeWaller,cosines);
        REQUIRE( ranges::equal( cosines, correctCosinesVec[i], equal ) );
      }
    } // THEN
  } // GIVEN
} // TEST CASE



/*
TEST_CASE( "incoherent elastic average scattering xs" ){
  GIVEN( "Material constants and neutron energy" ){
    double debyeWaller = 10.0, sigma_b = 2.5, xs;
    int numAtoms = 1;

    std::vector<double> 
      energyVec { 1e-5, 1e-4, 1e-3, 1e-2, 1e-1 },
      correctXS { 2.4995001, 2.4950067, 2.4506601, 2.0604997, 0.6135527 };

    THEN( "incoherent elastic scattering xs (integrated over all mu) "
          "is correctly generated" ){
      for ( size_t i = 0; i < energyVec.size(); ++i ){
        REQUIRE( correctXS[i] == Approx(getIncohElasticXS(sigma_b,debyeWaller,
                                        energyVec[i],numAtoms)).epsilon(1e-6) );
      }
    } // THEN

    WHEN( "bound scattering cross section is scaled by a constant" ){
      std::vector<double> scalingVec = ranges::view::iota(2,10);
      THEN( "scattering cross section is scaled by that same value" ){
        for ( size_t i = 0; i < energyVec.size(); ++i ){
          for ( const auto& scaling : scalingVec ){
            REQUIRE( scaling*correctXS[i] == 
                     Approx(getIncohElasticXS(scaling*sigma_b,debyeWaller,
                                    energyVec[i], numAtoms)).epsilon(1e-6) );
          }
        }
      } // THEN
    } // WHEN

    WHEN( "polyethylene values used" ){
      debyeWaller = 34.957;
      sigma_b = 162.88;
      numAtoms = 1;
  
      energyVec = { 0.00001, 0.0000178, 0.000025, 0.000035, 0.00005, 0.00007, 
        0.0001, 0.000126, 0.00016, 0.0002, 0.000253, 0.000297, 0.00035, 0.00042, 
      0.000506, 0.000615, 0.00075, 0.00087, 0.001012, 0.00123, 0.0015, 0.0018, 
      0.00203, 0.002277, 0.0026, 0.003, 0.0035, 0.004048, 0.0045, 0.005, 0.0056, 
      0.006325, 0.0072, 0.0081, 0.009108, 0.01, 0.01063, 0.0115, 0.012397, 
      0.0133, 0.01417, 0.015, 0.016192, 0.0182, 0.0199, 0.020493, 0.0215, 0.0228, 
      0.0253, 0.028, 0.030613, 0.0338, 0.0365, 0.0395, 0.042757, 0.0465, 0.05, 
      0.056925, 0.0625, 0.069, 0.075, 0.081972, 0.09, 0.096, 0.1035, 0.111573, 
      0.12, 0.128, 0.1355, 0.145728, 0.16, 0.172, 0.184437, 0.2, 0.2277, 
      0.2510392, 0.2705304, 0.2907501, 0.3011332, 0.3206421, 0.3576813, 0.39, 
      0.4170351, 0.45, 0.5032575, 0.56, 0.625, 0.7 };

      correctXS = { 162.766177, 162.677469, 162.595642, 162.482084, 162.311945, 
      162.085463, 161.746530, 161.453553, 161.071497, 160.623565, 160.032615, 
      159.544222, 158.958574, 158.189480, 157.251404, 156.073138, 154.630197, 
      153.362606, 151.880641, 149.642888, 146.932757, 143.999042, 141.803628, 
      139.496539, 136.556492, 133.032434, 128.801762, 124.376802, 120.885466, 
      117.181960, 112.946452, 108.114287, 102.669377, 97.475420, 92.104811, 
      87.711433, 84.795662, 81.005013, 77.362184, 73.945169, 70.871378, 
      68.123109, 64.464103, 58.980092, 54.913544, 53.604411, 51.499027, 
      48.982806, 44.702911, 40.772715, 37.524707, 34.157953, 31.720115, 
      29.372364, 27.174729, 25.013159, 23.275766, 20.455913, 18.634770, 
      16.880935, 15.531029, 14.210309, 12.942841, 12.133937, 11.254677, 
      10.440335, 9.707164, 9.100466, 8.596750, 7.993383, 7.280373, 6.772440, 
      6.315759, 5.824298, 5.115765, 4.640151, 4.305837, 4.006395, 3.868254, 
      3.632897, 3.256697, 2.986820, 2.793193, 2.588577, 2.314639, 2.080107, 
      1.863775, 1.664085 };

      THEN( "incoherent elastic scattering xs (integrated over all mu) "
            "is correctly generated" ){
        for ( size_t i = 0; i < energyVec.size(); ++i ){
          REQUIRE( correctXS[i] == Approx(getIncohElasticXS(sigma_b,debyeWaller,
                                          energyVec[i],numAtoms)).epsilon(1e-6) );
        }
      } // THEN

    } // WHEN
  } // GIVEN
} // TEST CASE

*/



/*
TEST_CASE( "convert some X, Y grid to be on a new X grid" ){
  std::vector<double> currentX(88), currentY(88), desiredX(153), desiredY(153), correctY(153);
    currentX = { 1E-5, 1.78E-5, 2.5E-5, 3.5E-5, 5E-5, 7E-5, 1E-4, 1.26E-4, 
    1.6E-4, 2.0E-4, 2.53E-4, 2.97E-4, 3.5E-4, 4.2E-4, 5.06E-4, 6.15E-4, 7.5E-4, 
    8.7E-4, 1.012E-3, 1.23E-3, 1.5E-3, 1.8E-3, 2.03E-3, 2.277E-3, 2.6E-3, 3.0E-3, 
    3.5E-3, 4.048E-3, 4.5E-3, 5.0E-3, 5.6E-3, 6.325E-3, 7.2E-3, 8.1E-3, 9.108E-3, 
    .01, .01063, .0115, .012397, .0133, .01417, .015, .016192, .0182, .0199, 
    .020493, .0215, .0228, .0253, .028, .030613, .0338, .0365, .0395, .042757, 
    .0465, .05, .056925, .0625, .069, .075, .081972, .09, .096, .1035, .111573, 
    .12, .128, .1355, .145728, .16, .172, .184437, .2, .2277, .2510392, .2705304, 
    .2907501, .3011332, .3206421, .3576813, .39, .4170351, .45, .5032575, .56, 
    .625, .7 };
    desiredX = {1E-5, 1.0625E-5, 1.125E-5, 1.1875E-5, 1.25E-5, 1.375E-5, 1.5E-5, 
    1.625E-5, 1.75E-5, 1.875E-5, 2.E-5, 2.1875E-5, 2.375E-5, 2.5625E-5, 2.75E-5, 
    2.9375E-5, 3.125E-5, 3.3125E-5, 3.5E-5, 3.875E-5, 4.25E-5, 4.625E-5, 5.E-5, 
    5.3125E-5, 5.625E-5, 5.9375E-5, 6.25E-5, 6.875E-5, 7.5E-5, 8.125E-5,8.75E-5, 
    9.375E-5, 1.E-4, 1.0625E-4, 1.125E-4, 1.1875E-4, 1.25E-4, 1.375E-4, 1.5E-4, 
    1.625E-4, 1.75E-4, 1.875E-4, 2.E-4, 2.1875E-4, 2.375E-4, 2.5625E-4, 2.75E-4, 
    2.9375E-4, 3.125E-4, 3.3125E-4, 3.5E-4, 3.875E-4, 4.25E-4, 4.625E-4, 5.E-4, 
    5.3125E-4, 5.625E-4, 5.9375E-4, 6.25E-4, 6.875E-4, 7.5E-4, 8.125E-4,8.75E-4, 
    9.375E-4, 1.E-3, 1.0625E-3, 1.125E-3, 1.1875E-3, 1.25E-3, 1.375E-3, 1.5E-3, 
    1.625E-3, 1.75E-3, 1.875E-3, 2.E-3, 2.1875E-3, 2.375E-3, 2.5625E-3, 2.75E-3, 
    2.9375E-3, 3.125E-3, 3.3125E-3, 3.5E-3, 3.875E-3, 4.25E-3, 4.625E-3, 5.E-3, 
    5.3125E-3, 5.625E-3, 5.9375E-3, 6.25E-3, 6.875E-3, 7.5E-3,8.125E-3, 8.75E-3, 
    9.375E-3, 1.E-2, 1.0625E-2, 1.125E-2, 1.1875E-2, 1.25E-2, 1.375E-2, 1.5E-2, 
    1.625E-2, 1.75E-2, 1.875E-2, 2.E-2, 2.1325E-2, 2.265E-2, 2.3975E-2, 2.53E-2, 
    2.684375E-2, 2.83875E-2,3.1475E-2, 3.45625E-2,3.765E-2,4.07375E-2,4.3825E-2, 
    4.69125E-2, 5.E-2, 5.3125E-2, 5.625E-2, 5.9375E-2, 6.25E-2, 6.875E-2,7.5E-2, 
    8.125E-2, 8.75E-2, 9.375E-2, 1.E-1, 1.0625E-1, 1.125E-1, 1.1875E-1, 1.25E-1, 
    1.375E-1, 1.5E-1, 1.625E-1, 1.75E-1, 1.875E-1, 2.E-1, 2.1875E-1, 2.375E-1, 
    2.5625E-1, 2.75E-1, 2.9375E-1, 3.125E-1, 3.3125E-1, 3.5E-1,3.875E-1,4.25E-1, 
    5.E-1, 6.25E-1, 6.250062E-1 };
    correctY = { 162.7661, 162.7590, 162.7519, 162.7448, 162.7377, 162.7235, 
    162.7093, 162.6950, 162.6808, 162.6666, 162.6524, 162.6311, 162.6098, 
    162.5885, 162.5672, 162.5459, 162.5246, 162.5033, 162.4820, 162.4395, 
    162.3969, 162.3544, 162.3119, 162.2765, 162.2411, 162.2057, 162.1703, 
    162.0996, 162.0289, 161.9582, 161.8876, 161.8170, 161.7465, 161.6760, 
    161.6055, 161.5351, 161.4648, 161.3241, 161.1837, 161.0434, 160.9033, 
    160.7633, 160.6235, 160.4141, 160.2051, 159.9964, 159.7881, 159.5802, 
    159.3726, 159.1654, 158.9585, 158.5459, 158.1347, 157.7249, 157.3166, 
    156.9774, 156.6391, 156.3019, 155.9656, 155.2960, 154.6301, 153.9682, 
    153.3100, 152.6557, 152.0051, 151.3583, 150.7151, 150.0756, 149.4398, 
    148.1792, 146.9327, 145.7007, 144.4824, 143.2781, 142.0874, 140.3266, 
    138.5956, 136.8935, 135.2204, 133.5748, 131.9574, 130.3666, 128.8017, 
    125.7509, 122.7996, 119.9447, 117.1819, 114.9492, 112.7749, 110.6608, 
    108.6006, 104.6462, 100.8962, 97.33701, 93.96291, 90.75717, 87.71143, 
    84.81823, 82.06832, 79.45196, 76.96076, 72.33060, 68.12310, 64.29536, 
    60.80814, 57.61356, 54.68917, 51.85511, 49.26298, 46.89451, 44.70291, 
    42.38435, 40.26329, 36.56183, 33.43830, 30.78347, 28.50265, 26.52597, 
    24.79749, 23.27576, 21.92675, 20.70359, 19.61954, 18.63477, 16.94278, 
    15.53102, 14.33758, 13.31482, 12.42628, 11.65028, 10.96497, 10.35495, 
    9.809819, 9.319582, 8.472700, 7.767551, 7.169103, 6.657084, 6.213723, 
    5.824298, 5.327597, 4.906062, 4.546385, 4.236182, 3.965564, 3.728002, 
    3.517853, 3.328842, 3.006203, 2.741299, 2.329899, 1.863775, 0.0 };


  WHEN( "directly specify the current Y vector" ){
    currentY = { 162.7661, 162.6774, 162.5956, 162.4820, 162.3119, 162.0854, 
    161.7465, 161.4535, 161.0714, 160.6235, 160.0326, 159.5442, 158.9585, 
    158.1894, 157.2514, 156.0731, 154.6301, 153.3626, 151.8806, 149.6428, 
    146.9327, 143.9990, 141.8036, 139.4965, 136.5564, 133.0324, 128.8017, 
    124.3768, 120.8854, 117.1819, 112.9464, 108.1142, 102.6693, 97.47542, 
    92.10481, 87.71143, 84.79566, 81.00501, 77.36218, 73.94516, 70.87137, 
    68.12310, 64.46410, 58.98009, 54.91354, 53.60441, 51.49902, 48.98280, 
    44.70291, 40.77271, 37.52470, 34.15795, 31.72011, 29.37236, 27.17472, 
    25.01315, 23.27576, 20.45591, 18.63477, 16.88093, 15.53102, 14.21030, 
    12.94284, 12.13393, 11.25467, 10.44033, 9.707164, 9.100466, 8.596750, 
    7.993383, 7.280373, 6.772440, 6.315759, 5.824298, 5.115765, 4.640151, 
    4.305837, 4.006395, 3.868254, 3.632897, 3.256697, 2.986820, 2.793193, 
    2.588577, 2.314639, 2.080107, 1.863775, 1.664085 };

    THEN( "current Y gets correctly mapped from current X to desiredX" ){
      desiredY = getOnRightGrid( currentX, currentY, desiredX );
      REQUIRE( ranges::equal(desiredY, correctY, equal) );
    } // WHEN
  } // WHEN 

  WHEN( "generate the inputs on our own" ){
    double sigma_b = 162.88, debyeWaller = 34.957;
    int numAtoms = 1;
    desiredY = getIncohElasticXSGrid_withInterpolation( sigma_b, debyeWaller, 
                                                numAtoms, currentX, desiredX );
    REQUIRE( ranges::equal(desiredY, correctY, equal) );
  } // WHEN 

  WHEN( "generate the currentY vector using incoh elastic xs function" ){
    AND_WHEN( "two scattering atoms are considered" ){
      double sigma_b = 10.0, debyeWaller = 20.0;
      int numAtoms = 2;

      correctY = {4.998001, 4.997876, 4.997751, 4.997626, 4.997501, 4.997251, 
      4.997001, 4.996751, 4.996502, 4.996252, 4.996002, 4.995628, 4.995253, 
      4.994879, 4.994504, 4.994130, 4.993755, 4.993381, 4.993007, 4.992258, 
      4.991510, 4.990761, 4.990013, 4.989390, 4.988767, 4.988144, 4.987521, 
      4.986275, 4.985030, 4.983785, 4.982541, 4.981297, 4.980053, 4.978810, 
      4.977567, 4.976325, 4.975083, 4.972601, 4.970120, 4.967640, 4.965163, 
      4.962687, 4.960212, 4.956504, 4.952799, 4.949098, 4.945401, 4.941708, 
      4.938018, 4.934331, 4.930649, 4.923295, 4.915955, 4.908630, 4.901320, 
      4.895239, 4.889169, 4.883108, 4.877058, 4.864987, 4.852956, 4.840964, 
      4.829013, 4.817101, 4.805228, 4.793395, 4.781601, 4.769846, 4.758129, 
      4.734813, 4.711648, 4.688638, 4.665778, 4.643068, 4.620507, 4.586943, 
      4.553708, 4.520798, 4.488212, 4.455940, 4.423988, 4.392344, 4.361005, 
      4.299246, 4.238675, 4.179269, 4.120999, 4.073303, 4.026356, 3.980169, 
      3.934697, 3.845926, 3.759922, 3.676574, 3.595837, 3.517568, 3.441694, 
      3.368148, 3.296844, 3.227693, 3.160609, 3.032415, 2.911691, 2.797969, 
      2.690802, 2.589601, 2.494076, 2.398626, 2.308712, 2.224073, 2.143958, 
      2.056533, 1.974520, 1.825758, 1.694526, 1.578476, 1.475381, 1.383444, 
      1.301102, 1.227105, 1.160164, 1.098876, 1.043754, 0.993262, 0.905398, 
      0.831268, 0.768123, 0.713744, 0.666356, 0.624878, 0.588199, 0.555521, 
      0.526302, 0.500015, 0.454591, 0.416762, 0.384654, 0.357182, 0.333394, 
      0.312500, 0.285850, 0.263232, 0.243934, 0.227290, 0.212770, 0.200024, 
      0.188749, 0.178607, 0.161296, 0.147083, 0.125010, 0.100000, 0.000000};
      desiredY = getIncohElasticXSGrid_withInterpolation( sigma_b, debyeWaller, 
                 numAtoms, currentX, desiredX );
      REQUIRE( ranges::equal(desiredY, correctY, equal) );

    } // AND WHEN
  } // WHEN
} // TEST CASE




*/












