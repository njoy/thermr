//#include <Eigen/Dense>
#include <iostream> 
#include <vector>
#include <cmath>
#include <tuple>
#include <string>
#include <time.h>
#include "ENDFtk.hpp"


int thermr( ){

  //clock_t t;
  //t = clock();
  using namespace njoy::ENDFtk;

  // convenience typedefs
  using ReactionProduct = 
  section::Type< 6 >::ReactionProduct;
//  using Multiplicity = 
//  section::Type< 6 >::Multiplicity;
//  using Distribution = 
//  section::Type< 6 >::Distribution;
  using ContinuumEnergyAngle = 
  section::Type< 6 >::ContinuumEnergyAngle;
  using LegendreCoefficients = 
  section::Type< 6 >::ContinuumEnergyAngle::LegendreCoefficients;
  
  //std::string chunk();
  /*
  std::string invalidLAW;
  std::string validSEND;
  std::string invalidSEND;
  */

    int zaid = 92235;
    double awr = 2.330248e+2;
    int jp = 0;
    int lct = 6; // Temporary flag for incoherent inelastic data
    int mt = 5;  // This is the number of product subsections
    std::vector< ReactionProduct > products = {
      ReactionProduct(
        // multiplicity
	// ZAP     AWP    LIP  LAW  NP   Interpolants
        { 1001., 0.9986234, 0, 1, { 4 }, { 2 },
	// Energies Vector
          { 1e-5, 1.1e+7, 1.147e+7, 2e+7 },
	// Multiplicities Vector
          { 0., 8.45368e-11, 6.622950e-8, 2.149790e-1 } },

        // distribution
        { ContinuumEnergyAngle(
	// LAW LEP  NE    NR
            1, 2, { 2 }, { 1 },
            { LegendreCoefficients(
	//  1st Energy ND LANG NEP    (note that NW=12 b/c length of next vec)
                  1e-5, 0, 1, 4,
	// These are ordered Energy, Coeff1, Coeff2, Energy, Coeff1, Coeff2 so 
	// that E1=1.0, E2=4.0, etc.
                  { 1., 2., 3., 4., 5., 6., 7., 8., 9., 10., 11., 12. } ),
              LegendreCoefficients(
        // 2nd Energy  ND LANG NEP   (NW=6  b/c length of vec)
                  2e+7, 0, 1, 2, { 1., 2., 3., 4., 5., 6.} ) } ) } ),
      ReactionProduct( 
        // multiplicity
        { 1., 1., 0, 1, { 2 }, { 2 },
          { 1.858639e+7, 2.e+7 },
          { 4., 4. } },
        // distribution
        { ContinuumEnergyAngle(
            1, 2, { 2 }, { 22 },
            { LegendreCoefficients(
                  1.858639e+7, 0, 0, 3, { 0., 0., 0.5, 2., 1., 0. } ),
              LegendreCoefficients(
                  2e+7, 0, 0, 3, { 0., 0., 0.5, 2., 1., 0. } ) } ) } ),
      ReactionProduct( 
        // multiplicity
        { 0., 0., 0, 1, { 3 }, { 2 },
          { 1.858639e+7, 1.9e+7, 2.e+7 },
          { 1., 2., 3. } },
        // distribution
        { ContinuumEnergyAngle(
          1, 2, { 2 }, { 5 },
          { LegendreCoefficients(
                 1.858639e+7, 0, 0, 3, { 0., 0., 1., 1., 2., 0. } ),
            LegendreCoefficients(
                 2e+7, 0, 0, 3, { 0., 0., 1., 1., 2., 0. } ) } ) } ) };



      section::Type< 6 > chunk( mt, zaid, awr, jp, lct, std::move( products ) );
      std::cout << chunk.ZA() << std::endl;
      std::cout << chunk.AWR() << std::endl;
      std::cout << chunk.JP() << std::endl;
      std::cout << chunk.LCT() << std::endl;
      std::cout << chunk.NK() << std::endl;
      std::cout << chunk.MT() << std::endl;
  
 std::cout << "hello, world" << std::endl;

  return 0;


}

