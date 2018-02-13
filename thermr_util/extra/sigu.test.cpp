#define CATCH_CONFIG_MAIN
#include "../../catch.hpp"
#include "sigu.h"


TEST_CASE( "computeCrossSections" ){
  GIVEN( "" ){
    std::vector<double> alpha { 2.25, 2.30, 2.35, 2.40 }, beta { 0.02, 0.2, 1.2 };
    double enow = 1.5, tev = 0.035, tevz = 0.0253, tolin = 5e-3, az = 11.898, u = -1.0;
    std::vector<std::vector<double>> sab { { 0.1, 0.2, 0.3 }, { 1.4, 1.5, 1.6 }, { 2.1, 2.2, 2.3, 2.4 }, { 3.4, 3.5, 3.6 } };
    int nemax = 5, lasym = 0, lat = 0;
    
    WHEN( "temperature is reasonable" ){ 
    // ( 111 113 116 150 160 ) repeated a few times, then 170
    {
    auto out = sigu( enow, u, tev, tevz, alpha, beta, sab, tolin, az, nemax, lasym, lat );
    std::vector<double> s = std::get<2>(out), x = std::get<0>(out), y = std::get<1>(out);

    std::vector<double> outS { 0, 4, 0, 0, 1.0708783, 0, 1.493, 0, 1.4993 };
    std::vector<double> outX { 1.5007, 1.4993, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };

    for ( size_t i = 0; i < outS.size(); ++i ){ REQUIRE( outS[i] == Approx( s[i] ).epsilon(1e-6) ); }
    for ( size_t i = 0; i < outX.size(); ++i ){ REQUIRE( 0.0     == Approx( y[i] ).epsilon(1e-6) ); }
    for ( size_t i = 0; i < outX.size(); ++i ){ REQUIRE( 0.0     == Approx( y[i] ).epsilon(1e-6) ); }
    }

    // ( 111 113 116 150 160 ) repeated a few times, then 170
    {
    enow = 3.0; u = -1.25; tev = 0.0055;

    auto out = sigu( enow, u, tev, tevz, alpha, beta, sab, tolin, az, nemax, lasym, lat );
    std::vector<double> s = std::get<2>(out), x = std::get<0>(out), y = std::get<1>(out);

    std::vector<double> outS { 0, 4, 0, 0, 2.0536991, 0, 2.9989, 0, 2.99989, 0 };
    std::vector<double> outX {   3.00011, 2.99989, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };

    for ( size_t i = 0; i < outS.size(); ++i ){ REQUIRE( outS[i] == Approx( s[i] ).epsilon(1e-6) ); }
    for ( size_t i = 0; i < outX.size(); ++i ){ REQUIRE( outX[i] == Approx( x[i] ).epsilon(1e-6) ); }
    for ( size_t i = 0; i < outX.size(); ++i ){ REQUIRE( 0.0     == Approx( y[i] ).epsilon(1e-6) ); }
    }

    // ( 111 113 116 150 160 ) repeated a few times, then 170
    {
    enow = 2.0; u = -5.95; tev = 0.0155; az = 18.43;

    auto out = sigu( enow, u, tev, tevz, alpha, beta, sab, tolin, az, nemax, lasym, lat );
    std::vector<double> s = std::get<2>(out), x = std::get<0>(out), y = std::get<1>(out);

    std::vector<double> outS { 0, 4, 0, 0, 0.94994594, 0, 1.9969, 0, 1.99969, 0 };
    std::vector<double> outX {  2.00031, 1.99969, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };

    for ( size_t i = 0; i < outS.size(); ++i ){ REQUIRE( outS[i] == Approx( s[i] ).epsilon(1e-6) ); }
    for ( size_t i = 0; i < outX.size(); ++i ){ REQUIRE( outX[i] == Approx( x[i] ).epsilon(1e-6) ); }
    for ( size_t i = 0; i < outX.size(); ++i ){ REQUIRE( 0.0     == Approx( y[i] ).epsilon(1e-6) ); }
    }

    // ( 111 113 116 150 160 ) repeated a few times, then 170
    {
     enow = 2.0; u = 5.95; tev = 0.0155; az = 18.43;

    auto out = sigu( enow, u, tev, tevz, alpha, beta, sab, tolin, az, nemax, lasym, lat );
    std::vector<double> s = std::get<2>(out), x = std::get<0>(out), y = std::get<1>(out);

    std::vector<double> outS { 0, 4, 0, 0, 1.9814, 0, 1.9969, 0, 1.99969, 0 };
    std::vector<double> outX {  2.00031, 1.99969, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };

    for ( size_t i = 0; i < outS.size(); ++i ){ REQUIRE( outS[i] == Approx( s[i] ).epsilon(1e-6) ); }
    for ( size_t i = 0; i < outX.size(); ++i ){ REQUIRE( outX[i] == Approx( x[i] ).epsilon(1e-6) ); }
    for ( size_t i = 0; i < outX.size(); ++i ){ REQUIRE( 0.0     == Approx( y[i] ).epsilon(1e-6) ); }


    }
    } // THEN

    // ( 111 113 113 116 150 160 ), then ( 111 113 116 150 160 ) repeated a few times, then 170
    WHEN( "temperature is pretty high" ){
     tev = 2.5; 
     THEN( "113 gets repeated once" ){

    auto out = sigu( enow, u, tev, tevz, alpha, beta, sab, tolin, az, nemax, lasym, lat );
    std::vector<double> s = std::get<2>(out), x = std::get<0>(out), y = std::get<1>(out);

    std::vector<double> outS { 0, 4, 0, 0, 0.999999, 0, 1.0708783, 0, 1.55, 0 };
    std::vector<double> outX {  2, 1.55, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };


    for ( size_t i = 0; i < outS.size(); ++i ){ REQUIRE( outS[i] == Approx( s[i] ).epsilon(1e-6) ); }
    for ( size_t i = 0; i < outX.size(); ++i ){ REQUIRE( outX[i] == Approx( x[i] ).epsilon(1e-6) ); }
    for ( size_t i = 0; i < y.size(); ++i ){ REQUIRE( 0.0     == Approx( y[i] ).epsilon(1e-6) ); }

     } // THEN

    } // WHEN

    // ( 111 113 113 116 150 160 ) a few times, then ( 111 113 116 150 160 ) repeated a few times, then 170
    WHEN( "temperature is really high" ){
     tev = 12.5; 
     THEN( "113 gets repeated a few times" ){

    auto out = sigu( enow, u, tev, tevz, alpha, beta, sab, tolin, az, nemax, lasym, lat );
    std::vector<double> s = std::get<2>(out), x = std::get<0>(out), y = std::get<1>(out);

    std::vector<double> outS { 0, 4, 0, 0, 1.0708783, 0, 1.75, 0, 4, 0 };
    std::vector<double> outX {16.500001, 4, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };


    for ( size_t i = 0; i < outS.size(); ++i ){ REQUIRE( outS[i] == Approx( s[i] ).epsilon(1e-6) ); }
    for ( size_t i = 0; i < outX.size(); ++i ){ REQUIRE( outX[i] == Approx( x[i] ).epsilon(1e-6) ); }
    for ( size_t i = 0; i < y.size(); ++i ){ REQUIRE( 0.0     == Approx( y[i] ).epsilon(1e-6) ); }

     } // THEN

    } // WHEN



      
    THEN( "" ){
      //REQUIRE( true );
    } // THEN
  } // GIVEN
} // TEST CASE
