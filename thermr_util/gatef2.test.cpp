#define CATCH_CONFIG_MAIN
#include "../catch.hpp"
#include "gatef2.h"

TEST_CASE( "gatef2" ){
  int mat;
  GIVEN( "only one material temperature is considered" ){
    std::vector<double> temp(1), eftemp(1), tempVals(8), effTemp(8);

    WHEN( "material considered is Benzine" ){
      mat = 1095;
      tempVals = { 200, 250, 296, 300, 315, 350, 600, 743 };
      effTemp = { 200, 250, 685.54, 685.54, 315, 712.02, 866.63, 743 };

      THEN( "the effective temperature is correctly identified in table" ){
        for ( size_t i = 0; i < tempVals.size(); ++i ){
          temp[0] = tempVals[i];
          gatef2( temp, eftemp, mat );
          REQUIRE( effTemp[i] == Approx( eftemp[0] ).epsilon(1e-6) );
          eftemp = { 0 };
        } // for
      } // THEN
    } // WHEN

    WHEN( "material considered is BeO" ){
      mat = 1099;
      tempVals = { 200, 250, 296, 300, 315, 350, 600, 743 };
      effTemp = { 200, 250, 427.8, 427.8, 315, 350, 671.3, 743 };

      THEN( "the effective temperature is correctly identified in table" ){
        for ( size_t i = 0; i < tempVals.size(); ++i ){
          temp[0] = tempVals[i];
          gatef2( temp, eftemp, mat );
          REQUIRE( effTemp[i] == Approx( eftemp[0] ).epsilon(1e-6) );
          eftemp = { 0 };
        } // for
      } // THEN
    } // WHEN
  } // GIVEN

  GIVEN( "multiple material temperatures are considered" ){
    std::vector<double> temp(4), eftemp(4);
    std::vector<std::tuple<double,double,double,double>> tempVals(4), effTemp(4);

    WHEN( "material considered is Benzine" ){
      mat = 1095;
      tempVals = { { 200.00, 220.00, 250.00, 260.00 }, 
        { 296.00, 300.00, 301.00, 405.00 }, { 300.00, 315.00, 320.00, 345.00 },
        { 400.00, 455.00, 520.00, 645.00 } };
      effTemp  = { { 200.00, 220.00, 250.00, 260.00 }, 
        { 685.54, 685.54, 685.54, 738.97 }, { 685.54, 315.00, 320.00, 712.02 },
        { 738.97, 768.10, 520.00, 645.00 } };

      THEN( "the effective temperature is correctly identified in table" ){
        for ( size_t i = 0; i < tempVals.size(); ++i ){
          temp = { std::get<0>(tempVals[i]), std::get<1>(tempVals[i]), 
                   std::get<2>(tempVals[i]), std::get<3>(tempVals[i]) };
          gatef2( temp, eftemp, mat );
          REQUIRE( std::get<0>(effTemp[i]) == Approx( eftemp[0] ).epsilon(1e-6) );
          REQUIRE( std::get<1>(effTemp[i]) == Approx( eftemp[1] ).epsilon(1e-6) );
          REQUIRE( std::get<2>(effTemp[i]) == Approx( eftemp[2] ).epsilon(1e-6) );
          REQUIRE( std::get<3>(effTemp[i]) == Approx( eftemp[3] ).epsilon(1e-6) );
          eftemp = { 0, 0, 0, 0 };
        } // end for
      } // THEN
    } // WHEN

    WHEN( "material considered is BeO" ){
      mat = 1099;
      tempVals = { { 200.00, 220.00, 250.00, 260.00 }, 
        { 296.00, 300.00, 301.00, 405.00 }, { 300.00, 315.00, 320.00, 345.00 },
        { 400.00, 455.00, 520.00, 645.00 } };
      effTemp  = { { 200.00, 220.00, 250.00, 260.00 }, 
        { 427.80, 427.80, 427.80, 502.80 }, { 427.80, 315.00, 320.00, 345.00 },
        { 502.80, 455.00, 520.00, 645.00 } };

      THEN( "the effective temperature is correctly identified in table" ){
        for ( size_t i = 0; i < tempVals.size(); ++i ){
          temp = { std::get<0>(tempVals[i]), std::get<1>(tempVals[i]), 
                   std::get<2>(tempVals[i]), std::get<3>(tempVals[i]) };
          gatef2( temp, eftemp, mat );
          REQUIRE( std::get<0>(effTemp[i]) == Approx( eftemp[0] ).epsilon(1e-6) );
          REQUIRE( std::get<1>(effTemp[i]) == Approx( eftemp[1] ).epsilon(1e-6) );
          REQUIRE( std::get<2>(effTemp[i]) == Approx( eftemp[2] ).epsilon(1e-6) );
          REQUIRE( std::get<3>(effTemp[i]) == Approx( eftemp[3] ).epsilon(1e-6) );
          eftemp = { 0, 0, 0, 0 };
        } // end for
      } // THEN
    } // WHEN
  } // GIVEN
} // TEST CASE
