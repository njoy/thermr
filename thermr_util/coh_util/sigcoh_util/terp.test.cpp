#define CATCH_CONFIG_MAIN
#include "../../../catch.hpp"
#include "terp.h"


TEST_CASE( "do220" ){
  std::vector<double> x, y;

  GIVEN( "inputs" ){
    REQUIRE( true );


  } // GIVEN
} // TEST CASE


TEST_CASE( "do230" ){
  std::vector<double> x, y, out;
  std::vector<double> inputs (4);
  int il = 2;
  double arg;

  GIVEN( "x, y vectors of equal length" ){
    x = { 296, 400, 500, 600, 700, 800, 1000, 1200, 1600, 2000 },
    y = { 2.1997, 2.7448, 3.2912, 3.851, 4.421, 4.9969, 6.1624, 7.3387, 
          9.6287, 11.992 };
    WHEN( "arg value is a value in x" ){
      arg = 296;

      inputs = { 5, 8, 2, 4 };
      out = { 2.094364, 2.1633, 2.176544, 2.1182 };

      for ( size_t i = 0; i < inputs.size(); ++i ){
        REQUIRE( out[i] == Approx( do230( x, y, arg, il, inputs[i] ) ).epsilon( 1e-6 ) );
      }

      arg = 590;
      out = { 3.78751, 3.84645, 3.78296, 3.794 };

      for ( size_t i = 0; i < inputs.size(); ++i ){
        REQUIRE( out[i] == Approx( do230( x, y, arg, il, inputs[i] ) ).epsilon( 1e-6 ) );
      }

    } // WHEN
    WHEN( "arg value is not a value in x" ){
      AND_WHEN( "arg value is within range" ){
        arg = 350;

        inputs = { 5, 8, 2, 4 };
        out = { 2.40535,  2.47245, 2.4716, 2.426 };

        for ( size_t i = 0; i < inputs.size(); ++i ){
          REQUIRE( out[i] == Approx( do230( x, y, arg, il, inputs[i] ) ).epsilon( 1e-6 ) );
        }

        arg = 1500;
        out = { 9.0282, 9.0562, 8.7552, 8.981 };

        for ( size_t i = 0; i < inputs.size(); ++i ){
          REQUIRE( out[i] == Approx( do230( x, y, arg, il, inputs[i] ) ).epsilon( 1e-6 ) );
        }

      } // AND WHEN
      AND_WHEN( "arg value is not within range" ){
        arg = 200;

        inputs = { 5, 8, 2, 4 };
        out = { 1.5415, 1.6137, 1.652, 1.571 };

        for ( size_t i = 0; i < inputs.size(); ++i ){
          REQUIRE( out[i] == Approx( do230( x, y, arg, il, inputs[i] ) ).epsilon( 1e-6 ) );
        }

        arg = 3000;
        out = { 17.6667,  17.6437, 16.9512, 17.531 };

        for ( size_t i = 0; i < inputs.size(); ++i ){
          REQUIRE( out[i] == Approx( do230( x, y, arg, il, inputs[i] ) ).epsilon( 1e-6 ) );
        }

      } // AND WHEN

    } // WHEN
  } // GIVEN
} // TEST CASE



TEST_CASE( "do250" ){
  std::vector<double> x, y, out;
  std::vector<std::tuple<int,int,int>> inputs (4);
  int il = 2;
  double arg;

  GIVEN( "x, y vectors of equal length" ){
    x = { 296, 400, 500, 600, 700, 800, 1000, 1200, 1600, 2000 },
    y = { 2.1997, 2.7448, 3.2912, 3.851, 4.421, 4.9969, 6.1624, 7.3387, 
          9.6287, 11.992 };
    WHEN( "arg value is a value in x" ){
      arg = 296;

      inputs = { { 4, 5, 6 }, { 9, 3, 2 }, { 4, 3, 1 }, { -1, -5, 0 } };
      out = { 2.094364, 2.1633, 2.176544, 2.1182 };

      for ( size_t i = 0; i < inputs.size(); ++i ){
        REQUIRE( out[i] == Approx( do250( x, y, arg, il, 
          std::get<0>(inputs[i]), std::get<1>(inputs[i]), 
          std::get<2>(inputs[i]) ) ).epsilon( 1e-6 ) );
      }

      arg = 590;
      out = { 3.78751, 3.84645, 3.78296, 3.794 };

      for ( size_t i = 0; i < inputs.size(); ++i ){
        REQUIRE( out[i] == Approx( do250( x, y, arg, il, 
          std::get<0>(inputs[i]), std::get<1>(inputs[i]), 
          std::get<2>(inputs[i]) ) ).epsilon( 1e-6 ) );
      }

    } // WHEN
    WHEN( "arg value is not a value in x" ){
      AND_WHEN( "arg value is within range" ){
        arg = 350;

        inputs = { { 4, 5, 6 }, { 9, 3, 2 }, { 4, 3, 1 }, { -1, -5, 0 } };
        out = { 2.40535,  2.47245, 2.4716, 2.426 };

        for ( size_t i = 0; i < inputs.size(); ++i ){
          REQUIRE( out[i] == Approx( do250( x, y, arg, il, 
            std::get<0>(inputs[i]), std::get<1>(inputs[i]), 
            std::get<2>(inputs[i]) ) ).epsilon( 1e-6 ) );
        }

        arg = 1500;
        out = { 9.0282, 9.0562, 8.7552, 8.981 };

        for ( size_t i = 0; i < inputs.size(); ++i ){
          REQUIRE( out[i] == Approx( do250( x, y, arg, il, 
            std::get<0>(inputs[i]), std::get<1>(inputs[i]), 
            std::get<2>(inputs[i]) ) ).epsilon( 1e-6 ) );
        }

      } // AND WHEN
      AND_WHEN( "arg value is not within range" ){
        arg = 200;

        inputs = { { 4, 5, 6 }, { 9, 3, 2 }, { 4, 3, 1 }, { -1, -5, 0 } };
        out = { 1.5415, 1.6137, 1.652, 1.571 };

        for ( size_t i = 0; i < inputs.size(); ++i ){
          REQUIRE( out[i] == Approx( do250( x, y, arg, il, 
            std::get<0>(inputs[i]), std::get<1>(inputs[i]), 
            std::get<2>(inputs[i]) ) ).epsilon( 1e-6 ) );
        }

        arg = 3000;
        out = { 17.6667,  17.6437, 16.9512, 17.531 };

        for ( size_t i = 0; i < inputs.size(); ++i ){
          REQUIRE( out[i] == Approx( do250( x, y, arg, il, 
            std::get<0>(inputs[i]), std::get<1>(inputs[i]), 
            std::get<2>(inputs[i]) ) ).epsilon( 1e-6 ) );
        }

      } // AND WHEN

    } // WHEN
  } // GIVEN
} // TEST CASE



TEST_CASE( "do260" ){
  std::vector<double> x, y, arg, out;
  int l = 1, il = 2;

  GIVEN( "x, y vectors of equal length" ){
    x = { 296, 400, 500, 600, 700, 800, 1000, 1200, 1600, 2000 },
    y = { 2.1997, 2.7448, 3.2912, 3.851, 4.421, 4.9969, 6.1624, 7.3387, 
          9.6287, 11.992 };

    WHEN( "arg inputs that are not in x" ){
      AND_WHEN( "inputs are in range" ){
        arg = { 312, 385, 423, 945, 1345, 1995 };
        out = { 2.283561, 2.666179, 2.865350, 5.601333, 7.697872, 11.10474 };
        THEN( "outputs" ){
          for ( size_t i = 0; i < arg.size(); ++i ){
            REQUIRE( out[i] == 
                     Approx( do260( x, y, arg[i], l, il ) ).epsilon(1e-6) );
          }
        } // THEN
      } // AND WHEN

      AND_WHEN( "inputs are not in range" ){
        arg = { 0, 5, 50, 185, 204, 2001, 2450, 3000, 10000 };
        out = { 0.6482615, 0.6744682, 0.9103288, 1.617910, 1.717496, 11.13619, 
                13.48955, 16.3723, 53.06172 };
        THEN( "outputs" ){
          for ( size_t i = 0; i < arg.size(); ++i ){
            REQUIRE( out[i] == 
                     Approx( do260( x, y, arg[i], l, il ) ).epsilon(1e-6) );
          }
        } // THEN
      } // AND WHEN
    } //  WHEN 

    WHEN( "inputs that are in x" ){
      arg = x;
      out = { 2.1997, 2.7448, 3.268934, 3.793069, 4.317203, 4.841338, 5.889607, 
              6.937876, 9.034415, 11.13095 };
      THEN( "outputs" ){
        for ( size_t i = 0; i < arg.size(); ++i ){
          REQUIRE( out[i] ==
                   Approx( do260( x, y, arg[i], l, il ) ).epsilon(1e-6) );
        }
      } // THEN
    } // GIVEN
  } // GIVEN
} // TEST CASE


TEST_CASE( "terp" ){
  std::vector<double> x, y, out;
  int nl = 10, il1 = 2;
  double arg;


  GIVEN( "" ){
    x = { 296, 400, 500, 600, 700, 800, 1000, 1200, 1600, 2000 },
    y = { 2.1997, 2.7448, 3.2912, 3.851, 4.421, 4.9969, 6.1624, 7.3387, 
          9.6287, 11.992 };

    WHEN( "120 260 300" ){ 
      arg = 50;
      REQUIRE( 0.910328 == Approx( terp( x, y, nl, arg, il1 ) ).epsilon(1e-6) );

      arg = 100;
      REQUIRE( 1.17239615 == Approx( terp( x, y, nl, arg, il1 ) ).epsilon(1e-6) );

      arg = 200;
      REQUIRE( 1.69653076 == Approx( terp( x, y, nl, arg, il1 ) ).epsilon(1e-6) );
    
      arg = 300;
      REQUIRE( 2.2206653 == Approx( terp( x, y, nl, arg, il1 ) ).epsilon(1e-6) );



    } // WHEN

    WHEN( "120 170 200 220 240 300" ){ 

      arg = 500;
      REQUIRE( 3.2912 == Approx( terp( x, y, nl, arg, il1 ) ).epsilon(1e-6) );

      arg = 700;
      REQUIRE( 4.421 == Approx( terp( x, y, nl, arg, il1 ) ).epsilon(1e-6) );

      arg = 800;
      REQUIRE( 4.9969 == Approx( terp( x, y, nl, arg, il1 ) ).epsilon(1e-6) );

      arg = 1000;
      REQUIRE( 6.1624 == Approx( terp( x, y, nl, arg, il1 ) ).epsilon(1e-6) );


    } // WHEN

    WHEN( "120 170 200 220 230 260 300" ){

      arg = 1300;
      REQUIRE( 7.9112 == Approx( terp( x, y, nl, arg, il1 ) ).epsilon(1e-6) );

      arg = 1500;
      REQUIRE( 9.0562 == Approx( terp( x, y, nl, arg, il1 ) ).epsilon(1e-6) );

    } // WHEN

    WHEN( "120 170 260 300" ){

      arg = 1800;
      REQUIRE( 10.81035 == Approx( terp( x, y, nl, arg, il1 ) ).epsilon(1e-6) );

      arg = 2800;
      REQUIRE( 16.7186 == Approx( terp( x, y, nl, arg, il1 ) ).epsilon(1e-6) );

    } // WHEN

    WHEN( "120 160 300" ){

      arg = 400;
      REQUIRE( 2.7448 == Approx( terp( x, y, nl, arg, il1 ) ).epsilon(1e-6) );


    } // WHEN

  } // GIVEN

} // TEST CASE
