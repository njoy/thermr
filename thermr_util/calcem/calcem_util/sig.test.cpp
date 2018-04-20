#define CATCH_CONFIG_MAIN
#include "../../../catch.hpp"
#include "sig.h"


TEST_CASE( "sig" ){
  GIVEN( "inputs" ){
    WHEN( "invalid option selected (iinc != 1,2)" ){
      double e = 1.0e-5, ep = 0.0, u = -1.0, tev = 2.5e-2, bbm = 0.0, az = 11.9,
        tevz = 2.53e-2, az2 = 0.0, teff2 = 0.0, cliq = 0.0, sb = 5.53, sb2 = 0.0,
        teff = 6.14e-2;
      int nalpha = 5, nbeta = 7, lasym = 0, lat = 1, iinc = 0;
      std::vector<double> alpha { 1.1, 2.2, 3.3, 4.5, 5.8 },
        beta { 0.1, 0.2, 1.3, 1.4, 2.5, 2.6, 3.7 };
      std::vector<std::vector<double>> sab(nalpha, std::vector<double>(nbeta,0));
      REQUIRE_THROWS( sig( e, ep, u, tev, nalpha, alpha, nbeta, beta, sab, 
        bbm, az, tevz, lasym, az2, teff2, lat, cliq, sb, sb2, teff, iinc ) );
    } // WHEN
    WHEN( "free gas option selected (iinc = 2)" ){
      AND_WHEN("ep = 0, so sigc = 0"){
        double e = 1.0e-5, ep = 0.0, u = -1.0, tev = 2.5e-2, bbm = 0.0, az = 11.9,
          tevz = 2.53e-2, az2 = 0.0, teff2 = 0.0, cliq = 0.0, sb = 5.53, sb2 = 0.0,
          teff = 6.14e-2;
        int nalpha = 5, nbeta = 7, lasym = 0, lat = 1, iinc = 1;
        std::vector<double> alpha { 1.1, 2.2, 3.3, 4.5, 5.8 },
          beta { 0.1, 0.2, 1.3, 1.4, 2.5, 2.6, 3.7 };
        std::vector<std::vector<double>> sab(nalpha, std::vector<double>(nbeta,0));
        auto sigVal = sig( e, ep, u, tev, nalpha, alpha, nbeta, beta, sab, 
          bbm, az, tevz, lasym, az2, teff2, lat, cliq, sb, sb2, teff, iinc );
        REQUIRE( 0.0 == Approx( sigVal ).epsilon(1e-6) );
      } // AND WHEN
      AND_WHEN("legitimate values are provided, so real sigVal is output" ){
        double e = 1.0e-5, ep = 1.2e-6, u = -1.0, tev = 2.5e-2, bbm = 0.0, az = 11.9,
          tevz = 2.53e-2, az2 = 0.0, teff2 = 0.0, cliq = 0.0, sb = 5.53, sb2 = 0.0,
          teff = 6.14e-2;
        int nalpha = 5, nbeta = 7, lasym = 0, lat = 1, iinc = 1;
        std::vector<double> alpha { 1.1, 2.2, 3.3, 4.5, 5.8 },
          beta { 0.1, 0.2, 1.3, 1.4, 2.5, 2.6, 3.7 };
        std::vector<std::vector<double>> sab(nalpha, std::vector<double>(nbeta,0));
        auto sigVal = sig( e, ep, u, tev, nalpha, alpha, nbeta, beta, sab, 
          bbm, az, tevz, lasym, az2, teff2, lat, cliq, sb, sb2, teff, iinc );
        REQUIRE( 1384.063418 == Approx( sigVal ).epsilon(1e-6) );
      } // AND WHEN
      AND_WHEN( "extremely small tev, so -arg < -225" ){
        double e = 1.0e-5, ep = 1.2e-6, u = -1.0, tev = 1.5e-8, bbm = 0.0, az = 11.9,
          tevz = 2.53e-2, az2 = 0.0, teff2 = 0.0, cliq = 0.0, sb = 5.53, sb2 = 0.0,
          teff = 6.14e-2;
        int nalpha = 5, nbeta = 7, lasym = 0, lat = 1, iinc = 1;
        std::vector<double> alpha { 1.1, 2.2, 3.3, 4.5, 5.8 },
          beta { 0.1, 0.2, 1.3, 1.4, 2.5, 2.6, 3.7 };
        std::vector<std::vector<double>> sab(nalpha, std::vector<double>(nbeta,0));
        auto sigVal = sig( e, ep, u, tev, nalpha, alpha, nbeta, beta, sab, 
          bbm, az, tevz, lasym, az2, teff2, lat, cliq, sb, sb2, teff, iinc );
        REQUIRE( 0.0 == Approx( sigVal ).epsilon(1e-6) );
      } // AND WHEN
    } // WHEN
    WHEN( "we have to use sct approximation (170)" ){
      AND_WHEN( "this is because of a > alpha(nalpha)" ){
        { 
        double e = 1.0e-2, ep = 1.2e-2, u = -1.0, tev = 1.5e-1, bbm = 0.0, az = 11.9,
          tevz = 2.2e-4, az2 = 0.0, teff2 = 0.0, cliq = 0.0, sb = 5.53, sb2 = 0.0,
          teff = 6.14e-2;
        int nalpha = 5, nbeta = 7, lasym = 0, lat = 1, iinc = 2;
        std::vector<double> alpha { 1.1, 2.2, 3.3, 4.5, 5.8 },
          beta { 0.1, 0.2, 1.3, 1.4, 2.5, 2.6, 3.7 };
        std::vector<std::vector<double>> sab(nalpha, std::vector<double>(nbeta,0));
        auto sigVal = sig( e, ep, u, tev, nalpha, alpha, nbeta, beta, sab, 
          bbm, az, tevz, lasym, az2, teff2, lat, cliq, sb, sb2, teff, iinc );
        REQUIRE( 55.838635 == Approx( sigVal ).epsilon(1e-6) );
      }
      { 
        double e = 1.0e-3, ep = 1.2e-3, u = -1.0, tev = 1.5e-1, bbm = 0.0, az = 11.9,
          tevz = 2.2e-4, az2 = 0.0, teff2 = 0.0, cliq = 0.0, sb = 5.53, sb2 = 0.0,
          teff = 6.14e-2;
        int nalpha = 5, nbeta = 7, lasym = 0, lat = 1, iinc = 2;
        std::vector<double> alpha { 0.1, 0.2, 0.3, 0.5, 0.8 },
          beta { 0.1, 0.2, 1.3, 1.4, 2.5, 2.6, 3.7 };
        std::vector<std::vector<double>> sab(nalpha, std::vector<double>(nbeta,0));
        auto sigVal = sig( e, ep, u, tev, nalpha, alpha, nbeta, beta, sab, 
          bbm, az, tevz, lasym, az2, teff2, lat, cliq, sb, sb2, teff, iinc );
        REQUIRE( 179.2164755 == Approx( sigVal ).epsilon(1e-6) );
      }

      } // AND WHEN
      AND_WHEN( "this is because of being out of range with beta" ){
        AND_WHEN( "bbm < beta[0]" ){ 
          double e = 1.0e-5, ep = 1.2e-6, u = -1.0, tev = 1.5e-8, bbm = 0.0, az = 11.9,
            tevz = 2.2e-4, az2 = 0.0, teff2 = 0.0, cliq = 0.0, sb = 5.53, sb2 = 0.0,
            teff = 6.14e-2;
          int nalpha = 5, nbeta = 7, lasym = 1, lat = 1, iinc = 2;
          std::vector<double> alpha { 1.1, 2.2, 3.3, 4.5, 5.8 },
            beta { 0.1, 0.2, 1.3, 1.4, 2.5, 2.6, 3.7 };
          std::vector<std::vector<double>> sab(nalpha, std::vector<double>(nbeta,0));
          auto sigVal = sig( e, ep, u, tev, nalpha, alpha, nbeta, beta, sab, 
            bbm, az, tevz, lasym, az2, teff2, lat, cliq, sb, sb2, teff, iinc );
          REQUIRE( 883.34692414 == Approx( sigVal ).epsilon(1e-6) );
        } // AND WHEN
        AND_WHEN( "bbm > beta[nbeta-1]" ){ 
          AND_WHEN( "-arg < sabflg" ){ 
            double e = 1.0e-5, ep = 1.2e-1, u = -1.0, tev = 1.5e-8, bbm = 0.0, az = 11.9,
              tevz = 2.2e-4, az2 = 0.0, teff2 = 0.0, cliq = 0.0, sb = 5.53, sb2 = 0.0,
              teff = 6.14e-2;
            int nalpha = 5, nbeta = 7, lasym = 1, lat = 1, iinc = 2;
            std::vector<double> alpha { 10.1, 20.2, 30.3, 40.5, 50.8 },
              beta { 0.1, 0.2, 1.3, 1.4, 2.5, 2.6, 3.7 };
            std::vector<std::vector<double>> sab(nalpha, std::vector<double>(nbeta,0));
            auto sigVal = sig( e, ep, u, tev, nalpha, alpha, nbeta, beta, sab, 
              bbm, az, tevz, lasym, az2, teff2, lat, cliq, sb, sb2, teff, iinc );
            REQUIRE( 0.0 == Approx( sigVal ).epsilon(1e-6) );
          } // AND WHEN
          AND_WHEN( "-arg > sabflg" ){ 
            {
            double e = 1.0e-5, ep = 1.2e-1, u = -1.0, tev = 1.5e-1, bbm = 0.0, az = 11.9,
              tevz = 2.2e-4, az2 = 0.0, teff2 = 0.0, cliq = 0.0, sb = 5.53, sb2 = 0.0,
              teff = 6.14e-2;
            int nalpha = 5, nbeta = 7, lasym = 1, lat = 1, iinc = 2;
            std::vector<double> alpha { 10.1, 20.2, 30.3, 40.5, 50.8 },
              beta { 0.1, 0.2, 1.3, 1.4, 2.5, 2.6, 3.7 };
            std::vector<std::vector<double>> sab(nalpha, std::vector<double>(nbeta,0));
            auto sigVal = sig( e, ep, u, tev, nalpha, alpha, nbeta, beta, sab, 
              bbm, az, tevz, lasym, az2, teff2, lat, cliq, sb, sb2, teff, iinc );
            REQUIRE( 12.9236152 == Approx( sigVal ).epsilon(1e-6) );
            }
            {
            double e = 1.0e-4, ep = 1.2e-1, u = -1.0, tev = 1.5e-1, bbm = 0.0, az = 11.9,
              tevz = 2.2e-4, az2 = 0.0, teff2 = 0.0, cliq = 0.0, sb = 5.53, sb2 = 0.0,
              teff = 6.14e-2;
            int nalpha = 5, nbeta = 7, lasym = 1, lat = 1, iinc = 2;
            std::vector<double> alpha { 10.1, 20.2, 30.3, 40.5, 50.8 },
              beta { 0.1, 0.2, 1.3, 1.4, 2.5, 2.6, 3.7 };
            std::vector<std::vector<double>> sab(nalpha, std::vector<double>(nbeta,0));
            auto sigVal = sig( e, ep, u, tev, nalpha, alpha, nbeta, beta, sab, 
              bbm, az, tevz, lasym, az2, teff2, lat, cliq, sb, sb2, teff, iinc );
            REQUIRE( 5.011738726 == Approx( sigVal ).epsilon(1e-6) );
            }

          } // AND WHEN
        } // AND WHEN
      } // AND WHEN
      AND_WHEN( "this is because of being out of range with beta" ){
        AND_WHEN( "lasym != 1" ){ 
          double e = 1.0e-4, ep = 5.2e-3, u = -1.0, tev = 1.5e-1, bbm = 0.0, az = 11.9,
            tevz = 2.2e-4, az2 = 0.0, teff2 = 0.0, cliq = 0.0, sb = 5.53, sb2 = 0.0,
            teff = 6.14e-2;
          int nalpha = 5, nbeta = 7, lasym = 0, lat = 1, iinc = 2;
          std::vector<double> alpha { 1.1, 2.2, 3.3, 4.5, 5.8 },
            beta { 0.1, 0.2, 1.3, 1.4, 2.5, 2.6, 3.7 };
          std::vector<std::vector<double>> sab(nalpha, std::vector<double>(nbeta,0));
          auto sigVal = sig( e, ep, u, tev, nalpha, alpha, nbeta, beta, sab, 
            bbm, az, tevz, lasym, az2, teff2, lat, cliq, sb, sb2, teff, iinc );
          REQUIRE( 795.18855771477308 == Approx( sigVal ).epsilon(1e-6) );
        } // AND WHEN
      } // AND WHEN
    } // WHEN
  } // GIVEN
} // TEST CASE

