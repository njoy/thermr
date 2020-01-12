#define CATCH_CONFIG_MAIN
#include "catch.hpp"
#include "calcem/calcem_util/e_ep_mu_util/sigl.h"
#include <range/v3/all.hpp>
#include "generalTools/testing.h"





/*


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

  addMidpointsRight(i,x,y,e,ep,tev,alpha,beta,sab,az,lasym,lat,sigma_b,sigma_b2,teff,iinc,tol);
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


*/



TEST_CASE( "Get pdf value" ){
  double e = 1e-2, ep = 1e-3, tev = 0.025, tolin = 5e-2*0.5, az = 0.99917, 
         sigma_b = 4.0, sigma_b2 = 0.0, teff = 0.12;
  int nL = -9, lat = 1, iinc = 2, lasym = 0, nbin = 8;
  bool equiprobableBins = true;
  std::vector<double> alphas { 1.1, 2.2, 3.3, 4.5, 5.8 },
                       betas { 0.1, 0.2, 1.3, 1.4, 2.5, 2.6, 3.7 };
  std::vector<double> sab(alphas.size()*betas.size());
  for ( size_t i = 0; i < alphas.size(); ++i ){
    for ( size_t j = 0; j < betas.size(); ++j ){
      sab[i*betas.size()+j] = 0.01*((j+1) + 0.1*(i+1));
    } 
  } 

  double pdfVal, correctPDF;


  e = 1.0;

  ep = 2.0;
  pdfVal = getPDF(ep,e,tev,tolin,lat,iinc,alphas,betas,sab,az,lasym,sigma_b,sigma_b2,teff);
  correctPDF = 0.0;
  REQUIRE( pdfVal == Approx(correctPDF).epsilon(1e-6) );



  //ep = 1.0;
  //correctPDF = 3.0989839243186652;
  //REQUIRE( pdfVal == Approx(correctPDF).epsilon(1e-6) );








  e = 1e-2;

  std::vector<double> correctPDFVals { 71.84429, 220.1616, 61.17360, 19.66197, 6.227846 };
  std::vector<double> finalEnergies  { 1e-1, 1e-2, 1e-3, 1e-4, 1e-5 };
  for ( size_t i = 0; i < finalEnergies.size(); ++i ){
    ep = finalEnergies[i];
    pdfVal = getPDF(ep,e,tev,tolin,lat,iinc,alphas,betas,sab,az,lasym,sigma_b,sigma_b2,teff);
    REQUIRE( pdfVal == Approx(correctPDFVals[i]).epsilon(1e-6) );
  }


  e  = 1e-1; 
  correctPDFVals = { 262.937411, 1.02706201, 0.32551817, 0.10296110 };
  finalEnergies = { 1e-2, 1e-3, 1e-4, 1e-5 };
  for ( size_t i = 0; i < finalEnergies.size(); ++i ){
    ep = finalEnergies[i];
    pdfVal = getPDF(ep,e,tev,tolin,lat,iinc,alphas,betas,sab,az,lasym,sigma_b,sigma_b2,teff);
    REQUIRE( pdfVal == Approx(correctPDFVals[i]).epsilon(1e-6) );
  }





} // TEST CASE
<<<<<<< HEAD
=======






/*
>>>>>>> ba824f1f8c460d0589cef6155b0f0a6c2ef7e894
*/



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
  ep = 1e-1;
  s = sigl(ep,e,tev,tolin,lat,iinc,alphas,betas,sab,az,lasym,sigma_b,sigma_b2,teff,nbin,equiprobableBins);


<<<<<<< HEAD
  /*
  ep = 1e-1;
  s = sigl(ep,e,tev,tolin,lat,iinc,alphas,betas,sab,az,lasym,sigma_b,sigma_b2,teff,nbin,equiprobableBins);
  correct_s = {-0.48647654336295154, -0.28577865744470676, -8.8254959343803965E-2, 0.10934278668830079, 0.30704325144135147, 0.50477104242836623, 0.70270321344353548, 0.90176859041415269};
  REQUIRE(ranges::equal(correct_s, s, equal));
=======
>>>>>>> ba824f1f8c460d0589cef6155b0f0a6c2ef7e894

  WHEN ( "2 bins requested" ){
      nbin = 2;

      e = 1.0;

      // E = 1.0 eV
      ep = 2.0;
      s = sigl(ep,e,tev,tolin,lat,iinc,alphas,betas,sab,az,lasym,sigma_b,
               sigma_b2,teff,nbin,equiprobableBins);
      correct_s = { 0.0, 0.0 };
      REQUIRE(ranges::equal(correct_s, s, equal));


      //ep = 1.0;
      //s = sigl(ep,e,tev,tolin,lat,iinc,alphas,betas,sab,az,lasym,sigma_b,
      //         sigma_b2,teff,nbin,equiprobableBins);
      //correct_s = { 0.91492184570190904, 0.99902962413315932 };
      //REQUIRE(ranges::equal(correct_s, s, equal));



      
      
      
      e = 1e-2;
      // E = 1e-2 eV

      std::vector<std::vector<double>> correctMuVecs {
        {  0.0000000, 0.0000000 },
        { -0.1877993, 0.6040790 },
        {  0.2550645, 0.8945141 },
        { -0.5002791, 0.4996727 },
        { -0.5000937, 0.4999014 },
        { -0.5000301, 0.4999693 } };
      std::vector<double> finalEnergies { 1e0, 1e-1, 1e-2, 1e-3, 1e-4, 1e-5 };
      for ( size_t i = 0; i < correctMuVecs.size(); ++i ){
        ep = finalEnergies[i];
        s = sigl(ep,e,tev,tolin,lat,iinc,alphas,betas,sab,az,lasym,sigma_b,
                 sigma_b2,teff,nbin,equiprobableBins);
        REQUIRE(ranges::equal(correctMuVecs[i], s, equal));
      }





  } // WHEN



  WHEN( "8 bins are requested" ){ 

      nbin = 8;
  ep = 1e-1;
      std::vector<std::vector<double>> correctMuVecs {
{-0.4864765, -0.285778, -0.0882549, 0.10934278, 0.30704325, 0.50477104, 0.70270321, 0.90176859 },
{-0.1447954,  0.165989,  0.4005134, 0.598550924, 0.75607791, 0.87463820, 0.95401954, 0.99332099 },
{-0.8750342, -0.625260, -0.3753526,-0.125469544, 0.12447684, 0.37451189, 0.62447532, 0.87522679 },
{-0.8750080, -0.625091, -0.3751193,-0.125156114, 0.12483210, 0.37484817, 0.62485204, 0.87507351 },
{-0.8750021, -0.625029, -0.3750385,-0.125050067, 0.12494679, 0.37495231, 0.62495468, 0.87502358 } };

      std::vector<double> finalEnergies { 1e-1, 1e-2, 1e-3, 1e-4, 1e-5 };
      for ( size_t i = 0; i < correctMuVecs.size(); ++i ){
        ep = finalEnergies[i];
        s = sigl(ep,e,tev,tolin,lat,iinc,alphas,betas,sab,az,lasym,sigma_b,
                 sigma_b2,teff,nbin,equiprobableBins);
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
  
      std::vector<double> correct_s, s;
      ep = 1e-1;
      s = sigl(ep,e,tev,tolin,lat,iinc,alphas,betas,sab,az,lasym,sigma_b,sigma_b2,teff,nbin,equiprobableBins);
     correct_s = { 0.20813984053303053, -0.12610223534110626, 3.1542005728193809E-2, 9.9854434619668919E-3, -5.0951456351513429E-2, 3.5256512100910638E-3, -2.3306274801414548E-2, -4.7364309932466368E-2 };
      REQUIRE(ranges::equal(correct_s, s, equal));


      ep = 1e-2;
      s = sigl(ep,e,tev,tolin,lat,iinc,alphas,betas,sab,az,lasym,sigma_b,sigma_b2,teff,nbin,equiprobableBins);
      correct_s = {0.57478936634776701, 0.21262714685127976, 0.14716901781442446, 0.16253311095511172, 0.10777361840677388, 7.4051505015408664E-2, 8.8236720971749408E-2, 7.0371173451007532E-2 };
      REQUIRE(ranges::equal(correct_s, s, equal));

      ep = 1e-3;
      correct_s = {-3.0320534205725835E-4, -7.8101962585028079E-3, 1.8046088909718511E-4, -2.5322851957153197E-2, 2.7573825028156737E-4, -4.8626772815391177E-2, 1.8907683649725671E-5, -6.3091000845413761E-2};
      s = sigl(ep,e,tev,tolin,lat,iinc,alphas,betas,sab,az,lasym,sigma_b,sigma_b2,teff,nbin,equiprobableBins);
      REQUIRE(ranges::equal(correct_s, s, equal));


      ep = 1e-4;
      s = sigl(ep,e,tev,tolin,lat,iinc,alphas,betas,sab,az,lasym,sigma_b,sigma_b2,teff,nbin,equiprobableBins);
      correct_s = {-9.6131247445313539E-5, -7.8041006589044393E-3, 6.3081720496883931E-5, -2.5478679834518559E-2, 8.6190495380771004E-5, -4.8743549347926809E-2, 5.1933656887576118E-6, -6.2982367634516187E-2 };
      REQUIRE(ranges::equal(correct_s, s, equal));


      ep = 1e-5;
      s = sigl(ep,e,tev,tolin,lat,iinc,alphas,betas,sab,az,lasym,sigma_b,sigma_b2,teff,nbin,equiprobableBins);
      correct_s = {-3.0413003964419039E-5, -7.8090977480863855E-3, 2.0593541745807242E-5, -2.5523322397632636E-2, 2.7494078227038865E-5, -4.8774081688794373E-2, 1.5976592289690394E-6, -6.2953246075833497E-2 };
      REQUIRE(ranges::equal(correct_s, s, equal));

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



<<<<<<< HEAD
  */
    }
=======
  } // GIVEN 
>>>>>>> ba824f1f8c460d0589cef6155b0f0a6c2ef7e894

}

/*
*/

/*

*/
