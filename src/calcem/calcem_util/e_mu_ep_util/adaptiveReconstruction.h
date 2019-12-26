#include "general_util/sigfig.h"
#include "calcem/calcem_util/e_mu_ep_util/sigu.h"

template <typename A, typename F>
auto do575(int& i, A& x, A& yy, const A& yu, const F& xm ){
  // test fails. add to stack and continue
  // 575 continue
  std::cout << 575 << std::endl;
  ++i;
  x[i-1] = x[i-2];
  x[i-2] = xm;
  yy[i-1]= yy[i-2];
  yy[i-2]= yu[0];

}


auto adaptiveReconstruction( const double& teff, 
  const int& iinc, const double& tevz, const int& lat, const int& lasym, 
  std::vector<double>& yy, std::vector<double>& yu, const double& sb, 
  const double& sb2, std::vector<double>& x, const std::vector<double>& alpha, 
  const std::vector<double>& beta, const std::vector<double>& sab, 
  const double& az, std::vector<double>& uj, std::vector<double>& sj, 
  const double& tol, const double& tolmin, const double& mumax, int& i, 
  double& sum, const int& imax, const double& enow, 
  const double& tev, int& j ){
    double xl = x[1];
    double yl = yy[1];

    // adaptive reconstruction
    do { 
      // 530 continue
      std::cout << 530 << std::endl;
      if (i != imax) {
        double xm = sigfig(0.5*(x[i-2]+x[i-1]),7,0);
        if (xm > x[i-1] and xm < x[i-2]) {
          sigu( int(yu.size()), enow, xm, tev, alpha, beta, sab, yu, tol, az, 
            tevz, iinc, lat, lasym, sb, sb2, teff );
          double ym = yy[i-1]+(xm-x[i-1])*(yy[i-2]-yy[i-1])/(x[i-2]-x[i-1]);
          if ( (x[i-2]-x[i-1]) > 0.25 or (std::abs(yu[0]-ym) > 2*tol*ym+tolmin) ){ 
            do575(i, x, yy, yu, xm );
            continue;
          }
        }
      } 
      // point passes.  save top point in stack and continue.
      std::cout << 560 << std::endl;
      ++j;
      if (j > mumax-1) { 
        throw std::exception(); // error('calcem','too many angles','see mumax')
      }
      uj[j-1] = x[i-1];
      sj[j-1] = yy[i-1];
      if (j > 1) {
         sum += 0.5 * (yy[i-1]+yl) * (x[i-1]-xl);
         xl = x[i-1];
         yl = yy[i-1];
      }
      --i;
    }  while ( i >= 2 );
    return std::make_tuple(xl,yl);
}


