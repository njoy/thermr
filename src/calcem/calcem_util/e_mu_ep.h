#include "iel/iel_util/terp1.h"
#include "coh/coh_util/sigcoh_util/legndr.h"
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


auto adaptiveReconstruction( const double& teff, const double& cliq, 
  const int& iinc, const double& tevz, const int& lat, const int& lasym, 
  std::vector<double>& yy, std::vector<double>& yu, const double& sb, 
  const double& sb2, std::vector<double>& x, const std::vector<double>& alpha, 
  const std::vector<double>& beta, const std::vector<std::vector<double>>& sab, 
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
            tevz, iinc, lat, lasym, cliq, sb, sb2, teff );
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


auto e_mu_ep( int matdp, int mtref, double t, double& teff, 
  double& teff2, std::vector<double>& scr, double za, double awr, int ncds,
  int nw, int nne, double cliq, int iinc, double emax, std::vector<double> egrid,
  double temp, double breakVal, std::vector<double>& esi, double tevz, int lat,
  int lasym, std::vector<double>& yy, std::vector<double>& yu, double sb,
  double sb2, std::vector<double>& x, std::vector<double> alpha, 
  std::vector<double> beta, const std::vector<std::vector<double>>&sab, 
  double az, std::vector<double>& uj, std::vector<double>& sj,
  double tol, double tolmin, double mumax, int imax, std::vector<double>& ubar,
  std::vector<double>& xsi, int nep, double yumin, int iend, int npage, int ib, int nb ){
  double bk = 8.6173845e-5;

  // compute kernel and write in mf6/law7 angle-energy format
  // 510 continue
  int ltt=4;
  int math=matdp;
  int mfh=6;
  int mth=mtref;
  double tev=t*bk;
  teff=teff*bk;
  teff2=teff2*bk;
  scr[1-1]=za;
  scr[2-1]=awr;
  scr[3-1]=0;
  scr[4-1]=ltt; // temporary flag for this format
  scr[5-1]=1;
  scr[6-1]=0;
  //call contio(0,0,nscr,scr,nb,nw);
  ncds=ncds+1;
  scr[1-1]=1;
  scr[2-1]=1;
  scr[3-1]=0;
  scr[4-1]=7;
  scr[5-1]=1;
  scr[6-1]=2;
  scr[7-1]=2;
  scr[8-1]=2;
  scr[9-1]=1.e-5;
  scr[10-1]=1;
  scr[11-1]=emax;
  scr[12-1]=1;
  nw=12;
  //call tab1io(0,0,nscr,scr,nb,nw)
  ncds=ncds+2;
  scr[1-1]=0;
  scr[2-1]=0;
  scr[3-1]=0;
  scr[4-1]=0;
  scr[5-1]=1;
  scr[6-1]=nne;
  scr[7-1]=nne;
  scr[8-1]=2;
  nw=8;
  //call tab2io(0,0,nscr,scr,nb,nw)
  ncds=ncds+2;
  cliq=0;
  if (iinc == 1) std::cout << "go to 515" << std::endl;
  //if (sab(1,1) > sab(2,1)) cliq=(sab(1,1)-sab(1,2))*alpha[1-1]/(beta[2-1]*beta[2-1]);

  // loop over given incident energy grid.
  // the first pass computes the cross section and mu grid
 
  //515 continue
  int ie=0;



  //520 continue
  while (true){
    std::cout << 520 << std::endl;
    ie=ie+1;
    double enow=egrid[ie-1];
    if (ie > 1 and temp > breakVal) enow=enow*temp/breakVal;
    //enow=sigfig(enow,8,0)
    esi[ie-1]=enow;
    int j=0;
    double sum=0;
 
    // adaptive reconstruction of angular cross section
    int u = -1;
    x[1] = u;
    sigu( int(yu.size()), enow, u, tev, alpha, beta, sab, yu, tol, az, tevz, 
      iinc, lat, lasym, cliq, sb, sb2, teff );
    yy[1] = yu[0];
    u = 1;
    x[0] = u;
    //call sigu(enow,u,tev,nalpha,alpha,nbeta,beta,sab,yu,nemax,tol)
    sigu( int(yu.size()), enow, u, tev, alpha, beta, sab, yu, tol, az, tevz, iinc, lat, 
        lasym, cliq, sb, sb2, teff );
    yy[0] = yu[0];
    int i = 2;

    // This does 530, 560, and 575
    auto out = adaptiveReconstruction( teff, cliq, iinc, tevz, lat, lasym, yy, 
      yu, sb, sb2, x, alpha, beta, sab, az, uj, sj, tol, tolmin, mumax, i, sum, 
      imax, enow, tev, j  );

    double xl = std::get<0>(out), yl = std::get<1>(out);

    std::cout << 580 << std::endl;
    // linearization complete.  write out result.
    // 580 continue
    ++j;
    uj[j-1]=x[0];
    sj[j-1]=yy[0];
    int nmu=j;
    ubar[ie-1]=0;
    sum=sum+0.5*(yy[1-1]+yl)*(x[1-1]-xl);
    xsi[ie-1]=sum/2;
    for ( int i = 1; i < nmu; ++i ){
       ubar[ie-1]=ubar[ie-1]+0.5*(uj[i]-uj[i-1])*(sj[i]+sj[i-1])*(uj[i]+uj[i-1]);
    }
    ubar[ie-1]=0.5*ubar[ie-1]/sum;
 
    return;

    // now loop through the mu grid to write out the distributions
    mth=mtref;
    scr[1-1]=0;
    scr[2-1]=enow;
    scr[3-1]=0;
    scr[4-1]=0;
    scr[5-1]=1;
    scr[6-1]=nmu;
    scr[7-1]=nmu;
    scr[8-1]=2;
    nw=8;
    //call tab2io(0,0,nscr,scr,nb,nw)
    ncds=ncds+1;

    for ( int il = 0; il < nmu; ++il ){
    //do il=1,nmu
       u = uj[il];
       sigu( int(yu.size()), enow, u, tev, alpha, beta, sab, yu, tol, az, tevz, 
         iinc, lat, lasym, cliq, sb, sb2, teff );

       nep = int(yu[1]);
       j=0;
       for ( int i = 1; i <= nep; ++i ){
       //do i=1,nep
         j = nep - i;
         if (yu[2*(nep-i)+4-1]/sum > yumin){ 
           std::cout << " exit " << std::endl; 
           throw std::exception(); 
         }
       }
       nep = j;

       scr[1-1]=0;
       scr[2-1]=u;
       scr[3-1]=0;
       scr[4-1]=0;
       scr[5-1]=1;
       scr[6-1]=nep;
       scr[7-1]=nep;
       scr[8-1]=2;
       int k=8;
       int istart=1;
       while (true) {
         //595 continue
         std::cout << 595 << std::endl;
         iend=nep;
         if ((iend-istart) >= npage/2) iend=istart+npage/2-1;
         j=k-1;
         ib=istart-1;

        //596 continue
         while (true){
         j=j+2;
         ib=ib+1;
         scr[j-1]=yu[1+2*ib-1];
         scr[j+1-1]=yu[2+2*ib-1]*2/sum;
         if (ib < iend) { continue; }// go to 596
         else { break; }
         }
         nw=j+1;
         if (k == 0) {
           std::cout << "go to 597" << std::endl;
           //597 continue
             //call moreio(0,0,nscr,scr,nb,nw)
             if (nb == 0) {
               std::cout << "go to 598" << std::endl;
               ncds=ncds+1+(j*(nep+1)+5)/6;
               continue;
             }
             istart=iend+1;
             std::cout << "go to 595" << std::endl;
             continue;
         }
         k=0;
         //call tab1io(0,0,nscr,scr,nb,nw)
         if (nb == 0) {
           std::cout << "go to 598" << std::endl;
           continue;
         } 
         istart=iend+1;
         std::cout << "go to 595" << std::endl;
         continue;

        //598 continue
        // ncds=ncds+1+(j*(nep+1)+5)/6;
       }
       /*
       */
    }
 
    //if (ie < nne) go to 520
  
       /*

    */
    if ( ie >= nne ){ 
      std::cout << "go to 610" << std::endl;
      return; 
    }
  } // while we want to do 520 
} 
