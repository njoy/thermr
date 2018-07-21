#include "iel/iel_util/terp1.h"
#include "coh/coh_util/sigcoh_util/legndr.h"
#include "general_util/sigfig.h"
#include "calcem/calcem_util/e_mu_ep_util/sigu.h"
#include "calcem/calcem_util/e_mu_ep_util/adaptiveReconstruction.h"
#include "calcem/calcem_util/e_mu_ep_util/mainLoop.h"


auto e_mu_ep( int matdp, int mtref, double t, double& teff, 
  double& teff2, std::vector<double>& scr, double za, double awr, int ncds,
  int nw, int nne, double cliq, int iinc, double emax, std::vector<double> egrid,
  double temp, double breakVal, std::vector<double>& esi, double tevz, int lat,
  int lasym, std::vector<double>& yy, std::vector<double>& yu, double sb,
  double sb2, std::vector<double>& x, std::vector<double> alpha, 
  std::vector<double> beta, const std::vector<std::vector<double>>&sab, 
  double az, std::vector<double>& uj, std::vector<double>& sj,
  double tol, double tolmin, double mumax, int imax, std::vector<double>& ubar,
  std::vector<double>& xsi, int nep, double yumin, int npage, int ib, int nb ){
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


       mainLoop(nep, npage, j, k, ib, scr, yu, sum, nw, ncds, nb );




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
