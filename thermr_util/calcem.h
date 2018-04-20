
auto calcem( double temp, int itemp, std::fstream& iold, std::fstream& inew,
    int ne, int next ){
  /*-------------------------------------------------------------------
   * Calculate incoherent inelastic scattering kernels from
   * s(alpha,beta) in endf mf7 format or from analytic models
   * of s(alpha,beta). The incident energy grid is fixed (see egrid).
   * For iform=0, the secondary energy and scattering cosine grids
   * are determined adaptively for linear interpolations.  A compact
   * angle representation using equally probable cosines is generated.
   * For iform=1, a fixed grid of mu values is used, and the secondary
   * distribution is determined adaptively for linear interpolation
   * for each mu.  In both cases, the results are written to a
   * scratch tape for tpend.
   *-------------------------------------------------------------------
   */
  


  //  character(len=60)::strng
   int nr,np,nwtab,nl,nlt,nlp,nlp1,nnl,jmax,nne;
   int i,ia,nl1,ltt,loc,l,jscr,ilog;
   int matd,itprnt,nb,nw,ni,nbeta,lt,it,nalpha;
   int itrunc,ib,ip,ir,idis,ie,nmu,nep,istart,iend;
   int jbeta,j,iskip,il,k,jnz,ll,jj,ii,nexn,ien,isave;
   int nll;

   double smz,t1,t,tmax,test,tempt,tt1,ttn,tnxt,xm;
   double uu,uum,ym,test2,xlast,ulast,xs;
   double b,diff,enow,ep,sabmin,tev,ylast;
   double u,xl,yl,sum;
   double tone,elo;
   const int ngrid = 118, nlmax = 65, nemax = 5000, mumax = 300, imax = 20;

   std::vector<double> ex(imax), x(imax), yt(nlmax), yy(imax), yu(2*nemax),
     ubar(ngrid);
   std::vector<std::vector<double>> y( nlmax, std::vector<double>(imax) );
   double u2,u2last,u3,u3last;
   std::vector<double> uj(mumax),sj(mumax),p2(ngrid),p3(ngrid),p(4);
   std::vector<double> alpha, beta;
   std::vector<std::vector<double>> sab;

   std::vector<double> egrid { 1.e-5, 1.78e-5, 2.5e-5, 3.5e-5, 5.0e-5, 7.0e-5,
      1.e-4, 1.26e-4, 1.6e-4, 2.0e-4, 0.000253, 0.000297, 0.000350, 0.00042, 
      0.000506, 0.000615, 0.00075, 0.00087, 0.001012, 0.00123, 0.0015, 
      0.0018, 0.00203, 0.002277, 0.0026, 0.003, 0.0035, 0.004048, 0.0045, 
      0.005, 0.0056, 0.006325, 0.0072, 0.0081, 0.009108, 0.01, 0.01063, 
      0.0115, 0.012397, 0.0133, 0.01417, 0.015, 0.016192, 0.0182, 0.0199, 
      0.020493, 0.0215, 0.0228, 0.0253, 0.028, 0.030613, 0.0338, 0.0365, 
      0.0395, 0.042757, 0.0465, 0.050, 0.056925, 0.0625, 0.069, 0.075, 
      0.081972, 0.09, 0.096, 0.1035, 0.111573, 0.120, 0.128, 0.1355, 0.145728, 
      0.160, 0.172, 0.184437, 0.20, 0.2277, 0.2510392, 0.2705304, 0.2907501, 
      0.3011332, 0.3206421, 0.3576813, 0.39, 0.4170351, 0.45, 0.5032575, 0.56, 
      0.625, 0.70, 0.78, 0.86, 0.95, 1.05, 1.16, 1.28, 1.42, 1.55, 1.70, 1.855,
      2.02, 2.18, 2.36, 2.59, 2.855, 3.12, 3.42, 3.75, 4.07, 4.46, 4.90, 5.35,
      5.85, 6.40, 7.00, 7.65, 8.40, 9.15, 9.85, 10.00 };

   double unity = 1.0, sabflg = -225.0, eps = 1.0e-4, tolmin = 5.0e-7, half = 0.5,
     therm = 0.0253, breakVal = 3000.0, em9 = 1.0e-9, zero = 0.0, tenth = 0.1,
     up = 1.1, dn = 0.9, uumin = 0.00001, yumin = 2.0e-7;

   int nlpmx = 10;

   // save nwtab,sabmin,nl,nlt,nlp,nlp1,nl1,nnl,jmax,nne
   
   double tevz=therm;
  }

