

auto doSCTApproximation( int lat, double a, double b, double c,
  double tevz, double rtev, double tfff, double teff, double tfff2, 
  double teff2, double arg, double az, double az2, double sigc, double s2,
  double sigVal, double u, double sb2, double e, double tev, double sigmin,
  double sabflg, double bb, double s, double sb ){
  // short collision time for large beta or alpha
  std::cout << "170" << std::endl;
   if (lat == 1) b=b*tevz*rtev;
   if (lat == 1) a=a*tevz*rtev;
   tfff=teff;
   tfff2=teff2;
          // if (az < 3) {
          //    if (e > 10) {
          //       tfff=tev
          //       tfff2=tev
          //    else { if (e > 2) {
          //       rat=(e-2)/8
          //       tfff= (1-rat)*teff+rat*tev
          //       tfff2= (1-rat)*teff2+rat*tev
          //    } // end if
          // } // end if
   s=0;
   arg=(a-b)*(a-b)*tev/(4*a*tfff)+(b+bb)/2;
   if (-arg > sabflg) s=exp(-arg)/(c*sqrt(a*tfff*rtev));
   sigVal=sigc*sb*s;
   if (sb2 > 0.0) {
      a=a*az/az2;
      arg=(a-b)*(a-b)*tev/(4*a*tfff2)+(b+bb)/2;
      s2=0;
      if (-arg > sabflg) s2=exp(-arg)/(c*sqrt(a*tfff2*rtev));
      sigVal=sigVal+sigc*sb2*s2;
   } // end if
    if (abs(e-10) < .01 and abs(u-.99219) < .0001) {
    } // end if
   if (sigVal < sigmin) sigVal=0;
   return sigVal;



}


auto sig( double e, double ep, double u, double tev, int nalpha, 
    std::vector<double>& alpha, int nbeta, std::vector<double>& beta,
    std::vector<std::vector<double>>& sab, double bbm, double az, double tevz,
    int lasym, double az2, double teff2, int lat, double cliq, double sb,
    double sb2, double teff, int iinc ){

  /*-------------------------------------------------------------------
   * Compute the differential scattering cross section from e to
   * ep through the angle with the cosine u from endf tabulated
   * data or an analytic law.
   *-------------------------------------------------------------------
   */
   int nb1,na1,i,ib,ia;
   double rtev,bb,a,sigc,b,c,bbb,s,s1,s2,s3,arg;
   double tfff,tfff2,rat;
   double sigmin=1.e-10, sabflg=-225.e0, amin=1.e-6, test1=0.2e0,
          test2=30.e0;
   double sigVal;

   // common factors.
   rtev=1/tev;
   bb=(ep-e)*rtev;
   a=(e+ep-2*u*sqrt(e*ep))/(az*tev);
   if (a < amin) { a=amin; }
   sigc=sqrt(ep/e)*rtev/2;
   b=abs(bb);
   c=sqrt(4*M_PI);

   // tabulated s(alpha,beta).


   // free-gas scattering.
   if ( iinc != 2 ){ 
     std::cout << "200" << std::endl; 

     if (iinc != 1){
       // other options not yet implemented.
       std::cout << "300" << std::endl; 
       // std::cout << "call error('sig','illegal option.',' ')" << std::endl;
       sigVal = 0;
       throw std::exception();
     }

     s=0;
     arg=(a+bb)*(a+bb)/(4*a);
     if (-arg > sabflg) s=exp(-arg)/(c*sqrt(a));
     sigVal=sigc*sb*s;
     if (sigVal < sigmin) sigVal=0;
     return sigVal;
   } // end free gas scattering option


   if (lat == 1) b=b*tev/tevz;
   if (lat == 1) a=a*tev/tevz;

   if (a > alpha[nalpha-1]) {
     // go to 170
     return doSCTApproximation( lat, a, b, c, tevz, rtev, tfff, teff, tfff2, 
         teff2, arg, az, az2, sigc, s2, sigVal, u, sb2, e, tev, sigmin, 
         sabflg, bb, s, sb );
   }

   if (lasym == 1) {
      bbm=bb;
      if (lat == 1) bbm=bb*tev/tevz;
      if ( bbm > beta[nbeta-1] or bbm < beta[0] ){
        // go to 170
        return doSCTApproximation( lat, a, b, c, tevz, rtev, tfff, teff, 
            tfff2, teff2, arg, az, az2, sigc, s2, sigVal, u, sb2, e, tev, 
            sigmin, sabflg, bb, s, sb );
      } 
   } // end if 
   else {
     if ( b > beta[nbeta-1] ){
        return doSCTApproximation( lat, a, b, c, tevz, rtev, tfff, teff, 
            tfff2, teff2, arg, az, az2, sigc, s2, sigVal, u, sb2, e, tev, 
            sigmin, sabflg, bb, s, sb );
      }
   } // end if

   nb1=nbeta-1;
   na1=nalpha-1;
   bbb=b;
   if (lasym == 1 and bb < 0.0 ) bbb=-b;
   for ( int i = 1; i <= nb1; ++i ){
      ib=i;
      if (bbb < beta[i+1-1]) break;
   } // end do

   for ( int i = 1; i <= na1; ++i ){
      ia=i;
      if (a < alpha[i+1-1]) break;
   } // end do

   //if (cliq == 0.0 or a >= alpha(1)) go to 150
   //if (lasym == 1) go to 150
   //if (b > test1) go to 150
   s=sab[1-1][1-1]+log(alpha[1-1]/a)/2-cliq*b*b/a;
   if (s < sabflg) s=sabflg;
   // go to 160
  // 150 continue
   //if (a*az < test2 and b < test2) go to 155
   //if (sab(ia,ib) <= sabflg) go to 170
   //if (sab(ia+1,ib) <= sabflg) go to 170
   //if (sab(ia,ib+1) <= sabflg) go to 170
   //if (sab(ia+1,ib+1) <= sabflg) go to 170
  //155 continue
   if (ia+1 == nalpha) ia=ia-1;
   if (ib+1 == nbeta) ib=ib-1;
   //call terpq(alpha(ia),sab(ia,ib),alpha(ia+1),sab(ia+1,ib),&
   //  alpha(ia+2),sab(ia+2,ib),a,s1)
   //call terpq(alpha(ia),sab(ia,ib+1),alpha(ia+1),sab(ia+1,ib+1),&
   //  alpha(ia+2),sab(ia+2,ib+1),a,s2)
   //call terpq(alpha(ia),sab(ia,ib+2),alpha(ia+1),sab(ia+1,ib+2),&
   //  alpha(ia+2),sab(ia+2,ib+2),a,s3)
   //call terpq(beta(ib),s1,beta(ib+1),s2,beta(ib+2),s3,bbb,s)
  //160 continue
   sigVal=0;
   if (s-bb/2 > sabflg) sigVal=exp(s-bb/2);
   sigVal=sigc*sb*sigVal;
   if (sigVal < sigmin) sigVal=0;
   return sigVal;

}


