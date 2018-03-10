
auto terpa(double y, double x, double xnext, int idis,
    std::vector<double> a, int ip, int ir){
  /*--------------------------------------------------------------------
   * Interpolate for y(x) in the TAB1 record packed in a.  Return  0 
   * if x is outside the range of the table.   Here xnext is the next
   * data grid point greater than x.  On entry, ip and ir are starting
   * estimates for the first data point greater than x and the
   * corresponding interpolation range.  Initialize them to 2 and 1
   * before first call to routine.
   *--------------------------------------------------------------------
   */
   int nr,np,jr,jp,intVar,it;
   double shade=1.00001, xbig = 1e12;

   // set up limits and pointers
   nr=round(a[5-1]);
   np=round(a[6-1]);
   if (ir > nr) ir=nr;
   if (ip > np) ip=np;
   jr=5+2*ir;
   jp=5+2*nr+2*ip;
   idis=0;

   // locate interpolation interval and law for x
   //110 continue
   if (x < a(jp)) go to 120
   if (ip == np) {
    // go to 150
    // special branch for last point and above
    // 150 continue
    if (x < shade*a(jp)) {
       // go to 160
       // 160 continue
       y=a[jp+1-1];
       xnext=xbig;
       if (y >  0 ) xnext=shade*shade*a[jp-1];
       return;
     }
     y=0;
     xnext=xbig;
     return;

   }


   // move up
   jp=jp+2;
   ip=ip+1;
   it=round(a(jr))
   if (ip <= it) go to 110
   jr=jr+2;
   ir=ir+1;
   //go to 110
   //120 continue
   if (x > a[jp-2-1]){
     // go to 130
     // interpolate for y in this interval
     // 130 continue
     intVar=round(a[jr+1-1]);
     // call terp1(a(jp-2),a(jp-1),a(jp),a(jp+1),x,y,intVar);
     xnext=a[jp-1];
     if (intVar == 1) idis=1;
     if (ip == np) return;
     if (a[jp+2-1] == xnext) idis=1;
     return;
   }

   if (x == a[jp-2-1]) {
     // go to 140
     // 140 continue
     y=a[jp-1-1];
     intVar=round(a[jr+1-1]);
     xnext=a[jp-1];
     if (intVar == 1) idis=1;
     if (ip == np) return;
     if (a[jp+2-1] == xnext) idis=1;
     return;
   }


   if (ip == 2) {
     // go to 170
     // special branch for x below first point
     // 170 continue
     y=0;
     xnext=a[jp-2-1];
     idis=1;
     return;
   }

   //! move down
   jp=jp-2;
   ip=ip-1;
   if (ir == 1) go to 110
   it=round(a[jr-2-1]);
   if (ip > it) go to 110
   jr=jr-2;
   ir=ir-1;
   go to 110

   end subroutine terpa


