
auto terpa(double y, double x, int idis,
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
  double shade=1.00001, xbig = 1.0e12, xnext;

  // set up limits and pointers
  nr=round(a[4]);
  np=round(a[5]);
  if (ir > nr) ir=nr;
  if (ip > np) ip=np;
  jr=5+2*ir;
  jp=5+2*nr+2*ip;
  idis=0;

  while (true) {
    // locate interpolation interval and law for x
    std::cout << "110  " << ip << "    " << jp << "    " << a[jp-1] << std::endl;
    if (x < a[jp-1]){
      std::cout << "120" << std::endl;
      if (x > a[jp-3]){
        std::cout << "130" << std::endl;
        // interpolate for y in this interval
        intVar=round(a[jr+1-1]);
        // call terp1(a(jp-2),a(jp-1),a(jp),a(jp+1),x,y,intVar);
        xnext=a[jp-1];
        if (intVar == 1) idis=1;
        if (ip == np) return xnext;
        if (a[jp+1] == xnext) idis=1;
        return xnext;
      }
      if (x == a[jp-3]) {
        std::cout << "140" << std::endl;
        y=a[jp-1-1];
        intVar=round(a[jr+1-1]);
        xnext=a[jp-1];
        if (intVar == 1) idis=1;
        if (ip == np) return xnext;
        if (a[jp+1] == xnext) idis=1;
        return xnext;
      }
      if (ip == 2) {
        // special branch for x below first point
        std::cout << "170" << std::endl;
        y = 0;
        idis = 1;
        return a[jp-3];
      }
      // move down
      jp=jp-2;
      ip=ip-1;
      if ( ir != 1 ){
        it=round(a[jr-3]);
        if ( ip <= it ){
          jr = jr-2;
          ir = ir-1;
        }
      }
      //go to 110
    }

    if (ip == np) {
      std::cout << "150" << std::endl;
      // special branch for last point and above
      if (x < shade*a[jp-1]) {
         std::cout << "160" << std::endl;
         y = a[jp];
         return ( y > 0 ) ? shade * shade * a[jp-1] : xbig;
      }
      y=0;
      return xbig;
    }

    // move up
    jp=jp+2;
    ip=ip+1;
    it=round(a[jr-1]);
    if ( ip < it ){
      jr=jr+2;
      ir=ir+1;
    }
    //go to 110
  }
}


