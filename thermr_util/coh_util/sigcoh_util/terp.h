#include <iostream>
#include <vector>


auto do260( std::vector<double> x, std::vector<double> y, double arg, int l, 
  int il ){

  // interpolation section
  double sum = 0, p, pk, in, inp;
  for ( size_t i = 0; i < il; ++i ){
    p = 1;
    pk = 1;
    in = l + i - 1;
    for ( size_t ip = 0; ip < il; ++ip ){
      if ( ip == i ){ continue; }
      inp = l + ip - 1;
      p  *= arg   - x[inp];
      pk *= x[in] - x[inp];
    }
    sum += p * y[in] / pk;
  } 
  return sum;
}


auto do250( std::vector<double> x, std::vector<double> y, double arg, int il, 
  int m, int il2, int iadd  ){

  int l = m-il2 + iadd;
  return do260( x, y, arg, l, il );
   
}




auto terp( std::vector<double> x, std::vector<double> y, int nl, double arg, 
  int il1 ){

  /*-------------------------------------------------------------------
   * This function does Lagrangian interpolation or
   * extrapolation of ilth order on x and y for arg
   * when the tables are either increasing or decreasing
   *      x          x array, independent variable
   *      y          y array, dependent variable
   *      nl         number of entries in tables of x and y
   *      arg        independent variable value
   *      il         order of interpolation
   *-------------------------------------------------------------------
   */

   int il,l,iadd,il2,ilow,ihi,iusel,iuseh,ibeg,iend,last,n,m,i,in,ip,inp;
   double sum,p,pk;
   double small=1.e-10;

   il = il1;
   // if (nl > il) go to 120
   // if (nl == il) go to 110
   // not enough entries in tables for this order interpolation
   il = nl;
  // 110 continue
   l = 1;
  //  go to 260
  // 120 continue
   il2 = il / 2;
   iadd = il%2;
   // check if tables in increasing or decreasing sequence
   if (x[0] <= x[nl-1]) {
      // increasing sequence
      ilow = il2 + 1;
      ihi = nl - il2 - iadd;
      iusel = 1;
      iuseh = nl - il + 1;
      ibeg = ilow + 1;
      iend = ihi - 1;
      last = iend - il2 + 1;
      iadd = 0;
    }
    else {
      // decreasing sequence
      ilow = nl - il2;
      ihi = il2 + iadd + 1;
      iusel = nl - il + 1;
      iuseh = 1;
      ibeg = ihi + 1;
      iend = ilow - 1;
      last = 2;
      iadd = 1 - iadd;
    }

    /*

   // checks if arg is smaller than table values
  // if (std::abs(arg-x(ilow)) < small*arg) go to 160
  // if (arg > x(ilow)) go to 170
   // smaller than smallest table value
   l = iusel;
   //go to 260
   
  160 continue
   terp = y(ilow);
   return; // go to 300
   // checks if arg is greater than table values
  170 continue
    if (abs(x(ihi)-arg).lt.small*arg) go to 190
   if (x(ihi) > arg) go to 200
   // arg greater than table value
   l = iuseh;
   go to 260
  190 continue
   terp = y(ihi);
   return; // go to 300
   ! searches x array to bracket arg
  200 continue
   do 230 n=ibeg,iend
   if (iusel.gt.1) go to 210
   m=n
   go to 220
  210 continue
   m=nl-n+1
  220 continue
   if (abs(x(m)-arg).lt.small*arg) go to 240
   if (x(m).gt.arg) go to 250
  230 continue
   l=last
   go to 260
   // equals argument, return ok
  240 continue
   terp=y(m)
   return; // go to 300
   // eureka
  250 continue
   l=m-il2+iadd
  260 continue
   // interpolation section
   sum=0
   for ( size_t i = 0; i < il; ++i ){
      p = 1;
      pk = 1;
      in = l + i - 1;
      for ( size_t ip = 0; ip < il; ++ip ){
         if (ip != i){
            inp = l + ip - 1;
            p = p * (arg-x(inp));
            pk = pk * (x(in)-x(inp));
          }
      }
      sum=sum+p*y(in)/pk
   } 
   terp = sum;
   
  */

  }
