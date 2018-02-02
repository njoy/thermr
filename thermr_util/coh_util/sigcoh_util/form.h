#include <iostream>
#include <vector>


auto form( int lat, int l1, int l2, int l3 ){
  /*-------------------------------------------------------------------
   * Compute form factors for the specified lattice.
   *       lat=1  graphite
   *       lat=2  be
   *       lat=3  beo
   *-------------------------------------------------------------------
   */
   int i;
   double beo1 = 7.54, beo2 = 4.24, beo3 = 11.31, formVal;

   if (lat == 1){ 
      // graphite
      i = l3 / 2;
      if ((2*i) != l3){
         formVal = std::pow(sin(M_PI*(l1-l2)/3),2);
      }
      else {
         formVal = ( 6 + 10*cos(2*M_PI*(l1-l2)/3) ) / 4;
      } 
   }
   else if (lat == 2){ 
      // beryllium
      formVal = 1 + cos(2*M_PI*(2*l1+4*l2+3*l3)/6);
   }
   else if (lat == 3){
      // beryllium oxide
      formVal=(1+cos(2*M_PI*(2*l1+4*l2+3*l3)/6))*(beo1+beo2+beo3*
        cos(3*M_PI*l3/4));
   }
   return formVal;
}

