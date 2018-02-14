#include <iostream>
#include <vector>

auto terpq( double x1, double y1, double x2, double y2, double x3, double y3,
  double x, double y ){
  /*-------------------------------------------------------------------
   * Compute y(x) by quadratic interpolation,
   * except use log-lin if x < x1 and lin-lin if x > x3.
   * and use lin-lin if the function takes big steps (corners).
   *-------------------------------------------------------------------
   */
  double b, c, sabflg = -225, step = 2;

   if (x < x1) {
      if (y1 > y2) {
         y = y1;
      } 
      else {
         // call terp1(x1,y1,x2,y2,x,y,3);
      }
   }
   else if (x > x3) {
      if (y3 > y2) {
         y = y3;
      }
      else {
         // call terp1(x2,y2,x3,y3,x,y,2)
      }
   }
   else if (std::abs(y1-y2) > step or std::abs(y2-y3) > step) {
      if (x < x2) {
         // call terp1(x1,y1,x2,y2,x,y,2)
      }
      else {
         // call terp1(x2,y2,x3,y3,x,y,2)
      }
   }
   else {
      b = (y2-y1)*(x3-x1)/((x2-x1)*(x3-x2));
      b = b-(y3-y1)*(x2-x1)/((x3-x1)*(x3-x2));
      c = (y3-y1)/((x3-x1)*(x3-x2)); 
      c = c-(y2-y1)/((x2-x1)*(x3-x2));
      y = y1+b*(x-x1)+c*(x-x1)*(x-x1);
   }
   if (y < sabflg) y = sabflg;
}


