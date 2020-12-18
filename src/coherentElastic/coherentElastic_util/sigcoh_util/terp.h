

template <typename Range, typename Float>
auto do260( const Range& x, const Range& y, const Float& arg, int l, int il ){
  // interpolation section
  Float sum = 0, p, pk;
  for ( int i = 0; i < il; ++i ){
    p  = 1.0; pk = 1.0;
    for ( int ip = 0; ip < il; ++ip ){
      if ( ip != i ){ 
        p  *= arg      - x[l+ip-1];
        pk *= x[l+i-1] - x[l+ip-1];
      }
    }
    sum += p * y[l+i-1] / pk;
  } 
  return sum;
}





template <typename Range, typename Float>
auto terp( Range x, Range y, Float arg, int il1, int nl = 0 ){
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
  int il=il1, l, iadd, il2, ilow, ihi, iusel, iuseh, ibeg, iend, last, m; 
  if (nl == 0){ nl = x.size(); }

  // If order of interpolation is too large, set it to the size of input vector
  if ( nl <= il ){
    l = 1; il = nl;              // Make sure nl is not less than il. 
    return do260(x,y,arg,l,il);  // Continue with interpolation
  } 
   
  il2 = il / 2;
  // check if tables in increasing or decreasing sequence
  if (x[0] <= x[nl-1]) {        // increasing sequence
    ilow  = il2 + 1;            ihi   = nl - il2 - il%2;
    iusel = 1;                  iuseh = nl - il + 1;
    ibeg  = ilow + 1;           iend  = ihi - 1;
    last  = iend - il2 + 1;     iadd  = 0;
  }
  else {                        // decreasing sequence
    ilow  = nl - il2;           ihi   = il2 + il%2 + 1;
    iusel = nl - il + 1;        iuseh = 1;
    ibeg  = ihi + 1;            iend  = ilow - 1;
    last  = 2;                  iadd  = 1 - il%2;
  }

  //if (arg > 0.17599 and arg < 0.1759987){ 
  //    std::cout << nl << "   " << iadd << "   " << ilow << "   " << ihi << std::endl; 
 // }
  // If arg is approximately equal to lowest value, return lowest known value
  // If arg is approximately equal to highest value, return highest value
  if ( std::fabs(arg-x[ilow-1]) < arg*1e-10 ) { return y[ilow-1]; }
  if ( std::fabs(arg-x[ihi-1])  < arg*1e-10 ) { return y[ihi -1]; }

  // If arg is lower than x[ilow-1], then set l = iusel and interpolate
  // If arg is higher than x[ihi-1], then set l = iuseh and interpolate
  if ( arg < x[ilow-1] ){ return do260( x, y, arg, iusel, il ); }
  if ( arg > x[ihi-1]  ){ return do260( x, y, arg, iuseh, il ); }

  // If arg is a reasonable size, continue
  for ( int n = ibeg; n <= iend; ++n ){
    // If we have a decreasing seqeuence, and the order of interpolation is less
    // than the size of the vector
    m = iusel > 1 ? nl - n + 1 : n;

    // If arg is approximately equal to a value in my x vector, just return its
    // corresponding value in the y vector
    if ( std::fabs(x[m-1] - arg) < 1e-10*arg ){ return y[m-1]; }

    // If while iterating through my x vector, I pass over arg, then I should 
    // calculate l and interpolate
    if ( x[m-1] > arg ) { return do260( x, y, arg, m-il2+iadd, il ); }
  }
  
  // If after iterating through x I haven't found my arg value to interpolate,
  // I assume the latest possible value and interpolate
  l = last;
  return do260( x, y, arg, l, il );
   
}
