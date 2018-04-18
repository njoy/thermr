


auto coh2( int lat, std::fstream& inew, int ne, int nex, double temp, 
  std::fstream& iold, double emax, int natom, std::vector<double>& fl, 
  std::vector<double>& bufo, std::vector<double>& bufn ){
 /*-------------------------------------------------------------------
  * Compute the coherent scattering cross sections for a crystalline
  * material. The cross section is computed on an energy grid
  * chosen adaptively to represent the actual function within a
  * given tolerance with linear-linear interpolation.
  *
  * The energy grid is stored on loada/finda scratch files that were 
  * used for the elastic cross section. The elastic cross section is 
  * converted to the coherent grid using Lagrangian interpolation (terp).
  *-------------------------------------------------------------------
  */
  int nl, imax, nx, nj, i, nlt, nlt1, nb, nw, ix, j, iex, il, isave, ltt;
  int nlmax = 6;
  double enext = 0, en, xm, ym, test;
  std::vector<double> s(nlmax), ej(20), ex(20), x(5), y(5), z(5), b(12);
  double half = 0.5, small = 1.0e-10, tolmin = 1.0e-6, eps = 3.0e-5; 


  // -------------- Begin ----------------
  // initialize.
  nl = 1;

  imax = 20;
  if (nl > nlmax) { 
    std::cout << "too many legendre orders" << std::endl; 
    return; 
  }

  nx = nl + 1;
  nj = nex + nl;
  for ( int i = 0; i < nl; ++i ){
     ex[nex+i] = 0;
  }
  nlt = 5;
  if (ne < nlt) nlt = ne;
  nlt1 = nlt - 1;

  //allocate(stk(nx,imax)) --> Use eigen? 
  // Oh well, we'll use this for now. 

  std::vector<std::vector<double>> stk (nx, std::vector<double> (imax) );


  // determine the energy grid adaptively and
  // store the cross sections in a scratch file.
  
  int nbragg = nl, k = 0;
  double e = 0, scon = 0;
  std::vector<double> p;

  auto wrk = sigcoh( e, enext, s, nl, lat, temp, emax, natom, fl, p, k, scon );

  ix = 1;
  j  = 0;
  iex = 0;



  // --------------  100  ----------------
  bool do100 = true, do100Inner = true;

  // Remove this once you're convinced your code is competent 
  int counter = 0;

  while ( do100 ){
    while ( do100Inner ){

      iex = iex + 1;
      // Read in the first nex values from iold file and put them into ex
      finda( iex, nex, iold, ex, bufo, bufo.size() );

      if (ex[0] > enext*(1+small)) { do100Inner = false; }
      else {
        x[0] = ex[0];
        y[0] = ex[1];
        z[0] = ex[nex-1];
        j = j + 1;
        if ( counter > 10 ){ break; }
        ++counter;
        // Take the elements of ex and put them in file inew, also in buffer bufn
        loada( j, nj, inew, bufn.size(), ex, bufn );
      }

    }

    // --------------  105  ----------------
    ix = ix + 1;
    x[ix] = ex[1];
    y[ix] = ex[2];
    z[ix] = ex[nex];

    // 105 --> 100
    if (ix >= nlt) { do100 = false; }

  }

  // (continuing with 105)
  
  // prime stack with first bragg edge
  e = enext;
  wrk = sigcoh( e, enext, s, nl, lat, temp, emax, natom, fl, p, k, scon );
   
  stk[0][0] = e;
   
  for ( int il = 0; il < nl; ++il ){
    stk[1+il][0] = s[il];
  } 

  i = 1;


  bool do120 = true;
  // --------------  120  ----------------
  // add next bragg edge to stack
  while (do120){
    e = enext;
    wrk = sigcoh( e, enext, s, nl, lat, temp, emax, natom, fl, p, k, scon );
    upstk(e, s, stk, nl, nx, i );
    // make sure input grid points are included
   
    // --------------  125  ----------------
   
    for ( int ix = 0; ix < nlt; ++ix ){
      if (x[ix] > stk[0][i]*(1+small)){  // go to 135
        // 135 continue
        if (x[ix] >= stk[0][i-1]*(1-small)){ // go to 140
          e = x[ix];
          wrk = sigcoh( e, enext, s, nl, lat, temp, emax, natom, fl, p, k, scon );
          upstk(e, s, stk, nl, nx, i );
          break; // Don't continue with 125's loop from ix = 0 --> nlt
        }
        else {
      }
    }
  
 

    // if 190 tells us we don't have to go back to 120, do120 = false
  } 
 


}



