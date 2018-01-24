#include <iostream>
#include <vector>


int main(){

  int nendf, nin, nout, nscr, nscr2, matdde, nbin, iprint, ncds, matdpp, 
    natom, ntemp, iinc, iform, mtref, ncdse;
  double za, awr, tol, emax, sb, az, tevz, teff, sb2, az2, teff2, cliq;
  int lasym, lat;

  // array for user temperatures
  std::vector<double> tempr;

  // array for thermal elastic data
  std::vector<double> fl;

  // arrays for thermal inelastic cross section
  std::vector<double> esi, xsi;

  const int nbuf = 1000, nwscr = 500000;

  /* Generate neutron scattering cross sections and point-to-point
   *  scattering kernels in the thermal range.  The coding can
   *  generate incoherent inelastic cross sections and distributions
   *  for a free gas, incoherent inelastic cross sections and
   *  distributions from read-in S(alpha,beta,T) data, coherent
   *  elastic scattering cross sections for crystalline materials,
   *  and incoherent elastic cross sections and angular distributions
   *  for hydrogenous solids.
   *
   *  The pointwise scattering cross sections and distributions are
   *  added to an existing PENDF tape.  Cross sections are added in
   *  MF3 and distributions are written in MF6, both using mtref for
   *  inelastic and mtref+1 for elastic (if any).
   *
   *  Multiple scattering types (i.e., H free and H in H2O) can be
   *  written on one PENDF tape by using different values of mtref
   *  for each thermr run.  If data for one mtref is already on the
   *  tape, it will be replaced with the new cross sections.
   *
   *  The energy grid for coherent scattering is determined
   *  adaptively so as to represent the sharp Bragg edges to a
   *  specified tolerance using linear interpolation.  No angular
   *  information is provided for coherent scattering.  It can be
   *  deduced by subsequent processing from the Bragg edges in the
   *  cross section.
   *
   *  The incident energy grid for incoherent inelastic scattering
   *  is fixed in the code, although it is scaled to higher values
   *  for very large incident energies.  There are two options for
   *  representing the energy and angle distributions for secondary
   *  inelastic neutrons.
   *
   *  The standard approach is to generate an energy grid for the
   *  secondary energy distribution integrated over angle
   *  adaptively.  The angular distribution for each secondary
   *  energy is given using a set of equally probable emission
   *  cosines.  This is E-E'-mu ordering, and it is represented
   *  using a special format in MF6.
   *
   *  An alternate approach introduced in NJOY2010 produces E-mu-E'
   *  ordering using MF6/Law7.  A grid of mu values is generated
   *  adaptively for linear interpolation.  In addition, secondary
   *  energy grids are determined adaptively for each mu value.
   *
   *  The sections of File 6 produced by thermr have several
   *  special features.  The LIP flag is used to identify the
   *  various kinds of thermal data:  LIP=-1 means incoherent
   *  inelastic data, LIP=-2 means incoherent inelastic data, and
   *  LIP=-nbrag means coherent elastic data with nbrag Bragg
   *  edges.  For incoherent elastic or inelastic data using
   *  LAW=1, the angular data are equally probable cosines.
   *
   *  For ENDF 3 to 5 formats, the constants used for coherent
   *  elastic, incoherent elastic, and short-collision-time calcu-
   *  lations are obtained from internal data statements based on
   *  the original general atomic report on the evaluations
   *  (GA-8774 revised, ENDF-269, July 1978).
   *
   *  For ENDF-6 format libraries, these constants are included
   *  in the format.
   *
   *---input specifications (free format)--------------------------
   *
   *  card 1
   *     nendf      endf tape for mf7 data
   *     nin        old pendf tape
   *     nout       new pendf tape
   *  card 2
   *     matde      material desired on endf tape
   *     matdp      material desired on pendf tape
   *     nbin       number of equi-probable angles
   *     ntemp      number of temperatures (default = 1)
   *     iinc       inelastic options
   *                   0     none
   *                   1     compute as free gas
   *                   2     read s(a,b) and compute matrix
   *     icoh       elastic options
   *                   0     none
   *                   1     compute using ENDF6 format data
   *                   --------or for earlier formats
   *                   1     graphite
   *                   2     beryllium
   *                   3     beryllium oxide
   *                  11     polyethylene
   *                  12     h(zrh)
   *                  13     zr(zrh)
   *     iform      output format for inelastic distributions
   *                  0      E-E'-mu ordering (MF6 special)
   *                  1      E-mu-E' ordering (MF6/Law7)
   *     natom      number of principal atoms
   *     mtref      mt for inelastic reaction (221-250 only)
   *     iprint     print option (0=minimum, 1=maximum,
   *                2=max. normal + intermediate results)
   *                (default=0)
   *  card 3
   *     tempr      temperatures (kelvin)
   *  card 4
   *     tol        tolerance
   *     emax       maximum energy for thermal treatment
   *                (for temperatures greater than 3000,
   *                emax and the energy grid are scaled by
   *                temp/3000.  free gas only.)
   *
   *       nendf can be ENDF6 format (e.g., from leapr) while
   *       nin and nout are ENDF4 or 5 format, if desired.
   *
   * ------------------------------------------------------------------
   */
  int icoh, i, icopy, idis, nb, nw, iold, inew;
  int iverp, itemp, nwb, it, ntape, nex, ne, np, isave, lthr;
  double e, enext, emaxs, time, sz2, t, templ, s, temp;
  std::vector<double> ex(2);
  std::vector<double> eftemp, eftmp2;
  double s1099 = 3.76e0, a1099 = 15.86e0, s1095 = 4.74e0, a1095 = 11.90e0,
    eps = 1.e-4, breakVal = 3000.e0, small = 1.e-10, up = 1.00001e0, zero = 0;

  std::vector<double> scr(nwscr,0.0), bufo(nbuf), bufn(nbuf);

  nendf = 24;
  nin = -23;
  nout = -25;

  if (nin == 0){ 
    std::cout << "OH NO! nin equals zero" << std::endl;
    return 0;
  }

  if ( (nin < 0 and nout > 0) or (nin > 0 and nout < 0) ){
    std::cout << "OH NO! Mode conversion not allowed" << std::endl;
  }
  
return 1;
}

