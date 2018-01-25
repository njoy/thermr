#include <iostream>
#include <vector>
//#include "thermr_util/openz.h"
#include "thermr_util/tpidio.h"
#include <fstream>


int main(){

  int nendf, nin, nout, nscr, nscr2, matde, nbin, iprint, ncds, matdp, 
    natom, ntemp, iinc, iform, mtref, ncdse;
  double za, awr, tol, emax, sb, az, tevz, teff, sb2, az2, teff2, cliq;
  int lasym, lat;

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
  
  iprint=0;
  ntemp=1;
  std::vector<double> tempr(ntemp), eftemp(ntemp,0.0), eftmp2(ntemp,0.0);
  if ( matde == 0 ){ 
    for ( int i = 0; i < ntemp; ++i ){ eftemp[i] = tempr[i]; }
  }

  matde = 101, matdp = 600, nbin = 8, ntemp = 1, iinc = 2, icoh = 1, iform = 0,
    natom = 1, mtref = 222, iprint = 2;
  

  if ( mtref < 221 or mtref > 250 ){ 
    std::cout << "OH NO! Illegal reference mt" << std::endl;
  }

  nendf = 24;
  nin = 23;
  nout = 25;
  std::ofstream nendfFile;
  std::ofstream ninFile;
  std::ofstream noutFile;
  nendfFile.open( "tape" + std::to_string(nendf) );
  ninFile.open  ( "tape" + std::to_string(nin) );
  noutFile.open ( "tape" + std::to_string(nout) );
  nendfFile.close();
  ninFile.close();
 
  
  /*
   !--check for endf-6 format data
   call tpidio(nendf,0,0,scr,nb,nw)
   call contio(nendf,0,0,scr,nb,nw)
   call contio(nendf,0,0,scr,nb,nw)
   if (n1h.ne.0) then
      iverf=4
   else if (n2h.eq.0) then
      iverf=5
   else
      iverf=6
   endif
   call skiprz(nendf,-1)

   !--set up internal data for older endf versions
   if (iverf.lt.6) then

      !--check for mixed s(a,b) cases (beo, benzine)
      !--supply appropriate parameters
      sz2=0
      az2=0
      if (matde.eq.1099) sz2=s1099
      if (matde.eq.1099) az2=a1099
      if (matde.eq.1095) sz2=s1095
      if (matde.eq.1095) az2=a1095
      sb2=sz2
      if (az2.ne.zero) sb2=sz2*((az2+1)/az2)**2

         !--default effective temperatures to standard values if
         !--available, otherwise set them to the material temperature
         if (matde.ne.0) then
            call gateff(tempr,eftemp,ntemp,matde)
            if (sz2.ne.zero) then
               call gatef2(tempr,eftmp2,ntemp,matde)
            endif
         endif
      endif

   !--write parameters.
   write(nsyso,'(/&
     &'' unit for endf tape ................... '',i10/&
     &'' unit for input pendf tape ............ '',i10/&
     &'' unit for output pendf tape ........... '',i10)')&
     nendf,nin,nout
   write(nsyso,'(/&
     &'' material to be processed (endf) ...... '',i10/&
     &'' material to be processed (pendf) ..... '',i10/&
     &'' number of angle bins ................. '',i10/&
     &'' number of temperatures ............... '',i10/&
     &'' inelastic option ..................... '',i10/&
     &'' elastic option ....................... '',i10/&
     &'' MF6 format option .................... '',i10/&
     &'' number of principal atoms ............ '',i10/&
     &'' reference mt ......................... '',i10/&
     &'' print option (0 min, 1 max) .......... '',i10)')&
     matde,matdp,nbin,ntemp,iinc,icoh,iform,natom,mtref,iprint
   write(nsyso,'(&
     &'' temperatures (kelvin) ................ '',1p,e10.4)')&
     tempr(1)
   if (ntemp.gt.1) write(nsyso,'(40x,1p,e10.4)')&
     (tempr(i),i=2,ntemp)
   write(nsyso,'(&
     &'' tolerance ............................ '',1p,e10.4/&
     &'' max energy for thermal treatment ..... '',e10.4)')&
     tol,emax
   if (iinc.eq.2.and.iverf.lt.6) then
      write(nsyso,'(/&
        &'' parameters for sct app.''/&
        &'' -----------------------'')')
      write(nsyso,'(&
        &'' effective temperatures ............... '',1p,e10.4)')&
        eftemp(1)
      if (ntemp.gt.1) write(nsyso,'(40x,1p,e10.4)')&
        (eftemp(i),i=2,ntemp)
      if (sz2.ne.zero) then
         write(nsyso,'(&
           &'' free xsec for second atom ............ '',1p,e10.4/&
           &'' mass ratio for second atom ........... '',1p,e10.4/&
           &'' eff. temps for second atom ........... '',1p,e10.4)')&
           sz2,az2,eftmp2(1)
         if (ntemp.gt.1) write(nsyso,'(40x,1p,e10.4)')&
           (eftmp2(i),i=2,ntemp)
      endif
   endif
   write(nsyso,'(/'' endf uses endf-'',i1,'' format'')') iverf

   !--initialize i/o units
   iold=10
   inew=11
   nscr=12
   if (nout.lt.0) nscr=-nscr
   nscr2=13
   if (nout.lt.0) nscr2=-nscr2
   call openz(-iold,1)
   call openz(-inew,1)
   call openz(nscr,1)
   call openz(nscr2,1)
   call repoz(nendf)
   call repoz(nout)
   call repoz(nin)
   call repoz(nscr2)

   !--copy through to desired material
   nsh=0
   call tpidio(nin,nout,0,scr,nb,nw)
  110 continue
   call contio(nin,0,0,scr,nb,nw)
   if (math.eq.matdp) go to 120
   if (math.lt.0)&
     call error('thermr','desired material not on pendf tape.',' ')
   call contio(0,nout,0,scr,nb,nw)
   call tomend(nin,nout,0,scr)
   go to 110
  120 continue
   call contio(nin,0,0,scr(7),nb,nw)
   if (n1h.ne.0) then
      iverp=4
   else if (n2h.eq.0) then
      iverp=5
   else
      iverp=6
   endif
   call skiprz(nin,-1)
   write(nsyso,'(/'' pendf uses endf-'',i1,'' format'')') iverp

   !--loop over desired temperatures
   itemp=1
   icopy=0
  130 continue
   call contio(0,0,nscr2,scr,nb,nw)
   if (math.ne.matdp)&
     call error('thermr','desired material not on pendf tape.',' ')
   nwb=0
   za=scr(1)
   awr=scr(2)
   az=scr(2)
   if (iverp.ge.5) call contio(nin,0,nscr2,scr,nb,nw)
   if (iverp.ge.6) call contio(nin,0,nscr2,scr,nb,nw)
   call hdatio(nin,0,nscr2,scr,nb,nw)
   t=scr(1)
   if (abs(t-tempr(itemp)).le.eps*tempr(itemp)) go to 180
   if (t.gt.tempr(itemp))&
     call error('thermr','desired temperature not on tape.',' ')

   !--check for skipped temperatures to be copied
   if (itemp.eq.1) icopy=icopy+1
   if (itemp.eq.1) go to 160
   templ=tempr(itemp-1)
   it=0
  150 continue
   it=it+1
   if (it.ge.itemp) go to 160
   if (abs(t-tempr(it)).le.eps*tempr(it)) go to 150
   if (abs(t-templ).le.eps*templ) go to 150
   if (t.lt.templ) go to 150
   icopy=icopy+1
   go to 150
  160 continue
   if (icopy.eq.0) call repoz(nscr2)
   ntape=nscr2
   if (icopy.eq.0) ntape=0
   call tomend(nin,0,ntape,scr)
   call contio(nin,0,0,scr,nb,nw)
   go to 130
  180 continue
   call findf(matdp,3,2,nin)

   !--save elastic cross section on scratch file.
   call contio(nin,0,0,scr,nb,nw)
   nex=2
   ne=0
   e=0
   call gety1(e,enext,idis,s,nin,scr)
   emaxs=emax
   if (t.gt.break) emax=emax*t/break
  200 continue
   e=enext
   call gety1(e,enext,idis,s,nin,scr)
   ex(1)=e
   ex(2)=s
   if (e.gt.emax*(1+small)) ex(2)=0
   ne=ne+1
   if (e.gt.emax*(1+small)) ne=-ne
   call loada(ne,ex,nex,inew,bufn,nbuf)
   if (ne.lt.0) go to 210
   if (enext.gt.emax*(1+small)) enext=emax
   if (abs(e-emax).lt.small*emax) enext=up*emax
   go to 200
  210 continue
   np=-ne
   emax=emaxs
   isave=iold
   iold=inew
   inew=isave
   call skiprz(nin,-1)

   !--set up for elastic calculation
   if (iverf.ge.6.and.nendf.ne.0) then
      call findf(matde,7,0,nendf)
      call contio(nendf,0,0,scr,nb,nw)
      if (mth.eq.2) then
         lthr=l1h
         icoh=10*lthr
         temp=tempr(itemp)
         call rdelas(temp,lthr,nendf,nwb)
      endif
   endif
   write(*,*) "HERE"
   return
*/
return 1;
}

