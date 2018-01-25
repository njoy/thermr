#include <iostream>
#include <string>
#include <fstream>

int tpidio( int nin, int nout, int nscr, int nb, int nw ){
 /*--------------------------------------------------------------------
  * Utility routine for ENDF bcd and binary tapes.  Read, write,
  * and/or convert the tape identification record to/from a.  If any
  * unit is zero, it is not used.  Positive units are bcd, and
  * negative units are binary.
  *--------------------------------------------------------------------
  */
  std::vector<double> a;
  std::vector<double> rb(17);
  std::string hb;

  // equivalence (rb(1),hb(1))
  
  int inin, inout, inscr, i;

  if (nin < 0){
    inin = abs(nin);
    //read(inin) math,mfh,mth,nb,nw,(a(i),i=1,17)
    std::ofstream ninFile;
    ninFile.open( "tape" + std::to_string(nin) );
    ninFile.close();
    
  }
  else if (nin > 0) {
    //read(nin,'(16a4,a2,i4,i2,i3,i5)') (hb(i),i=1,17),math,mfh,mth
    nw=17;
    for ( int i = 0; i < 17; ++i ){
      a[i] = rb[i];
    } 
  }

  /*
   !--output
   nb=0
   inout=iabs(nout)
   inscr=iabs(nscr)
   if (nout.lt.0) then
      write(inout) math,mfh,mth,nb,nw,(a(i),i=1,17)
   endif
   if (nscr.lt.0) then
      write(inscr) math,mfh,mth,nb,nw,(a(i),i=1,17)
   endif
   if (nout.le.0.and.nscr.le.0) return
   if (nout.gt.0) then
      do i=1,17
         rb(i)=a(i)
      enddo
      write(nout,'(16a4,a2,i4,i2,i3,i5)') (hb(i),i=1,17),math,mfh,mth,nsh
      nsh=nsh+1
   endif
   if (nscr.gt.0) then
      do i=1,17
         rb(i)=a(i)
      enddo
      write(nscr,'(16a4,a2,i4,i2,i3,i5)') (hb(i),i=1,17),math,mfh,mth,nsc
      nsc=nsc+1
   endif
   return
   end subroutine tpidio


   */
  return 0;
}

