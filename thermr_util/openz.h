#include <iostream>
#include <vector>
#include <string>


auto openz( int lun, int newInt ){
 /*--------------------------------------------------------------------
  * System-dependent routine to open files.
  * No action if iabs(lun) ge 5 and le 7.
  * Mode--coded (formatted) if lun gt 0
  *       binary (unformatted) if lun lt 0
  * Destroy--on close or job termination if iabs(lun) ge 10 and
  *                                      iabs(lun) lt 20
  * Status--if new=1, destroy lun if it already exists, then
  *         open a new version.
  *--------------------------------------------------------------------
  */

  int nun;
  bool there;
  std::string forString, age;
  // character::for*15,age*7,fn*6
  nun = std::abs(lun);

  if ((nun >= 5 and nun <= 7) or nun == 0 ) return;
  if (nun > 99) { std::cout << "illegal unit number" << std::endl; }
  
  // construct file name
  // set format based on sign of unit number
  forString = lun < 0 ? "unformatted" : "formatted";

  if (nun > 10 and nun <= 19){
    // scratch
    age = "scratch";
    
  } 
  return;
}
/*

   ! set format based on sign of unit number
      ! scratch units
      age='scratch'
      open(nun,form=for,status=age)
   else
      ! regular units
      if (new.ne.1) then
         ! existing units
         age='old'
      else
         ! new units
         age='new'
         inquire(file=fn,exist=there)
         if (there) then
            open(nun,file=fn,status='old')
            close(nun,status='delete')
         endif
      endif
      ! open the connection
      open(nun,file=fn,form=for,status=age)
   endif
   return
   end subroutine openz

   */


