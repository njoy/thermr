#include <iostream>
#include <vector>
#include <cmath>

auto gatef2( const std::vector<double>& temp, std::vector<double>& eftmp2,
  int mat ){ 
 /* 
  *-------------------------------------------------------------------
  * Default effective temperatures for second atom in mixed
  * moderators from GA report:
  *       1095      benzine (c6h6)
  *       1099      beo
  *-------------------------------------------------------------------
  */
  int i, j, jmat, ntabl = 16;
  double diff, test;
  std::vector<double> tabl = { 1095., 296.0, 685.54, 1095., 350., 
     712.02, 1095., 400., 738.97, 1095., 
     450., 768.10, 1095., 500., 799.22,
     1095., 600., 866.63, 1095., 800., 
     1017.3,  1095., 1000., 1182.3, 
     1099., 296.0, 427.8, 1099., 400., 
     502.8, 1099., 500., 584.3, 1099., 
     600., 671.3, 1099., 700., 761.6, 
     1099., 800., 854.2, 1099., 1000., 
     1043.7, 1099., 1200., 1236.6 };

   for ( size_t i = 0; i < temp.size(); ++i ){
      if (eftmp2[i] == 0) {
        for ( auto j = 0; j < ntabl; ++j ){
            jmat = tabl[j];
            //jmat=nint(tabl(1,j));
            if (jmat == mat) {
               test = 5;
               diff = tabl[j] - temp[i];
               //diff=tabl(2,j)-temp[i];
               if (std::abs(diff) <= test) {
                  eftmp2[i]=tabl[j];
                  //eftmp2[i]=tabl(3,j);
              } // end if 
           } // end if 
         } // end for
         if (eftmp2[i] == 0) eftmp2[i]=temp[i];
      } // end if 
   } // end for 
}


