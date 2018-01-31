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
  int i, j, jmat;
  double diff, test;
  std::vector<double> tabl1, tabl2, tabl3;
  tabl1 = { 1095, 1095, 1095, 1095, 1095, 1095, 1095, 1095, 1099, 1099, 1099, 
    1099, 1099, 1099, 1099, 1099 };
  tabl2 = { 296, 350, 400, 450, 500, 600, 800, 1000, 296, 400, 500, 600, 700, 
    800, 1000, 1200 };
  tabl3 = { 685.54, 712.02, 738.97, 768.10, 799.22, 866.63, 1017.3,  1182.3, 
    427.8, 502.8, 584.3, 671.3, 761.6, 854.2, 1043.7, 1236.6 };
 
   for ( size_t i = 0; i < temp.size(); ++i ){
      if (eftmp2[i] == 0) {
        for ( auto j = 0; j < tabl2.size(); ++j ){
            jmat = tabl1[j];
            if (jmat == mat) {
               test = 5;
               diff = tabl2[j] - temp[i];
               if (std::abs(diff) <= test) {
                  eftmp2[i]=tabl3[j];
              } // end if 
           } // end if 
         } // end for
         if (eftmp2[i] == 0) eftmp2[i] = temp[i];
      } // end if 
   } // end for 
}


