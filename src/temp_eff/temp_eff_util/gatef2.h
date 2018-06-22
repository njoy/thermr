
#include <map>
inline auto gatef2( const std::vector<double>& temp, 
  std::vector<double>& eftmp2, int mat ){ 
  /* 
   *-------------------------------------------------------------------
   * Default effective temperatures for second atom in mixed
   * moderators from GA report:
   *       1095      benzine (c6h6)
   *       1099      beo
   *-------------------------------------------------------------------
   */

  std::vector<double> tabl1, tabl2, tabl3;
  tabl1 = { 1095, 1095, 1095, 1095, 1095, 1095, 1095, 1095, 1099, 1099, 1099, 
            1099, 1099, 1099, 1099, 1099 };
  tabl2 = { 296, 350, 400, 450, 500, 600, 800, 1000, 296, 400, 500, 600, 700, 
            800, 1000, 1200 };
  tabl3 = { 685.54, 712.02, 738.97, 768.10, 799.22, 866.63, 1017.3,  1182.3, 
            427.8, 502.8, 584.3, 671.3, 761.6, 854.2, 1043.7, 1236.6 };

  std::map<int, std::vector<double>> c { 
    { 1095, { 0.00000000000, 2.71841971e-4, 3.58624831e-1, 5.53282468e2 } },
    { 1099, { 0.00000000000, 1.32622103e-4, 7.04607881e-1, 2.02821917e2 } }
  };

  for ( size_t i = 0; i < temp.size(); ++i ){
    if (eftmp2[i] == 0) {
      for ( size_t j = 0; j < tabl2.size(); ++j ){
        if (tabl1[j] == mat and std::abs(tabl2[j] - temp[i]) <= 5 ) {
          eftmp2[i] = tabl3[j];
        } // end if 
      } // end for
      if (eftmp2[i] == 0) { eftmp2[i] = temp[i]; }
    } // end if 
  } // end for 
}




