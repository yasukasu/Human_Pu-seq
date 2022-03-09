#include <Rcpp.h>
using namespace Rcpp;


// [[Rcpp::export]]


NumericVector rcpp_moving_ave(NumericVector data, int n){
  
  if(data.length()<n)stop("THe number of averaging bins is too many!");
  
  int bin_n = data.length();
  NumericVector r = rep(0.0 ,bin_n);
  
  for(int i=0;i<bin_n;i++){ 
    
    //IntegerVector st = IntegerVector::create(0,(i-n));
    //IntegerVector en = IntegerVector::create(data.length()-1, i+n);
    if(data[i]==NA || data[i]==NAN || data[i]==R_PosInf || data[i]==R_NegInf){
      r[i]=NA_REAL;
    }
    else {
      NumericVector d_data = data[seq(max(IntegerVector::create(0,(i-n))), min(IntegerVector::create(data.length()-1, i+n)))];
      d_data = d_data[!is_infinite(d_data)]; // exluding "Inf"
      d_data = d_data[!is_nan(d_data)]; // exluding "NaN"
      d_data = d_data[!is_na(d_data)];
      r[i]=mean(d_data);
      
    }
    //std::cout << mean(d_data);
  }
  return(r);
}
  
  
// NumericVector rcpp_moving_ave_v2(NumericVector data, int n){
// 
//     if(data.length()<n)stop("THe number of averaging bins is too many!");
// 
//     int dataN = data.length();
//     NumericVector r = rep(0.0 , dataN);
//   
//     //IntegerVector s (0,n);
//     //IntegerVector e (dataN-1,n);
//     int s = 0;
//     int e = 1;
//     IntegerVector start = IntegerVector::create(s, seq(0,(dataN-1-n)));
//     IntegerVector end   = IntegerVector::create(seq(n,dataN-1), e);
//     IntegerMatrix se = cbind(start, end);
// 
//    
//     for(int i=0;i<(dataN-1);i++)
//     {
//       IntegerVector vecs = seq(start[i], end[i])  ;
//       //r[i] = mean(data[vecs]);
//     }
// 
//     return(r);
// }

  