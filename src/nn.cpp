#include <Rcpp.h>
using namespace Rcpp;

//' @name nnMean
//' @title nnMean is C++ code for calculating nearest neighbour means in a 
//' 2D matrix
//' @description nnMean is C++ code for calculating nearest neighbour means in a 
//' 2D matrix
//' @param y NumericVector
//' @param nrows int
//' @param ncols int
//' @return eY 
// [[Rcpp::export]]
NumericVector nnMean(NumericVector y, int nrows, int ncols)
{
  //printf("nearestNeighbourMean\n");

  int i,j,ij,im1j,ijp1,ip1j,ijm1,num,nr=nrows,nc=ncols;
  double tot;
  NumericVector eY(nr * nc);
  // assume matrix is surrounded by rows and columns of NAs
  // so the indices here can run from 2 to nr-1 etc

  for(i=2;i<=(nr-1);i++) {
    for(j=2;j<=(nc-1);j++) {
      
      //printf("%d %d \n",i,j);

      ij = (j-1)*nr + (i-1);   //i,j
      im1j = (j-1)*nr + (i-2); //i-1,j
      ijp1 = (j)*nr + (i-1);   //i,j+1
      ip1j = (j-1)*nr + (i);   //i+1,j
      ijm1 = (j-2)*nr + (i-1);  //i,j-1
       
      tot = 0;
      num = 0;
      if(!ISNA(y[im1j])) {
	tot = tot + y[im1j];
	num++;
      }
      if(!ISNA(y[ijp1])) {
	tot = tot + y[ijp1];
	num++;
      }
      if(!ISNA(y[ip1j])) {
	tot = tot + y[ip1j];
	num++;
      }
      if(!ISNA(y[ijm1])) {
	tot = tot + y[ijm1];
	num++;
      }
      
      //printf("%f %d\n",tot,num);

      eY[ij] = (tot / num);
      
      //printf("%f\n",eY[ij]);

    }
  }
  return(eY);

}

      
