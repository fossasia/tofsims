#include <RcppArmadillo.h>
using namespace arma;
using namespace Rcpp;

// Below is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp 
// function (or via the Source button on the editor toolbar)

//' @name EigenDecompose
//' @title EigenDecompose for the MNF analysis
//' @description EigenDecompose for the MNF analysis
//' @param A NumericMatrix
//' @param B NumericMatric
//' @param startIndex int
//' @param endIndex int
//' @return eigval eigvec mA mB
// [[Rcpp::export]]
List EigenDecompose(NumericMatrix A, 
                    NumericMatrix B, 
                    int startIndex, 
                    int endIndex){
  List eig;
  cx_vec eigval;
  cx_mat eigvec;
  // Convert Rcpp matrix to RcppArmadillo matrix
  int nA = A.nrow(), kA = A.ncol(), nB = B.nrow(), kB = B.ncol();
  mat mA(A.begin(), nA, kA, false);
  mat mB(B.begin(), nB, kB, false);
  
  eig_pair( eigval, eigvec, mA, mB );
  eig["eigval"] = eigval;
  eig["eigvec"] = eigvec;
  return eig;
}






