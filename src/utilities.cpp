#include <Rcpp.h>
#include <vector>
#include <algorithm>

using namespace std;
using namespace Rcpp;

#define STATE_START "start"
#define STATE_END "end"

//// version 2 of table function
////[[Rcpp::export]]
//NumericVector cTable(NumericVector vect);
//
//NumericVector cTable(NumericVector vect) {
//  int dataSize = vect.size();
//  cout << "Convert to C++ vector type." << endl;
//  vector<float> dataVect = Rcpp::as<std::vector<float> >(vect);
//  vector<float> vars;
//  vector<int> freqs;
//  int position;
//  float last;
//  
//  cout << "Sort vector" << endl;
//  sort (dataVect.begin(), dataVect.end());
//  
//  cout << "Start statistic..." << endl;
//  last = dataVect[0];
//  vars.push_back(last);
//  freqs.push_back(1);
//  for(int i = 1; i < dataSize; i++) {
//    if(last == dataVect[i]) {
//      position = vars.size() - 1;
//      freqs[position] = freqs[position] + 1;
//    } else {
//      last = dataVect[i];
//      vars.push_back(last);
//      freqs.push_back(1);
//    }
//  }
//  
//  cout << "Convert back to R numeric vector ..." << endl;
//  NumericVector rfreqs( freqs.begin(), freqs.end() );
//  rfreqs.names() = vars;
//  
//  return(rfreqs);
//}

// version 1 of table function
//' @name ctable
//' @title ctable is a C++ implementation to make contingency tables
//' @description ctable is a C++ implementation to make contingency tables
//' @param vect NumericVector
//' @return vars freqs
//[[Rcpp::export]]
List cTable(NumericVector vect);
IntegerVector cParIndicesSearch(NumericVector rawVector, 
                                NumericVector mz, 
                                IntegerVector mzsOrder, 
                                std::string startOrEnd);

List cTable(NumericVector vect) {
    int dataSize = vect.size();
    List cTable;
//  cout << "Convert to C++ vector type." << endl;
    vector<float> dataVect = Rcpp::as<std::vector<float> >(vect);
    vector<float> vars;
    vector<int> freqs;
    int position;
    float last;
  
//  cout << "Sort vector" << endl;
    sort (dataVect.begin(), dataVect.end());
  
//  cout << "Start statistic..." << endl;
    last = dataVect[0];
    vars.push_back(last);
    freqs.push_back(1);
    for(int i = 1; i < dataSize; i++) {
        if(last == dataVect[i]) {
            position = vars.size() - 1;
            freqs[position] = freqs[position] + 1;
        } else {
            last = dataVect[i];
            vars.push_back(last);
            freqs.push_back(1);
        }
    }
  
    cTable["vars"] = vars;
    cTable["freqs"] = freqs;
    return(cTable);
}

NumericVector rcppRev(NumericVector x) {
    NumericVector revX = clone<NumericVector>(x);
    std::reverse(revX.begin(), revX.end());
    ::Rf_copyMostAttrib(x, revX); 
    return revX;
}

//[[Rcpp::export]]
IntegerVector cParIndicesSearch(NumericVector rawVector, 
                                NumericVector mzs, 
                                IntegerVector mzsOrder, 
                                std::string startOrEnd) {
    int mzSize = mzs.size();
    int rawVectSize = rawVector.size();
    IntegerVector peakIndices(mzSize);
    int counter = 0;
    double raw;
  
    if(startOrEnd.compare(STATE_START) == 0) {
        for(int i = 0; i < rawVectSize; i++) {
            raw = rawVector[i];
            if(raw >= mzs[counter]) {
                peakIndices[mzsOrder[counter] - 1] = i + 1;
                counter++;
                if(counter == mzSize){
                    break;
                }
            }
        }
    } else if(startOrEnd.compare(STATE_END) == 0) {
        for(int i = rawVectSize; i > 0; i--) {
            raw = rawVector[i-1];
            if(raw <= mzs[counter]){
                peakIndices[mzsOrder[counter] - 1] = i;
                counter++;
                if(counter == mzSize){
                    break;
                }
            }
        }
    }
  
    return(peakIndices);
}













