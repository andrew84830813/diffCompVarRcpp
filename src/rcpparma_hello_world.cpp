// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

// we only include RcppArmadillo.h which pulls Rcpp.h in for us
#include "RcppArmadillo.h"

// via the depends attribute we tell Rcpp to create hooks for
// RcppArmadillo so that the build process will know what to do
//
// [[Rcpp::depends(RcppArmadillo)]]


// ks function
double KS(arma::colvec x, arma::colvec y) {
    int n = 100;
    arma::colvec cc = join_cols(x, y);
    double minx = cc.min();
    double maxx = cc.max();
    arma::colvec r = arma::linspace(minx, maxx, n);
    int N = r.n_rows;
    
    arma::colvec fx = arma::linspace(minx, maxx, n);
    arma::colvec v = arma::sort(x);
    arma::colvec s = v;
    for (int i=0; i<n;i++) {
        double X = r(i);
        s.fill(0.0); s.elem( find(v < X) ).ones();
        fx(i) = mean(s);
    }
    
    arma::colvec fx1 = arma::linspace(minx, maxx, n);
    arma::colvec v1 = arma::sort(y);
    arma::colvec s1 = v1;
    for (int i=0; i<n;i++) {
        double X = r(i);
        s1.fill(0.0); s1.elem( find(v1 < X) ).ones();
        fx1(i) = mean(s1);
    }
    
    arma::colvec diff = arma::abs(fx-fx1);
    return arma::max(diff);
}


//' Computes columnwise K-S statistic
//'
//' @param mt n-sample by p-feature matrix for class1
//' @param mt2 n-sample by p-feature matrix for class2
//' @return p-vector of K-S statistics
// [[Rcpp::export]]
Rcpp::NumericVector K_S(arma::mat mt,arma::mat mt2) {
    int n = mt.n_cols;
    Rcpp::NumericVector results(n);
    for (int i=0; i<n;i++) {
        arma::colvec x=mt.col(i);
        arma::colvec y=mt2.col(i);
        results[i] = KS(x, y);
    }
    return results;
}

//' Computes column median using rcpp
//'
//' @param x n-sample by p-feature matrix
//' @return p-vector of columnwise medians
// [[Rcpp::export]]
Rcpp::NumericVector column_median(arma::mat x) {
    int  ncol =  x.n_cols;
    Rcpp::NumericVector out(ncol);
    
    arma::mat B   = arma::median(x);
    out = B  ;
    return out;
}

//' Computes column variance based on median
//'
//' @param x n-sample by p-feature matrix
//' @param m p-vector of column medians to compute variance from
//' @return p-vecor of columnwise variance based on supplied median value
// [[Rcpp::export]]
Rcpp::NumericVector column_medianVar(arma::mat x, arma::mat m) {
    int  ncol =  x.n_cols, nrow = x.n_rows;
    Rcpp::NumericVector out(ncol);
    
    for (int i=0; i<ncol;i++) {
        double total = 0;
        for (int j=0; j<nrow;j++) {
            total += ( x(j, i)- m(0,i) ) * ( x(j,i)- m(0,i) );
        }
        
        out[i] = total;
    }
    return out;
}



