#include <Rcpp.h>
using namespace Rcpp;

//genes.raw <- apply(rawmat, 2, function(x)(sum(x>0)))
// [[Rcpp::export]]
double rcpp_add(NumericVector v){
  double add = 0;
  for(int i=0; i<v.length(); ++i){
    if(v[i]!=0){
      add += 1;
    }
  }
  return(add);
}

// [[Rcpp::export]]
NumericVector apply_cpp_col(NumericMatrix x ) {
  NumericVector output(x.ncol());
  for ( int i=0; i<x.ncol();i++){
    int temp=rcpp_add(x(_,i));
    if ( temp >=1){
      output[i]=temp;
    }
    else {
      output[i]=0;
    }
  }
  return output;
}

// der<- apply(rawmat,1,function(x)(sum(x>0)))/ncol(rawmat)
// [[Rcpp::export]]
NumericVector apply_cpp_row_der(NumericMatrix x ) {
  NumericVector output(x.nrow());
  for ( int i=0; i<x.nrow();i++){
    double temp=rcpp_add(x(i,_));
      output[i]=temp/x.ncol();
  }
  return output;
}


//norm.mat<- apply(norm.mat,2,function(x)(x <- x-mean(x)))

// [[Rcpp::export]]
double rcpp_sum(NumericVector v){
  double sum = 0;
  for(int i=0; i<v.length(); ++i){
      sum += v[i];
  }
  return(sum);
}

// [[Rcpp::export]]
DataFrame apply_cpp_col_norm(NumericMatrix x ) {
  int num_vars = x.ncol();
  List long_list(num_vars);
  for ( int i=0; i<x.ncol();i++){
    double temp=rcpp_sum(x(_,i));
    NumericVector t=rep(temp/x.rows(),x.rows());
    NumericVector temp_vec=NumericVector(x(_,i))-NumericVector(t);
    long_list[i]=temp_vec;
  }
  return long_list;
}

//basel <- apply(norm.mat.smooth[, which(colnames(norm.mat.smooth) %in% norm.cell.names)],1,median)
// [[Rcpp::export]]
double cpp_med(NumericVector x){
  std::size_t size = x.size();
  std::sort(x.begin(), x.end());
  if (size  % 2 == 0) return (x[size / 2 - 1] + x[size / 2]) / 2.0;
  return x[size / 2];
}

// [[Rcpp::export]]
NumericVector apply_cpp_row_basel(NumericMatrix x ) {
  NumericVector output(x.nrow());
  for ( int i=0; i<x.nrow();i++){
    double temp=cpp_med(x(i,_));
    output[i]=temp;
  }
  return output;
}




//DR2 <- apply(rawmat3,1,function(x)(sum(x>0)))/ncol(rawmat3)
// [[Rcpp::export]]
NumericVector apply_cpp_row_DR2(NumericMatrix x ) {
  NumericVector output(x.nrow());
  for ( int i=0; i<x.nrow();i++){
    double temp=rcpp_add(x(i,_))/x.ncol();
    output[i]=temp;
  }
  return output;
}


//temp=apply(uber.mat.adj[,which(com.pred=="diploid")], 1, mean)
// [[Rcpp::export]]
NumericVector apply_cpp_row_temp(NumericMatrix x ) {
  NumericVector output(x.nrow());
  for ( int i=0; i<x.nrow();i++){
    double temp=rcpp_sum(x(i,_))/x.ncol();
    output[i]=temp;
  }
  return output;
}


//cf.h <- apply(results.com.rat.norm, 1, sd)
// [[Rcpp::export]]
NumericVector apply_cpp_row_cfh(NumericMatrix x ) {
  NumericVector output(x.nrow());
  for ( int i=0; i<x.nrow();i++){
    double temp=sd(x(i,_));
    output[i]=temp;
  }
  return output;
}

//base <- apply(results.com.rat.norm, 1, mean)
// [[Rcpp::export]]
NumericVector apply_cpp_row_base(NumericMatrix x ) {
  NumericVector output(x.nrow());
  for ( int i=0; i<x.nrow();i++){
    double temp=mean(x(i,_));
    output[i]=temp;
  }
  return output;
}

//temp_a=apply(adj.results,2,mean))
// [[Rcpp::export]]
NumericVector apply_cpp_col_tempa(NumericMatrix x ) {
  NumericVector output(x.ncol());
  for ( int i=0; i<x.ncol();i++){
    double temp=mean(x(_,i));
    output[i]=temp;
  }
  return output;
}


//Aj <- apply(RNA[which(RNA.mat$hgnc_symbol %in% shr), ], 2, median)
// [[Rcpp::export]]
NumericVector apply_cpp_col_aj(NumericMatrix x ) {
  NumericVector output(x.ncol());
  for ( int i=0; i<x.ncol();i++){
    double temp=cpp_med(x(_,i));
    output[i]=temp;
  }
  return output;
}

