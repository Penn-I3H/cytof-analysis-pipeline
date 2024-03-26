#include "RcppArmadillo.h"
using namespace Rcpp;

// utility function: get intersection length
// for two sorted vectors A, B of length k
int get_intersection_length(arma::irowvec& A, arma::irowvec& B, int k) {
  int i1=0, i2=0, val1, val2;
  int l=0;
  while (i1 < k && i2 < k) {
    val1 = A(i1);
    val2 = B(i2);

    if (val1 == val2)
      ++l;

    if (val1<=val2)
      ++i1;

    if (val1>=val2)
      ++i2;
  }

  return l;
}


// [[Rcpp::export]]
List get_edgelist(arma::imat& nbrs) {
  int n = nbrs.n_rows, k = nbrs.n_cols, edge_cnt=0;
  arma::imat edgelist = arma::imat(n*k,2);
  arma::vec jac = arma::vec(n*k);
  arma::irowvec A,B;
  std::unordered_set<double> uniq;

  for (int i=0; i<n; i++) {
    A = sort(nbrs.row(i));

    for (int j=0; j<k; j++) {
      int nbr=A(j);
      double idx;

      if (i<nbr-1)
        idx = (i+1)*(n+0.0)+(nbr);
      else
        idx = nbr*(n+0.0)+(i+1);

      if (uniq.find(idx) != uniq.end())
        continue;
      uniq.insert(idx);

      B = sort(nbrs.row(nbr-1));
      int l = get_intersection_length(A, B, k);

      if (l==0)
        continue;

      edgelist(edge_cnt,0) = i+1;
      edgelist(edge_cnt,1) = nbr;
      jac(edge_cnt) = (double)l / (2*k-l);
      edge_cnt++;
    }
  }

  edgelist.resize(edge_cnt,2);
  jac.resize(edge_cnt);
  jac=jac/2;

  List L = List::create(Named("edgelist") = edgelist , _["jac"] = jac);
  return L;
}


// [[Rcpp::export]]
NumericMatrix bhatt(arma::mat& kde_mat, int d, int b) {
  int k = kde_mat.n_cols;
  arma::vec vec1, vec2, vec_all, vec_ch;
  NumericMatrix bhatt_prod(k,k);
  std::fill(bhatt_prod.begin(), bhatt_prod.end(), 1);

  for (int i=0; i<k-1; i++) {
    vec1 = kde_mat.col(i);

    for (int j=i+1; j<k; j++) {
      vec2 = kde_mat.col(j);

      vec_all = sqrt(vec1 % vec2);

      for (int l=0; l<d; l++) {
        vec_ch = vec_all.subvec(l*b, (l+1)*b-1);
        bhatt_prod(i,j) *= sum(vec_ch);
      }
    }
  }

  return bhatt_prod;
}

