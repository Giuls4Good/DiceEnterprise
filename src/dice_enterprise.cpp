// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
arma::mat construct_discrete_simplex(int d, int m) {
  //Construct a matrix with m columns where the members are
  //integers (possibly null) and the sum of the rows is d.
  arma::mat discrete_simplex = zeros<mat>(R::choose(d+m-1,m-1),m);
  int i;
  int counter = 0;
  if(m == 1) { //There is only one element
    discrete_simplex[0,0] = d;
    return(discrete_simplex);
  }
  for(i=0; i<=d; i++) {
    //We set the first element equal to i and all the other elements
    //we get them recursively
    arma::mat aux = construct_discrete_simplex(d-i,m-1);
    discrete_simplex.submat(counter,0,counter+aux.n_rows-1,0) = i*ones<mat>(aux.n_rows,1);
    discrete_simplex.submat(counter,1,counter+aux.n_rows-1,m-1) = aux;
    counter = counter+aux.n_rows;
  }
  return(discrete_simplex);
}

// [[Rcpp::export]]
double expected_tosses_bound_cpp(arma::mat A) {
  //Perform eigen decomposition
  arma::cx_vec eigval;
  arma::cx_mat eigvec;
  arma::eig_gen(eigval, eigvec, A);

  arma::mat Q = real(eigvec);
  arma::colvec D_vec = real(eigval);
  arma::mat Q_inv = inv(Q);

  //Compute expected number of tosses
  double expected_tosses = 0.0, new_expected_tosses = 1.0;
  int t = 1;
  while(new_expected_tosses - expected_tosses > 1e-4) {
    if(t > 1) {expected_tosses = new_expected_tosses;}
    new_expected_tosses = expected_tosses + accu(Q*diagmat(pow(D_vec,t))*Q_inv);
    t++;
  }
  return(new_expected_tosses);
}

// // NOT FINISHED
// List generate_ladder_initial_Cpp(List G_poly, int degree, int m, int v) {
//   //This function makes sure that:
//   //1) The polynomials are homogeneous (same degree)
//   //2) The coefficients of the polynomials are all positive
//   //3) All the polynomials have the same degree
//
//   //G_poly is a list of v elements. Each element is itself a list of 3 elements
//   //where the first elements are the coefficients R, the second element
//   //is a matrix of powers of the p_i, the 3rd element is the degree associated
//   //to each coefficients (it's just the sum of the rows of the second element)
//
//   int i, j, k;
//
//   for(i = 0; i<G_poly.length(); i++) {
//     colvec R = G_poly[i][0]; //Coefficients of the polynomial
//     mat M = G_poly[i][1]; //Powers of the polynomial
//     ivec degree_R = G_poly[i][2]; //Degree of each coefficients
//     for(j=0; j<R.length(); j++) {
//       //1 STEP: Convert to homogeneous polynomial
//       if(degree_R[j] < degree) { //The degree of the coefficient is smaller than the degree of the polynomial
//         //We multiply the coefficient by (p_1+...+p_m)^(degree-degree_R) to get
//         //new coefficients and new degrees
//         mat new_M = construct_discrete_simplex(degree-degree_R[j],m); //Construct the discrete simplex
//         colvec new_R = zeros<colvec>(new_M.n_rows); //This will store the new coefficients
//         for(k = 0; k<new_M.n_rows; k++) { //Discrete simplex is now the new M
//           new_R[k] = R[j] % multinomial_coeff(degree-degree_R[j],new_M.row(k));
//           new_M.row(k) = new_M.row(k) + M.row(j); //Update new_M to have the right powers
//         }
//       }
//
//     }
//
//   }
// }
