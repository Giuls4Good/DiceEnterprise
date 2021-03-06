// [[Rcpp::depends(RcppArmadillo)]]
//#include <Rcpp.h>
#include <RcppArmadilloExtensions/sample.h>
using namespace Rcpp;
using namespace std;

// [[Rcpp::export]]
int findIntervalSingle(double U, NumericVector breaks) {
  int out;
  for(out = 0; out < breaks.length(); out++) {
    if(U < breaks[out]) {
      return(out);
    }
  }
  return(out);
}


// [[Rcpp::export]]
int updateFunCpp(int currentState, const arma::uvec& B,
                 const arma::vec& U, bool connected, bool fine,
                 const List& P_cumsum, const List& P_moves_list) {
  //Check if it is a valid ladder and a valid request
  if(!connected || !fine || B.n_elem!= U.n_elem) {
    stop("Not a valid call. Check if the ladder is fine and connected and that B and U have the same length.");
  }

  int nextState_index;
  List aux;
  List aux2;
  NumericVector coeff_cumsum;
  NumericVector nextState_possibilities;

  for(int c=0; c<B.n_elem; c++) {
    aux = P_cumsum[currentState-1];
    coeff_cumsum = aux[B[c]-1];
    if(coeff_cumsum.length() > 0) {
      //If there are no possible moves -> stays still
      nextState_index = findIntervalSingle(U[c], coeff_cumsum) + 1;
      aux2 = P_moves_list[currentState-1];
      nextState_possibilities = aux2[B[c]-1];
      if(nextState_index <= nextState_possibilities.length()) {
        //Otherwise U > coeff and it stays still
        currentState = nextState_possibilities[nextState_index-1];
      }
    }
  }

  return(currentState);
}

// [[Rcpp::export]]
arma::rowvec updateFunVecCpp(const arma::vec& states, const arma::uvec& B,
                             const arma::vec& U, bool connected, bool fine,
                          const List& P_cumsum, const List& P_moves_list,
                          const arma::vec& mapped_states, int k, int t_mapped_states,
                          int current_time) {
  // We reuse some previous results. mapped_states is a k-long vector where the first entry
  // tells where the first state is mapped at time 0 starting from time -t_mapped_states;
  // the second entry tells where the second state is mapped at time 0 starting from time -t_mapped_states
  // and so on. Some entries may not be available, in this case we need to compute it.

  //Check if it is a valid ladder and a valid request
  if(!connected || !fine || B.n_elem!= U.n_elem) {
    stop("Not a valid call. Check if the ladder is fine and connected and that B and U have the same length.");
  }

  arma::rowvec current_states(states.n_elem);
  int new_state;

  for(int j=0; j<states.n_elem; j++) {
    // Move state[j] from current_time to t_mapped_states
    arma::vec U1 = U.head(current_time-t_mapped_states);
    arma::uvec B1 = B.head(current_time-t_mapped_states);
    new_state = updateFunCpp(states[j], B1, U1, connected, fine, P_cumsum, P_moves_list);
    //Check if this state has already been mapped
    if(mapped_states[new_state-1] >= 1 && mapped_states[new_state-1] <= k) {
      //Rcout << "Already in! " << states[j] << " -> " << new_state << " -> " << mapped_states[new_state-1] << endl;
      new_state = mapped_states[new_state-1]; //-1 cause indices start from 0
    } else {
      //It hasn't been mapped yet
      arma::vec U2 = U.tail(t_mapped_states);
      arma::uvec B2 = B.tail(t_mapped_states);
      new_state = updateFunCpp(new_state, B2, U2, connected, fine, P_cumsum, P_moves_list);
      //Rcout << "Not there yet! " << states[j] << " -> " << new_state << endl;
    }
    current_states[j] = new_state;
  }

  return(current_states);
}

arma::uvec rollDieGivenProbs(int n, NumericVector probs) {
  //Roll a die given the probability of each face
  NumericVector n_range(probs.length());
  for(int i=0; i<probs.length(); i++) {
    n_range[i] = i+1;
  }
  arma::vec aux = RcppArmadillo::sample(n_range,n,true,probs);
  arma::uvec res = arma::conv_to<arma::uvec>::from(aux);
  return(res);
}

bool isUnique(arma::uvec x) {
  //Check if there is a unique element in x
  int el = x[0];
  for(int i=1; i<x.n_elem; i++) {
    if(x[i] != el) return(false);
  }
  return(true);
}

// [[Rcpp::export]]
List CFTPCpp(int k, NumericVector probs, bool connected, bool fine,
             List P_cumsum, List P_moves_list, bool monotonic, int min, int max, bool verbose) {
  //This function is for debug/performance use only and it is used when the probabilities of the given die
  //are known
  Rcout << "ciao";
  arma::uvec X(k);
  if(monotonic) {
    X.set_size(2);
    X[0] = min;
    X[1] = max;
  } else {
    //Populate with all starting values
    for(int i=0; i<k; i++) {
      X[i] = i+1;
    }
  }

  //Initialize as vector (to adapt size)
  arma::uvec B = rollDieGivenProbs(1, probs);
  arma::vec U = arma::randu<arma::vec>(1);
  arma::uvec B_roll;
  arma::vec U_roll;
  int T = 2;
  while(!isUnique(X)) {
    //Roll the die and the uniform
    B_roll.set_size(T/2);
    U_roll.set_size(T/2);
    B_roll = rollDieGivenProbs(T/2, probs);
    U_roll = arma::randu<arma::vec>(T/2);
    //Reshape B and U
    B.reshape(T,1);
    U.reshape(T,1);
    //Populate B and U with the new rolls
    for(int i = 0; i < B_roll.n_elem; i++) {
      B[T/2+i] = B_roll[i];
      U[T/2+i] = U_roll[i];
    }
    //Watch out: they need to be reversed to be used by the update function!
    if(monotonic) {
      //Start from the minimum and maximum state
      X[0] = updateFunCpp(min, arma::flipud(B), arma::flipud(U), connected, fine, P_cumsum, P_moves_list);
      X[1] = updateFunCpp(max, arma::flipud(B), arma::flipud(U), connected, fine, P_cumsum, P_moves_list);
    } else {
      //Start from all the states
      for(int i=1; i<=k; i++) {
        X[i-1] = updateFunCpp(i, arma::flipud(B), arma::flipud(U), connected, fine, P_cumsum, P_moves_list);
      }
    }

    T = 2*T;

  }

  List res;
  res["res"] = X[0];
  if(verbose) res["rolls"] = T;
  return(res);
}
