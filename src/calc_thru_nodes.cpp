#include <vector>


#include <Rcpp.h>
using namespace Rcpp;



Rcpp::List calThruNodes_cpp(ances
                            NumericMatrix states
                               loglik
                               forTime
                            List parameter
                               methode
                               phy) {

  NumericMatrix lambdas <- parameter[1];
  NumericMatrix mus <- parameter[2]
  NumericMatrix Q = parameter[3]

  /* R:
   nb_node <- phy$Nnode
   reltol <- 1e-12
   abstol <- 1e-16
   ly <- ncol(states)
   d <- ncol(states) / 2
   focal <- ances
   desRows <- which(phy$edge[, 1] == focal)
   desNodes <- phy$edge[desRows, 2]
   nodeM <- numeric()
   nodeN <- numeric()
  */
  double reltol = 1e-12;
  double abstol = 1e-16;
  double ly = states.ncol();
  double d = stats.ncol() * 0.5;
  NumericVector focal = ances;

  NumericVector desRows = 0; // which(phy$edge[, 1] == focal)
  NumericVector desNodes = 0; // desNodes <- phy$edge[desRows, 2]




}