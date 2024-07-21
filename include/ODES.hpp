#include<armadillo>

#ifndef __ODES_hpp__  
#define __ODES_hpp__


arma::mat::fixed<1,5> ode(double t, arma::mat::fixed<1,5> y, double m);

arma::mat rk4(arma::mat::fixed<1,5> y0 ,arma::vec::fixed<2> tspan, int N, double m);


#endif