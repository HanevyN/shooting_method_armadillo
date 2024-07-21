#include "ODES.hpp"
#include<iostream>


arma::mat::fixed<1,5> ode(double t, arma::mat::fixed<1,5> y, double m ){
  
   arma::mat::fixed<1,5> sol;

   double const Pr = 0.72;

   sol(0,0) = y(1);
   sol(0,1) = y(2);
   sol(0,2) = (y(1)*y(1) - y(0)*y(2) + m/( (1 + m*y(3))*(1+m*y(3)) ) *y(2)*y(4)  )*(1 + m*y(3));
   sol(0,3) = y(4);
   sol(0,4) = -Pr*( y(0)*y(4));
   
    return sol;

}
arma::mat rk4(arma::mat::fixed<1,5> y0 ,arma::vec::fixed<2> tspan, int N,double m){

        // y0      -  initial conditions
        // tspan   - length of integrtion domain
        // N       - number of steps 
        

        arma::vec t = arma::linspace(tspan(0), tspan(1), N);

        arma::mat sol(N,5) ;
        arma::mat increment(1,5);
        // Set initial condition
        sol.row(0) = y0;

        arma::mat::fixed<1,5>  k1;
        arma::mat::fixed<1,5>  k2;
        arma::mat::fixed<1,5>  k3;
        arma::mat::fixed<1,5>  k4;

        
        double dt = (tspan[1] - tspan[0])/N;

        for(int i = 0 ; i < N-1; i++){            
            k1 = dt*ode(arma::as_scalar(t(i)),        sol.row(i)       ,m);
            k2 = dt*ode(arma::as_scalar(t(i)) + dt/2, sol.row(i) + k1/2,m);
            k3 = dt*ode(arma::as_scalar(t(i)) + dt/2, sol.row(i) + k2/2,m);
            k4 = dt*ode(arma::as_scalar(t(i+1)),      sol.row(i) + k3  ,m);
           
            sol.row(i+1) = sol.row(i) + (k1 + 2*k2 + 2*k3 + k4)/6;
    
        }
        
        return sol;


}






