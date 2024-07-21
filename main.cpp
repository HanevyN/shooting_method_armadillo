#include<iostream>
#include"ODES.hpp"
#include<cmath>
// To complie code use the command below
// g++ main.cpp src/*.cpp -I include -o main -std=c++11 -O2 -larmadillo]
// Run with ./main

int main(){

    arma::vec tspan = {0,40};
    // arma::mat y0 = {{0.0,1.0,-0.7296,1.0,-0.4982}};

    double m = -0.4;
    int N = 40*100;

    // parameters for secant shooting method
    double dp = 1e-6;
    double err = 1;
    arma::mat K(2,2);

    arma::vec ICS = {-0.7, -0.5};
    arma::vec H(2);

    while( err > 1e-10){ 

    arma::mat y1 = rk4(arma::mat{{0, 1 ,ICS(0) ,     1.0, ICS(1)     }}, tspan, N, m);
    arma::mat y2 = rk4(arma::mat{{0, 1 ,ICS(0) - dp ,1.0, ICS(1)     }}, tspan, N, m);
    arma::mat y3 = rk4(arma::mat{{0, 1 ,ICS(0) ,     1.0, ICS(1) - dp}}, tspan, N, m);
    K(0,0) = (y2(N-1,1) - y1(N-1,1))/dp;
    K(1,0) = (y2(N-1,3) - y1(N-1,3))/dp;
    K(0,1) = (y3(N-1,1) - y1(N-1,1))/dp;
    K(1,1) = (y3(N-1,3) - y1(N-1,3))/dp;


    arma::vec current_guess = {y1(N-1, 1), y1(N-1,3) };

    H = arma::conv_to< arma::vec >::from( K.i()*current_guess );

    ICS = ICS  + H;
    err = std::abs( arma::as_scalar(arma::max(H) ) );
    
    }
    // recompute solution
    arma::mat y1 = rk4(arma::mat{{0, 1 ,ICS(0) ,     1.0, ICS(1)     }}, tspan, N, m);

    std::cout << "Target boundary conditions " << std::endl;
    std::cout << "0 \t 0" << std::endl ;
    
    std::cout << "Computed boundary conditions "  << std::endl;
    std::cout <<  y1(N-1, 2) << "\t" << y1(N-1, 4) << std::endl ;


    


    return 0;
}