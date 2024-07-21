# shooting_method_armadillo

Code to solve a system of boundary layer ODEs which arise due to the linear strething of a sheet coupled 
to an energy equation via a temperature dependent viscosity function.
``` math
\displaylines{
   \frac{\partial u}{\partial x}  + \frac{\partial v}{\partial y} = 0,  \\
  u \frac{\partial u}{\partial x} + v \frac{\partial u}{\partial y} =  \frac{\partial}{\partial y} \left( \mu (T)   \frac{\partial u}{\partial y} \right),\\
    \frac{\partial^2 T}{\partial y^2}  - \text{Pr} v \frac{\partial T}{\partial y} = 0.
}
```
Subject to the boundary conditions
```math
  \displaylines{
  u - x = v = T - 1 = 0, \text{ at } y = 0, \\
  u = T = 0 , \text{ as } y \to \infty.
  }
```
The system is solved by making a similarity tranformation
```math
  v = -f(y) , \quad u = x f'(y), \quad T = T(y). 
```
and solving as a system of first order equations using a fourth order Rungga-Kutta method
see https://en.wikipedia.org/wiki/Runge%E2%80%93Kutta_methods.

The code was ran using Armadillo on ubuntu see https://arma.sourceforge.net/ 
Provided Armadillo is on your path (see the link above for details), the code can be 
compiled using
```
  	 g++ main.cpp src/*.cpp -I include -o main -std=c++11 -O2 -larmadillo
```
