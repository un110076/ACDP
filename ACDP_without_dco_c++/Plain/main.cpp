/*
Uwe Naumann hereby disclaims all copyright interest in the
Adjoint Code Design Patterns (ACDP) software.

Uwe Naumann, Aachen, Germany, 11 March 2019
*/

/*
ACDP is free software: you can redistribute it and/or modify
it under the terms of the GNU Lesser General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

ACDP is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License
along with this program.  If not, see <https://www.gnu.org/licenses/>.
*/

#include<stack>
#include<cmath>
#include<iostream>
#include<memory>

#include "ACDP_LateRecording.hpp"

void primal_euler(double &x, const int nsteps, const double& dt) {
  for (int i=0;i<nsteps;i++) x+=dt*sin(x);
}

void augmented_primal_euler(double &x, const int nsteps, const double& dt, std::stack<double>& tbr) {
  for (int j=0;j<nsteps;j++) {
    tbr.push(1+dt*cos(x)); // recording partial derivatives
    x+=dt*sin(x);
  }
}

void adjoint_euler(double &x_a, const int nsteps, std::stack<double>& tbr) {
  for (int j=nsteps;j>=1;j--) {
    x_a*=tbr.top(); tbr.pop();
  }
}

int main(int c, char* v[]) {
  if (c!=2) throw;
  int nsteps=std::stoi(v[1]);
  AT x(AT(1,1));
  double dt=1.0/(3*nsteps);
  std::stack<double> tbr; 

  // augmented primal
  augmented_primal_euler(x.first,3*nsteps,dt,tbr);
  std::cout << "x=" <<  x.first << std::endl;

  // adjoint
  adjoint_euler(x.second,3*nsteps,tbr);
  std::cout << "x_a=" <<  x.second << std::endl;
  return 0;
}
