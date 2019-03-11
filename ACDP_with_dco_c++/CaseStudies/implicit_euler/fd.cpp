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

#include <iostream>
#include <vector>
#include<cfloat>

#include "implicit_euler.hpp"

int main(int, char* v[]) {
  int n=std::stoi(v[1]);
  const double x0=1, T=1, eps=1e-15;            
  std::vector<double> p(n,1);
  double x=x0; 
  implicit_euler(x,p,T,eps);
  for (int i=0;i<n;i++) {
    double h=sqrt(DBL_EPSILON);
    double xp=x0;
    p[i]+=h;
    implicit_euler(xp,p,T,eps);
    std::cout << "dx/dp[" << i << "]=" << (xp-x)/h << std::endl;
    p[i]-=h;
  }
  return 0;
}

