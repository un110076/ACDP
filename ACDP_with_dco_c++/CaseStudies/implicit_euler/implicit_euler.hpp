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

#pragma once

#include<cmath>
#include<vector>

template<typename A_FLT, typename P_FLT>
A_FLT g(const A_FLT &x, const A_FLT &p, const P_FLT &t) {
  return p*sin(x*t);
}

template<typename A_FLT, typename P_FLT>
A_FLT dgdx(const A_FLT &x, const A_FLT &p, const P_FLT &t) {
  return p*t*cos(x*t);
}

template<typename A_FLT, typename P_FLT>
A_FLT f(const A_FLT &x, const A_FLT &x_prev, const A_FLT &p, const int &i, const P_FLT &dt) {
  return x-x_prev-dt*g(x,p,i*dt);
}

template<typename A_FLT, typename P_FLT>
void newton(A_FLT &x, const A_FLT &x_prev, const A_FLT& p, const int &i, const P_FLT& dt, const P_FLT &eps) {
  A_FLT r=f(x,x_prev,p,i,dt);
  do {
    x=x-r/(1-dt*dgdx(x,p,i*dt));
    r=f(x,x_prev,p,i,dt);
  } while (fabs(r)>eps); 
}

template<typename A_FLT, typename P_FLT>
void implicit_euler(A_FLT &x, const std::vector<A_FLT> &p, const P_FLT &T, const P_FLT &eps) {
  int n=p.size();
  P_FLT dt=T/n;
  for (int i=0;i<n;i++) {
    A_FLT x_prev=x;
    newton(x,x_prev,p[i],i+1,dt,eps);
  }
}
