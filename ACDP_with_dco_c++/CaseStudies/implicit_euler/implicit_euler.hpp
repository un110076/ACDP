/*
Adjoint Code Design Patterns with dco/c++
case study: implicit Euler
author: Uwe Naumann (2018)
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
  while (fabs(r)>eps) {
    x=x-r/(1-dt*dgdx(x,p,i*dt));
    r=f(x,x_prev,p,i,dt);
  }
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
