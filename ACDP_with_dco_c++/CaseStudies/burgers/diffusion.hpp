/*
Adjoint Code Design Patterns with dco/c++
case study: diffusion
author: Uwe Naumann (2018)
*/

#pragma once

#include <vector>
#include <cmath>
using namespace std;

#include "Eigen/SparseLU"

template <typename T, int N=Eigen::Dynamic>
using VT=Eigen::Matrix<T,N,1>;

template <typename T>
using MT=Eigen::SparseMatrix<T>;

// rhs of ode
template <typename T, int N=Eigen::Dynamic, typename PT=double>
inline void g(const PT& d, const VT<T,N>& y, VT<T,N>& r) {
  int n=y.size();
  r[0]=y[0];
  for (int i=1;i<n-1;i++) {
    T diffusion = d*(n-1)*(n-1)*(y(i-1)-2*y(i)+y(i+1));
    T advection = -y(i)*(n-1);
    if (advection < 0) { advection *= y(i)-y(i-1); }
    else               { advection *= y(i+1)-y(i); }
    r[i] = diffusion + advection;
  }
  r[n-1]=y[n-1];
}

// tangent of rhs of ode
template <typename T, int N=Eigen::Dynamic, typename PT=double>
inline void g_t(const PT& d, const VT<T,N>& y, const VT<T,N>& y_t, VT<T,N>& r_t) {
  int n=y.size();
  r_t[0]=y_t[0];
  for (int i=1;i<n-1;i++)  {
    T diffusion_t = d*(n-1)*(n-1)*(y_t[i-1]-2*y_t[i]+y_t[i+1]);
    T advection = -y[i]*(n-1);
    T advection_t = -y_t[i]*(n-1);
    if (advection < 0) { advection_t = advection * (y_t[i]-y_t[i-1]) + advection_t * (y[i]-y[i-1]); }
    else               { advection_t = advection * (y_t[i+1]-y_t[i]) + advection_t * (y[i+1]-y[i]); }
    r_t[i] = diffusion_t + advection_t;
  }
  r_t[n-1]=y_t[n-1];
}

// Jacobian of rhs of ode
template <typename T, int N=Eigen::Dynamic, typename PT=double>
inline void dgdy(const PT& d, const VT<T,N>& y, MT<T>& A) {
  int n=y.size();
  VT<T,N> y_t(n),r_t(n);
  std::vector<Eigen::Triplet<T>> entries;
  entries.reserve(3*n-2);
  for (int i=0;i<n;i++) {
    y_t=VT<T,N>::Unit(n,i);
    g_t(d,y,y_t,r_t);
    for (int j=max(i-1,0);j<min(i+2,n);j++) 
      entries.push_back(Eigen::Triplet<T>(j,i,r_t(j)));
  }
  A.setFromTriplets(entries.begin(),entries.end());
}

// residual of nls
template <typename T, int N=Eigen::Dynamic, typename PT=double>
inline void f(const int m, const PT& d,
    const VT<T,N>& y_prev, const VT<T,N>& y, VT<T,N>& r) {
  g(d,y,r);
  r=y-y_prev-r/m;
}

// Jacobian of residual of nls wrt. state
template <typename T, int N=Eigen::Dynamic, typename PT=double>
inline void dfdy(const int& m, const PT& d, const VT<T,N>& y, MT<T>& A) {
  int n=y.size();
  dgdy(d,y,A);
  A/=-m; 
  MT<T> B(n,n); B.setIdentity(); A+=B;
}

// Newton solver for nls
template <typename T, int N=Eigen::Dynamic, typename PT=double>
inline void newton(
    const int& m, const PT& d, const VT<T,N>& y_prev, VT<T,N>& y
) {
  int n=y.size();
  const double eps=1e-15;
  MT<T> A(n,n); A.reserve(3*n-2);
  VT<T,N> r(VT<T,N>::Zero(n));
  f(m,d,y_prev,y,r);
  Eigen::SparseLU<MT<T>> solver;
  while (r.norm()>eps) {
    dfdy(m,d,y,A);
    solver.analyzePattern(A);
    solver.factorize(A);
    r=solver.solve(r);
    y-=r;
    f(m,d,y_prev,y,r);
  }
}

// implicit Euler integration
template <typename T, int N=Eigen::Dynamic, typename PT=double>
inline void euler( const int& m, const PT& d, VT<T,N>& y) {
  for (int j=0;j<m;j++) {
    VT<T,N> y_prev=y;
    newton(m,d,y_prev,y);
  }
}
