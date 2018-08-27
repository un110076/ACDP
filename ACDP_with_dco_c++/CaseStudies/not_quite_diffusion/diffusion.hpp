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

// diffusivity as function of time, space and state
template <typename T>
inline T c(const T& y) { return sin(y); }

// derivative of diffusivity wrt. state
template <typename T>
inline T dcdy(const T& y) { return cos(y); }

template <typename T, int N=Eigen::Dynamic>
using VT=Eigen::Matrix<T,N,1>;

template <typename T>
using MT=Eigen::SparseMatrix<T>;

// rhs of ode
template <typename T, int N=Eigen::Dynamic, typename PT=double>
inline void g(
    const VT<T,N>& y, const PT& yl, const PT& yr,
    VT<T,N>& r
) {
  int n=y.size();
  double dt_rec=(n+1)*(n+1);
  r(0)=c(y(0))*dt_rec*(yl-2*y(0)+y(1));
  for (int i=1;i<n-1;i++)
    r(i)=c(y(i))*dt_rec*(y(i-1)-2*y(i)+y(i+1));
  r(n-1)=c(y(n-1))*dt_rec*(y(n-2)-2*y(n-1)+yr);
}

// tangent of rhs of ode
template <typename T, int N=Eigen::Dynamic, typename PT=double>
inline void g_t(
    const VT<T,N>& y, const VT<T,N>& y_t, const PT& yl, const PT& yr,
    VT<T,N>& r_t
) {
  int n=y.size(), dt_rec=(n+1)*(n+1);
  r_t(0)=dcdy(y(0))*dt_rec*(yl-2*y(0)+y(1))*y_t(0)
        +c(y(0))*dt_rec*(-2*y_t(0)+y_t(1));
  for (int i=1;i<n-1;i++)
    r_t(i)=dcdy(y(i))*dt_rec*(y(i-1)-2*y(i)+y(i+1))*y_t(i)
          +c(y(i))*dt_rec*(y_t(i-1)-2*y_t(i)+y_t(i+1));
  r_t(n-1)=dcdy(y(n-1))*dt_rec*(y(n-2)-2*y(n-1)+yr)*y_t(n-1)
        +c(y(n-1))*dt_rec*(y_t(n-2)-2*y_t(n-1));
}

// Jacobian of rhs of ode
template <typename T, int N=Eigen::Dynamic, typename PT=double>
inline void dgdy(
    const VT<T,N>& y, const PT& yl, const PT& yr,
    MT<T>& A
) {
  int n=y.size();
  VT<T,N> y_t(n),r_t(n);
  std::vector<Eigen::Triplet<T>> entries;
  entries.reserve(3*n-2);
  for (int i=0;i<n;i++) {
    y_t=VT<T,N>::Unit(n,i);
    g_t(y,y_t,yl,yr,r_t);
    for (int j=max(i-1,0);j<min(i+2,n);j++)
      entries.push_back(Eigen::Triplet<T>(j,i,r_t(j)));
  }
  A.setFromTriplets(entries.begin(),entries.end());
}

// residual of nls
template <typename T, int N=Eigen::Dynamic, typename PT=double>
inline void f(
    const int& m, const VT<T,N>& y, const  PT& yl, const PT& yr, const VT<T,N>& y_prev,
    VT<T,N>& r
) {
  g(y,yl,yr,r);
  r=y-y_prev-r/m;
}

// Jacobian of residual of nls wrt. state
template <typename T, int N=Eigen::Dynamic, typename PT=double>
inline void dfdy(
    const int& m, const VT<T,N>& y, const PT& yl, const PT& yr,
    MT<T>& A
) {
  int n=y.size();
  dgdy(y,yl,yr,A);
  A/=-m; 
  MT<T> B(n,n); B.setIdentity(); A+=B;
}

// Newton solver for nls
template <typename T, int N=Eigen::Dynamic, typename PT=double>
inline void newton(
    const int& m, const VT<T,N>& y_prev, const PT& yl, const PT& yr,
    VT<T,N>& y
) {
  int n=y.size();
  const double eps=1e-15;
  MT<T> A(n,n); A.reserve(3*n-2);
  VT<T,N> r(VT<T,N>::Zero(n));
  f(m,y,yl,yr,y_prev,r);
  Eigen::SparseLU<MT<T>> solver;
  while (r.norm()>eps) {
    dfdy(m,y,yl,yr,A);
    solver.analyzePattern(A);
    solver.factorize(A);
    r=solver.solve(r);
    y-=r;
    f(m,y,yl,yr,y_prev,r);
  }
}

// implicit Euler integration
template <typename T, int N=Eigen::Dynamic, typename PT=double>
inline void euler(
    const int& m, const int& ncs, const PT& yl, const PT& yr,
    VT<T,N>& y
) {
  for (int j=0;j<m;j+=ncs) 
    for (int i=j;i<min(j+ncs,m);i++) {
      VT<T,N> y_prev=y;
      newton(m,y_prev,yl,yr,y);
    }
}
