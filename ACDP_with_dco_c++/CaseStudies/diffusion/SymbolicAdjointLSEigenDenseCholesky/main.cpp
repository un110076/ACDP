/*
Adjoint Code Design Pattern with dco/c++
case study: diffusion
author: Uwe Naumann (2018)
*/

#include "ACDP_SymbolicAdjointLSEigenDenseCholesky.hpp"

#include "../diffusion.hpp"

inline void newton(
    DCO_TT* context_tape,
    const int& m, const VT<DCO_A>& y_prev, const double& yl, const double& yr,
    VT<DCO_A>& y
) {
  int n=y.size();
  const double eps=1e-15;
  MT<DCO_A> A(n,n); A.reserve(2*n-1);
  VT<DCO_A> r(VT<DCO_A>::Zero(n));
  f(m,y,yl,yr,y_prev,r);
  while (r.norm()>eps) {
    dfdy(m,y,yl,yr,A);
    ACDP_SymbolicAdjointLSEigenDenseCholesky::DMT<DCO_A> A_dense(A); 
    ACDP_SymbolicAdjointLSEigenDenseCholesky *e= new ACDP_SymbolicAdjointLSEigenDenseCholesky(context_tape,n);
    e->register_A(A_dense);
    e->register_b(r);
    context_tape->register_acdp(e);
    e->evaluate_augmented_primal();
    y-=r;
    f(m,y,yl,yr,y_prev,r);
  }
}

inline void euler(
    DCO_TT* context_tape,
    const int& m, const int& ncs, const double& yl, const double& yr,
    VT<DCO_A>& y
) {
  for (int j=0;j<m;j+=ncs) 
    for (int i=j;i<min(j+ncs,m);i++) {
      VT<DCO_A> y_prev=y;
      newton(context_tape,m,y_prev,yl,yr,y);
    }
}

int main(int, char* v[]){
  int n=stoi(v[1]), nsteps=stoi(v[2]), ncs=stoi(v[3]);
  VT<DCO_A> y(n), y_indep(VT<DCO_A>::Ones(n));
  double yl=0, yr=0;
  DCO_TT* context_tape=DCO_TT::create();
  for (int i=0;i<n;i++) context_tape->register_variable(y_indep(i));
  y=y_indep;
  euler(context_tape,nsteps,ncs,yl,yr,y);  
  context_tape->register_output_variable(y((n-1)/2));
  dco::derivative(y((n-1)/2))=1.;
  context_tape->interpret_adjoint();
  for(int i=0;i<n;i++)
    cout << dco::derivative(y_indep[i]) << endl;
  DCO_TT::remove(context_tape);
  return 0;
}
