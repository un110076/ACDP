/*
Adjoint Code Design Pattern with dco/c++
pattern: EarlyTangentPreaccumulation
author: Uwe Naumann (2018)
*/

#pragma once

#include "../Base/ACDP_Base.hpp"

struct ACDP_EarlyTangentPreaccumulation : ACDP_AdjointBase {
  // passive context
  int n,m;
  // local Jacobian
  std::vector<DCO_B> jac;
  // target
  std::unique_ptr<ACDP_PrimalBase<DCO_T>> tgt;
  // target i/o
  std::vector<DCO_T> x, y; 

  ACDP_EarlyTangentPreaccumulation(
    DCO_TT *context_tape, 
    ACDP_PrimalBase<DCO_T>* tgt, // EarlyTangentPreaccumulation_Target
    int n, // number of active inputs
    int m // number of active outputs
  ) : ACDP_AdjointBase(context_tape), n(n), m(m), 
      jac(std::vector<DCO_B>(m*n,0)), tgt(tgt), x(std::vector<DCO_T>(n,0)), y(std::vector<DCO_T>(m,0)) {}

  void link_target() {
    if (!tgt->linked) {
      for (auto& v : x) tgt->register_input(v);
      for (auto& v : y) tgt->register_output(v);
      tgt->linked=true;
    }
  }

  void evaluate_augmented_primal() {
    for (int i=0;i<n;i++) {
      for (int j=0;j<n;j++) x[j]=input_value(j);
      dco::derivative(x[i])=1;
      tgt->evaluate_primal();
      for (int j=0;j<m;j++) 
        jac[j*n+i]=dco::derivative(y[j]);
    }
    for (int j=0;j<m;j++) 
      output_value(j)=dco::value(y[j]);
  }

  void evaluate_adjoint() {
    std::vector<DCO_B> tmp(m);
    for (int j=0;j<m;j++) {
      tmp[j]=output_adjoint(j); output_adjoint(j)=0;   
    }
    for (int i=0;i<n;i++) 
      for (int j=0;j<m;j++) {
        input_adjoint(i)+=jac[j*n+i]*tmp[j];
      }
  }

};
