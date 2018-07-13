/*
Adjoint Code Design Pattern with dco/c++
pattern: SymbolicAdjointLSEigenDenseCholesky
author: Uwe Naumann (2018)
*/

#pragma once

#include "Eigen/Cholesky"
#include "ACDP_Base.hpp"

struct ACDP_SymbolicAdjointLSEigenDenseCholesky : ACDP_AdjointBase {

  template <typename T, int N=Eigen::Dynamic>
  using VT=Eigen::Matrix<T,N,1>;

  template <typename T, int N=Eigen::Dynamic>
  using DMT=Eigen::Matrix<T,N,N>;

  Eigen::LLT<DMT<DCO_B>> solver;

  int n;
  DMT<DCO_B> A;
  VT<DCO_B> x;

  ACDP_SymbolicAdjointLSEigenDenseCholesky(DCO_TT *owning_tape, const int &n) : ACDP_AdjointBase(owning_tape), n(n), A(DMT<DCO_B>::Zero(n,n)), x(VT<DCO_B>::Zero(n)) { }

  void register_A(DMT<DCO_A>& A_in) { 
    for (int i=0;i<n;i++) 
      for (int j=0;j<n;j++) { 
        register_input(A_in(i,j)); 
        A(i,j)=dco::value(A_in(i,j));
      }
  }
  void register_b(VT<DCO_A>& x_in) { 
    for (int i=0;i<n;i++) { 
      register_input(x_in(i)); 
      register_output(x_in(i)); 
      x(i)=dco::value(x_in(i)); 
    }
  }

  void evaluate_augmented_primal() {
    solver.compute(A);
    x=solver.solve(x);
    for (int i=0;i<n;i++) output_value(i)=x(i);
  }

  void evaluate_adjoint() {
    VT<DCO_B> x_a=VT<DCO_B>::Zero(n);
    for (int i=0;i<n;i++) { x_a(i)=output_adjoint(i); output_adjoint(i)=0; }
    x_a=solver.solve(x_a);
    for (int i=0;i<n;i++)
      input_adjoint(n*n+i)+=x_a(i);
    for (int j=0;j<n;j++) 
      for (int i=0;i<n;i++)
        input_adjoint(j*n+i)-=x_a(j)*x(i);
  }
};

