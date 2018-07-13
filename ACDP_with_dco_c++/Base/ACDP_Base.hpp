/*
Adjoint Code Design Pattern with dco/c++
pattern: Context
author: Uwe Naumann (2018)
*/

#pragma once
#include "dco.hpp"

typedef double DCO_B;
typedef dco::gt1s<DCO_B>::type DCO_T;
typedef dco::ga1sm<DCO_B> DCO_AM;
typedef DCO_AM::type DCO_A;
typedef DCO_AM::tape_t DCO_TT;
typedef DCO_TT::position_t DCO_TPT;

template<typename T>
struct ACDP_PrimalBase : dco::ACDP::PrimalBase<T> {
  bool linked;
  ACDP_PrimalBase() : linked(false) {}
};

struct ACDP_AdjointBase : dco::ACDP::AdjointBase<DCO_TT,DCO_A> {
  bool linked;
  ACDP_AdjointBase(DCO_TT* tape) 
    : dco::ACDP::AdjointBase<DCO_TT,DCO_A>(tape), linked(false) {}
};
