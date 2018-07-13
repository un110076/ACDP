/*
Adjoint Code Design Pattern with dco/c++
pattern: ArgCP 
author: Uwe Naumann (2018)
*/

#pragma once
#include "../Base/ACDP_Base.hpp"
#include<stack>

struct ACDP_ArgCP : ACDP_AdjointBase {
  std::stack<std::vector<DCO_B>> args;
  ACDP_ArgCP(DCO_TT* t) : ACDP_AdjointBase(t) {}

  // standard stack treatment assumes checkpoint to
  // consist of active input arguments
  virtual void push_args() {
    std::vector<DCO_B> cp;
    for (unsigned int i=0;i<input_count();i++)
      cp.push_back(input_value(i));
    args.push(cp);
  }

  virtual void top_args() {
    std::vector<DCO_B> cp=args.top();
    for (unsigned int i=0;i<input_count();i++)
      input_value(i)=cp[i];
  }

  virtual void pop_args() { args.pop(); }
};
