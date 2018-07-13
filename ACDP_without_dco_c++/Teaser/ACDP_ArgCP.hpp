#pragma once

#include<vector>
#include<stack>

#include "ACDP_Base.hpp"

struct ACDP_ArgCP : ACDP_AdjointBase {
  std::stack<std::vector<double>> args;

  // default argument checkpoint 
  // consist of active input arguments
  virtual void push_args() {
    std::vector<double> cp;
    for (int i=0;i<input_count();i++)
      cp.push_back(input_value(i));
    args.push(cp);
  }

  virtual void top_args() {
    std::vector<double> cp=args.top();
    for (int i=0;i<input_count();i++)
      input_value(i)=cp[i];
  }

  virtual void pop_args() { args.pop(); }

};

