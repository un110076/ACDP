#include<vector>
#include<stack>
#include<cmath>
#include<iostream>
#include<memory>
using namespace std;

#include "ACDP_ArgCP.hpp"

struct ACDP_LateRecording : ACDP_AdjointBase {
  unique_ptr<ACDP_ArgCP> tgt;
  ACDP_LateRecording(ACDP_ArgCP* tgt) : tgt(tgt) {}

  void link_target() {
    if (!tgt->linked) {
      for (int i=0;i<input_count();i++)
        tgt->register_input(input(i));
      for (int i=0;i<output_count();i++)
        tgt->register_output(output(i));
      tgt->linked=true;
    }
  }

  void evaluate_primal() {
    tgt->evaluate_primal();
  }
  void evaluate_augmented_primal() {
    tgt->push_args();
    tgt->evaluate_primal();
  }
  void evaluate_adjoint() {
    tgt->top_args(); tgt->pop_args();
    tgt->evaluate_augmented_primal();
    tgt->evaluate_adjoint();
  }
};

