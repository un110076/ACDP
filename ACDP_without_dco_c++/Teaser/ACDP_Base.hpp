#pragma once

#include<vector>

// association of values and adjoints by address
typedef pair<double,double> AT;

struct ACDP_AdjointBase {
  bool linked;
  std::vector<AT*> inputs, outputs;
  ACDP_AdjointBase() : linked(false) {}
  void register_input(AT& x) { inputs.push_back(&x); }
  void register_output(AT& x) { outputs.push_back(&x); }
  AT& input(int i) { return *(inputs[i]); }
  double& input_value(int i) { return inputs[i]->first; }
  double& input_adjoint(int i) { return inputs[i]->second; }
  double& output_value(int i) { return outputs[i]->first; }
  double& output_adjoint(int i) { return outputs[i]->second; }
  AT& output(int i) { return *(outputs[i]); }
  int input_count() { return inputs.size(); }
  int output_count() { return outputs.size(); }
  virtual void evaluate_primal()=0;
  virtual void evaluate_augmented_primal()=0;
  virtual void evaluate_adjoint()=0;
};

