/*
Uwe Naumann hereby disclaims all copyright interest in the
Adjoint Code Design Patterns (ACDP) software.

Uwe Naumann, Aachen, Germany, 11 March 2019
*/

/*
ACDP is free software: you can redistribute it and/or modify
it under the terms of the GNU Lesser General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

ACDP is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License
along with this program.  If not, see <https://www.gnu.org/licenses/>.
*/

#pragma once

#include<vector>

// association of values and adjoints by address
typedef std::pair<double,double> AT;

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

