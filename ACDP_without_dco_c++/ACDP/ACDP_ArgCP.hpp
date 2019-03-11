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

