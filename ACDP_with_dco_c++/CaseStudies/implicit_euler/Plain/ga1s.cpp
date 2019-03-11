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

#include "dco.hpp"
typedef dco::ga1sm<double> DCO_AM;
typedef DCO_AM::type DCO_A;
typedef DCO_AM::tape_t DCO_TT;

#include "implicit_euler.hpp"

int main(int, char* v[]) {
  int n=std::stoi(v[1]);
  const double T=1, eps=1e-15;            
  std::vector<DCO_A> p(n,1); DCO_A x=1;  
  DCO_TT* tape=DCO_TT::create();
  tape->register_variable(p);
  implicit_euler(x,p,T,eps);
  dco::derivative(x)=1;
  tape->interpret_adjoint();
  for (int i=0;i<n;i++)
    std::cout << "dx/dp[" << i << "]="  
              << dco::derivative(p[i]) << std::endl;
  DCO_TT::remove(tape);
  return 0;
}
