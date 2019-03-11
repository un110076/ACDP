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

#include "burgers.hpp"

int main(int c, char* v[]){
  typedef dco::ga1sm<double> DCO_AM;
  typedef DCO_AM::type DCO_A;
  typedef DCO_AM::tape_t DCO_TT;
  if (c!=3) throw;
  int n=stoi(v[1]), m=stoi(v[2]);
  VT<DCO_A> y(n), y_indep(VT<DCO_A>::Zero(n));
  for (int i=1;i<n-1;i++) y_indep[i]=sin((2*pi*i)/n);
  DCO_TT* tape=DCO_TT::create();
  for (int i=0;i<n;i++) tape->register_variable(y_indep[i]);
  y=y_indep;
  euler(m,d,y);  
  tape->register_output_variable(y[(n-1)/2]);
  dco::derivative(y[(n-1)/2])=1.;
  tape->interpret_adjoint();
  for(int i=0;i<n;i++)
    cout << dco::derivative(y_indep[i]) << endl;
  DCO_TT::remove(tape);
  return 0;
}

