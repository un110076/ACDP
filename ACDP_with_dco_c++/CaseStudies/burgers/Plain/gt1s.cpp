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
  typedef dco::gt1s<double>::type DCO_T;
  if (c!=3) throw;
  int n=stoi(v[1]), m=stoi(v[2]);
  for(int i=0;i<n;i++) { 
    VT<DCO_T> y(VT<DCO_T>::Zero(n));
    for (int i=1;i<n-1;i++) y[i]=sin((2*pi*i)/n);
    dco::derivative(y[i])=1;
    euler(m,d,y);
    cout << dco::derivative(y[(n-1)/2]) << endl;
  }
  return 0;
}
