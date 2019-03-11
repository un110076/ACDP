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

#include <iostream>
#include <vector>

#include "burgers.hpp"

int main(int c, char* v[]){
  if (c!=3) throw;
  int n=std::stoi(v[1]), m=std::stoi(v[2]);
  VT<double> y(VT<double>::Zero(n)); 
  for (int i=1;i<n-1;i++) y[i]=sin((2*pi*i)/n);
  euler(m,d,y);
  std::cout << y << std::endl;
  return 0;
}

