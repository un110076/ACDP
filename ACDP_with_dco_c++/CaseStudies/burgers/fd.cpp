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
#include <cfloat>
#include <cmath>
using namespace std;

#include "burgers.hpp"

int main(int c, char* v[]){
  if (c!=3) throw;
  int n=stoi(v[1]), m=stoi(v[2]);
  VT<double> y(VT<double>::Zero(n));
  for (int i=1;i<n-1;i++) y[i]=sin((2*pi*i)/n);
  for(int j=0;j<n;j++) {
    VT<double> yph(y), ymh(y);
    double h=(y[j]<10) ? sqrt(DBL_EPSILON) : sqrt(DBL_EPSILON)*abs(y[j]); 
    ymh[j]-=h;
    euler(m,d,ymh);
    yph[j]+=h;
    euler(m,d,yph);
    cout << j << " " << (yph[(n-1)/2]-ymh[(n-1)/2])/(2*h) << endl;
  }
  return 0;
}
