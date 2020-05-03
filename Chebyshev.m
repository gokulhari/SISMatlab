function R = Chebyshev(n,y)
% CHEBYSHEV
% Part of SISMatlab
% Compute n-th Chebyshev polynomial 
% n -- index that tells what polynomial we want to compute
% y -- grid

% SISMatlab is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% 
% SISMatlab is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with SISMatlab.  If not, see <https://www.gnu.org/licenses/>.
%
% Written by Mihailo Jovanovic


[row,col] = size(y);
T0 = ones(row,1);
T1 = y;

if (n == 0)
    R = T0;
  elseif(n == 1)
    R = T1;
  else 
    for i = 1:n-1
      T2 = 2*y.*T1 - T0;
      T0 = T1;
      T1 = T2;
    end
    R = T2;
end
