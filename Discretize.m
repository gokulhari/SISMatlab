% Discretize: 
% Part of SISMatlab
% This produces respective operators for linear operators using
% spectral integration.
% 
% INPUTS: 
% For simple operator:
% n: Highest order of the variable in the linear block-matrix operator.
% N: Number of basis functions
% L: A linear operator in the block-matrix operator using a cell (therefore
% can be constant or spatially varying coefficients).
%
% For block matrix operators:
% n: vector of size number of columns of L, storing highest order of each
% variable.
% N: Number of basis functions
% L: The block matrix operator as cell of cells

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
% Written by Gokul Hariharan
function out = Discretize(n,N,L)

if (~iscell(L{1,1})) % If simple operator
    [~,~,D] = Matgen(n,N);
    out = zeros(N + 1, N + 1 + n);
    diffn = n - length(L) + 1;

    for k = 0:length(L) - 1
          if (isscalar(L{k+1}))
              out = out + L{k + 1} * D{k + diffn + 1}(1:N+1,1:N + 1 + n);
          else 
              out = out + MultMat(L{k + 1})*D{k + diffn + 1}(1:N+1,1:N + 1 + n);
          end
    end
else % use only if all internal operators have same order
    [r,c] = size(L);
    out = zeros(r*(N+1) ,c*(N+1) + sum(n));
    for i = 1:r
        for j = 1:c
            out((i-1)*(N+1)+1 : i*(N+1) , ...
                (j-1)*(N+1)+1 + sum(n(1:j-1)) : j*(N+1)+ sum(n(1:j))) = Discretize(n(j),N,L{i,j});
        end
    end
end
      
          
          
