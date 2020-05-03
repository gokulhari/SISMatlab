function out = ChebMat2CellMat(in,N)
% CHEBMAT2CELLMAT
% Part of SISMatlab
% Converts matrix of Chebyshev basis of size m N x n to cells of size m x n
% INPUTS:
% V -- The matrix of size m N x n
% N -- Number of Chebyshev basis functions
% 
% OUTPUT:
% V -- Cell of size m x n

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
[r,c] = size(in);
if (rem(r,N+1) ~= 0)
    error('Wrong dimensions');
end

rows = r/(N +1);
cols = c;

out = cell(rows,cols);

for i = 1:rows
    for j = 1:cols
        out{i,j} = in((i-1)*(N + 1) + 1 : i*(N + 1),j);
    end
end

