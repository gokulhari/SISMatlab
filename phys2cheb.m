function out = phys2cheb(in)
% PHYS2CHEB
% Converts from values in physical space to Chebyshev coefficients. 
% Input can be 
% (a) a simple vector
% (b) an m by n cell of vectors, 
% where each vector is transformed to Chebyshev coefficients

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
if ~iscell(in)
    n = length(in)-1;
    out = sqrt(2/(n+1))*dct(in);
    out(1) = out(1)*sqrt(2);
else
    [r,c] = size(in);
    out = cell(r,c);
    for i = 1:r
        for j = 1:c
            out{i,j} = phys2cheb(in{i,j});
        end
    end
end