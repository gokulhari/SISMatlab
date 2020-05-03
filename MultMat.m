function out = MultMat(func)
% MultMat
% Part of SISMatlab
% Produces respective matrix for nonconstant coefficients
%
% INPUT:
% func: A function in physical space.

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
N = length(func) - 1;
teop = zeros(N + 1);
hank = zeros(N + 1);
func = phys2cheb(func);
for i = 0:N
      teop(i+ 1, i+ 1 : N + 1) = func(1: N + 1 - i).';
end

teop = teop + teop.';
teop(1:N+2:end) = diag(teop)/2;
for i = 1:N
      hank(i+1,1:N - i) = func(i + 1:N).';
end
out = 0.5*(teop + hank);
out = out.';
