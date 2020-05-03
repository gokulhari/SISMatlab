function out = cheb_mul(a,b)
% CHEB_MUL
% Part of SISMatlab
% a and b input functions in cheb space. a and b are taken to physical space,
% multiplied there, and brought back to cheb space

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
a = cheb2phys(a);
b = cheb2phys(b);
out = a.*b;
out = phys2cheb(out);