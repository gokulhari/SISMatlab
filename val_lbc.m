% VAL_LBC
% Part of SISMatlab
% Evaluates the value of a function at y = -1. 
% Input function in phys-space.

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

function out = val_lbc(in)
% Evaluates the value of a function at lbc 
% input function is in phys-space
in = phys2cheb(in);
n = length(in)-1;
in(1) = in(1)/2.0;
for i = 1:2:n
    in (i+1) = -in(i+1);
end
out = sum(in);

