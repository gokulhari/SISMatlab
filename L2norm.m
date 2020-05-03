function out = L2norm(f)
% L2NORM
% Part of SISMatlab
% Computes L2norm of a function or a cell of functions in phys-space

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
if ~iscell(f)
    f = f.*f;
    f = integ(f);
    %f = phys2cheb(f);
    out = val_rbc(f) - val_lbc(f);
else 
    out = zeros(size(f));
    [r,c] = size(f);
    for i = 1:r
        for j = 1:c
            out(i,j) = L2norm(f{i,j});
        end
    end
end