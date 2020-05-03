function du = chebdiff(u)
% CHEBDIFF
% Part of SISMatlab
% The input has to be a vector
% The input has to be in phys-space, otherwise output will be garbage

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
 u = phys2cheb(u);
if ~isvector(u)
   disp('Error: Input has to be a vector');
   exit(1);
else
   n = length(u)-1;
   du = zeros(size(u));
   du(n) = 2.0*n*u(n+1);
   for i = 1:(n-1)/2
      du(n-(2*i+1)+1) = du(n-(2*i-1)+1) + 2.0*(n-2*i)*u(n-2*i+1);
   end
   du(n+1) = 0.0;

   for i = 1:(n-1)/2
      du(n-2*i+1) = du(n-2*(i-1)+1) + 2.0*(n-(2*i-1))*u(n+1-2*i+1);
   end
   
      
end

