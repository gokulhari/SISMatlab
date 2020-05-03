function [Df,C,D] = Matgen(n,N)
% Matgen
% Part of SISMatlab
% This function generates discretization matrices, and constants of
% integration matrices needed for spectral integration.
%
% INPUTS: 
% n -- The highest order of all variables, and 
% N -- no. of basis functions
%
% OUTPUT: The Discretization matrices as cells: like D{1} etc for both
% constants of integration and the coeficients.

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

D = cell(n+1,1);
Df = cell(n+1,1);
C = cell(n,1);
for i = 1:n+1
    Df{i,1} = zeros(N + 1 + 2 * n, N + 1 + 2 * n);
    D{i,1} = zeros(N + 1, N + 1 + n);
end

for i = 1:n
    C{i} = zeros(n,n);
end

% Now fill them in:

Df{1}(1:N+1, 1:N+1) = eye(N + 1);
for i = 1 : n
      Df{i+1}(1,:) = 0.5 * Df{i}(2,:);
      for  j = 1 : (N + 2 * n - 1)
        Df{i+1}(j+1, :) = (Df{i}(j,:) - Df{i}(j+2,:)) / (2.0*(j));
      end
end
 
if (n > 0) 
    C = cell(n,1);
    for i = 0: n -1 
      C{i + 1} = zeros(n, n);
    end
    % Next we set the matrices for the constants
    C{n} = eye(n, n);
    for  i = n - 2:-1:0
      for  j = 0: n-1
        C{i+1}(:,j+1) = chebdiff(C{i+2}(:,j+1));
      end
    end
end

for i = 0:n
      D{i+1}(1:N+1,1:N+1) = Df{i+1}(1: N + 1, 1:N + 1);
end
for i = 0:n-1
      D{i+2}(1:n, N + 2: N + 1 + n) = C{i+1};
end
end

    
function du = chebdiff(u)
% The input has to be a vector
% The input has to be in phys-space, otherwise output will be garbage
%u = phys2cheb(u);
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
end
