function out = BcMat(n,varargin)
% BCMAT
% Part of SISMatlab
% Generates Boundary condition matrices for spectral integration
%
% mat = BcMat(n,N,eval,L); 
% mat = BcMat(n,N,bcs);
%
% INPUTS for simple operators:
%
% n -- Overall order of Differential system
% N -- number of basis functions
% eval -- where to evaluate
% L -- cell representing operator acting at the boundary.
% bcs -- Boundary conditions via BCs class
%
% Inputs for block matrix operators:
%
% n -- vector of orders of variables
% N -- number of basis functions
% eval -- A matrix with evaluation points
% L -- a block matrix operator, cell of cells, with operators corresponding
% to eval.
% bcs -- Boundary conditions via BCs class

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


if (length(varargin) == 3) 
    N = varargin{1};
    eval = varargin{2};
    L = varargin{3};
elseif (length(varargin) == 2)
    N = varargin{1};
    eval = varargin{2}.Points;
    L = varargin{3}.Operator;
end

if (~iscell(L{1,1}))
    [Df,C,D] = Matgen(n,N);
    out = zeros(1, N + 1 + n);
    vals = zeros(1, N + 1 + 2 * n);
    valcons = zeros(1,n);
    for k = 0:N + 2 * n
          vals(k+1) = cos( k * acos(eval));
    end
    for k = 0:n-1
          valcons(k+1) = cos(k * acos(eval));
    end

    vals(1) = 0.5*vals(1);
    if n > 0
       valcons(1) = 0.5*valcons(1);
    end
    %vals
    %valcons
    diffn = n - length(L) + 1;

    for k = 0:length(L) - 1
          if (isscalar(L{k+1}))
            matemp = L{k+1} * vals* Df{k + diffn + 1};
            out(1:N+1) = out(1:N+1) + matemp(1:N + 1);
            if (k + diffn - 1 > -1) 
              out(N + 2: N + 1+ n) = out(N + 2: N + 1+ n) + ...
                  (L{k+1} * valcons * C{k + diffn});
            end
          else
            matemp = cheb_eval(L{k+1},eval) * vals * Df{k + diffn + 1};
            out(1:N + 1) = out(1:N + 1) + matemp(1:N + 1);
            if (k + diffn - 1 > -1) 
              out(N + 2:N+1+n) = out(N + 2:N+1+n) + (cheb_eval(L{k+1},eval) * ...
                                            valcons * C{k + diffn});
            end
          end
    end
else
    [r,c] = size(L);
    out = zeros(r,c*(N+1) + sum(n));
    for i = 1:r
        for j = 1:c
            out(i, ...
                (j-1)*(N+1)+1 + sum(n(1:j-1)) : j*(N+1)+ sum(n(1:j))) = BcMat(n(j),N,eval(i,j),L{i,j});
        end
    end
end

