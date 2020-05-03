%function sol = sisSolves(n,N,A,bc,b)
function sol = sisSolves(n,varargin)
% SISSOLVES
% Solves for a forcing to a linear operator, uses sparse matrices, and
% sparse solver
%
% INPUTS:
%
% sol = sisSolves(n,N,A,bc,b)
% sol = sisSolves(n,N,A,bc,b,'Full');
%
% n -- A vector specifying the differential order of each variable in the
% block matrix operator. If simple linear operator, then n is scalar.
% N -- The number of Chebyshev polynomials in the approximation
% A -- Linear operator or block-matrix operator
% b -- input, in phys-space.
% bcs -- BCs class to specify boundary conditions. Look up help BCs
%
% Optional argument: 'Full' - uses the full matrix representation.
% 
% OUTPUT:
%
% The solution.

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
N = varargin{1};
A = varargin{2};
b = varargin{4};
bc = varargin{3};
if (length(varargin) == 5)
    full = true;
else
    full = false;
end

if (sum(n)~=(bc.r))
    error(['The system is ill-posed. Number of boundary conditions must' ...
     ' equal sum of the highest derivatives of all variables']);
end



if (iscell(b))
    iscellb = true;
    b = phys2cheb(b);
    [r,c] = size(b);
    temp = zeros((N+1)*r,c);
    for i = 1:r
        for j = 1:c
            temp((N+1)*(i-1) + 1: (N+1)*i,j) = b{i,j};
        end
    end
    b = temp;
else
    iscellb = false;
    b = ChebMat2CellMat(b,N);
    b = phys2cheb(b);
    [r,c] = size(b);
    temp = zeros((N+1)*r,c);
    for i = 1:r
        for j = 1:c
            temp((N+1)*(i-1) + 1: (N+1)*i,j) = b{i,j};
        end
    end
    b = temp;
end

[~,c] = size(b);
A = Discretize(n,N,A);
bcs = BcMat(n,N,bc.Points,bc.Operator);

if (~full)
    AA = sparse([A;bcs]);
    sol = AA\[b;repmat(bc.Values,[1,c])];
else
    sol = [A;bcs]\[b;repmat(bc.Values,[1,c])];
end

% Pad boundary conditions between spectral coefficients of highest
% derivatives
specCoef = sol(1:(N+1)*length(n),:);
intConsts = sol((N+1)*length(n) + 1 : end,:);

col_counter =  1;
cc = [1,cumsum(n)];
for i = 1:length(n)
    sol(col_counter: col_counter+(N+1) -1,:) = specCoef((N+1)*(i-1)+1:(N+1)*i,:); 
    sol(col_counter + N + 1 :col_counter + N + 1 + n(i) - 1,:) = intConsts(cc(i):cc(i+1),:);
    col_counter = col_counter + N+1 + n(i) - 1;
end
% Find the lowest derivatives:
Js = cell(length(n),1);
for i = 1:length(n)
    [~,~,J] = Matgen(n(i),N);
    Js{i} = J{end};
end
Jbig = blkdiag(Js{:});
sol = Jbig*sol;
sol = ChebMat2CellMat(sol,N);
sol = cheb2phys(sol);
if (~iscellb)
    [r,c] = size(sol);
    mat = zeros(r*(N+1),c);
    for i = 1:r
        for j = 1:c
            mat((N+1)*(i-1) + 1: (N+1)*i,j) = sol{i,j};
        end
    end
    sol = mat;
end