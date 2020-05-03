function varargout = sisEig(n,N,F,E,bc)
% SISEIG
% Computes eigenvalues of block matrix operators. 
%
% INPUTS:
%
%  [V,lam] = sisEig(n,N,F,E,bc)
% n -- A vector specifying the differential order of each variable in the
% block matrix operator. If simple linear operator, then n is scalar.
% N -- The number of Chebyshev polynomials in the approximation
% F -- The linear (block-matrix) operator in generalized eigenvalue problem
% with (F,E)
% E -- The linear (block-matrix) operator in generalized eigenvalue problem
% with (F,E)
% bc -- Boundary conditions using the class BCs(), look up "help BCs"
%
% OUTPUTS:
%
% V -- Eigenvectors as cells of eigenfunctions.
% lam -- Eigenvalues

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

Fmat = Discretize(n,N,F);
Emat = Discretize(n,N,E);
bcs = BcMat(n,N,bc.Points,bc.Op);
bcs_null = null(bcs);
% Project 
Fmat = Fmat*bcs_null;
Emat = Emat*bcs_null;
[V,lam] = eig(Fmat,Emat);
V = bcs_null*V;
% Find the lowest derivatives:
Js = cell(length(n),1);
for i = 1:length(n)
    [~,~,J] = Matgen(n(i),N);
    Js{i} = J{end};
end
Jbig = blkdiag(Js{:});
V = Jbig*V;
V = ChebMat2CellMat(V,N);
V = cheb2phys(V);
lam = diag(lam);
[lam, ind] = sort(lam,'descend');
V = V(:,ind);
if ( nargout < 2 )  % Return the eigenvalues only
    varargout = { lam };
else
    varargout = {V,diag(lam)};
end
