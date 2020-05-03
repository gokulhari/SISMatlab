function varargout = sisEigs(n,varargin)
% SISEIGS
% Computes eigenvalues of block matrix operators, uses the sparse 
% eigenvalue solver, eigs by default. 
%
% INPUTS:
%
%  [V,lam] = sisEigs(n,N,F,E,bc)
%  [V,lam] = sisEigs(n,N,F,E,bc,K)
%  [V,lam] = sisEigs(n,N,F,E,bc,K,SIGMA)
%
% n -- A vector specifying the differential order of each variable in the
% block matrix operator. If simple linear operator, then n is scalar.
% N -- The number of Chebyshev polynomials in the approximation
% F -- The linear (block-matrix) operator in generalized eigenvalue problem
% with (F,E)
% E -- The linear (block-matrix) operator in generalized eigenvalue problem
% with (F,E)
% bc -- Boundary conditions using the class BCs(), look up "help BCs"
%
% Optional inputs that mimic eigs() in Matlab:
%
% K -- By default it computes the 6 eigenvalues with largest real part, 
% use K to specify number of eigenvalues
%
% SIGMA -- finds K eigenvalues. If SIGMA is a scalar, the eigenvalues 
% found are the ones closest to SIGMA. Similar to Matlab's eigs, use 
% 'LR' and 'SR' for the eigenvalues of largest and smallest
% real part, and 'LM' and 'SM' for largest and smallest magnitude.
% Use 'Full' to use eig() instead of eigs(), in this special case,
% K is not used and all singular values from the finite-dimensional 
% approximation are the output.
%
% OUTPUTS:
%
% V -- Eigenvectors as cells of eigenfunction.
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


N = varargin{1};
F = varargin{2};
E = varargin{3};
bc = varargin{4};
if (numel(varargin) == 4)
    numval = 6;
    sigflag = true;
elseif (numel(varargin) == 5)
    if (varargin{5} > (N+1)*length(n))
        error(['Cannot compute more eigenvalues than the' ...
            'finite-dimensional approximation']);
    end
    numval = varargin{6};
    sigflag = true;  
elseif (numel(varargin)>5 && numel(varargin)<7)
    if (varargin{5} > (N+1)*length(n))
        error(['Cannot compute more eigenvalues than the' ...
            'finite-dimensional approximation']);
    end
    numval = varargin{5};
    sigflag = false;
    sigma = varargin{6};
else
    error('Wrong number of arguments for svdfr');
end



Fmat = Discretize(n,N,F);
Emat = Discretize(n,N,E);
bcs = BcMat(n,N,bc.Points,bc.Operator);
bcs_null = null(bcs);
% Project 
Fmat = Fmat*bcs_null;
Emat = Emat*bcs_null;
% Make feedback and calculate eigenvalues:
if (sigflag)
    [V,lam] = eigs(Fmat,Emat,numval,'LR');
else
    if strcmp(sigma,'Full')
        [V,lam] = eig(Fmat,Emat);
    else
        [V,lam] = eigs(Fmat,Emat,numval,sigma);
    end
end
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

if (strcmp(sigma,'Full'))
    lam = diag(lam);
    [lam, ind] = sort(lam,'descend','ComparisonMethod','real');
    V = V(:,ind);
    lam = diag(lam);
end
if ( nargout < 2 )  % Return the eigenvalues only
    varargout = { diag(lam) };
else
    varargout = {V,lam};
end
