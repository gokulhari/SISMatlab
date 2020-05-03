function varargout = sisSvdfrs(n,varargin)
% SISVDFRS
% Computes singular values of the frequency response operator in the system
% representation:
% A \phi = B d,
% \xi = C \phi,
% uses Matlab's EIGS by default. 
%
% INPUTS:
%
% Specify adjoint boundary conditions explicitly:
%  [V,lam] = sisSvdfrs(n,nad,N,A,B,C,bcReg,bcAdj);
%  [V,lam] = sisSvdfrs(n,nad,N,A,B,C,bcReg,bcAdj,K);
%  [V,lam] = sisSvdfrs(n,nad,N,A,B,C,bcReg,bcAdj,K,SIGMA);
%
% n -- A vector specifying the differential order of each variable in the
% block matrix operator. If simple linear operator, then n is scalar.
% nad -- A vector specifying the differential order of each adjiont 
% variable in the block matrix operator. If simple linear operator, then 
% nad is scalar.
% N -- The number of Chebyshev polynomials in the approximation
% A -- The linear (block-matrix) operator, A in the system representation
% B -- The linear (block-matrix) operator, B in the system representation
% C -- The linear (block-matrix) operator, C in the system representation
% bcReg -- Boundary conditions on regular variable, specified using class
% BCs(), look up "help BCs"
% bcAdj -- Boundary conditions on adjoint variable, specified using class
% BCs(), look up "help BCs"
%
% Optional inputs that mimic eigs() in Matlab:
%
% K -- By default it computes the 6 largest singular values, use K to
% specify number of singular values
%
% SIGMA -- finds K singular values. If SIGMA is a scalar, the singular 
% values found are the ones closest to SIGMA. Similar to Matlab's eigs, use 
% 'LR' and 'SR' for the eigenvalues of largest and smallest
% real part, and 'LM' and 'SM' for largest and smallest magnitude.
% Use 'Full' to use eig() instead of eigs(), in this special case,
% K is not used and all singular values from the finite-dimensional 
% approximation are the output.
%
% OUTPUTS:
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

nad = varargin{1};
N = varargin{2};
A = varargin{3};
B = varargin{4};
C = varargin{5};
bcReg = varargin{6};
bcAdj = varargin{7};

if (numel(varargin) == 7)
    numval = 6;
    sigflag = true;
elseif (numel(varargin) == 8)
    if (varargin{8} > (N+1)*length(n))
        error(['Cannot compute more eigenvalues than the' ...
            'finite-dimensional approximation']);
    end
    numval = varargin{8};
    sigflag = true;  
elseif (numel(varargin)>8 && numel(varargin)<10)
    if (varargin{8} > (N+1)*length(n))
        error(['Cannot compute more eigenvalues than the' ...
            'finite-dimensional approximation']);
    end
    numval = varargin{8};
    sigflag = false;
    sigma = varargin{9};
else
    error('Wrong number of arguments for svdfr');
end

Aad = AdjointFormal(A);
Bad = AdjointFormal(B);
Cad = AdjointFormal(C);

BBad = MultOps(B,Bad);
CadC = MultOps(Cad,C);

A = Discretize(n,N,A);
Aad = Discretize(n,N,Aad);
BBad = Discretize(n,N,BBad);
CadC = Discretize(n,N,CadC);

bcsReg = BcMat(n,N,bcReg.Points,bcReg.Operator);
bcsAdj = BcMat(n,N,bcAdj.Points,bcAdj.Operator);

bcsReg_null = null(bcsReg);
bcsAdj_null = null(bcsAdj);

A = A*bcsReg_null;
CadC = CadC*bcsReg_null;

Aad = Aad*bcsAdj_null;
BBad = BBad*bcsAdj_null;

Zs = zeros(size(Aad));


% Make feedback and calculate eigenvalues:
if (sigflag)
    [V,gamma] = eigs([Zs, BBad; CadC, Zs],[A, Zs; Zs, Aad],numval,'LR');
else
    if strcmp(sigma,'Full')
        [V,gamma] = eig([Zs, BBad; CadC, Zs],[A, Zs; Zs, Aad]);
    else
        [V,gamma] = eigs([Zs, BBad; CadC, Zs],[A, Zs; Zs, Aad],numval,sigma);
    end
end

zs = zeros(size(bcsReg_null));
V = [bcsReg_null, zs; zs, bcsAdj_null]*V;

% Integrate to get to lower derivatives:
Js = cell(length(n)+length(nad),1);
for i = 1:length(n)
    [~,~,J] = Matgen(n(i),N);
    Js{i} = J{end};
end
for i = 1:length(nad)
    [~,~,J] = Matgen(n(i),N);
    Js{length(n) + i} = J{end};
end
Jbig = blkdiag(Js{:});


V = Jbig*V;
V = ChebMat2CellMat(V,N);
V = cheb2phys(V);
if (~sigflag)
    if (strcmp(sigma,'Full'))
        gamma = diag(gamma);
        [gamma, ind] = sort(gamma,'descend','ComparisonMethod','real');
        V = V(:,ind);
        gamma = diag(gamma);
    end
end
if ( nargout < 2 )  % Return the eigenvalues only
    varargout = { diag(gamma) };
else
    varargout = {V,gamma};
end