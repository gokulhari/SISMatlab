function varargout = sisSvdfr(n,varargin)
% SISVDFR
% Computes singular values of the frequency response operator in the system
% representation:
% A \phi = B d,
% \xi = C \phi,
%
% INPUTS:
%
% Need to specify adjoint boundary conditions explicitly:
%  [V,lam] = sisSvdfr(n,nad,N,A,B,C,bcReg,bcAdj)
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

if (length(varargin) ~= 9)
    error('Wrong number of inputs. Look up help sisSvfdfr.');
end
nad = varargin{1};
N = varargin{2};
A = varargin{3};
B = varargin{4};
C = varargin{5};
bcPointsReg = varargin{6};
bcOpReg = varargin{7};
bcPointsAd = varargin{8};
bcOpAd = varargin{9};
Aad = AdjointFormal(A);
Bad = AdjointFormal(B);
Cad = AdjointFormal(C);

BBad = MultOps(B,Bad);
CadC = MultOps(Cad,C);

A = Discretize(n,N,A);
Aad = Discretize(n,N,Aad);
BBad = Discretize(n,N,BBad);
CadC = Discretize(n,N,CadC);

bcsReg = BcMat(n,N,bcPointsReg,bcOpReg);
bcsAd = BcMat(n,N,bcPointsAd,bcOpAd);

bcsReg_null = null(bcsReg);
bcsAd_null = null(bcsAd);

A = A*bcsReg_null;
CadC = CadC*bcsReg_null;

Aad = Aad*bcsAd_null;
BBad = BBad*bcsAd_null;

Zs = zeros(size(Aad));

% Make feedback and calculate eigenvalues:
[V,gamma] = eig([Zs, BBad; CadC, Zs],[A, Zs; Zs, Aad]);
%Aad = AdjointFormal(varargin{3});
zs = zeros(size(bcsReg_null));
V = [bcsReg_null, zs; zs, bcsAd_null]*V;
% OrderMat = zeros(size(Aad));

% if (iscell(Aad{1,1})) % If block-matrix operator
%     [r,c] = size(Aad);
%     for i = 1:r
%         for j = 1:c
%             OrderMat(i,j) = length(Aad{i,j}) - 1;
%         end
%     end
% else % Else 
%     OrderMat = length(Aad) - 1;
% end
% nad = max(OrderMat);
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
gamma = diag(gamma);
[gamma, ind] = sort(gamma,'descend','ComparisonMethod','real');
V = V(:,ind);
if ( nargout < 2 )  % Return the eigenvalues only
    varargout = { gamma };
else
    varargout = {V,diag(gamma)};
end