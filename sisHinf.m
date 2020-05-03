function [omega_opt,gamma_opt] = sisHinf(n,varargin)
% SISHINF
% Computes H-infinity or L-infinity norm
%
% INPUTS:
%
% n -- Orders of variables as an array
% N -- Number of basis functions
% E -- E in the statespace
% F -- F in the statespace
% B -- B in the statespace
% C -- C in the statespace
% kerM: Null-space that system has to lie in, if none, supply Identity.
%
% or,
% n -- Orders of variables as an array
% N -- Number of basis functions
% E -- E in the statespace
% F -- F in the statespace
% B -- B in the statespace
% C -- C in the statespace
% bcReg -- Boundary conditions on regular variable, specified using class
% BCs(), look up "help BCs"
% bcAdj -- Boundary conditions on adjoint variable, specified using class
% BCs(), look up "help BCs"

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



% Precision of the H-infinity norm, adjust if desired:
Hinf_eps = 1e-4;
% iteration counter
count = 1;

if (length(varargin) == 7)
    N = varargin{1};
    E = varargin{2};
    F = varargin{3};
    B = varargin{4};
    C = varargin{5};
    bcReg = varargin{6};
    bcAdj = varargin{7};
    % Compute the null-space:
    bcs_reg = BcMat(n,N,bcReg.Points,bcReg.Operator);
    bcs_ad = BcMat(n,N,bcAdj.Points,bcAdj.Operator);
    zs = zeros(size(bcs_reg));
    M_null = null([bcs_reg,zs; zs, bcs_ad]);
elseif (length(varargin) == 6)
    N = varargin{1};
    E = varargin{2};
    F = varargin{3};
    B = varargin{4};
    C = varargin{5};
    M_null = varargin{6};
else
    error('Wrong number of inputs, look up help sisSvdfr')
end
[~,~,J] = Matgen(2,N);
Ebig = [J{3},zeros(size(J{3})),zeros(size(J{3})),zeros(size(J{3}));...
    zeros(size(J{3})),J{3},zeros(size(J{3})),zeros(size(J{3}));...
    zeros(size(J{3})),zeros(size(J{3})),J{3},zeros(size(J{3}));...
    zeros(size(J{3})),zeros(size(J{3})),zeros(size(J{3})),J{3}];
Ebig = [Ebig,zeros(size(Ebig));...
        zeros(size(Ebig)),Ebig];

Ead = AdjointFormal(E);
Fad = AdjointFormal(F);
Bad = AdjointFormal(B);
Cad = AdjointFormal(C);
BBad = MultOps(B,Bad);
CadC = MultOps(Cad,C);

E = Discretize(n,N,E);
Ead = Discretize(n,N,Ead);
F = Discretize(n,N,F);
Fad = Discretize(n,N,Fad);
CadC = Discretize(n,N,CadC);
BBad = Discretize(n,N,BBad);

% Estimate lower bound at omega = 0
omega_opt = 0;
A2  = -F;
A2s = -Fad;

M = [A2, zeros(size(A2));
zeros(size(A2)), A2s];
L = [zeros(size(BBad)),BBad;
    CadC,zeros(size(BBad))];


L = L*M_null;
M = M*M_null;

[~,evals] = eigs(M,L,1,'SM');
evals = diag(evals);
gamma_lb = abs(1/min(abs(evals)));

while (true)
      disp(['In iteration ' int2str(count)] );
      count = count + 1;
      gamma = (1.0 + 2.0 * Hinf_eps) * gamma_lb;
      M = [E, zeros(size(E));
            zeros(size(E)), Ead];
      L = [F,BBad/(gamma);
          -CadC/gamma,-Fad];
      M = M*M_null;
      L = L*M_null;

      [V,evals] = eig(L,M);
      
      % Unpivot and go to the lowest derivative
      V = Ebig*M_null*V;
      V = ChebMat2CellMat(V,N);
      % Keep only eigenvalues that have converged to machine precision
      [~,evals] = keepConverged(V,evals,N);
      
      % find purely imaginary evals
      evals = diag(evals);
      omegas = evals(abs(real(evals)) < 1e-11)/1i;
      if (isempty(omegas))
          gamma_ub = gamma;
          break;
      else
        ms = (omegas(1:length(omegas)-1) + omegas(2:end))/2;
        gammas = zeros(size(ms));
      
        for i = 1:length(ms)
          A2  =  ((1i*ms(i)*E) - F);
          A2s = ((-1i*ms(i)*Ead) - Fad);
          
          M = [A2, zeros(size(A2));
            zeros(size(A2)), A2s];
          L = [zeros(size(BBad)),BBad;
                CadC,zeros(size(BBad))];
          L = L*M_null;
          M = M*M_null;
          
        [~,evals] = eigs(M,L,1,'SM');
        evals = diag(evals);
        gammas(i) = abs(1/min(abs(evals)));
        end
        [gamma_lb,in] = max(gammas);
        omega_opt = ms(in);
      end
end
gamma_opt = (gamma_lb + gamma_ub) / 2.0; 