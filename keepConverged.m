function [V,evals] = keepConverged(V,evals,N)
% KEEPCONVERGED 
% This function only keeps those eigenvalues that have 
% converged to machine precision and deletes the rest
%
% Inputs: 
% V -- Eigenvectors as a cell of functions in cheb-space
% evals -- Eigenvalues

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

% 

% TODO: Need to have variable outputs for numMp, MP, and minVal.
%
% numMp stores how many basis functions were needed to reach machine
% precision. For well-resolved calculations, this is << N + 1.
numMp = zeros(length(evals),1);

% minVal stores what is the minimum value of the spectral coefficients, if
% the solution is not good up to machine precision. If the solution is good
% up to machine precision, then minVal(i) = 0.
minVal = zeros(length(evals),1);

% Mp is a logical array that indicates Machine precision  or not. If a
% solution is in machine precision, it stores 1, else 0.
Mp = zeros(length(evals),1);

% Largest machine epsilon is typically in the order of 1e-11 after all
% matrix computations etc.
epsilon = 1e-11;
epsilon = abs(epsilon);

% Convert to cells for functions

V1 = V;
[r,c] = size(V1);
%l2norms = zeros(size(V1));
% Function with largest L2norm ascertains sorting criteria for 
% block matrices:

% Compute the l2 norm of every function
%for i = 1:r
%    for j = 1:c
%        l2norms(i,j) = L2norm(cheb2phys(V1{i,j}));
%    end
%end
l2norms = L2norm(cheb2phys(V1));

% Find the maximum columnwise
[~,ins] = max(l2norms,[],1);

% Find average columnwise
averageL2norm = mean(l2norms,1);

% Store functions in each column with the max L2 norm:
V2 = cell(1,c);
for i = 1:c
    V2{i} = V1{ins(i),i};
end

if (~isvector(evals))
    evals = diag(evals);
    isDiag = true;
else 
    isDiag = false;
end
% Find if machine precision by checking the last spectral coefficient
% Note that I take the average of the last two coefficients as functions
% can be either even or odd sometimes.
for i = 1:c
    % Eigenvector has to be non-zero, if zero I consider that it is not
    % machine precision and set it for deletion
    temp = abs(V2{i});
    if(abs(averageL2norm(i)) > 1e-11)
        
        numMp(i) = 1;
        while((numMp(i) < N) && ((temp(numMp(i) + 1) + temp(numMp(i) + 2))/2 > epsilon))
            numMp(i) = numMp(i) + 2;
        end
    else
        numMp(i) = N + 1;
    end
    
    % There should be non nans, if there, delete this solution
    if isnan(averageL2norm(i))
        numMp(i) = N + 1;
    end
    
    % If infinite, remove
    if (abs(evals(i)) > 1e10)
        numMp(i) = N + 1;
    end
    
    % If not in machine precision make note of how good is the
    % approximation, judging by how close is the average of the last two
    % spectral coefficients to machine precision.
    % Else mark that the solution has reached machine precision
    if (numMp(i) >= N)
        minVal(i) = (temp(N) + temp(N+1))/2;
        Mp(i) = 1;
    else
        Mp(i) = 0;
    end
end
Mp = logical(Mp);

evals(Mp==1) = [];
V(:,Mp==1) = [];
if (isDiag)
    evals = diag(evals);
end


