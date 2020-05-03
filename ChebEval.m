function out = ChebEval(A,pts)
% ChebEval(A,pts)
% Part of SISMatlab
% Evaluating a chebpolynomial at a desired point/points
% If A is a vector, then we assume that A is in cheb-space
% If A is a 2D matrix, then the columns are in cheb-space, that is A(i,:)
% is in cheb-space
% If A is a nD matrix, then the last dimension is assumed to be in
% chebspace
% pts has to be a 1D vector

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
dim = ndims(A);
if ~isvector(pts)
    fprintf('Error: pts has to be a 1D vector\n');
    exit;
end
[r,~] = size(pts);
if r ==1
    pts = pts';
end
if isvector(A)
    N = length(A)-1;
    [R,~] = size(A);
else
    N = size(A,dim)-1;
end

Ty = zeros(length(pts),N+1);
for k = 1:N+1
Ty(:,k) = cos((k-1)*acos(pts'));
end
Ty(:,1) = 0.5*Ty(:,1);

if ~ismatrix(A) %if A is not a 2D matrix
    out = zeros(size(A,1),size(A,2),size(pts,1));
    for i = 1:size(A,1)
        out(i,:,:) = squeeze(A(i,:,:))*Ty.';
    end
else
    if ~isvector(A) %if A is not a vector
        out = A*Ty.';
    else
        if (R ==1)
            out = A*Ty.';
        else
            out = Ty*A;
        end
    end
end
    




