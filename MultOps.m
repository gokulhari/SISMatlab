function out = MultOps(L1,L2)
% MULTOPS
% Part of SISMatlab
% Computes the composition of two operators, L1 and L2
% For block matrix operators, specify as cell of cells.

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
if (~iscell(L1{1,1}))
    out = cell(length(L1) + length(L2) - 1,1); % Length of composition
    for  i =1: length(out)
        out{i} = 0.0;
    end
    for  i = 0: length(L1)-1
      out = addops(out, MultFunOp(L1{i+1},diff(L2, length(L1)- 1 - i)));
    end
else
    [lr,lc] = size(L1);
    [rr,rc] = size(L2);
    out = cell(lr,rc);
    if (lc~=rr)
        error('Incompatible sizes for operator multiplication')
    end
    for i = 1:lr
        for j = 1:rc
            out{i,j} = {0};
        end
    end
    
    for i= 1:lr
        for j = 1:rc
            for k = 1:lc
                out{i,j} = addops(out{i,j},MultOps(L1{i,k},L2{k,j}) );
            end
        end
    end
    
end

end

function out = addops(L,R)

% if addition of simple linear operators then:
if (~iscell(L{1,1}))
    if (length(L) > length(R))
        diffn = length(L) - length(R);
        out = cell(size(L));
        for i = length(R):-1:1
            out{i+diffn} = L{i + diffn} + R{i};
        end
        for i = 1:diffn
            out{i} = L{i};
        end
    else
        diffn = length(R) - length(L);
        for i = length(L):-1:1
            out{i+diffn} = R{i + diffn} + L{i};
        end
        for i = 1:diffn
            out{i} = R{i};
        end
    end
else % if block matrix operators then:
    if (size(L)~=size(R) || size(L,2)~=size(R,2))
        error('wrong dimensions')
    end
    out = cell(size(L));
    [r,c] = size(L);
    for i = 1:r
        for j = 1:c
            out{i,j} = addops(L{i,j},R{i,j});
        end
    end
end

end

function out = diff(L,times)
if (times == 0)
    out = L;
else
    out = cell(length(L)+1,1);
    out{1} = L{1};
    out{length(out)} = diffFun(L{length(L)},1);
    for i = 1:length(out)-2
        out{i+1} = diffFun(L{i},1) + L{i+1};
    end
    out = diff(out,times-1);
end
end

function out = diffFun(fun,times)
    if (times == 0)
        out = fun;
    else
        out = chebdiff(fun);
        out = diff(out,times-1);
    end
end

function du = chebdiff(u)
    if(isscalar(u) || isempty(u))
                du = 0;
    else
        u = phys2cheb(u);
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
        du = cheb2phys(du);
    end
end

function out = MultFunOp(fun,L)
out = cell(size(L));
for i = 1:length(L)
    out{i} = fun.*L{i};
end
end