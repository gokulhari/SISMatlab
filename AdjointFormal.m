% ADJOINTFORMAL
% Part of SISMatlab.
% Computes the formal adjoint of a block matrix operator (can be non square)
% 
% Input: either a simple linear operator in cell format, see Discretize.m
% or a block matrix operator as cell of cells.
%
% For example, for the reaction-diffusion operator (D^2 - k^2), 
% input a cell as
%
% L = cell(3,1);
% L{1} = 1;
% L{2} = 0;
% L{3} = -k^2;
%
% and call AdjointFormal(L).
% 
% For block matrix operator, specify elements of block matrix operator as
% above and the make a cell of cells, for example, for a 2 x 2 block matrix
% operator,
%
% L = cell(2,2);
% L{1,1} = L11;
% L{1,2} = L12;
% L{2,1} = L21;
% L{2,2} = L22;
%
% where L11, L12 etc are simple operators like the reaction-diffusion
% operator.


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
% Written by Gokul Hariharan

function out = AdjointFormal(Lop)

% Find if it is a simple operator, if yes, then make a 1 x 1 block-matrix
% operator and call adjointFormal
if (~iscell(Lop{1,1}))
    L = cell(1,1);
    L{1,1} = Lop;
    out = AdjointFormal(L);
    out = out{1};
else 
    % Find the size of the linear operator:
    [lr,lc] = size(Lop);
    
    % Find the highest derivative of each variable, columnwise.
    highest_each_col = zeros(lc,1);
    for j = 1:lc
        temp = zeros(lr,1);
        for i = 1:lr
            temp(i) = length(Lop{i,j}) - 1;
        end
        highest_each_col(j) = max(temp);
    end
    
    % Overall order:
    n = max(highest_each_col);
    
    % Set alphais, spatially varying coefficients of regular operator
    % derivative wise:
    alphais = cell(n+1,1);
    for i = 1:length(alphais)
        alphais{i} = cell(lr,lc);
    end
    for i = 1:length(alphais)
        alphais{i} = cell(lr,lc);
        for j = 1:lr
            for k = 1:lc
                alphais{i}{j,k} = 0;
            end
        end
    end
    
    for k = 0:n
        for i = 1:lr
            for j = 1:lc
                if (length(Lop{i,j})-1-k >= 0)
                    alphais{k+1}{i,j} =  Lop{i,j}{length(Lop{i,j}) - k}; 
                end       
            end
        end
    end
    
    % Calculate betais, spatially varying coefficients of adjoint operators
    % derivative wise:
    betais = cell(n+1,1);
    for i = 1:length(betais)
        betais{i} = cell(lr,lc);
        for j = 1:lr
            for k = 1:lc
                betais{i}{j,k} = 0;
            end
        end
    end
    
    for i = 0:n
        for j = i:n
            betais{i+1} = addop(betais{i+1}, ...
                multop((-1)^j*nchoosek(j,i),diffop(alphais{j+1},j-i)) );
        end
    end
    
    % take conjugate transpose:
    for i = 1:length(betais)
        betais{i} = ctransposeop(betais{i});
    end
    
    out = cell(lc,lr);
    % pack output:
    for i = 1:lc
        for j = 1:lr
            for k = 0:n
                out{i,j}{n - k + 1} = betais{k+1}{i,j};
            end
        end
    end
    
    temp = cell(size(out));
    
    % Reduce order so that orders of Linops are transpose of orders
    for i = 1:lc
        for j= 1:lr
            temp{i,j} = cell(size(Lop{j,i}));
            for k = 0:length(Lop{j,i})-1
                temp{i,j}{length(Lop{j,i}) -k} = out{i,j}{length(out{i,j})-k};
            end
        end
    end
    
    out = temp;        
end
end

function out = addop(Lop1,Lop2)
    [lr1,lc1] = size(Lop1);
    [lr2,lc2] = size(Lop2);
    
     % Throw error if not same size
     if (lr1 ~= lr2 || lc1 ~= lc2)
         error('Cannot add two block matrix operators of different sizes')
     end
     out = Lop1;
     for i = 1:lr1
         for j = 1:lc1
            out{i,j} = Lop1{i,j} + Lop2{i,j};
         end
     end
end
%             if (length(Lop1{i,j}) > length(Lop2{i,j}))
%                 diffn = length(Lop1{i,j}) - length(Lop2{i,j});
%                 out1 = cell(length(Lop1{i,j}), 1);
%                 for k = length(Lop2{i,j})-1:-1:0
%                     out1{k+diffn + 1} = Lop1{i,j}{k+diffn + 1} + Lop2{i,j}{k+1};
%                 end
%                 for k = 0:diffn-1
%                     out1{k+1} = Lop1{i,j}{k+1};
%                 end
%             out{i,j} = out1;
%             else
%                diffn = length(Lop2{i,j}) - length(Lop1{i,j}); 
%                out1 = cell(length(Lop2{i,j}), 1);
%                for k = length(Lop1)-1:-1:0
%                    out1{k+diffn+1} = Lop2{i,j}{k+diffn+1} + Lop1{i,j}{k+1};
%                end
%                for k = 0:diffn-1
%                    out1{k+1} = Lop2{i,j}{k+1};
%                end
%                out{i,j} = out1;
%             end
%         end
%     end
% end

function out = multop(num,Lop)
    out = Lop;
    [r,c] = size(Lop);
    for i = 1:r
        for j = 1:c
            out{i,j} = num*Lop{i,j};
        end
    end
end

function out = diffop(Lop,times)
    [r,c] = size(Lop);
    out = Lop;
    for i = 1:r
        for j = 1:c
            out{i,j} = diff(Lop{i,j},times); 
        end
    end
end

function out = ctransposeop(Lop)
    [r,c] = size(Lop);
    out = cell(c,r);
    for i = 1:r
        for j = 1:c
            out{j,i} = conj(Lop{i,j});
        end
    end
end

function out = diff(fun,times)
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

