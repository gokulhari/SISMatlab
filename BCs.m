classdef BCs 
% BCS Boundary conditions class to specify boundary conditions in
% SISMatlab.
%
% Specify boundary conditions using this class. For example, to specify
% Robin boundary conditions,
% 2 Dy u(-1) + 3 u(-1) = 4,
% 5 Dy u(+1) + 6 u(+1) = 7,
% Where Dy is the first derivative operator, and u is the variable on which
% the boundary conditions apply, then specify boundary conditions in the
% following manner:
%
% bc = BCs(2,1); % Two boundary conditions, on one variable
%
% OpLeft = cell(2,1); OpLeft{1} = 2; OpLeft{2} = 3; % for 2 Dy + 3 I
% OpRight = cell(2,1); OpRight{1} = 5; OpRight{2} = 6; % for 5 Dy + 6 I
%
% bc.Operator = {OpLeft;OpRight};
% bc.Points = [-1;1];
% bc.Values = [4;7];

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

properties (SetAccess = private)
    r; % Number of constraints (rows)
    c; % Number of variables (columns)
end
properties (SetAccess = public)
  Operator; % Linear operators that act at boundaries
  Points; % Points in the domain where the constraints hold
  Values; % Value at the boundary
end
methods
    function out = BCs(r_,c_)
        out.r = r_;
        out.c = c_;
        out.Operator = cell(r_,c_);
    end
    
    function out = adjoint(obj,A,N)
        out = BCs(obj.r,obj.c);
        % Check if A is simple operator or block matrix operator:
        if (~iscell(A{1,1}))
            out = obj.adjoint({A}); % 1 x 1 linear operator
        else
            % Find the size of the linear operator:
            [lr,lc] = size(A);

            % Find the highest derivative of each variable, columnwise.
            highest_each_col = zeros(lc,1);
            for j = 1:lc
                temp = zeros(lr,1);
                for i = 1:lr
                    temp(i) = length(A{i,j}) - 1;
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

            % Compute Amat based on alphais, based on Cij in theorem. 
            for k = 0:n
                for i = 1:lr
                    for j = 1:lc
                        if (length(A{i,j})-1-k >= 0)
                            alphais{k+1}{i,j} =  A{i,j}{length(A{i,j}) - k}; 
                        end       
                    end
                end
            end
            Amat = cell(n,n);
            for i = 1:n
                for j = 1:n
                    Amat{i,j} = cell(lr,lc);
                    for k = 1:lr
                        for l = 1:lc
                            Amat{i,j}{k,l} = 0;
                        end
                    end
                end
            end
            for i = 0:n-1
                for j = 0:n-i-1
                    for k = i:n-j-1
                       Amat{i + 1,j + 1} =  addop(Amat{i + 1,j + 1}, multop((-1)^k * nchoosek(k, i), ...
                           diffop(alphais{k + j + 2},k-i)));
                       %addop(betais{i+1}, ...
                %multop((-1)^j*nchoosek(j,i),diffop(alphais{j+1},j-i)) );
                    end
                end
            end
            % Evaluate Amat at the endpoints to get Aplus and Amin:
            Aplus = eval(Amat,1);
            Amin = eval(Amat,-1);
%             tplus = ones(N+1,1);
%             tminus = ((-1).^(0:N))';
%             tplus = repmat(tplus,1,obj.c);
%             tminus = repmat(tminus,1,obj.c);
%            % Aplus = tplus*Aplus*tplus';
%            % Amin = tminus*Amin*tminus';
%             [~,~,J] = Matgen(n,N);
%             J1 = J(end:-1:2);
%             
%             matP = 0;
%             matM = 0;
%             
%             for i = 0:n-1
%                 for j = 0:(n-1-i)
%                     matP = matP + J1{i+1}'*tplus*Aplus(obj.c*i + 1: obj.c*(i+1), ...
%                         obj.c*j + 1: obj.c*(j+1))*tplus'*J1{j+1};
%                     matM = matM + J1{i+1}'*tminus*Amin(obj.c*i + 1: obj.c*(i+1), ...
%                         obj.c*j + 1: obj.c*(j+1))*tminus'*J1{j+1};
%                 end
%             end
            
            % Next find the null-space of the boundary conditions, and
            % multiply and find the boundary conditions that satisfy
            % adjoints.
            Blbc = zeros(obj.r,n*obj.c);
            for  k = 0:n
                for i = 0:obj.r-1
                    for j = 0:obj.c-1
                        if (length(obj.Op{i+1,j+1})-1- k >= 0)
                            Blbc(i+1, obj.c * k + j + 1) = ...
                                obj.Op{i+1, j+1}{length(obj.Op{i+1, j+1}) - k};
                      %Blbc(i, Lbc.n * k + j) = Lbc.L(i, j).coef[Lbc.L(i, j).n - k];
                        end
                    end
                end
            end
    
           f= 0; 
           L = null(Blbc);
           Blbcad = Aplus*L;
            
%             for i = 0:n-1
%                 for j = 0:n-1
%                     for l = 0:lr-1
%                         for m=0:lc-1
%                             Aplus(i*lr + l + 1, j*lc +m+1) = ...
%                                 eval(Amat,1);
%                             Amin(i*lr + l + 1, j*lc +m+1) = ...
%                                 eval(Amat,-1);
%                         end
%                     end
%                 end
%             end
            
        end
    end
    
        
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

function out = eval(A,pt)
    % Finds the evaluation of a function at a particular point, used for
    % the cells with chebfuns or single entry.
    [r,c] = size(A);
    [r1,c1] = size(A{1,1});
    out = zeros(r*r1,c*c1);
    for i = 1:r
        for j = 1:c
            for l = 1:r1
                for m = 1:c1
                    out((i-1)*r1+l, (j-1)*c1+m) = chebeval(A{i,j}{l,m},pt);
                end
            end
        end
    end
end

function out = chebeval(fun,pt)
    if isvector(fun)
        out = ChebEval(phys2cheb(fun),pt);
    else
       out = fun;
    end

end