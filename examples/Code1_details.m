% Set problem data
% number of basis functions, differential order, spatial variable, parameters
N = 63; % N has to be an odd number
m = 2; % order of differential operator
y = sety(N); % spatially-independent variable

% parameter eps in TPBVP
eps = 1; eps2 = eps*eps;


% Represent operator in TPBVP
Delta = cell(m+1,1); % cell array with coefficients
Delta{1} = 1.0; Delta{2} = 1./(y.^2 + 1); Delta{3} = -eps2;

% Use Discretize to represent Delta
Delta = Discretize(m,N,Delta); % inputs: differential order m, N, and the operator


% Dirichlet BCs
% Identity operator
I = cell(1,1); % cell array with coefficients
I{1} = 1.0;

% Dirichlet BCs at y = \pm 1
bcOp = {I; I};
bcPts = [1; -1];
% inputs to BcMat: differential order m, N, boundary points, and boundary operator
bc1 = BcMat(m,N,bcPts,bcOp);

% Neumann BCs
% 1st derivative operator: 1.0 Dy + 0.0 I
Dy = cell(2,1); % cell array with coefficients
Dy{1} = 1.0; Dy{2} = 0.0;

% Neumann BCs at y = \pm 1
bcOp = {Dy; Dy};
bcPts = [1; -1];
bc2 = BcMat(m,N,bcPts,bcOp);

% Robin BCs
% Operator 4 Dy + 3 I
Op = cell(2,1); % cell array with coefficients
Op{1} = 4; Op{2} = 3;

% Robin BCs at y = \pm 1
bcOp = {Op; Op};
bcPts = [1; -1];
bc3 = BcMat(m,N,bcPts,bcOp);

% Operator in the equation for spectral coefficients
A1 = [Delta; bc1]; % Dirichlet BCs
A2 = [Delta; bc2]; % Neumann BCs
A3 = [Delta; bc3]; % Robin BCs

% Forcing in physical space
d = 1 + y + y.^2;

% Forcing in the basis of Chebyshev polynomials
d = phys2cheb(d);

% Solve differential equation
Phi2IntConst1 = A1\[d; 1; -1]; % Dirichlet BCs
Phi2IntConst2 = A2\[d; 2; -2]; % Neumann BCs
Phi2IntConst3 = A3\[d; 3; -3]; % Robin BCs

% Use Matgen to represent J
[~,~,J] = Matgen(m,N); % inputs: the order of the differential equation and N
% Compute spectral coefficients in the basis representation of \phi (y)
Phi01 = J{end}*Phi2IntConst1; % Dirichlet BCs
Phi02 = J{end}*Phi2IntConst2; % Neumann BCs
Phi03 = J{end}*Phi2IntConst3; % Robin BCs

% Solutions in physical space
solution1 = cheb2phys(Phi01); % Dirichlet BCs
solution2 = cheb2phys(Phi02); % Neumann BCs
solution3 = cheb2phys(Phi03); % Robin BCs

% Visualize solution for Dirichlet BCs
plot(y,solution1);

% Compare with Chebfun
A = chebop([-1 1]);
yc = chebfun('y');
A.op = @(y,u) diff(u,2) + (1/(yc^2 + 1))*diff(u) - eps2*u;
A.lbc = @(u) u + 1;
A.rbc = @(u) u - 1;
solcDir = A\(1+yc+yc^2);

% Visualize solution for Dirichlet BCs
plot(y,solution1,'-b',yc,solcDir,'--r');
xlabel('$y$');
ylabel('$\phi (y)$');
legend('SISMatlab','Chebfun','location','northwest');
print('-painters','-dsvg','../docs/pics/cod1Dir');

% Visualize solution for Neumann BCs
A.lbc = @(u) diff(u) + 2;
A.rbc = @(u) diff(u) - 2;
solcNeu = A\(1+yc+yc^2);
plot(y,solution2,'-b',yc,solcNeu,'--r');
xlabel('$y$');
ylabel('$\phi (y)$');
legend('SISMatlab','Chebfun','location','north');
print('-painters','-dsvg','../docs/pics/cod1Neu');


% Visualize solution for Robin BCs
A.lbc = @(u) 4*diff(u) + 3*u + 3;
A.rbc = @(u) 4*diff(u) + 3*u - 3;
solcRob = A\(1+yc+yc^2);
plot(y,solution3,'-b',yc,solcRob,'--r');
xlabel('$y$');
ylabel('$\phi (y)$');
legend('SISMatlab','Chebfun','location','northwest');
print('-painters','-dsvg','../docs/pics/cod1Rob');

