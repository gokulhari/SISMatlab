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

% Input in physical space
d = 1 + y + y.^2;

% Dirichlet BCs
bc1 = BCs(2,1); % two constraints on one variable
% Identity operator
I = cell(1,1); % cell array with coefficients
I{1} = 1.0;

% Dirichlet BCs at y = \pm 1
bc1.Operator = {I; I}; % boundary operator
bc1.Points = [1; -1];  % boundary points
bc1.Values = [1; -1];  % BCs at the boundaries ([0; 0] for homogeneous BCs)

% Neumann BCs
bc2 = BCs(2,1); % two constraints on one variable
% 1st derivative operator: 1.0 Dy + 0.0 I
Dy = cell(2,1); % cell array with coefficients
Dy{1} = 1.0; Dy{2} = 0.0;

% Neumann BCs at y = \pm 1
bc2.Operator = {Dy; Dy};
bc2.Points = [1; -1];
bc2.Values = [2; -2];

% Robin BCs
bc3 = BCs(2,1); % two constraints on one variable
% Operator 4 Dy + 3 I
Op = cell(2,1); % cell array with coefficients
Op{1} = 4; Op{2} = 3;

% Robin BCs at y = \pm 1
bc3.Operator = {Op; Op};
bc3.Points = [1; -1];
bc3.Values = [3; -3];

% Use sisSolves to solve TPBVP
% inputs to sisSolve: differential order m, N, operator, BCs, input forcing d
solution1 = sisSolves(m,N,Delta,bc1,d); % Dirichlet BCs
solution2 = sisSolves(m,N,Delta,bc2,d); % Neumann BCs
solution3 = sisSolves(m,N,Delta,bc3,d); % Robin BCs

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

