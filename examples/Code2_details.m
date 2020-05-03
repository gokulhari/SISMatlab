% Code2: Reaction-diffusion equation

%% Eigenvalues
clear;
clc;
close all;

% Set problem data
% number of basis functions, differential order, spatial variable, parameters
N = 63; % N has to be an odd number
m = 2; % order of differential operator
y = sety(N); % spatially-independent variable

% parameter eps Eq. (1)
eps = 1; eps2 = eps*eps;

% Represent operator in Eq. (1)
Delta = cell(m+1,1); % cell array with coefficients
Delta{1} = 1.0; Delta{2} = 0; Delta{3} = -eps2;
% Identity operator:
I = cell(1,1);
I{1} = 1;

% Use Discretize to represent Delta
Delta = Discretize(m,N,Delta); % inputs: differential order m, N, and the operator
Imat = Discretize(m,N,I);

% The first derivative operator 1.0 Dy + 0.0:
Dy = cell(2,1); Dy{1} = 1.0; Dy{2} = 0.0;

% Neumann boundary conditions at y = \pm 1
bcOp = {Dy; Dy};
bcPoints = [1; -1];
bcs = BcMat(2,N,bcPoints,bcOp);

% Null-space of the boundary conditions:
bcs_null = null(bcs);
% Effective operators for Dirichlet boundary conditions
Delta = Delta*bcs_null; Imat = Imat*bcs_null;

% The eigen solution
[V,lam] = eig(Delta,Imat);
% Account for null-space parametrization:
Phi2IntConst = bcs_null*V;
lam = diag(lam);

% Generate J:
[~,~,J] = Matgen(2,N);
% Compute spectral coefficients in the basis representation of \phi (y)
Phi0 = J{end}*Phi2IntConst;

% Convert to a cell of functions:
Phi0 = ChebMat2CellMat(Phi0,N);

% Go to physical space:
Phi0 = cheb2phys(Phi0);

% Sort them:
[lam, ind] = sort(lam,'descend');
Phi0 = Phi0(:,ind);


plot(y,Phi0{1,2});

% Compare with Chebfun eigenvalue:
L = chebop([-1 1]);
yc = chebfun('y');
L.op = @(y,u) diff(u,2) - eps2*u;
L.lbc = @(u) diff(u);
L.rbc = @(u) diff(u);

[Vc,lamc] = eigs(L,6);
plot(y,Phi0{1,2}/val_rbc(Phi0{1,2}),'-b',yc,Vc.blocks{5},'--r');
xlabel('$y$');
ylabel('$\phi (y)$');
legend('SISMatlab','Chebfun','location','North');
print('-painters','-dsvg','../docs/pics/Code2_1');

%% Frequency responses
clear;
close all;
clc;

% Set problem data
% number of basis functions, differential order, spatial variable, parameters
N = 63; % N has to be an odd number
m = 2; % order of differential operator
y = sety(N); % spatially-independent variable

% parameter eps Eq. (2)
eps = 1; eps2 = eps*eps;
omega = 0;

% Represent operator in Eq. (3)
A = cell(m+1,1); % cell array with coefficients
A{1} = -1.0; A{2} = 0; A{3} = 1i*omega + eps2;
% Identity operator:
I = cell(1,1);
I{1} = 1;

B = I;
C = I;

%sisSvdfr(2,N,A,B,C,)

% Find adjoint operators:
Aad = AdjointFormal(A);
Bad = AdjointFormal(B);
Cad = AdjointFormal(C);

% Form composite operators:
BBad = MultOps(B,Bad);
CadC = MultOps(Cad,C);

% Use Discretize to represent Delta
A = Discretize(m,N,A); % inputs: differential order m, N, and the operator
Aad = Discretize(m,N,Aad);
BBad = Discretize(m,N,BBad);
CadC = Discretize(m,N,CadC);

% The first derivative operator 1.0 Dy + 0.0:
Dy = cell(2,1); Dy{1} = 1.0; Dy{2} = 0.0;

% Zero operator:
Z = cell(1,1); Z{1} = 0;
% Neumann boundary conditions at y = \pm 1
bcOp = {Dy; Dy};
bcPoints = [1; -1];

bcs = BcMat(2,N,bcPoints,bcOp);

bcs_null = null(bcs);
A = A*bcs_null;
Aad = Aad*bcs_null;
BBad = BBad*bcs_null;
CadC = CadC*bcs_null;

% Compute the eigenvalues of the feedback interconnected system:
Z = zeros(size(A));
[V, gamma] =  eig([Z, BBad; CadC, Z], [A, Z; Z, Aad]);
gamma = diag(gamma);

% Retrive the eigenvectors from null-space parametrization
z = zeros(size(bcs_null));
Phi2Psi2IntConst = [bcs_null, z ; z, bcs_null]*V;

% Go to the lowest derivatives using the integration operator
[~,~,J] = Matgen(m,N);
Z = zeros(size(J{end}));
Phi0Psi0 = [J{end}, Z; Z, J{end}]*Phi2Psi2IntConst;

% Convert to a cell of functions:
Phi0Psi0 = ChebMat2CellMat(Phi0Psi0,N);

% Go to physical space:
Phi0Psi0 = cheb2phys(Phi0Psi0 );

% Sort them:
[gamma, ind] = sort(gamma,'descend');
Phi0Psi0 = Phi0Psi0(:,ind);

plot(y,Phi0Psi0{1,2}/val_rbc(Phi0Psi0{1,2}));
print('-painters','-dsvg','../docs/pics/Code2_2');


