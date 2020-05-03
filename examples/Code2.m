% Code 2: Eigenvalues and frequency responses of the reaction-diffusion
% equation.

%% Eigenvalues
% Set problem data
% number of basis functions, differential order, spatial variable, parameters
N = 63; % N has to be an odd number
m = 2; % order of differential operator
y = sety(N); % spatially-independent variable

% parameter eps in reaction-diffusion equation
eps = 1; eps2 = eps*eps;

% Represent operators Delta = D^2 - eps^2 I and I in e-value problem
Delta = cell(m+1,1); % cell array with coefficients
Delta{1} = 1.0; Delta{2} = 0; Delta{3} = -eps2;
I = cell(1,1); % Identity operator
I{1} = 1;

% The first derivative operator 1.0 Dy + 0.0
Dy = cell(2,1); Dy{1} = 1.0; Dy{2} = 0.0;

% Neumann boundary conditions at y = \pm 1
bc = BCs(2,1); % two constraints on one variable
bc.Operator = {Dy; Dy};
bc.Points = [1; -1];

% Use sisEig to conduct the eigenvalue decomposition
[V,lambda] = sisEigs(m,N,Delta,I,bc,10,'Full');

%% Frequency responses
clear;
close all;
clc;

% Set problem data
% number of basis functions, differential order, spatial variable, parameters
N = 63; % N has to be an odd number
m = 2; % order of differential operator
y = sety(N); % spatially-independent variable

% parameter eps in TPBVP
eps = 1; eps2 = eps*eps;
% temporal frequency
omega = 0;

% Represent operators A, B, and C in the frequency response operator
A = cell(m+1,1); % cell array with coefficients
A{1} = -1.0; A{2} = 0; A{3} = 1i*omega + eps2;

% Identity operator
I = cell(1,1); I{1} = 1;

% Input and output operators B and C
B = I; C = I;

% The first derivative operator 1.0 Dy + 0.0:
Dy = cell(2,1); Dy{1} = 1.0; Dy{2} = 0.0;

% Neumann boundary conditions at y = \pm 1
bc = BCs(2,1); % two constraints on one variable
bc.Operator = {Dy; Dy};
bc.Points = [1; -1];

% Solve
[Phi0Psi0,gamma] = sisSvdfrs(2,2,N,A,B,C,bc,bc);
plot(y,Phi0Psi0{1,2}/val_rbc(Phi0Psi0{1,2}));