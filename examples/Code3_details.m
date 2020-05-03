% Code3 Navier-Stokes equations

clear;
close all;
clc;

%% Eigenvalues

N = 127;
[~,~,J] = Matgen(2,N);
y = sety(N);

% Set parameters:
Re = 2000;
kx = 1;
kz = 1;
U = 1-y.^2; % Poiseuille flow.
Uy = -2*y;

% Make operators: identity, zero, and first derivative operator:
I = cell(1,1); I{1} = 1; Z = cell(1,1); Z{1} = 0;
Dy = cell(2,1); Dy{1} = 1; Dy{2} = 0;

% Make the diagonal of F:
k2 = kx*kx + kz *kz;
F11 = cell(3,1);
F11{1} = 1/Re;
F11{2} = 0;
F11{3} = (-k2/Re) - 1i*kx*U;

% Make other operators for F12, F14 etc.
F12 = cell(1,1) ; F12{1} = -Uy;
F14 = cell(1,1); F14{1} = -1i*kx;
F24 = cell(2,1); F24{1} = -1; F24{2} = 0;
F34 = cell(1,1); F34{1} = -1i*kz;
F41 = cell(1,1); F41{1} = 1i*kx;
F42 = cell(2,1); F42{1} = 1; F42{2} = 0;
F43 = cell(1,1); F43{1} = 1i*kz;

% Make the operator:

F = {F11, F12, Z, F14;...
     Z, F11, Z, F24;...
     Z, Z, F11, F34;...
     F41, F42, F43, Z};

% Make E:
E = {I, Z, Z, Z;...
    Z, I, Z, Z;...
    Z, Z, I, Z;...
    Z, Z, Z, Z};


% Make boundary conditions matrix:
bcOp = {I,Z,Z,Z;...
	Z,I,Z,Z;...
	Z,Z,I,Z;...
	Z,Dy,Z,Z;
	I,Z,Z,Z;...
	Z,I,Z,Z;...
	Z,Z,I,Z;...
	Z,Dy,Z,Z};
bcPoints = [ones(4,4);  % as first 4 rows for bc at y  = 1,
	-ones(4,4)]; %  next four rows for y = -1.

% Differential orders of u, v, w, and p:
n = [2,2,2,2];


% Make boundary conditions matrix
M = BcMat(n,N,bcPoints,bcOp);

% Find null space of M:
M_null = null(M);

% Find eigenvalues
Ff = Discretize(n,N,F)*M_null;
Ef = Discretize(n,N,E)*M_null;
[~,evals] = eig(Ff,Ef);

plot(real(evals),imag(evals),'xk');
ylabel('$\mathrm{Im}(\lambda)$');
xlabel('$\mathrm{Re}(\lambda)$');
ax = gca;
ax.XLim = [-2,0];
ax.YLim = [-1,-0.2];
ax.YTick = [-1 -0.8 -0.6 -0.4 -0.2];
ax.XTick = [-2 -1.5 -1 -0.5 0];
print('-painters','-dsvg','../docs/pics/Code3_1');

%% Frequency responses

% Specify the frequency:
omegas = linspace(-2,0,20);

% Store first two singular values:
eval1 = zeros(20,1);
eval2 = zeros(20,1);

% Forms feedback boundary conditions:
Mf = [M, zeros(size(M));...
zeros(size(M)), M];

% Find null space:
Mf_null = null(Mf);

for i = 1:length(omegas)
	omega = omegas(i);
	% Make operator A:
	A = 1i*omega*Discretize(n,N,E)  - Discretize(n,N,F);
	% Make operator A-adjoint:
	Ead = AdjointFormal(E);
	Fad = AdjointFormal(F);
	Aad = -1i*omega*Discretize(n,N,Ead)  - Discretize(n,N,Fad);

	% BBadjoint and CadjointC is the same, given by
	BBad = {I,Z,Z,Z;...
		Z,I,Z,Z;...
		Z,Z,I,Z;...
		Z,Z,Z,Z};
	CadC = BBad;

	BBad = Discretize(n,N,BBad);
	CadC = Discretize(n,N,CadC);


	%Form feedback operator and find eigenvalues:
	Zs = zeros(size(A));
	lams = eigs([Zs,BBad; CadC,Zs]*Mf_null, [A,Zs; Zs, Aad]*Mf_null, 2, 'LR');
	eval1(i) = lams(1);
	eval2(i) = lams(2);
	disp(i);
end

semilogy(omegas, real(eval1),'-ok',omegas, real(eval2),'-xk');
ylabel('$\sigma_i$');
xlabel('$\omega$');
hh = legend('$\sigma_0$','$\sigma_1$','location','northwest');
hh.Interpreter = 'latex';
ax = gca;
ax.YTick = [1 10 100];
print('-painters','-dsvg','../docs/pics/Code3_2');


%% H-infinity norm

kz = 1;
kx = 1;
F11 = cell(3,1);
F11{1} = 1/Re;
F11{2} = 0;
F11{3} = (-k2/Re) - 1i*kx*U;

% Make other operators for F12, F14 etc.
F12 = cell(1,1) ; F12{1} = -Uy;
F14 = cell(1,1); F14{1} = -1i*kx;
F24 = cell(2,1); F24{1} = -1; F24{2} = 0;
F34 = cell(1,1); F34{1} = -1i*kz;
F41 = cell(1,1); F41{1} = 1i*kx;
F42 = cell(2,1); F42{1} = 1; F42{2} = 0;
F43 = cell(1,1); F43{1} = 1i*kz;

% Make the operator:

F = {F11, F12, Z, F14;...
     Z, F11, Z, F24;...
     Z, Z, F11, F34;...
     F41, F42, F43, Z};

% Make E:
E = {I, Z, Z, Z;...
    Z, I, Z, Z;...
    Z, Z, I, Z;...
    Z, Z, Z, Z};

% Make B:
B = {I,Z,Z;...
    Z,I,Z;...
    Z,Z,I;...
    Z,Z,Z};

% Make C:
C = {I,Z,Z,Z;...
     Z,I,Z,Z;...
     Z,Z,I,Z};
n = [2,2,2,2];
[omega_opt,Hinf] = HinfNorm(n,N,E,F,B,C,Mf_null);
disp(['Omega_opt: ' num2str(omega_opt)]);
disp(['Hinfinity norm: ' num2str(Hinf)]);