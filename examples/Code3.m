% Code3 Navier-Stokes equations

clear;
close all;
clc;

%% Eigenvalues

% Set problem data
% number of basis functions, differential order, spatial variable, parameters
N = 127; % number of basis functions, has to be an odd number
m = 2;   % highest differential order
y = sety(N); % spatially-independent variable

% Set parameters and base velocity
Re = 2000; % Reynolds number
kx = 1; kz = 1; % wall-parallel wavenumbers
k2 = kx*kx + kz*kz;
U = 1 - y.^2; % base flow (Poiseuille; U = y for Couette)
Uy = -2*y;    % derivative of U (Poiseuille; Uy = 1 for Couette)

% Represent zero, identity, and first derivative operators
Z = cell(1,1); Z{1} = 0; % zero operator
I = cell(1,1); I{1} = 1; % identity operator
Dy = cell(2,1); Dy{1} = 1; Dy{2} = 0; % first derivative operator

% Represent diagonal entries of the operator F
F11 = cell(3,1);
F11{1} = 1/Re;
F11{2} = 0;
F11{3} = (-k2/Re) - 1i*kx*U;

% Represent off-diagonal entries of the operator F
F12 = cell(1,1); F12{1} = -Uy;
F14 = cell(1,1); F14{1} = -1i*kx;
F24 = cell(2,1); F24{1} = -1; F24{2} = 0;
F34 = cell(1,1); F34{1} = -1i*kz;
F41 = cell(1,1); F41{1} = 1i*kx;
F42 = cell(2,1); F42{1} = 1; F42{2} = 0;
F43 = cell(1,1); F43{1} = 1i*kz;

% Represent the operator F
F = {F11, F12,   Z, F14;
       Z, F11,   Z, F24;
       Z,   Z, F11, F34;
     F41, F42, F43, Z};

% Represent the operator E
E = {I, Z, Z, Z;
     Z, I, Z, Z;
     Z, Z, I, Z;
     Z, Z, Z, Z};


% Boundary conditions
bc = BCs(8,4); % eight constraints on four variables

% Boundary conditions at y = \pm 1
% Boundary operator
bc.Operator = {I,  Z, Z, Z;  % u(+1) = 0
	       Z,  I, Z, Z;  % v(+1) = 0
	       Z,  Z, I, Z;  % w(+1) = 0
	       Z, Dy, Z, Z;  % D v(+1) = 0
	       I,  Z, Z, Z;  % u(-1) = 0
	       Z,  I, Z, Z;  % v(-1) = 0
	       Z,  Z, I, Z;  % w(-1) = 0
	       Z, Dy, Z, Z}; % D v(-1) = 0

% Boundary points
bc.Points = [ones(4,4);  % first 4 rows for BCs at y = +1
	    -ones(4,4)]; %  next 4 rows for BCs at y = -1

% Differential orders of u, v, w, and p
nd = [2, 2, 2, 2];

% Use sisEig to conduct the eigenvalue decomposition
[V,lambda] = sisEigs(nd,N,F,E,bc,100,'Full');

plot(real(lambda),imag(lambda),'xk');
ylabel('$\mathrm{Im}(\lambda)$');
xlabel('$\mathrm{Re}(\lambda)$');
ax = gca;
ax.XLim = [-2,0];
ax.YLim = [-1,-0.2];
ax.YTick = [-1 -0.8 -0.6 -0.4 -0.2];
ax.XTick = [-2 -1.5 -1 -0.5 0];
print('-painters','-dsvg','../docs/pics/Code3_1');

%% Frequency responses

% Set problem data
% number of basis functions, differential order, spatial variable, parameters
N = 127; % number of basis functions, has to be an odd number
nd = [2, 2, 2, 2]; % differential orders of u, v, w, and p
y = sety(N); % spatially-independent variable

% Set parameters and base velocity
Re = 2000; % Reynolds number
kx = 1; kz = 1; % wall-parallel wavenumbers
k2 = kx*kx + kz*kz;
U = 1 - y.^2; % base flow (Poiseuille; U = y for Couette)
Uy = -2*y;    % derivative of U (Poiseuille; Uy = 1 for Couette)

% Temporal frequency
omval = linspace(-2,0,20);

% Store first two singular values
sval1 = zeros(length(omval),1);
sval2 = zeros(length(omval),1);

% Start computations
for ind = 1:length(omval)
    
    % temporal frequency
	om = omval(ind);

	% Represent zero, identity, and first derivative operators
	Z = cell(1,1); Z{1} = 0; % zero operator
	I = cell(1,1); I{1} = 1; % identity operator
	Dy = cell(2,1); Dy{1} = 1; Dy{2} = 0; % first derivative operator

    
    % Represent diagonal entries of the operator A
    A11 = cell(3,1); A11{1} = -1/Re; A11{2} = 0; A11{3} = 1i*om + k2/Re + 1i*kx*U;

    % Represent off-diagonal entries of the operator A
    A12 = cell(1,1) ; A12{1} = Uy;
    A14 = cell(1,1); A14{1} = 1i*kx;
    A24 = cell(2,1); A24{1} = 1; A24{2} = 0;
    A34 = cell(1,1); A34{1} = 1i*kz;
    A41 = cell(1,1); A41{1} = -1i*kx;
    A42 = cell(2,1); A42{1} = -1; A42{2} = 0;
    A43 = cell(1,1); A43{1} = -1i*kz;

    % Represent the operator A
    A = {A11, A12, Z, A14;...
         Z, A11, Z, A24;...
         Z, Z, A11, A34;...
         A41, A42, A43, Z};
	
	% Represent operators B and C
	B = {I, Z, Z;
	     Z, I, Z;
	     Z, Z, I;
	     Z, Z, Z};

	C = {I, Z, Z, Z;
	     Z, I, Z, Z;
	     Z, Z, I, Z};
    
	% Boundary conditions
	bc = BCs(8,4); % eight constraints on four variables

	% Boundary conditions at y = \pm 1
	% Boundary operator
	bc.Operator = {I,  Z, Z, Z;  % u(+1) = 0
		       Z,  I, Z, Z;  % v(+1) = 0
		       Z,  Z, I, Z;  % w(+1) = 0
		       Z, Dy, Z, Z;  % D v(+1) = 0
		       I,  Z, Z, Z;  % u(-1) = 0
		       Z,  I, Z, Z;  % v(-1) = 0
		       Z,  Z, I, Z;  % w(-1) = 0
		       Z, Dy, Z, Z}; % D v(-1) = 0

	% Boundary points
	bc.Points = [ones(4,4);  % first 4 rows for BCs at y = +1
	            -ones(4,4)]; %  next 4 rows for BCs at y = -1

	% Use function sisSvdfrs to compute singular values
	svals = sisSvdfrs(nd,nd,N,A,B,C,bc,bc,2);

	% Store first two singular values
	sval1(ind) = svals(1);
	sval2(ind) = svals(2);
    disp(ind)
end


% Compute Hinf norm

% Represent diagonal entries of the operator F
F11 = cell(3,1); F11{1} = 1/Re; F11{2} = 0; F11{3} = (-k2/Re) - 1i*kx*U;

% Represent off-diagonal entries of the operator F
F12 = cell(1,1); F12{1} = -Uy;
F14 = cell(1,1); F14{1} = -1i*kx;
F24 = cell(2,1); F24{1} = -1; F24{2} = 0;
F34 = cell(1,1); F34{1} = -1i*kz;
F41 = cell(1,1); F41{1} = 1i*kx;
F42 = cell(2,1); F42{1} = 1; F42{2} = 0;
F43 = cell(1,1); F43{1} = 1i*kz;

% Represent the operator F
F = {F11, F12,   Z, F14;
       Z, F11,   Z, F24;
       Z,   Z, F11, F34;
     F41, F42, F43, Z};

% Represent the operator E
E = {I, Z, Z, Z;
     Z, I, Z, Z;
     Z, Z, I, Z;
     Z, Z, Z, Z};

[omega_opt,Hinf] = sisHinf(nd,N,E,F,B,C,bc,bc);
disp(['Omega_opt: ' num2str(omega_opt)]);
disp(['Hinfinity norm: ' num2str(Hinf)]);

% Visualize results
semilogy(omval, real(sval1),'-ok',omval, real(sval2),'-xk');
ylabel('$\sigma_i$');
xlabel('$\omega$');
hh = legend('$\sigma_0$','$\sigma_1$','location','northwest');
hh.Interpreter = 'latex';
ax = gca;
ax.YTick = [1 10 100];
print('-painters','-dsvg','../docs/pics/Code3_2');