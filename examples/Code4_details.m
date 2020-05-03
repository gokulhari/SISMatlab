% Code 4: 2D inertialess channel flow of a Viscoelastic fluid

clear;
clc;
close all;

m = 4; % Differential order of the system

% Flow parameters:
Wes = 1:20; % Vector of Weissenberg numbers
beta = 0.5; % beta
omega = 0;  % Temporal frequency
kx = 1;     % Spatial wavenumber

% Couette flow:
type = 'Couette';

svalsN255 = zeros(20,1);
svalsN383 = zeros(20,1);

% Calculations are performed with two values of N to demonstrate
% grid-independence.

%% N = 255

% Set parameters
N = 255; % Number of basis functions
y = sety(N); % Spatially independent variable
% Make the boundary conditions matrix:
I = cell(1,1); I{1} = 1; % Identity operator
Z = cell(1,1); Z{1} = 0; % Zero operator
Dy = cell(2,1); Dy{1} = 1; Dy{2} = 0; % First derivative operator
bcOp = {I;Dy;I;Dy};
bcPoints = [1;1;-1;-1];
bcs = BcMat(m,N,bcPoints, bcOp);
bcs_null = null(bcs);

for i = 1:20
    We = Wes(i);
    ii = sqrt(-1);
    Re = 0;
    beta = 0.5;
    kx = 1;
    k2 = kx*kx;

    if strcmp(type,'Poiseuille')
        U = 1 - y.*y;
        Uy = -2*y;
        Uyy = -2*ones(N+1,1);
    else 
        U = y;
        Uy = ones(N+1,1);
        Uyy = zeros(N+1,1);
    end
    
    omega = 0;
    c = (ii * omega + 1.0 + (ii * kx * We * U));
    cy = ii * We * kx * Uy;
    cyy = ii * We * kx * Uyy;
    
    a4 = -beta + (-1 + beta)./c;
    
    a3 = (-2*(-1 + beta)*(cy - (ii)*c.*kx.*Uy.*We))./(c.^2);
    
    a2 = 2*beta*k2 - ((-1 + beta)*(-2*cy.^2 - 4*kx*Uy.*We.*(ii*cy + kx*Uy*We) + ...
         c.^2.*kx.*(2*kx - 3*ii*Uyy.*We + 2*kx*Uy.^2*We^2) + c.*(cyy + 2*ii*kx*We*(Uyy + Uy.*(cy + ii*kx*Uy*We)))))./...
     c.^3;
    
    a1 = ((-2*ii)*(-1 + beta)*kx*(6*cy.*Uy.*We.*(cy - ii*kx*Uy*We) + c.^3.*kx.*Uy.*We.*(kx - 2*ii*Uyy*We) + ...
        c.^2.*(Uyy.*We.*(cy + 3*ii*kx*Uy*We) + kx*(ii*cy - 2*kx*(Uy.^3)*We^3)) - ...
        2*c.*We.*(2*cy.*(Uyy + ii*kx*Uy.^2*We) + Uy.*(cyy + kx*We*((-2*ii)*Uyy + kx*Uy.^2*We)))))./c.^4;
    
    a0 = (kx*(-(beta*c.^4*kx^3) + 12*(-1 + beta)*cy.*kx.*Uy.^2.*We^2.*(cy - ii*kx.*Uy.*We) + ...
        (-1 + beta)*c.^3*k2.*(kx - ii*Uyy*We + 2*kx*Uy.^2*We^2) + ...
        (-1 + beta)*c.^2.*(-(cyy*kx) + ii*(-cyy + 2*k2).*Uyy*We + 2*kx*Uyy.^2.*We^2 + ...
           2*ii*cy.*kx.*Uy.*We.*(kx + 2*ii*Uyy*We) + 2*kx*Uy.^2.*We^2.*(-cyy + k2 + (6*ii)*kx*Uyy.*We)) + ...
        2*(-1 + beta)*c.*(-2*cyy.*kx.*Uy.^2*We^2 + cy.*kx.*(cy - 2*ii*kx*Uy*We).*(1 + 2*Uy.^2.*We^2) + ...
           ii*Uyy.*We.*(cy.^2 + (6*ii)*cy.*kx.*Uy*We + 6*k2*Uy.^2.*We^2))))./c.^4;
       a3 = a3;
       a2 = a2;
       a1 = a1;
       a0 = a0;
    b0 = 1;
    b1 = -ii*kx;
    
    % Make operators:
    A = cell(5,1);
    A{1} = a4; A{2} = a3; A{3} = a2; A{4} = a1; A{5} = a0;
    
    B = cell(1,2);
    B{1,1} = Dy; B{1,2} = {-1i*kx};
    
    C = cell(2,1);
    C{1,1} = Dy; C{2,1} = {-1i*kx};
    
    Aad = AdjointFormal(A);
    Bad = AdjointFormal(B);
    Cad = AdjointFormal(C);
    
    % Form composite operators:
    BBad = MultOps(B,Bad);
    CadC = MultOps(Cad,C);
    
    % Discretize:
    A = Discretize(m,N,A)*bcs_null;
    Aad = Discretize(m,N,Aad)*bcs_null;
    BBad = Discretize(m,N,BBad)*bcs_null;
    CadC = Discretize(m,N,CadC)*bcs_null;
    
    Zs = zeros(N+1);
    sval = eigs((A\BBad)*(Aad\CadC),1,'LR');
    
    svalsN255(i) = abs(sqrt(sval));
    disp(i)
end

%% N = 383
% Set parameters
N = 383; % Number of basis functions
y = sety(N); % Spatially independent variable
% Make the boundary conditions matrix:
I = cell(1,1); I{1} = 1; % Identity operator
Z = cell(1,1); Z{1} = 0; % Zero operator
Dy = cell(2,1); Dy{1} = 1; Dy{2} = 0; % First derivative operator
bcOp = {I;Dy;I;Dy};
bcPoints = [1;1;-1;-1];
bcs = BcMat(m,N,bcPoints, bcOp);
bcs_null = null(bcs);
for i = 1:20
    We = Wes(i);
    ii = sqrt(-1);
    Re = 0;
    beta = 0.5;
    kx = 1;
    k2 = kx*kx;

    if strcmp(type,'Poiseuille')
        U = 1 - y.*y;
        Uy = -2*y;
        Uyy = -2*ones(N+1,1);
    else 
        U = y;
        Uy = ones(N+1,1);
        Uyy = zeros(N+1,1);
    end
    
    omega = 0;
    c = (ii * omega + 1.0 + (ii * kx * We * U));
    cy = ii * We * kx * Uy;
    cyy = ii * We * kx * Uyy;
    
    a4 = -beta + (-1 + beta)./c;
    
    a3 = (-2*(-1 + beta)*(cy - (ii)*c.*kx.*Uy.*We))./(c.^2);
    
    a2 = 2*beta*k2 - ((-1 + beta)*(-2*cy.^2 - 4*kx*Uy.*We.*(ii*cy + kx*Uy*We) + ...
         c.^2.*kx.*(2*kx - 3*ii*Uyy.*We + 2*kx*Uy.^2*We^2) + c.*(cyy + 2*ii*kx*We*(Uyy + Uy.*(cy + ii*kx*Uy*We)))))./...
     c.^3;
    
    a1 = ((-2*ii)*(-1 + beta)*kx*(6*cy.*Uy.*We.*(cy - ii*kx*Uy*We) + c.^3.*kx.*Uy.*We.*(kx - 2*ii*Uyy*We) + ...
        c.^2.*(Uyy.*We.*(cy + 3*ii*kx*Uy*We) + kx*(ii*cy - 2*kx*(Uy.^3)*We^3)) - ...
        2*c.*We.*(2*cy.*(Uyy + ii*kx*Uy.^2*We) + Uy.*(cyy + kx*We*((-2*ii)*Uyy + kx*Uy.^2*We)))))./c.^4;
    
    a0 = (kx*(-(beta*c.^4*kx^3) + 12*(-1 + beta)*cy.*kx.*Uy.^2.*We^2.*(cy - ii*kx.*Uy.*We) + ...
        (-1 + beta)*c.^3*k2.*(kx - ii*Uyy*We + 2*kx*Uy.^2*We^2) + ...
        (-1 + beta)*c.^2.*(-(cyy*kx) + ii*(-cyy + 2*k2).*Uyy*We + 2*kx*Uyy.^2.*We^2 + ...
           2*ii*cy.*kx.*Uy.*We.*(kx + 2*ii*Uyy*We) + 2*kx*Uy.^2.*We^2.*(-cyy + k2 + (6*ii)*kx*Uyy.*We)) + ...
        2*(-1 + beta)*c.*(-2*cyy.*kx.*Uy.^2*We^2 + cy.*kx.*(cy - 2*ii*kx*Uy*We).*(1 + 2*Uy.^2.*We^2) + ...
           ii*Uyy.*We.*(cy.^2 + (6*ii)*cy.*kx.*Uy*We + 6*k2*Uy.^2.*We^2))))./c.^4;
       a3 = a3;
       a2 = a2;
       a1 = a1;
       a0 = a0;
    b0 = 1;
    b1 = -ii*kx;
    
    % Make operators:
    A = cell(5,1);
    A{1} = a4; A{2} = a3; A{3} = a2; A{4} = a1; A{5} = a0;
    
    B = cell(1,2);
    B{1,1} = Dy; B{1,2} = {-1i*kx};
    
    C = cell(2,1);
    C{1,1} = Dy; C{2,1} = {-1i*kx};
    
    Aad = AdjointFormal(A);
    Bad = AdjointFormal(B);
    Cad = AdjointFormal(C);
    
    % Form composite operators:
    BBad = MultOps(B,Bad);
    CadC = MultOps(Cad,C);
    
    % Discretize:
    A = Discretize(m,N,A)*bcs_null;
    Aad = Discretize(m,N,Aad)*bcs_null;
    BBad = Discretize(m,N,BBad)*bcs_null;
    CadC = Discretize(m,N,CadC)*bcs_null;
    
    Zs = zeros(N+1);
    sval = eigs((A\BBad)*(Aad\CadC),1,'LR');
    
    svalsN383(i) = abs(sqrt(sval));
    disp(i)
end

% Plot
plot(Wes,svalsN255,'ok',Wes,svalsN383,'xk');
axis tight
xlabel('$W\! e$');
ylabel('$\sigma_{\mathrm{max}}$');
hh = legend('$N = 255$','$N = 383$');
hh.Interpreter = 'latex';
print('-painters','-dsvg','../docs/pics/Code4_1');