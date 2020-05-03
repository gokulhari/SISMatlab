% Code 4: 2D inertialess channel flow of a Viscoelastic fluid

clear;
clc;
close all;

% Frequency response analysis for 2D inertialess channel flow of a viscoelastic fluid

% Set problem data
% number of basis functions, differential order, spatial variable, parameters
N = 255; % number of basis functions, has to be an odd number
m = 4; % differential order of the evolution model
y = sety(N); % spatially-independent variable

% Parameters
Weval = 1:20;  % Weissenberg numbers
beta = 0.5;    % viscosity ration
omega = 0;     % temporal frequency
kx = 1;        % streamwise wavenumber
k2 = kx*kx;    % kx^2

% Flow type
type = 'Couette'; % or 'Poiseuille'

% Store singular values
svalsN255 = zeros(length(Weval),1);

% Computations are done for two values of N to demonstrate grid-independence

% Represent identity and first derivative operators
I = cell(1,1); I{1} = 1; % Identity operator
Dy = cell(2,1); Dy{1} = 1; Dy{2} = 0; % First derivative operator

% Boundary conditions
bc = BCs(4,1); % four constraints on one variable

% Boundary conditions at y = \pm 1
% Boundary operator
bc.Operator = {I; Dy; I; Dy};

% Boundary points
bc.Points = [1; 1; -1; -1];

% Start computations
for ind = 1:length(Weval),

	% Weissenberg number
    We = Weval(ind);

    if strcmp(type,'Poiseuille')

        U = 1 - y.^2;
        Uy = -2*y;
        Uyy = -2*ones(N+1,1);

    else

        U = y;
        Uy = ones(N+1,1);
        Uyy = zeros(N+1,1);

    end

    % Specify coefficients in the frequency response
    c = (1i * omega + 1.0 + (1i * kx * We * U));
    cy = 1i * We * kx * Uy;
    cyy = 1i * We * kx * Uyy;

    % Operator A
    a4 = -beta + (-1 + beta)./c;

    a3 = (-2*(-1 + beta)*(cy - 1i*c.*kx.*Uy.*We))./(c.^2);

    a2 = 2*beta*k2 - ((-1 + beta)*(-2*cy.^2 - 4*kx*Uy.*We.*(1i*cy + kx*Uy*We) + ...
         c.^2.*kx.*(2*kx - 3*1i*Uyy.*We + 2*kx*Uy.^2*We^2) + ...
         c.*(cyy + 2*1i*kx*We*(Uyy + Uy.*(cy + 1i*kx*Uy*We)))))./c.^3;

    a1 = ((-2*1i)*(-1 + beta)*kx*(6*cy.*Uy.*We.*(cy - 1i*kx*Uy*We) + ...
         c.^3.*kx.*Uy.*We.*(kx - 2*1i*Uyy*We) + ...
         c.^2.*(Uyy.*We.*(cy + 3*1i*kx*Uy*We) + kx*(1i*cy - 2*kx*(Uy.^3)*We^3)) - ...
         2*c.*We.*(2*cy.*(Uyy + 1i*kx*Uy.^2*We) + ...
         Uy.*(cyy + kx*We*((-2*1i)*Uyy + kx*Uy.^2*We)))))./c.^4;

    a0 = (kx*(-(beta*c.^4*kx^3) + ...
         12*(-1 + beta)*cy.*kx.*Uy.^2.*We^2.*(cy - 1i*kx.*Uy.*We) + ...
         (-1 + beta)*c.^3*k2.*(kx - 1i*Uyy*We + 2*kx*Uy.^2*We^2) + ...
         (-1 + beta)*c.^2.*(-(cyy*kx) + 1i*(-cyy + 2*k2).*Uyy*We + 2*kx*Uyy.^2.*We^2 + ...
         2*1i*cy.*kx.*Uy.*We.*(kx + 2*1i*Uyy*We) + ...
         2*kx*Uy.^2.*We^2.*(-cyy + k2 + (6*1i)*kx*Uyy.*We)) + ...
         2*(-1 + beta)*c.*(-2*cyy.*kx.*Uy.^2*We^2 + ...
         cy.*kx.*(cy - 2*1i*kx*Uy*We).*(1 + 2*Uy.^2.*We^2) + ...
         1i*Uyy.*We.*(cy.^2 + (6*1i)*cy.*kx.*Uy*We + 6*k2*Uy.^2.*We^2))))./c.^4;

    % Represent operators A, B, and C
    % Operator A
    A = cell(5,1); A{1} = a4; A{2} = a3; A{3} = a2; A{4} = a1; A{5} = a0;

    % Operator B
    B = cell(1,2); B{1,1} = Dy; B{1,2} = {-1i*kx};

    % Operator C
    C = cell(2,1); C{1,1} = Dy; C{2,1} = {-1i*kx};

    % Use function sisSvdfrs to compute singular values
    svals = sisSvdfrs(m,m,N,A,B,C,bc,bc,1);
    svalsN255(ind) = real(svals(1));
    disp(ind);

end

%% N = 383
% Set parameters
N = 383; % Number of basis functions
y = sety(N); % Spatially independent variable

% Store singular values
svalsN383 = zeros(length(Weval),1);

% Start computations
for ind = 1:length(Weval),

	% Weissenberg number
    We = Weval(ind);

    if strcmp(type,'Poiseuille')

        U = 1 - y.^2;
        Uy = -2*y;
        Uyy = -2*ones(N+1,1);

    else

        U = y;
        Uy = ones(N+1,1);
        Uyy = zeros(N+1,1);

    end

    % Specify coefficients in the frequency response
    c = (1i * omega + 1.0 + (1i * kx * We * U));
    cy = 1i * We * kx * Uy;
    cyy = 1i * We * kx * Uyy;

    % Operator A
    a4 = -beta + (-1 + beta)./c;

    a3 = (-2*(-1 + beta)*(cy - 1i*c.*kx.*Uy.*We))./(c.^2);

    a2 = 2*beta*k2 - ((-1 + beta)*(-2*cy.^2 - 4*kx*Uy.*We.*(1i*cy + kx*Uy*We) + ...
         c.^2.*kx.*(2*kx - 3*1i*Uyy.*We + 2*kx*Uy.^2*We^2) + ...
         c.*(cyy + 2*1i*kx*We*(Uyy + Uy.*(cy + 1i*kx*Uy*We)))))./c.^3;

    a1 = ((-2*1i)*(-1 + beta)*kx*(6*cy.*Uy.*We.*(cy - 1i*kx*Uy*We) + ...
         c.^3.*kx.*Uy.*We.*(kx - 2*1i*Uyy*We) + ...
         c.^2.*(Uyy.*We.*(cy + 3*1i*kx*Uy*We) + kx*(1i*cy - 2*kx*(Uy.^3)*We^3)) - ...
         2*c.*We.*(2*cy.*(Uyy + 1i*kx*Uy.^2*We) + ...
         Uy.*(cyy + kx*We*((-2*1i)*Uyy + kx*Uy.^2*We)))))./c.^4;

    a0 = (kx*(-(beta*c.^4*kx^3) + ...
         12*(-1 + beta)*cy.*kx.*Uy.^2.*We^2.*(cy - 1i*kx.*Uy.*We) + ...
         (-1 + beta)*c.^3*k2.*(kx - 1i*Uyy*We + 2*kx*Uy.^2*We^2) + ...
         (-1 + beta)*c.^2.*(-(cyy*kx) + 1i*(-cyy + 2*k2).*Uyy*We + 2*kx*Uyy.^2.*We^2 + ...
         2*1i*cy.*kx.*Uy.*We.*(kx + 2*1i*Uyy*We) + ...
         2*kx*Uy.^2.*We^2.*(-cyy + k2 + (6*1i)*kx*Uyy.*We)) + ...
         2*(-1 + beta)*c.*(-2*cyy.*kx.*Uy.^2*We^2 + ...
         cy.*kx.*(cy - 2*1i*kx*Uy*We).*(1 + 2*Uy.^2.*We^2) + ...
         1i*Uyy.*We.*(cy.^2 + (6*1i)*cy.*kx.*Uy*We + 6*k2*Uy.^2.*We^2))))./c.^4;

    % Represent operators A, B, and C
    % Operator A
    A = cell(5,1); A{1} = a4; A{2} = a3; A{3} = a2; A{4} = a1; A{5} = a0;

    % Operator B
    B = cell(1,2); B{1,1} = Dy; B{1,2} = {-1i*kx};

    % Operator C
    C = cell(2,1); C{1,1} = Dy; C{2,1} = {-1i*kx};

    % Use function sisSvdfrs to compute singular values
    svals = sisSvdfrs(m,m,N,A,B,C,bc,bc,1);
    svalsN383(ind) = real(svals(1));
    disp(ind);

end

% Visualize results
plot(Weval,svalsN255,'ok',Weval,svalsN383,'xk');
axis tight
xlabel('$W\! e$');
ylabel('$\sigma_{\mathrm{max}}$');
hh = legend('$N = 255$','$N = 383$');
hh.Interpreter = 'latex';
print('-painters','-dsvg','../docs/pics/Code4_1');