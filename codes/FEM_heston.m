%  MAIN_HESTONCALL computes European call wit the Heston model using finite 
%          elements
%
function C = FEM_heston(Nx,Ny,T,K,kappa,mt,beta,rho,vol)
%clear all;
%close all;

% -------------------------------------------------------------------------
%  Set parameters
% -------------------------------------------------------------------------

%Nx = 51;               % number of nodes in x-ccordinate
%Ny = 51;               % number of nodes in y-ccordinate
m = 50;                 % number of time steps
R_1 = 4;                % domain (-R_1,R_1)
R_2 = 3.2;              % domain (0,R_2)
%T = 1/2;               % maturity
%K = 1;                 % strike
%rho = -0.5;            % correlation
%kappa = 2.5;           % rate of mean reversion
%mt = 0.06;             % level of mean reversion
%beta = 0.5;            % volatility of volatility

% -------------------------------------------------------------------------
%  Discretization
% -------------------------------------------------------------------------

hx = (2*R_1)/(Nx+1);          % mesh size in x-coordinate 
hy = (R_2)/(Ny+1);            % mesh size in y-coordinate
x = linspace(-R_1,R_1,Nx+2)'; % mesh nodes in x-coordinate 
y = linspace(0,R_2,Ny+2)';    % mesh nodes in y-coordinate 
dt = T/m;                     % time steps

% -------------------------------------------------------------------------
%  Compute Generator/ Source Term/ Initial Data
% -------------------------------------------------------------------------

% functions for the matrices
zero_Handle = @(x) 0;
one_Handle = @(x) 1;
y_Handle = @(x) x;

% non-weighted matrices
e = ones(Nx,1);
M1 = hx/6*spdiags([e,4*e,e], -1:1 , Nx, Nx);
B1 = 0.5*spdiags([-e, zeros(Nx,1),e], -1:1, Nx, Nx);
S1 = spdiags([-e, 2*e, -e], -1:1, Nx, Nx)/hx; 

e = ones(Ny+2,1);
M2 = hy/6*spdiags([e,4*e,e], -1:1 , Ny+2, Ny+2);
B2 = 0.5*spdiags([-e, zeros(Ny+2,1),e], -1:1, Ny+2, Ny+2);


% weighted matrices
Sy = stiff(y, y_Handle, zero_Handle, zero_Handle);
My = stiff(y, zero_Handle, zero_Handle, y_Handle);
By = stiff(y, zero_Handle, y_Handle, zero_Handle);

% incoorporate hom. Neumann boundary conditions
M2(1,1) = 0.5*M2(1,1); M2(Ny+2,Ny+2) = 0.5*M2(Ny+2,Ny+2);
B2(1,1) = -0.5; B2(Ny+2,Ny+2) = 0.5;
%Sy(1,1) = 0.5*Sy(1,1); Sy(end,end) = 0.5*Sy(end,end);
%My(1,1) = 0.5*My(1,1); My(end,end) = 0.5*My(end,end);
%By(1,1) = 0.5*By(1,1); By(end,end) = 0.5*By(end,end);

% define matrices Y1, Y2
Y1 = -beta*rho*By+0.5*My;
Y2 = beta^2/2*Sy+(0.5*beta^2-kappa*mt)*B2+kappa*By;    

% tensor product
%dofx = 1:Nx;
M = kron(M1,M2);                            
A = 0.5*kron(S1,My) + kron(B1,Y1) + kron(M1,Y2); 

% initial data
u0x = max(0,exp(x(2:end-1))-K); 
u0y = e;

u0 = kron(u0x,u0y);  

% -------------------------------------------------------------------------
%  Solver
% -------------------------------------------------------------------------

theta = 0.5;
B = M+dt*theta*A; C = M-(1-theta)*dt*A;

% lopp over time points
u = u0;
for i = 0:m-1
    u = B\(C*u);
end  

% -------------------------------------------------------------------------
%  Postprocessing
% -------------------------------------------------------------------------

% area of interest
idxd = find(x <= -1,1,'last');
idxu = find(x >= 1,1);
idyd = find(y <= 0.1,1,'last');
idyu = find(y >= 1.2,1);

% compute exact solution
S = exp(x(idxd:idxu)); y = y(idyd:idyu);
uex = heston_call(S,y,T,K,rho,kappa,mt,beta,0);
          
% plot option price
u = reshape(u,Ny+2,Nx);    
u = [zeros(Ny+2,1),u,zeros(Ny+2,1)]; 
u = u(idyd:idyu,idxd:idxu); 
u0 = reshape(u0,Ny+2,Nx);    
u0 = [zeros(Ny+2,1),u0,zeros(Ny+2,1)]; 
u0 = u0(idyd:idyu,idxd:idxu); 
[X,Y] = meshgrid(S,y);

[m,idx] = min(abs(y - vol))
%C = u
C = u(idx,:);
S

%figure(1)
%surf(X,Y,u), 
%hold on
%mesh(X,Y,u0)
%title('European Call option in Heston model')
%xlabel('S'), ylabel('y'), zlabel('u')

% plot error
%figure(2)
%mesh(exp(X),Y,abs(u-uex))
%xlabel('S'), ylabel('y'), zlabel('|e|')
end
