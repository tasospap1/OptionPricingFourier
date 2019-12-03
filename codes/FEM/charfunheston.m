clear all
%theta = long run variance
%sigma = Volatility of volatility
%v0 = initial Variance
%rho = correlation
%T = Time till maturity
%r = interest rate
%s0 = initial asset price
x0 = log(80);
kappa = 1.5768;
theta = 0.0398;
sigma = 0.5751; 
%strike = 100;
rho = -0.571;
r = 0.01;
T = 1; 
v0 = 0.0175;

alpha = 1.1;
%N= 4096;
%c = 600;
%eta = c/N;
%b =pi/eta;
%u = [0:N-1]*eta;
%lambda = 2*b/N;
%position = (log(strike) + b)/lambda + 1; %position of call


%value in FFT
%matrix
v = 1 - (alpha+1)*i;
zeta = -.5*(v.^2 +i*v);
gamma = kappa - rho*sigma*v*i;
PHI = sqrt(gamma.^2 - 2*sigma^2*zeta);
A = i*v*(x0 + r*T);
B = v0*((2*zeta.*(1-exp(-PHI.*T)))./(2*PHI - ...
(PHI-gamma).*(1-exp(-PHI*T))));
C = -kappa*theta/sigma^2*(2*log((2*PHI - ...
(PHI-gamma).*(1-exp(-PHI*T)))./ ...
(2*PHI)) + (PHI-gamma)*T);
charFunc = exp(A + B + C);
u = 1;
ModifiedCharFunc = charFunc*exp(-r*T)./(alpha^2 ...
+ alpha - u.^2 + i*(2*alpha +1)*u)