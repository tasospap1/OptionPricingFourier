clear all; 

mu = 0.0; % Interest rate
T = 1; % Maturity
kappa = 1.5768; 
sigma = 0.5751; 
theta = 0.0398; 
u0 = 0.0175; % initial volatility 
rho = -0.5711;
K = 100; % Strike Price


% T = 1/2;               % maturity
% K = 1;                 % strike
% rho = -0.5;            % correlation
% kappa = 2.5;           % rate of mean reversion
% theta = 0.06;             % level of mean reversion
% sigma = 0.5;    

S = 80:1:120; 
fft = PriceOption_fft(K,S,kappa,sigma,theta,u0,rho,mu,T,1);
cos = PriceOption_Heston_COS(mu, ...
           K, T, S, 1, kappa, sigma, theta, u0, rho);
%[fem, fem_S] = FEM_heston(100,100,T, K, kappa, theta, sigma, rho, u0); 
Settle = datenum('05-Jun-2005');
Maturity = datemnth(Settle, 12);
fem = optByHestonNI(mu,S,Settle,Maturity,'Put',K,u0,theta,kappa,sigma,rho)';

payoff = @(S) ( max(1*(K-S), 0) );

set(groot,'defaultLegendInterpreter','latex');
plot(S, fft, 'DisplayName', '$FFT$', 'LineWidth', 2); hold on; 
plot(S, cos', 'DisplayName', '$COS$', 'LineWidth', 2); hold on;
plot(S, fem, '--', 'DisplayName', '$FEM$', 'LineWidth', 2); hold on;
plot(S, payoff(S), 'DisplayName', '$(-S_T+K)^{+}$', 'LineWidth', 2); hold on;
xlabel('$S_0$', 'Interpreter', 'latex'); 
ylabel('$P_0$', 'Interpreter', 'latex'); 
title('European Put Option');
legend;
grid on;
xlim([85 110]); 
set(gca,'FontSize',14)