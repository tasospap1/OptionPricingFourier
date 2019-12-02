% Pricing a put/call options under the Black-Scholes model 
% using the Direct Integration. 
% ---------------------------------------------------
clear all; 

alpha = -1; % Price European Call or Put 
           % 1 for Put, -1 for Call
           
N = 100000; % Num of nodes in integration domain
K = 100; % Strike Price
r = 0.01; % Interest rate
T = 1; % Maturity
sigma = 0.25; % Volatility 

S = 80:1:120; % Stock-price
L = (log(K./S) - T*(r-0.5*sigma^2)) / sigma /sqrt(T); 

z=L; 
Pt = 0.001 * normpdf(L) .* (S.*exp((r-0.5*sigma^2)*T + sigma*sqrt(T)*L) - K);
for k=1:N
   z = z + 0.001; 
   Pt = Pt + 0.001* normpdf(z) .* (S.*exp((r-0.5*sigma^2)*T + sigma*sqrt(T)*z)-K);
end
Pt = Pt * exp(-r*T) ; 

% Exact Black-Scholes solution of Put
Pt_Exact = putExact_BS(K,r,sigma,S,0,T);

% Payoff of Put
Payoff = @(S) ( max(alpha*(K-S), 0) );

if alpha==1 % Pricing Put
    plot(S, Pt_Exact, 'DisplayName', 'Exact P_0'); hold on;
    plot(S, Payoff(S), 'DisplayName', 'Payoff'); hold on; 
    plot(S, Pt, 'DisplayName', 'P_0^{Integ}');
    xlabel('S_0'); ylabel('P_0'); title('European Put Option');
    legend; 
elseif alpha ==-1 % Pricing Call using Put-Call Parity
    Ct = Pt + S - K*exp(-r*T); 
    Ct_Exact = Pt_Exact + S - K*exp(-r*T); 
    plot(S, Ct_Exact, 'DisplayName', 'Exact C_0'); hold on;
    plot(S, Payoff(S), 'DisplayName', 'Payoff'); hold on; 
    plot(S, Ct , 'DisplayName', 'C_0^{Integ}');
    xlabel('S_0'); ylabel('C_0'); title('European Call Option');
    legend;
end
set(gca,'FontSize',14)
    