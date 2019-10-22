% Pricing a put/call option under the Black-Scholes model 
% using the Fourier Cosine Expansion method. 
% ---------------------------------------------------

alpha = 1; % Price European Call or Put 
           % 1 for Put, -1 for Call
N = 2^6; % Truncation of Fourier series
K = 60; % Strike Price
r = 0.1; % Interest rate
T = 1; % Maturity
sigma = 0.25; % Volatility 
c1 = r*T; c2 = sigma^2*T; c4 = 0; L = 10; % Truncated domain [a,b]
a = c1 - L*sqrt(c2 + sqrt(c4)); b = c1 + L*sqrt(c2 + sqrt(c4));

ksi_k = @(k) ( 1/(1 +(k*pi/(b-a))^2) * ( cos(-k*pi*a/(b-a)) ...
                - exp(a) + k*pi/(b-a)*sin(-k*pi*a/(b-a)) ) ); 

psi_k = @(k) (((b-a)/k/pi * sin(-k*pi*a/(b-a))).*(k>0) );
                    
VPut_k = @(k) ( 2*K/(b-a) * (-ksi_k(k) + psi_k(k)) ); 

phiLevy = @(w) ( exp(1i*r*T*w - 0.5*sigma^2*T*w.^2) ); 

S = 20; % Stock-price
x = log(S'/K); % Log-price = log(S/K)

Pt = 1/2 * real(phiLevy(0) .* exp(zeros(length(S),1))) ...
                *(2*K/(b-a)*(-ksi_k(0)-a));
for k=1:N-1
   Pt = Pt + ...
       real(phiLevy(k*pi/(b-a)) .* exp(1i*k*pi*(x-a)/(b-a))).*VPut_k(k);  
end
Pt = Pt * exp(-r*T) ; 

% Exact Black-Scholes solution of Put
Pt_Exact = putExact_BS(K,r,sigma,S,0,T)';
% Payoff of Put
Payoff = @(S) ( max(alpha*(K-S), 0) );

if alpha==1 % Pricing Put
    plot(S, Pt_Exact, 'DisplayName', 'Exact P_0'); hold on;
    plot(S, Payoff(S), 'DisplayName', 'Payoff'); hold on; 
    plot(S, Pt, 'DisplayName', 'P_0^{Fourier}');
    xlabel('S_0'); ylabel('P_0'); title('European Put Option');
    legend; 
elseif alpha ==-1 % Pricing Call using Put-Call Parity
    Ct = Pt + S' - K*exp(-r*T); 
    Ct_Exact = Pt_Exact + S' - K*exp(-r*T); 
    plot(S, Ct_Exact, 'DisplayName', 'Exact C_0'); hold on;
    plot(S, Payoff(S), 'DisplayName', 'Payoff'); hold on; 
    plot(S, Ct , 'DisplayName', 'C_0^{Fourier}');
    xlabel('S_0'); ylabel('C_0'); title('European Call Option');
    legend;
end
set(gca,'FontSize',14)
    