% Pricing a put option under the Black-Scholes model 
% using the Fourier Cosine Expansion method. 
% ---------------------------------------------------

N = 2^6; % Truncation of Fourier series
K = 100; % Strike Price
r = 0.1; % Interest rate
T = 0.1; % Maturity
sigma = 0.25; % Volatility 
c1 = r*T; c2 = sigma^2*T; c4 = 0; L = 10;% Truncated domain [a,b]
a = c1 - L*sqrt(c2 + sqrt(c4)); b = c1 + L*sqrt(c2 + sqrt(c4));

ksi_k = @(k) ( 1/(1 +(k*pi/(b-a))^2) * ( cos(-k*pi*a/(b-a)) ...
                - exp(a) + k*pi/(b-a)*sin(-k*pi*a/(b-a)) ) ); 

psi_k = @(k) ( -a*(k==0) + ( (b-a)/k/pi * ...
                    sin(-k*pi*a/(b-a)) )*(k>0) );
                    
VPut_k = @(k) ( 2*K/(b-a) * (-ksi_k(k) + psi_k(k)) ); 

phiLevy = @(w) ( exp(1i*r*T*w - 0.5*sigma^2*T*w.^2) ); 

S = 80:0.1:120; % Stock-price
x = log(S'/K); % Log-To-Maturity-price = log(S/K)

vPut = 0.5 * real(phiLevy(0) * exp(zeros(length(S),1)));
for k=1:1:N-1
   vPut = vPut + ...
       real(phiLevy(k*pi/(b-a)) * exp(1i*k*pi*(x-a)/(b-a)))*VPut_k(k);  
end
vPut = vPut * exp(-r*T); 

% Payoff of Put
PutPayoff = @(S) ( max(K-S, 0) );

% Exact Black-Scholes solution
vPut_Exact = putExact_BS(K,r,sigma,S,0,T)';

plot(S, vPut_Exact, 'DisplayName', 'Exact V_t'); hold on; 
plot(S, PutPayoff(S), 'DisplayName', 'Payoff'); hold on; 
plot(S, vPut, 'DisplayName', 'V_t^{Fourier}'); 
legend; 
