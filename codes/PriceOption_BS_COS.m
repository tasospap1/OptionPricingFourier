function [Price] = PriceOption_BS_COS(r, sigma, K, T, S, PutOrCall)

% Price European Call or Put with "PutOrCall"
% 1 for Put, -1 for Call
N = 2^6; % Truncation of Fourier series

c1 = r*T; c2 = sigma^2*T; c4 = 0; L = 10; % Truncated domain [a,b]
a = c1 - L*sqrt(c2 + sqrt(c4)); b = c1 + L*sqrt(c2 + sqrt(c4));

ksi_k = @(k) ( 1/(1 +(k*pi/(b-a))^2) * ( cos(-k*pi*a/(b-a)) ...
    - exp(a) + k*pi/(b-a)*sin(-k*pi*a/(b-a)) ) );

psi_k = @(k) (((b-a)/k/pi * sin(-k*pi*a/(b-a))).*(k>0) );

VPut_k = @(k) ( 2*K/(b-a) * (-ksi_k(k) + psi_k(k)) );

phiLevy = @(w) ( exp(1i*(r-sigma^2/2)*T*w - 0.5*sigma^2*T*w.^2) );

x = log(S'/K); % Log-price = log(S/K)

Pt = 1/2 * real(phiLevy(0) * exp(zeros(length(S),1))) ...
    *(2*K/(b-a)*(-ksi_k(0)-a));
for k=1:N-1
    Pt = Pt + ...
        real(phiLevy(k*pi/(b-a)) .* exp(1i*k*pi*(x-a)/(b-a))).*VPut_k(k);
end
Pt = Pt * exp(-r*T) ;

if PutOrCall==1
    Price = Pt;
elseif PutOrCall==-1
    Price = Pt + S' - K*exp(-r*T);
end

end

