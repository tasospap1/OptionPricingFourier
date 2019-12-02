function Ct = PriceOption_fft(K,S,kappa,sigma,theta,u0,rho,mu,T,PUTorCALL, alpha)

clear all; 
%PUTorCALL = -1; % Price European Call (-1) or Put (1) 
N = 2^12; % Truncation of Fourier series
%alpha = 1.25; % or 1.5


upperlimit = 500;
eta = upperlimit/N;
b =pi/eta;
v = [0:N-1] * eta;
lambda = 2*b/N;


%K = 100; % Strike Price
index = round((log(K) + b)/lambda)+ 1; %position of call

%mu = 0.01; % Interest rate
%T = 1; % Maturity
%kappa = 1.5768; 
%sigma = 0.5751; 
%theta = 0.0398; 
%u0 = 0.0175; % initial volatility 
%rho = -0.5711;
%S = 80:1:120;
x0 = log(S); 


D = @(w)  sqrt((kappa-1i*rho*sigma*w)^2 + (w^2 + 1i*w) * sigma^2) ; 
G = @(w)  (kappa - 1i*rho*sigma*w - D(w)) ...
                    / (kappa - 1i*rho*sigma*w + D(w)) ;

phiHes = @(w,s) exp(1i*w*mu*T + u0/sigma^2 *(1-exp(-D(w)*T))/(1-G(w) ...
                *exp(-D(w)*T))*(kappa-1i*rho*sigma*w-D(w))) * ...
                exp(kappa*theta/sigma^2*(T*(kappa-1i*rho*sigma*w-D(w))- ...
                2*log((1-G(w)*exp(-D(w)*T))/(1-G(w))))) * exp(1i * w * s); 
            
            
psi_T = @(w,s) exp(-mu*T) * phiHes(w-(alpha+1)*1i,s) / (alpha^2+alpha - w^2 ...
                + 1i*(2*alpha+1)*w);
            
Ct = zeros(1,size(x0,2));
for j = 1 : size(x0,2)
    charFunc = zeros(1,N);
    for i = 1 : N
        charFunc(i) = psi_T(v(i),x0(j));
    end
    SimpsonW = 1/3*(3 + (-1).^[1:N] - [1, zeros(1,N-1)]);
    FftFunc = exp(1i*b*v).*charFunc*eta.*SimpsonW;
    payoff = real(fft(FftFunc));
    C_values = exp(-log(K)*alpha)*payoff/pi;
    Ct(j) = C_values(index);
end
if PUTorCALL == 1
    Ct = Ct - S + K*exp(-mu*T);
end
end