function [Price] = PriceOption_Heston_COS(mu, ...
           K, T, S, PutOrCall, lambda, eta, u_bar, u0, rho)
       
       N = 2^12; % Truncation of Fourier series
       
       c1 = mu*T + (1-exp(-lambda*T)) * (u_bar-u0)/2/lambda - 1/2*u_bar*T ;
       c2 = 1/8/lambda^3*(eta*T*lambda*exp(-lambda*T)*(u0-u_bar)*(8*lambda*rho -4*eta)...
           + lambda*rho*eta*(1-exp(-lambda*T))*(16*u_bar - 8*u0) + ...
           2*u_bar*lambda*T*(-4*lambda*rho*eta + eta^2 + 4*lambda^2) + ...
           eta^2*((u_bar-2*u0)*exp(-2*lambda*T) + u_bar*(6*exp(-lambda*T) - 7) + 2*u0) + ...
           8*lambda^2*(u0-u_bar)*(1-exp(-lambda*T)));
       
       a = c1 - 12 * sqrt(abs(c2));
       b = c1 + 12 * sqrt(abs(c2));
       
       ksi_k = @(k) ( 1/(1 +(k*pi/(b-a))^2) * ( cos(-k*pi*a/(b-a)) ...
           - exp(a) + k*pi/(b-a)*sin(-k*pi*a/(b-a)) ) );
       
       psi_k = @(k) (((b-a)/k/pi * sin(-k*pi*a/(b-a))).*(k>0) );
       
       VPut_k = @(k) ( 2*K/(b-a) * (-ksi_k(k) + psi_k(k)) );
       
       D = @(w) ( sqrt((lambda-1i*rho*eta*w)^2 + (w^2 + 1i*w) * eta^2) );
       G = @(w) ( (lambda - 1i*rho*eta*w - D(w)) ...
           / (lambda - 1i*rho*eta*w + D(w)) );
       
       phiHes = @(w) ( exp(1i*w*mu*T + u0/eta^2 *((1-exp(-D(w)*T))/(1-G(w) ...
           *exp(-D(w)*T)))*(lambda-1i*rho*eta*w-D(w))) * ...
           exp(lambda*u_bar/eta^2*(T*(lambda-1i*rho*eta*w-D(w))- ...
           2*log((1-G(w)*exp(-D(w)*T))/(1-G(w))))) );
       
       x = log(S'/K); % Log-price = log(S/K)
       
       Pt = 1/2 * real(phiHes(0) .* exp(zeros(length(S),1))) ...
           *(2*K/(b-a)*(-ksi_k(0)-a));
%        Delta = zeros(length(S),1);
%        Gamma = zeros(length(S),1);
       for k=1:N-1
           Pt = Pt + ...
               real(phiHes(k*pi/(b-a)) .* exp(1i*k*pi*(x-a)/(b-a))).*VPut_k(k);
%            Delta = Delta + ...
%                real(phiHes(k*pi/(b-a)) .* exp(1i*k*pi*(x-a)/(b-a)) .* 1i*k*pi/(b-a)).*VPut_k(k)./S';
%            Gamma = Gamma + ...
%                real(phiHes(k*pi/(b-a)) .* exp(1i*k*pi*(x-a)/(b-a)) .*(-1i-1)*k*pi/(b-a) ).*VPut_k(k)./S.^2';
       end
       Pt = Pt * exp(-mu*T) ;
       
       if PutOrCall==1
           Price = Pt;
       elseif PutOrCall==-1
           Price = Pt + S' - K*exp(-mu*T);
       end

end

