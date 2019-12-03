% HESTON_CALL computes the call price for the Heston model
%
%     [P] = heston_call(S,y,T,K,rho,k,m,nu,lambda))
%     
%     Input:    S      ...  stock price at time 0
%               y      ...  volatility at time 0
%               T      ...  maturity
%               K      ...  strike
%               rho    ...  correlation
%               k      ...  rate of mean reversion
%               m      ...  level of mean reversion
%               nu     ...  volatility of volatility
%               lambda ...  price of volatility risk
%
%                                       
%      Output: P       ...  call price

function P = heston_call(S,y,T,K,rho,k,m,nu,lambda)

% parameters
n = 1000;
left = 1000; 

% define certain constants
x = log(S); 
[phi,w] = gauleg(n);  % the integration points and weights

% integrate
phi=left/2*(phi+1);   % the element mapping

% initialize
P = zeros(length(y),length(x));

for r = 1:length(x)
    for j = 1:length(y)
        [I1,I2] = hestonintegrand(phi,x(r),y(j),T,K,rho,k,m,...
                nu,lambda);
        P(j,r) = S(r)*(1/2+1/pi*left/2*sum(w.*I1))-...
                K*(1/2+1/pi*left/2*sum(w.*I2));
    end
end

return 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [I1,I2] = hestonintegrand(phi,x,y,t,K,rho,k,m,nu,lambda)

T = 1; tau = T-t;

% define certain constants
u1 = 1/2; a = k*m; b1 = k+lambda-rho*nu; 
u2 = -1/2; a = k*m; b2 = k+lambda; r = 0;

% define helper functions
d1 = sqrt((rho*nu*i*phi-b1).^2-nu^2*(2*u1*i*phi-phi.^2));
d2 = sqrt((rho*nu*i*phi-b2).^2-nu^2*(2*u2*i*phi-phi.^2));
g1 = (b1-rho*nu*phi*i+d1)./(b1-rho*nu*phi*i-d1);
g2 = (b2-rho*nu*phi*i+d2)./(b2-rho*nu*phi*i-d2);
C1 = r*i*phi*tau+a/nu^2*((b1-rho*nu*i*phi+d1)*tau-...
    2*log((1-g1.*exp(d1*tau))./(1-g1)));
C2 = r*i*phi*tau+a/nu^2*((b2-rho*nu*i*phi+d2)*tau-...
    2*log((1-g2.*exp(d2*tau))./(1-g2)));
D1 = 1/nu^2*(b1-rho*nu*i*phi+d1).*((1-exp(d1*tau))./...
    (1-g1.*exp(d1*tau)));
D2 = 1/nu^2*(b2-rho*nu*i*phi+d2).*((1-exp(d2*tau))./...
    (1-g2.*exp(d2*tau)));

% define the characteristic functions
f1 = exp(C1+D1*y+i*x*phi); f2 = exp(C2+D2*y+i*x*phi); 

% the integrands
I1 = real(exp(-i*phi*log(K)).*f1./(i*phi)); 
I2 = real(exp(-i*phi*log(K)).*f2./(i*phi));

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [x,w] = gauleg(n)
 
x = zeros(n,1); w=zeros(n,1);
m = (n+1)/2; xm = 0.0; xl = 1.0;
for i = 1:m
    z = cos(pi*(i-0.25)/(n+0.5));
    while 1
      p1 = 1.0; p2 = 0.0;
      for j = 1:n
        p3 = p2; p2 = p1; p1 = ((2.0*j-1.0)*z*p2-(j-1.0)*p3)/j;
      end
      pp = n*(z*p1-p2)/(z*z-1.0); z1 = z; z = z1-p1/pp;
      if (abs(z-z1)<eps), break, end
    end
    x(i) = xm-xl*z; x(n+1-i) = xm+xl*z;
    w(i) = 2.0*xl/((1.0-z*z)*pp*pp); w(n+1-i) = w(i);
end
x = x'; w = w';
return
