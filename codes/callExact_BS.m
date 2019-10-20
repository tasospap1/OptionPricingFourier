function [Ct] = callExact_BS(K,r,sigma,S,t,T)
% Documentation 
d1 = ( log(S/K) + (r+sigma^2/2)*(T-t) ) / sigma / sqrt(T-t); 
d2 = d1 - sigma * sqrt(T-t); 
Ct = -K*exp(-r*(T-t)) .* normcdf(d2) + S.*normcdf(d1); 
end

