function [res] = calibration(mu, lambda, eta, u_bar, rho) 

data = readtable('data.csv');
date_start = datenum(num2str(data.date, '%d'), 'yyyymmdd');
date_end = datenum(num2str(data.exdate, '%d'), 'yyyymmdd');
T = (date_end - date_start)/365;
K = data.strike_price/1000;
S = data.forward_price;
v0 = data.impl_volatility;
price = 0.5*abs(data.best_bid - data.best_offer);

N = size(price,1);
diff_price = zeros(1,N);

for i=1:1:N
    diff_price(i) = price(i) - PriceOption_Heston_COS(mu, K(i), T(i), S(i), -1, ...
        lambda, eta, u_bar, v0(i), rho);
    res = sum(diff_price.^2);
end

end