global price; global T; global S; global K; global mu; 

data = readtable('data.csv');
date_start = datenum(num2str(data.date, '%d'), 'yyyymmdd');
date_end = datenum(num2str(data.exdate, '%d'), 'yyyymmdd');
T = (date_end - date_start)/365;
K = data.strike_price/1000;
S = data.forward_price;
mu = 0.1; 
%v0 = data.impl_volatility;
price = 0.5*abs(data.best_bid - data.best_offer);


N = size(price,1);
diff_price = zeros(1,N);
x0 = [1.5768 0.5751 0.0898 -0.5711 0.15]; 
lb = [0.00001 0.00001 0.00001 -1 0]; 
ub = [6 1 1 0 0.5]; 
res = fmincon(@fun, x0, [], [], [], [], lb, ub);
minval = fun(res); 

function  [res] = fun(x)
global price; global T; global S; global K; global mu; 
N = 100; 
objective = zeros(1,N);
x0 = [3 0.5751 0.0998 -0.5711 0.15]; 
for i=1:1:N
    objective(i) = (PriceOption_Heston_COS(mu, K(i), T(i), S(i), -1, ...
        x(1), x(2), x(3), x(5), x(4)) - price(i)) / price(i) + ...
        norm(x-x0)^2;
end
res = sum(objective.^2);
end


