global price; 
global T; 
global S; 
global K; 
global v0; 

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
x0 = [0.01 1.5768 0.5751 0.0398 -0.5711]; 
lb = [0 0 0 0 -1]; 
ub = [1 4 1 1 0]; 
res = lsqnonlin(@Differences, x0, lb, ub);

    function  [res] = Differences(x)
    global price;
    global T;
    global S;
    global K;
    global v0; 
        for i=1:1:100
            diff_price(i) = price(i) - PriceOption_Heston_COS(x(1), K(i), T(i), S(i), -1, ...
                x(2), x(3), x(4), v0(i), x(5));
            i
        end
        res = diff_price'; 
    end


