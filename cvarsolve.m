function cvarvec = cvarsolve(returns, esg, tickers, esg_range, beta)
k = size(returns,1);
mu = mean(returns);
n = size(mu,2);
covar = corrcoef(returns);

%covar = findcovar_mean(returns); % if scenarios have uneven prob
b = 1;
rho = 0;
iter = 20;
rhovec = zeros(1, iter);
v = 1/n*ones(1, n)';

retvec = [];
cvarvec = [];
portfolio = [];
factor = 1/(1-beta);
prob = 1/k * ones(1,k);

for i=1:length(esg_range)
    cvx_begin quiet
        variable y(n) nonnegative;
        variable z(k) nonnegative;
        variable t nonnegative;
        minimise(t + prob*z/(1 - beta));  
        subject to
            for j = 1:k
                z(k) + t >= 0;
                z(k) + t + returns(j,:)*y >= 1.0005;
            end
            sum(y) == 1;
            mu*y >= 1.0005;
            esg*y <= esg_range(i);
            for j = 1:n
                y(j) <= 0.2;
            end
    cvx_end
    obj = t + prob*z/(1 - beta);
    retvec= [retvec mu*y];
    portfolio = [portfolio y]; 
    cvarvec = [cvarvec obj];
end

figure;
plot(esg_range, cvarvec, LineWidth=1)
xlabel('Permitted ESG risk');
ylabel('CV@R');
end