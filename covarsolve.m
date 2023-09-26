function covarvec = covarsolve(returns, esg, tickers, esg_range)

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
covarvec = [];
portfolio = [];

for i=1:length(esg_range)
    cvx_begin quiet
        variables x(n) 
        minimise (x'*covar*x)  
        subject to
            sum(x) == 1;
            mu*x >= 1.0005;
            esg*x <= esg_range(i);
            x >= 0;
            for j = 1:n
                x(j) <= 0.2;
            end
    cvx_end
    retvec= [retvec mu*x];
    portfolio = [portfolio x];
    covarvec = [covarvec x'*covar*x];
end

figure;
plot(esg_range, covarvec, LineWidth=1)
title('Effect of ESG Consideration on Risk Management: Portfolio Variance');
xlabel('Permitted ESG risk');
ylabel('Portfolio Variance');

end