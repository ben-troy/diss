function omegavec = omegasolve(returns, esg, tickers, esg_range)

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
omegavec = [];
portfolio = [];

for i=1:length(esg_range)
    cvx_begin quiet
        variables s(n) q(k) z  
        maximise (mu*s - 1.0005*z)  
        subject to
            for j = 1:k
                q(j) >= 0;
                q(j) >= 1.0005*z - returns(j,:)*s;
            end
            sum(q) == 1;
            sum(s) == z;
            %mu*s >= z*1.0005;
            esg*s <= z*esg_range(i);
            z >= 0;
            for j = 1:n
                s(j) <= z*0.2;
                s(j) >= 0;
            end
    cvx_end
    x = (1/z)*(s);
    retvec= [retvec mu*x];
    portfolio = [portfolio x];
    denom = 0;
    for i = 1:k
        att = 1.0005 - returns(i,:)*x;
        if att > 0
            new = att;
        else
            new = 0;
        end 
        denom = denom + new;
    end
    om = k*(mu*x - 1.0005)/(denom);
    omegavec = [omegavec om];
end

figure;
plot(esg_range, omegavec, LineWidth=1)
title('Effect of ESG Consideration on Risk Management: Omega Ratio');
xlabel('Permitted ESG risk');
ylabel('Omega Ratio');

end