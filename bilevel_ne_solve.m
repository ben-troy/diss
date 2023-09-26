returns = readmatrix("dow_returns.csv");
esg = readmatrix("dow_esg.csv");
vol = readmatrix('dow_volumes.csv');
tickers = ["AXP", "AMGN",  "AAPL", 'BA', 'CAT', 'CSCO', 'CVX', 'GS','HD', 'HON', 'IBM', 'INTC', 'JNJ', 'KO', 'JPM', 'MCD', 'MMM', 'MRK' 'MSFT', 'NKE', 'PG', 'TRV', 'UNH', 'CRM', 'VZ', 'V', 'WBA', 'WMT', 'DIS'];

m = 4; % number of accounts

k = size(returns,1); % number of scenarios
n = size(returns,2); % number of assets
dim = n*m;

mu = mean(returns); % mean return for each asset
covar = corrcoef(returns); % covariance matrix for the assest
vol_bar = mean(vol); % average volume of each asset traded

% Calculate market impact matrix
theta = zeros(n,n);
for i = 1:n
    for j = 1:n
        num = (vol(:,i)- vol_bar(i)).'*(vol(:,j) - vol_bar(j));
        den = sqrt((vol(:,i)- vol_bar(i)).'*(vol(:,i) - vol_bar(i))) * sqrt((vol(:,j)- vol_bar(j)).'*(vol(:,j) - vol_bar(j)));
        theta(i,j) = num/den;
    end
end

% Generate risk aversion parameters and budgets
rng('default')
rng(42)
s = rng;
B = 10*rand(m, 1);
rho = 50*rand(m,1);

% generate Q, d, and C
Q = zeros(dim,dim); d = zeros(dim,1); C = zeros(m, dim);

for i = 1:m
    for j = 1:m
        is = (i-1)*n + 1;
        ifin = i*n;
        js = (j-1)*n + 1;
        jfin = j*n;
        if i == j
            Q(is:ifin, js:jfin) = B(i)^2.*(rho(i).*covar + 2.*theta);
        else
            Q(is:ifin, js:jfin) = B(i).*B(j).*theta;
        end
    end
    d(is:ifin) = -B(i).*mu;
    C(i, is:ifin) = B(i).*esg';
end

taui = zeros(m,1);

beta = max(B).*esg*esg';
epsilon = 0.01;
iter = 0;
Fyvec = [];
while iter < 200
    wi = nash_ll(Q, d, C, taui, n, m, mu);
    [taui, y_i, Fy] = nash_ul(Q, d, C, taui, n, m, beta, wi, esg, B, mu, epsilon);
    if iter == 0
        Fy_init = Fy;
        y_init = y_i;
        tau_init = taui;
    end
    Fyvec = [Fyvec Fy];
    Fy_prev = Fy;
    iter = iter +1;
end

figure;
plot(Fyvec,LineWidth=1.5);
title('ESG Risk of the Firm`s investments');
xlabel('Iteration');
ylabel('F(y)');

risk_ne_0 = [y_init(1:29)'*covar*y_init(1:29), y_init(30:58)'*covar*y_init(30:58), y_init(59:87)'*covar*y_init(59:87),y_init(88:116)'*covar*y_init(88:116)];
risk_ne_star = [y_i(1:29)'*covar*y_i(1:29), y_i(30:58)'*covar*y_i(30:58), y_i(59:87)'*covar*y_i(59:87),y_i(88:116)'*covar*y_i(88:116)];
esg_ne_0 = [esg*y_init(1:29), esg*y_init(30:58), esg*y_init(59:87),esg*y_init(88:116)];
esg_ne_star = [esg*y_i(1:29), esg*y_i(30:58), esg*y_i(59:87),esg*y_i(88:116)];
sum_ne_0 = y_init(1:29)+ y_init(30:58) + y_init(59:87)+ y_init(88:116);
TC_ne_0 = [y_init(1:29)'*theta*sum_ne_0, y_init(30:58)'*theta*sum_ne_0, y_init(59:87)'*theta*sum_ne_0, y_init(88:116)'*theta*sum_ne_0];
sum_ne_star = y_i(1:29)+ y_i(30:58) + y_i(59:87)+ y_i(88:116);
TC_ne_star =[y_i(1:29)'*theta*sum_ne_star, y_i(30:58)'*theta*sum_ne_star, y_i(59:87)'*theta*sum_ne_star, y_i(88:116)'*theta*sum_ne_star];

ysw_0 = sw_ll(tau_init, n, m, mu, B, rho, covar, theta,esg);
ysw_i = sw_ll(taui, n, m, mu, B, rho, covar, theta,esg);

risk_sw_0 = [ysw_0(1:29)'*covar*ysw_0(1:29), ysw_0(30:58)'*covar*ysw_0(30:58), ysw_0(59:87)'*covar*ysw_0(59:87),ysw_0(88:116)'*covar*ysw_0(88:116)];
risk_sw_star = [ysw_i(1:29)'*covar*ysw_i(1:29), ysw_i(30:58)'*covar*ysw_i(30:58), ysw_i(59:87)'*covar*ysw_i(59:87),ysw_i(88:116)'*covar*ysw_i(88:116)];
esg_sw_0 = [esg*ysw_0(1:29), esg*ysw_0(30:58), esg*ysw_0(59:87),esg*ysw_0(88:116)];
esg_sw_star = [esg*ysw_i(1:29), esg*ysw_i(30:58), esg*ysw_i(59:87),esg*ysw_i(88:116)];
sum_sw_0 = ysw_0(1:29)+ ysw_0(30:58) + ysw_0(59:87)+ ysw_0(88:116);
TC_sw_0 = [ysw_0(1:29)'*theta*sum_sw_0, ysw_0(30:58)'*theta*sum_sw_0, ysw_0(59:87)'*theta*sum_sw_0, ysw_0(88:116)'*theta*sum_sw_0];
sum_sw_star = ysw_i(1:29)+ ysw_i(30:58) + ysw_i(59:87)+ ysw_i(88:116);
TC_sw_star =[ysw_i(1:29)'*theta*sum_sw_star, ysw_i(30:58)'*theta*sum_sw_star, ysw_i(59:87)'*theta*sum_sw_star, ysw_i(88:116)'*theta*sum_sw_star];

ESG_total =  esg*(B(1).*y_i(1:29) + B(2).*y_i(30:58) + B(3).*y_i(59:87) + B(4).*y_i(88:116));

y_hat = free_tau(ESG_total, taui, n, m, mu, B, rho, covar, theta, esg);

risk_sw_hat = [y_hat(1:29)'*covar*y_i(1:29), y_hat(30:58)'*covar*y_hat(30:58), y_hat(59:87)'*covar*y_hat(59:87),y_hat(88:116)'*covar*y_hat(88:116)];
esg_sw_hat = [esg*y_hat(1:29), esg*y_hat(30:58), esg*y_hat(59:87),esg*y_hat(88:116)];
sum_sw_hat = y_hat(1:29)+ y_hat(30:58) + y_hat(59:87)+ y_hat(88:116);
TC_hat =[y_hat(1:29)'*theta*sum_sw_hat, y_hat(30:58)'*theta*sum_sw_hat, y_hat(59:87)'*theta*sum_sw_hat, y_hat(88:116)'*theta*sum_sw_hat];

y_risk = [risk_ne_star(1), risk_sw_hat(1); risk_ne_star(2), risk_sw_hat(2); risk_ne_star(3), risk_sw_hat(3); risk_ne_star(4), risk_sw_hat(4)];
y_costs = [TC_ne_star(1), TC_hat(1); TC_ne_star(2), TC_hat(2); TC_ne_star(3), TC_hat(3); TC_ne_star(4), TC_hat(4)];
y_esg = [esg_ne_star(1), esg_sw_hat(1); esg_ne_star(2), esg_sw_hat(2); esg_ne_star(3), esg_sw_hat(3); esg_ne_star(4), esg_sw_hat(4)];


figure;
bar(y_risk);
xlabel('Account');
ylabel('Portfolio Variance');
legend('Nash', 'Social Welfare');

figure;
bar(y_costs);
xlabel('Account');
ylabel('Transaction Costs');
legend('Nash', 'Social Welfare');

figure;
bar(y_esg);
xlabel('Account');
ylabel('Portfolio ESG Risk');
legend('Nash', 'Social Welfare');

strt = []; fin = [];
for i = 1:m
    strt = [strt ((i-1)*n + 1) ];
    fin = [fin i*n] ;
end

U_ne = zeros(m,1);
U_sw = zeros(m,1);
for i = 1:m
    U_ne(i) = B(i)*1.0005 - B(i)^2.*rho(i).*risk_ne_star(i) - - B(i)*taui(i)*esg_ne_star(i) - B(i).*B(i).*TC_ne_star(i);
    U_sw(i) = B(i)*1.0005 - B(i)^2.*rho(i).*risk_sw_hat(i) - - B(i)*taui(i)*esg_sw_hat(i) - B(i).*B(i).*TC_hat(i);
end

y_utility = [U_ne(1), U_sw(1); U_ne(2), U_sw(2); U_ne(3), U_sw(3); U_ne(4), U_sw(4)];
figure;
bar(y_utility);
xlabel('Account');
ylabel('Account Utility');
legend('Nash', 'Social Welfare');

