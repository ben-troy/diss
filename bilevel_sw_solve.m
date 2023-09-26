%wi = sw_ll(taui, n, m, mu, B, rho, covar, theta, esg);
%[tau_i, y_i, Fy] = sw_ul(C, taui, n, m, wi, esg, B, mu, epsilon, rho, covar, theta);

taui = zeros(m,1);
epsilon = 0.1;
iter = 0;
Fyvec = [];
while iter < 3
    wi = sw_ll(taui, n, m, mu, B, rho, covar, theta, esg);
    [taui, y_i, Fy] = sw_ul(taui, wi, n, m, esg, B, mu, epsilon, rho, covar, theta);
    Fyvec = [Fyvec Fy];
    Fy_prev = Fy;
    iter = iter +1;
end

figure;
plot(Fyvec,LineWidth=1.5);
title('ESG Risk of the Firm`s investments');
xlabel('Iteration');
ylabel('F(y)');


risk_sw_hat = [y_i(1:29)'*covar*y_i(1:29), y_i(30:58)'*covar*y_i(30:58), y_i(59:87)'*covar*y_i(59:87),y_i(88:116)'*covar*y_i(88:116)];
esg_sw_hat = [esg*y_i(1:29), esg*y_i(30:58), esg*y_i(59:87),esg*y_i(88:116)];
sum_sw_hat = y_i(1:29)+ ysw_i(30:58) + y_i(59:87)+ y_i(88:116);
TC_star =[y_i(30:58)'*theta*sum_sw_hat, y_i(59:87)'*theta*sum_sw_hat, y_i(88:116)'*theta*sum_sw_hat];

