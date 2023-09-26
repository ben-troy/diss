phi_base = zeros(n,1);
for i = 1:n
    if esg(i) <= 17
        phi_base(i) = 1;
    end
end

phi = [phi_base; phi_base; phi_base; phi_base];
k=size(returns,1);
b = 2*k*max(B)^2;

[x_rec, z1] = recourse_solve(mu, covar, theta, rho, B, phi, b, m, n, returns);

failures17 = 0;
for i = 1:k
    if z1(i) ~= 0
        failures17 = failures17 +1;
    end
end

for i = 1:m
    strt = [strt ((i-1)*n + 1) ];
    fin = [fin i*n] ;
    y_rec(strt(i):fin(i)) = (1/B(i)).*x_rec(strt(i):fin(i)); 
end


risk_rec17 = [y_rec(1:29)'*covar*y_rec(1:29), y_rec(30:58)'*covar*y_rec(30:58), y_rec(59:87)'*covar*y_rec(59:87),y_rec(88:116)'*covar*y_rec(88:116)];
esg_rec17 = [esg*y_rec(1:29), esg*y_rec(30:58), esg*y_rec(59:87),esg*y_rec(88:116)];
sum_rec17 = y_rec(1:29)+ y_rec(30:58) + y_rec(59:87)+ y_rec(88:116);
TC_rec17 =[y_rec(1:29)'*theta*sum_rec, y_rec(30:58)'*theta*sum_rec, y_rec(59:87)'*theta*sum_rec, y_rec(88:116)'*theta*sum_rec];

phi_base = zeros(n,1);
for i = 1:n
    if esg(i) <= 16
        phi_base(i) = 1;
    end
end

phi = [phi_base; phi_base; phi_base; phi_base];
k=size(returns,1);
b = 2*k*max(B)^2;

[x_rec, z2] = recourse_solve(mu, covar, theta, rho, B, phi, b, m, n, returns);

for i = 1:m
    strt = [strt ((i-1)*n + 1) ];
    fin = [fin i*n] ;
    y_rec(strt(i):fin(i)) = (1/B(i)).*x_rec(strt(i):fin(i)); 
end


risk_rec16 = [y_rec(1:29)'*covar*y_rec(1:29), y_rec(30:58)'*covar*y_rec(30:58), y_rec(59:87)'*covar*y_rec(59:87),y_rec(88:116)'*covar*y_rec(88:116)];
esg_rec16 = [esg*y_rec(1:29), esg*y_rec(30:58), esg*y_rec(59:87),esg*y_rec(88:116)];
sum_rec16 = y_rec(1:29)+ y_rec(30:58) + y_rec(59:87)+ y_rec(88:116);
TC_rec16 =[y_rec(1:29)'*theta*sum_rec, y_rec(30:58)'*theta*sum_rec, y_rec(59:87)'*theta*sum_rec, y_rec(88:116)'*theta*sum_rec];

failures16 = 0;
for i = 1:k
    if z2(i) ~= 0
        failures16 = failures16 +1;
    end
end

phi_base = zeros(n,1);
for i = 1:n
    if esg(i) <= 15
        phi_base(i) = 1;
    end
end

[x_rec, z3] = recourse_solve(mu, covar, theta, rho, B, phi, b, m, n, returns);

for i = 1:m
    strt = [strt ((i-1)*n + 1) ];
    fin = [fin i*n] ;
    y_rec(strt(i):fin(i)) = (1/B(i)).*x_rec(strt(i):fin(i)); 
end


risk_rec15 = [y_rec(1:29)'*covar*y_rec(1:29), y_rec(30:58)'*covar*y_rec(30:58), y_rec(59:87)'*covar*y_rec(59:87),y_rec(88:116)'*covar*y_rec(88:116)];
esg_rec15 = [esg*y_rec(1:29), esg*y_rec(30:58), esg*y_rec(59:87),esg*y_rec(88:116)];
sum_rec15 = y_rec(1:29)+ y_rec(30:58) + y_rec(59:87)+ y_rec(88:116);
TC_rec15 =[y_rec(1:29)'*theta*sum_rec, y_rec(30:58)'*theta*sum_rec, y_rec(59:87)'*theta*sum_rec, y_rec(88:116)'*theta*sum_rec];

failures15 = 0;
for i = 1:k
    if z3(i) ~= 0
        failures15 = failures15 +1;
    end
end