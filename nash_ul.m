function [tau_i, y_i, Fy] = nash_ul(Q, d, C, taui, n, m, beta, wi, esg, B, mu, epsilon)

dim = n*m;
strt = []; fin = [];
for i = 1:m
    strt = [strt ((i-1)*n + 1) ];
    fin = [fin i*n] ;
end

cvx_begin quiet
    variable y(dim) nonnegative;
    variable t(m) nonnegative;
    minimise (esg*(B(1).*y(1:29) + B(2).*y(30:58) + B(3).*y(59:87) + B(4).*y(88:116)))
    subject to
        for i = 1:m
            mu*y(strt(i):fin(i)) >= 1.0005;
            sum(y(strt(i):fin(i))) == 1;
        end
        for j = 1:dim
            y(j) <= 0.2;
        end
        1/2*y'*Q*y + d'*y  + 0.5*beta.*t'*t <= 1/2*wi'*Q*wi + d'*wi + t'*C*wi + 0.5*beta.*taui'*taui + (C*wi)'*t + (beta.*taui)'*t - (C*wi)'*taui - (beta.*taui)'*taui + epsilon;
cvx_end

y_i = y;
tau_i = t;
Fy = esg*(B(1).*y(1:29) + B(2).*y(30:58) + B(3).*y(59:87) + B(4).*y(88:116));
