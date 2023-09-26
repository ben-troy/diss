function [tau_i, y_i, Fy] = sw_ul(taui, wi, n, m, esg, B, mu, epsilon, rho, covar, theta)

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
        -(B(1).*mu*y(1:29) - B(1)^2.*rho(1).*y(1:29)'*covar*y(1:29) - B(1).*taui(1)*esg*y(1:29)  - B(1).*B(1).*y(1:29)'*theta*y(1:29)) - (B(2).*mu*y(30:58) - B(2)^2.*rho(2).*y(30:58)'*covar*y(30:58) - B(2).*taui(2)*esg*y(30:58)  - B(2).*B(2).*y(30:58)'*theta*y(30:58)) - (B(3).*mu*y(59:87) - B(3)^2.*rho(3).*y(59:87)'*covar*y(59:87) - B(3).*taui(3)*esg*y(59:87)  - B(3).*B(3).*y(59:87)'*theta*y(59:87)) - (B(4).*mu*y(88:116) - B(4)^2.*rho(4).*y(88:116)'*covar*y(88:116) - B(4).*taui(4)*esg*y(88:116)  - B(4).*B(4).*y(88:116)'*theta*y(88:116)) <= -(B(1).*mu*wi(1:29) - B(1)^2.*rho(1).*wi(1:29)'*covar*wi(1:29) - B(1).*taui(1)*esg*wi(1:29)  - B(1).*B(1).*wi(1:29)'*theta*wi(1:29)) - (B(2).*mu*wi(30:58) - B(2)^2.*rho(2).*wi(30:58)'*covar*wi(30:58) - B(2).*taui(2)*esg*wi(30:58)  - B(2).*B(2).*wi(30:58)'*theta*wi(30:58)) - (B(3).*mu*wi(59:87) - B(3)^2.*rho(3).*y(59:87)'*covar*wi(59:87) - B(3).*taui(3)*esg*wi(59:87)  - B(3).*B(3).*wi(59:87)'*theta*wi(59:87)) - (B(4).*mu*wi(88:116) - B(4)^2.*rho(4).*wi(88:116)'*covar*wi(88:116) - B(4).*t(4)*esg*wi(88:116)  - B(4).*B(4).*wi(88:116)'*theta*wi(88:116));
        esg*(B(1).*y(1:29) + B(2).*y(30:58) + B(3).*y(59:87) + B(4).*y(88:116)) == 
cvx_end

y_i = y;
tau_i = t;
Fy = esg*(B(1).*y(1:29) + B(2).*y(30:58) + B(3).*y(59:87) + B(4).*y(88:116));