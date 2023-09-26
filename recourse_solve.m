function [y, z] = recourse_solve(mu, covar, theta, rho, B, phi, b, m, n, returns)

k = size(returns,1);
dim = n*m;
strt = []; fin = [];
for i = 1:m
    strt = [strt ((i-1)*n + 1) ];
    fin = [fin i*n] ;
end
cvx_begin quiet
    variable z(k) nonnegative;
    variable y(dim) nonnegative;
    minimise (B(1)^2.*rho(1).*y(1:29)'*covar*y(1:29)+  B(1).*B(1).*y(1:29)'*theta*y(1:29) + B(2)^2.*rho(2).*y(30:58)'*covar*y(30:58) + B(2).*B(2).*y(30:58)'*theta*y(30:58) +  B(3)^2.*rho(3).*y(59:87)'*covar*y(59:87) + B(3).*B(3).*y(59:87)'*theta*y(59:87) + B(4)^2.*rho(4).*y(88:116)'*covar*y(88:116) + B(4).*B(4).*y(88:116)'*theta*y(88:116) + b/k*sum(z));   
    subject to
        for i = 1:m
            sum(y(strt(i):fin(i))) == B(i);
            mu*y(strt(i):fin(i)) >= 1.0005*B(i);
            for j = strt(i):fin(i)
                y(j) <= 0.2*B(i);
            end
        end
        for i = 1:k
            returns_all = [returns(i,:)';returns(i,:)';returns(i,:)';returns(i,:)'];
            (phi.*returns_all)'*y + z(k) >= 0.7*sum(B);
        end
cvx_end
    
end