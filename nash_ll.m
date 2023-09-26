function wi = nash_ll(Q, d, C, taui, n, m, mu)

dim = n*m;
strt = []; fin = [];
for i = 1:m
    strt = [strt ((i-1)*n + 1) ];
    fin = [fin i*n] ;
end

cvx_begin quiet
    variable y(dim) nonnegative
    minimise (1/2*y.'*Q*y + d.'*y + taui.'*C*y);
    subject to
        for i = 1:m
            mu*y(strt(i):fin(i)) >= 1.0005;
            sum(y(strt(i):fin(i))) == 1;
        end
        for j = 1:dim
            y(j) <= 0.2;
        end
cvx_end
opt_prev = 1/2*y.'*Q*y + d.'*y + taui.'*C*y;
wi = y;