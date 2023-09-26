returns = readmatrix("dow_returns.csv");
esg = readmatrix("dow_esg.csv");
tickers = ["AXP", "AMGN",  "AAPL", 'BA', 'CAT', 'CSCO', 'CVX', 'GS','HD', 'HON', 'IBM', 'INTC', 'JNJ', 'KO', 'JPM', 'MCD', 'MMM', 'MRK' 'MSFT', 'NKE', 'PG', 'TRV', 'UNH', 'CRM', 'VZ', 'V', 'WBA', 'WMT', 'DIS'];
esg_range = [];
for i = 1:30
    esg_range = [ esg_range, 25 - 0.25*i];
end

%dowcov = covarsolve(returns, esg, tickers, esg_range);

%dowomega = omegasolve(returns, esg, tickers, esg_range);

betavec = [0.9, 0.95, 0.99];
for i = 1:3
   beta = betavec(i);
   cvarvec = cvarsolve(returns, esg, tickers, esg_range, beta);
end

