dowtickers = ["AXP", "AMGN",  "AAPL", 'BA', 'CAT', 'CSCO', 'CVX', 'GS','HD', 'HON', 'IBM', 'INTC', 'JNJ', 'KO', 'JPM', 'MCD', 'MMM', 'MRK' 'MSFT', 'NKE', 'PG', 'TRV', 'UNH', 'CRM', 'VZ', 'V', 'WBA', 'WMT', 'DIS'];

%% Initialise matrices to contain return and volume data
ret = zeros(1007,29);
vol = zeros(1007,29);

%% For each asset, add a column of historic returns and trading volume
for i = 1:length(dowtickers)
    data = getMarketDataViaYahoo(dowtickers{i}, '1-Jan-2019', '31-Dec-2022', '1d');
    ret(:,i) = get_returns(data);
    vol(:,i) = get_volume(data);
end

%% Write data to CSV files

%writematrix(ret,'dow_returns.csv')
%writematrix(vol, 'dow_volumes.csv')

%% ESG ratings were copied manually into their csv