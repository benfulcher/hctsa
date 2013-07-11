function out = CO_fmac(y)
% Outputs the first minimum of the linear autocorrelation function using CO_autocorr
% Really badly coded by Ben Fulcher, 2008

N = length(y);

a = zeros(N-1,1);
for i = 0:N-1
    a(i+1) = CO_autocorr(y,i);
    if (i > 1) && (a(i-1)-a(i) > 0) && (a(i)-a(i+1) < 0); % minimum
        out = i-1;
        return
    end
end

% If no minimum across all lags, output the length of the time series
out = N;

end