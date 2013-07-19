% EN_pdfent
% 
% Estimates the (log_2) Shannon entropy from the probability distribution of the time
% series.
% 
% INPUT:
% y, a time series

function out = EN_pdfent(y)
% Ben Fulcher, 2009

N = length(y); % time-series length
tmp = sort(y);
tmp = diff(tmp); % incremental differences of sorted set of values in y
if size(tmp,1) < size(tmp,2)
    tmp = tmp';
end

pdf = diff([0; find(tmp > 0); N])/N; % Empirical probability distribution
out = -sum(pdf.*log2(pdf));

end