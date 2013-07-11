function out = EN_pdfent(y)
% Input time series, y
% Calculates the (log_2) entropy of the probability distribution
% Ben Fulcher, 2009

N = length(y);
tmp = sort(y);
tmp = diff(tmp); % incremental differences of sorted set of values in y
if size(tmp,1) < size(tmp,2)
    tmp = tmp';
end
pdf = diff([0; find(tmp > 0); N])/N; % Empirical probability distribution
out = -sum(pdf.*log2(pdf));

end