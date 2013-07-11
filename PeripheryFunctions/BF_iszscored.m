function iszscored = BF_iszscored(x)
% Given input time series, x, does a quick heuristic to see whether
% (numerically, to within eps) z-scored
% Ben Fulcher 2013

iszscored = ((abs(mean(x)) > eps) || (abs(std(x)-1) > eps));

end