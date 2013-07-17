function z = BF_zscore(x)
% Avoids using a Statistics Toolbox licence to do zscores
% Ben Fulcher
	z = (x-mean(x)) / std(x);
end