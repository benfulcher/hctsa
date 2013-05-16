function z = benzscore(x)
% avoids using a Statistics Toolbox licence
	z = (x-mean(x))/std(x);
end