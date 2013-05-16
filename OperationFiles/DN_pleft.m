function out=DN_pleft(y,th)
p = quantile(abs(y-mean(y)),1-th);
% a proportion th of the data lie outside p of the mean
out = p/std(y);
end