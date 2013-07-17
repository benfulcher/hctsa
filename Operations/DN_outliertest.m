function out = DN_outliertest(y,p)
% Removes the pth percent of highest and lowest outliers from the z-scored 
% input time series, y (i.e., 2*p percent removed from the time series in total)
% Returns the ratio of the mean and std before and after doing this
% Ben Fulcher, 2009

% mean of the middle (100-2*p)% of the data
out.mean = mean(y(y > prctile(y,p) & y < prctile(y,100-p)));

% std of the middle (100-2*p)% ofthe data
out.std = std(y(y > prctile(y,p) & y < prctile(y,100-p))) / std(y); % [although std(y) should be 1]

end