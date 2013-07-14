function out = SY_drifting(y,howl,l)
% A function inspired by this bit of advice from a MATLAB forum:
%> http://www.mathworks.de/matlabcentral/newsreader/view_thread/136539
%> It seems to me that you are looking for a measure for a drifting mean.
%> If so, this is what I would try:
%> 
%> - Decide on a frame length N
%> - Split your signal in a number of frames of length N
%> - Compute the means of each frame
%> - Compute the variance for each frame
%> - Compare the ratio of maximum and minimum mean
%>   with the mean variance of the frames.
%> 
%> Rune
%
% Ben Fulcher, 2009

% INPUTS:
% howl can be either 'fix' (this given length; for l) or 'num' (this number
% of windows)

N = length(y); % length of the input time series

% Check inputs
if strcmp(howl,'num')
    l = floor(N/l);
else
    if ~strcmp(howl,'fix')
        return
    end
end

% ++BF 19/3/2010
if N < l % doesn't make sense to split into more windows than there are data points
    fprintf(1,'Time Series (N = %u < l = %u) is too short for this operation\n',N,l)
    out = NaN; return
end

% Get going
nfits = floor(N/l);
z = zeros(l,nfits);
for i = 1:nfits; % number of times l fits completely into N
    z(:,i) = y((i-1)*l+1:i*l);
end
zm = mean(z);
zv = var(z);
meanvar = mean(zv);
maxmean = max(zm);
minmean = min(zm);
meanmean = mean(zm);

out.max = maxmean/meanvar;
out.min = minmean/meanvar;
out.mean = meanmean/meanvar;
out.meanmaxmin = (out.max+out.min)/2;
out.meanabsmaxmin = (abs(out.max)+abs(out.min))/2;

end