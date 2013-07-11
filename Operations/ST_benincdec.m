function out = ST_benincdec(y,binarymeth)
% returns stats on the consecutive (increase/decrease properties) (or other 
% binary transformation) of the time series y.
% Ben Fulcher

switch binarymeth
    case 'diff' % 1 if 
        y = ((sign(diff(y)))+1)/2; % binary signal, equal to one for stepwise increases
    case 'mean' % 1 if above mean, zero otherwise
        y = (sign(y)+1)/2;
    case 'iqr' % 1 if inside interquartile range, 0 otherwise
        iqr = quantile(y,[0.25 0.75]);
        iniqr = find(y>iqr(1) & y<=iqr(2));
        y = zeros(length(y),1);
        y(iniqr) = 1;
end

N = length(y); % length of signal - 1 (difference operation)

pup = sum(y == 1)/N;
pdown = 1 - pup;
p = [pup pdown];

out.pup = pup;
out.pupstat2 = sum(y(floor(end/2)+1:end) == 1)/sum(y(1:floor(end/2)) == 1);

% Shannon entropy
out.h = - sum(p(p>0).*log(p(p>0)));

% longest consecutive string of ones / zeros (normalized by length)
difffy = diff(find([1;y;1]));
stretch0 = difffy(difffy~=1)-1;

difffy = diff(find([0;y;0] == 0));
stretch1 = difffy(difffy~=1)-1;

% pstretches
% number of different stretches as proportion of time series
out.pstretch1 = length(stretch1)/N;
out.pstretch0 = length(stretch0)/N;
out.pstretches = (length(stretch0)+length(stretch1))/N;

if isempty(stretch0) % all 1s (almost never happens)
    out.longstretch0 = 0;
    out.meanstretch0 = 0;
    out.stdstretch0 = NaN;
else
    out.longstretch0 = max(stretch0)/N; % longest consecutive stretch of zeros
    out.meanstretch0 = mean(stretch0)/N; % mean stretch of zeros
    out.stdstretch0 = std(stretch0); % standard deviation of stretch lengths of consecutive zeros
end

if isempty(stretch1) % all zeros (almost never happens)
    out.longstretch1 = 0;
    out.meanstretch1 = 0;
    out.stdstretch1 = NaN;
else
    out.longstretch1 = max(stretch1)/N; % longest consecutive stretch of ones
    out.meanstretch1 = mean(stretch1)/N;
    out.stdstretch1 = std(stretch1);
end

out.meanstretchrat = out.meanstretch1/out.meanstretch0;
out.stdstretchrat = out.stdstretch1/out.stdstretch0;

a=length(find(stretch1 == 1)); b = length(find(stretch1 == 2));
if b>0, out.rat21stretch1=a/b;
else out.rat21stretch1 = NaN; 
end

a=length(find(stretch0 == 1)); b = length(find(stretch0 == 2));
if b>0, out.rat21stretch0 = a/b;
else out.rat21stretch0 = NaN; 
end

% out.rat21stretch0=length(find(stretch0 == 1))/length(find(stretch0 == 2));


end