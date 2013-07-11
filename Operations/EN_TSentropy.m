function entr = TSentropy(x, q)
%   TSentropy estimates the entropy of signals
%   entr      : The entropy estimate
%   q           : input parameter, q >=1; 
%   x           : The time series to be analyzed
%   q           : Tsallis non-extensive parameter value, q >= 1;
%                 if q == 1 then Tsallis' entropy concides with Shannon's
%   http://www.TSResearchGroup.com
%   D. Tolstonogov
%   Copyright (c) by Trade Smart Research
%   08/04/2004 

[NRow, NCol] = size(x);

if NRow ~= 1
    x = x';
end

% The number of cells of the histogram   
minx = min(x);
maxx = max(x);
delta = (maxx-minx)/(length(x)-1);
ncell = ceil(sqrt(length(x)));

% Histogram
h = hist(x, ncell); 
entr = 0;
count = 0;

% Shannon's entropy
if q == 1
    for n = 1:ncell
        if h(n)~=0 
            hn = log(h(n));
        else
            hn = 0;
        end
        count = count+h(n);
        entr = entr-h(n)*hn;
    end
    entr = entr/count+log(count); 
else
% Tsallis entropy
    for n = 1:ncell
        hn = (h(n))^q;        
        count = count+h(n);
        entr = entr+hn;
    end
    entr = entr/(count)^q;
    entr = (1-entr)/(q-1);
end

end