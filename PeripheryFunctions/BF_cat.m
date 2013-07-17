function out = BF_cat(s,d,surr)
% d is the delimiter
% s is the list of strings to be catinated
% surrounds each string in optional string surr

if nargin < 2
    d = ', '; % default delimiter is a comma
else
    d = [d, ' '];
end

if nargin < 3 || isempty(surr)
    sumstr = [];
    if iscellstr(s)
        for i = 1:length(s), sumstr = [sumstr, s{i}, d]; end
    elseif isnumeric(s)
        for i = 1:length(s), sumstr = [sumstr, num2str(s(i)), d]; end
    end
else
    sumstr = [];
    if iscellstr(s)
        for i = 1:length(s), sumstr = [sumstr, surr, s{i}, surr, d]; end
    elseif isnumeric(s)
        for i = 1:length(s), sumstr = [sumstr, surr, num2str(s(i)), surr, d]; end
    end    
end

out = sumstr(1:end-length(d));

end