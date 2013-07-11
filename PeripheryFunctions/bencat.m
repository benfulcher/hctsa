function out = bencat(s,x,surr)
% x is the delimiter
% s is the list of strings to be catinated
% surrounds each string in optional string surr

if nargin < 2, x= ', '; % default delimiter is a comma
else x=[x ' '];
end

if nargin < 3 || isempty(surr)
    sumstr=[];
    if iscellstr(s)
        for i=1:length(s), sumstr=[sumstr s{i} x]; end
    elseif isnumeric(s)
        for i=1:length(s), sumstr=[sumstr num2str(s(i)) x]; end
    end
else
    sumstr=[];
    if iscellstr(s)
        for i=1:length(s), sumstr=[sumstr surr s{i} surr x]; end
    elseif isnumeric(s)
        for i=1:length(s), sumstr=[sumstr surr num2str(s(i)) surr x]; end
    end    
end

out=sumstr(1:end-length(x));

end