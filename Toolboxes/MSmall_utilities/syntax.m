function syntax(filename)
%function syntax(filename)
%
% return the syntax of the file filename
% (i.e. either the first line, or the lines preceeding the first blank
% line)
%
%Michael Small
%3/3/2005
%ensmall@polyu.edu.hk

file=which([filename,'.m']);
if isempty(file),
    warning(['File ',filename,' not found.']);
    return;
end;
ret=textread(file,'%s','delimiter','\n','whitespace','');
n=1;i=1;lret=length(ret);
while i<lret,
    if strcmp(deblank(ret{i}),'%'),
        n=i;
        i=lret;
    end;
    i=i+1;
end;

if n>1,
    for i=2:(n-1),
        disp(ret{i}(2:end));
    end;
else
    disp(ret{1});
end;
