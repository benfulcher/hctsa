function out = ST_benbin(x,zey)

N = length(x);
%convert to a binary series
x(x>0) = 1;
x(x<=0) = 0;


% THIS CODE DOESN'T DO WHAT IT SHOULD DO!

switch zey
    case 'lseq1'
        % longest stretch of 1s (this doesn't actually measure this!)
        out = max(diff(sgnchange(diff(find(x==1))-1.5)))/N;
    case 'lseq0'
        % longest stretch of 0s (this doesn't actualy measure this!)
        out = max(diff(sgnchange(diff(find(x==0))-1.5)))/N;
end

if isempty(out), out=0; end

end