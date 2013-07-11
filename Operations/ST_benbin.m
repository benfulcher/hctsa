function out = ST_benbin(x,stretchwhat)
% Converts the time series, x, to a binary series
% FYI: THIS CODE DOESN'T DO WHAT IT THINKS IT DOES!
% Ben Fulcher, 2009

N = length(x); % length of the time series
x(x > 0) = 1;
x(x <= 0) = 0;


switch stretchwhat
    case 'lseq1'
        % longest stretch of 1s (this doesn't actually measure this!)
        out = max(diff(sgnchange(diff(find(x == 1))-1.5)))/N;
    case 'lseq0'
        % longest stretch of 0s (this doesn't actualy measure this!)
        out = max(diff(sgnchange(diff(find(x == 0))-1.5)))/N;
    otherwise
        error('Unknown input %s',stretchwhat)
end

if isempty(out)
    out = 0;
end

end