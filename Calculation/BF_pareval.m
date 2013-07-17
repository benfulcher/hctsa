function [out, T] = BF_pareval(x,y,s,bevocal)
% x and y are possible elements of the string s to be evaluated
% Stores any text output as T

if bevocal
    % Any text output from operations is printed to screen
    out = eval(s); 
else
    % Text output from operations is suppressed and stored in T
    [T, out] = evalc(s);
end

end