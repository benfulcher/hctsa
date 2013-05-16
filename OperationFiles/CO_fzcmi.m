function out=CO_fzcmi(y)
% Ben Fulcher 19/3/2010 -- corrected for error and cleaned up a few things
% from this old, messy code. In fact, I think it's used only as CO_fmmi?
% The title is wrong, in any case (!) -- it doesn't look for a 'fzc' =
% first zero crossing.

N = length(y);

a = zeros(N-3,1);
a(1) = CO_information(y,y);
for i=2:N-3
    try
        a(i) = CO_information(y(1:end-i),y(1+i:end));
    catch emsg
        % some crazy thing happened
        out = NaN;
        return
    end
    if a(i-1)*a(i)<0
        out = i;
        return
    end
end

% if no minimum in entire spectrum
out = N;

end