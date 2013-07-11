function out = CO_fmmi(y)
% Finds the first minimum of the automutual information function
% Uses Rudy Moddemeijer's information.m code (RM_information.m here) that may or may not be great...
% Clumsily coded by Ben Fulcher, 2008

N = length(y);

a = zeros(N-1,1);
for i = 0:N-1
    a(i+1) = RM_information(y(1:end-i),y(1+i:end));
    if i > 1 && a(i-1)-a(i) > 0 && a(i)-a(i+1) < 0; % minimum
        out = i; % I found the first minimum!
        return
    end
end

% If there is no minimum in the automutual information
out = N;

end