% LA_permen
% 
% Originally logisticPE.m
% http://people.ece.cornell.edu/land/PROJECTS/Complexity/
% http://people.ece.cornell.edu/land/PROJECTS/Complexity/logisticPE.m
% Original code by Bruce Land and Damian Elias
%
% Modified slightly by Ben Fulcher, 2009
%

function permen = LA_permen(y,ord)

N = length(y); % time-series length, N

permlist = perms(1:ord);
% c(1:length(permlist))=0;
c = zeros(size(permlist,1),1);

for j = 1:N-ord
    [a, iv] = sort(y(j:j+ord-1));
    for jj = 1:length(permlist)
        if (abs(permlist(jj,:)-iv)) == 0
            c(jj) = c(jj) + 1 ;
        end
    end
end

p = max(1/N,c/(N-ord));
permen = -sum(p .* log(p))/(ord-1);

end