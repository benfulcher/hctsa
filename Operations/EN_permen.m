function permen = EN_permen(y,ord)
% Permutation entropy algorithm
% Adapted from logisticPE, downloaded from an online complexity article:
% http://people.ece.cornell.edu/land/PROJECTS/Complexity/
% Code adapted by Ben Fulcher, but adapted September 2009

% Ensure y is a column vector
if size(y,1) > size(y,2);
    y = y';
end
% length of the time series, N
N = length(y);

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