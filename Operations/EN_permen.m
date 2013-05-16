function e = EN_permen(y,ord)
% Adapted from logisticPE, downloaded from an online complexity article
% Ben Fulcher approx September 2009


if size(y,1)>size(y,2); y=y'; end
N = length(y);

permlist = perms(1:ord);
% c(1:length(permlist))=0;
c = zeros(size(permlist,1),1);

for j = 1:N-ord
    [a,iv] = sort(y(j:j+ord-1));
    for jj = 1:length(permlist)
        if (abs(permlist(jj,:)-iv)) == 0
            c(jj) = c(jj) + 1 ;
        end
    end
end

p = max(1/N,c/(N-ord));
e = -sum(p .* log(p))/(ord-1);

end