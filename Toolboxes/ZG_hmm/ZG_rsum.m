% row sum
% function Z = ZG_rsum(X)

function Z = ZG_rsum(X)

Z = zeros(size(X(:,1)));

for i = 1:length(X(1,:))
    Z = Z + X(:,i);
end

end