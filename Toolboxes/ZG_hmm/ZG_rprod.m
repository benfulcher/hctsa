% row product
% function Z = ZG_rprod(X,Y)

function Z = ZG_rprod(X,Y)

[n, m] = size(X);

if (length(Y(:,1)) ~=n) || (length(Y(1,:)) ~=1)
    error('Error in RPROD');
end

Z = X.*(Y*ones(1,m));

end