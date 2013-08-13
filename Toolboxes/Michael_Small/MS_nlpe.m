% function e = MS_nlpe(y,v);
%
% Compute the normalised "drop-one-out" constant interpolation nonlinear
% prediction error for embedding dimension de and lag tau or for embedding
% strategy v (v>0)
%
% Michael Small
% michael.small@uwa.edu.au, http://school.maths.uwa.edu.au/~small/
% 3/3/2005
% For further details, please see M. Small. Applied Nonlinear Time Series
% Analysis: Applications in Physics, Physiology and Finance. Nonlinear Science
% Series A, vol. 52. World Scientific, 2005. (ISBN 981-256-117-X) and the
% references therein.
% (Minor edits by Ben Fulcher, 2010)

function e = MS_nlpe(y,de,tau);

if min(size(y)) > 1
    x = y;
    y = y(1,2:end);
    x = x(:,1:(end-1));
elseif max(size(de))>1,
    v = de(de > 0);
    [x, y] = MS_embed(y,v-1);
else
    [x, y] = MS_embed(y,[-1, 0:tau:((de-1)*tau)]);
end

if isempty(x)
    error('Error embedding the time series.')
end

[de, n] = size(x);

% dd=sparse(n,n);
dd = zeros(n,n);
for i = 1:de %loop on de and compute the distance.^2
    dd = dd + (ones(n,1)*x(i,:)-x(i,:)'*ones(1,n)).^2;
end
% dd is the distance .^2 with inf on the diag
% speye = sparse(1:n,1:n,1);
warning('off','MATLAB:divideByZero')
dd = dd+1./(1-eye(n,n));
% dd=dd+1./(1-speye);
warning('on','MATLAB:divideByZero')
%near is the index of the nearest neighbour of each point
[dist, near] = min(dd);
%the prediction error is
e = y(near)-y;
% e = mean(e.^2);
% _________________________________________________________ %
end