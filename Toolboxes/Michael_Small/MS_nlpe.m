% _____________________________________________________________
function e = MS_nlpe(y,de,tau);
% function e = MS_nlpe(y,v);
%
%compute the normalised "drop-one-out" constant interpolation nonlinear
%prediction error for embedding dimension de and lag tau or for embedding
%strategy v (v>0)
%
%Michael Small
%3/3/2005
%ensmall@polyu.edu.hk
%
% Edited to work with HCTS package, Ben Fulcher, 2010

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

% ++BF: this bit can cause memory pains for large time series
% Let's do this dirty cheat
maxn = 3000;
if n > maxn
    x = x(:,1:maxn);
    y = y(1:maxn);
    n = maxn;
    fprintf(1,'Michael Small''s ''nlpe'' is only being evaluated on the first %u samples...\n',maxn);
end
if n < 20
    fprintf(1,'Time series (N = %u) is too short\n',length(y))
    out = NaN;
    return
end

% dd=sparse(n,n);
dd = zeros(n,n);
for i = 1:de %loop on de and compute the distance.^2
    dd = dd + (ones(n,1)*x(i,:)-x(i,:)'*ones(1,n)).^2;
end;
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