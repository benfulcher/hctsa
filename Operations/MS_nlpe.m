function out = MS_nlpe(y,de,tau)
% Wrapped up version of Michael Small's code:
% http://small.eie.polyu.edu.hk/matlab/
% Ben Fulcher 19/2/2010


% Do my own inputs
% Embedding dimension
if nargin < 2 || isempty(de)
    de=3;
end


% Time-delay, tau
if nargin < 3 || isempty(tau)
    tau=1;
end
if strcmp(tau,'ac')
    tau = CO_fzcac(y);
end
if strcmp(tau,'mi')
    tau = CO_fmmi(y);
end

% Do false nearest neighbours if needed
if strcmp(de,'fnn')
    de = MS_fnn(y,1:10,tau,5,1,1,0.05);
end

% normalized
% y=y-mean(y(:));
% y=y/std(y(:));


% _____________________________________________________________
%function e=nlpe(y,de,tau);
%function e=nlpe(y,v);
%
%compute the normalised "drop-one-out" constant interpolation nonlinear
%prediction error for embedding dimension de and lag tau or for embedding
%strategy v (v>0)
%
%Michael Small
%3/3/2005
%ensmall@polyu.edu.hk
%


if min(size(y))>1,
    x=y;
    y=y(1,2:end);
    x=x(:,1:(end-1));
elseif max(size(de))>1,
    v=de(de>0);
    [x,y]=MS_embed(y,v-1); % changed from embed to MS_embed ++BF
else
    [x,y]=MS_embed(y,[-1 0:tau:((de-1)*tau)]); % changed from embed to MS_embed ++BF
end;

if isempty(x)
    disp('forget about this -- can''t embed...')
    out = NaN;
    return
end

[de,n]=size(x);


% ++BF: this bit can cause memory pains for large time series
% Let's do this dirty cheat
maxn = 3000;
if n>3000
    x = x(:,1:3000);
    y = y(1:3000);
    n = maxn;
end
if n<20
    disp('time series too short')
    out = NaN;
    return
end

% dd=sparse(n,n);
dd = zeros(n,n);
for i=1:de, %loop on de and compute the distance.^2
    dd=dd+(ones(n,1)*x(i,:)-x(i,:)'*ones(1,n)).^2;
end;
% dd is the distance .^2 with inf on the diag
% speye = sparse(1:n,1:n,1);
warning off MATLAB:divideByZero
dd=dd+1./(1-eye(n,n));
% dd=dd+1./(1-speye);
warning on MATLAB:divideByZero
%near is the index of the nearest neighbour of each point
[dist,near] = min(dd);
%the prediction error is
e = y(near)-y;
% e = mean(e.^2);
% _________________________________________________________ %

% Thanks Michael, and now back to me

%% Get outputs
out.msqerr = mean(e.^2);
% out.rmserr = sqrt(mean(e.^2));
% out.mabserr = mean(abs(e));
% out.meanres = mean(e);

% Use MF_residanal on the residuals
% 1) Get statistics on residuals
residout = MF_residanal(e);

% convert these to local outputs in quick loop
fields = fieldnames(residout);
for k=1:length(fields);
    eval(['out.' fields{k} ' = residout.' fields{k} ';']);
end


end