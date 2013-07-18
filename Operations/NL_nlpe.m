function out = NL_nlpe(y,de,tau)
% Wrapped up version of Michael Small's code, 'nlpe':
% http://small.eie.polyu.edu.hk/matlab/
% Ben Fulcher 19/2/2010


% Do my own inputs
% Embedding dimension
if nargin < 2 || isempty(de)
    de = 3;
end

% Time-delay, tau
if nargin < 3 || isempty(tau)
    tau = 1;
end
if strcmp(tau,'ac')
    tau = CO_fzcac(y);
end
if strcmp(tau,'mi')
    tau = CO_firstmin(y,'mi');
end

% Do false nearest neighbours if needed
if strcmp(de,'fnn')
    de = NL_MS_fnn(y,1:10,tau,5,1,1,0.05);
end

% normalized
% y=y-mean(y(:));
% y=y/std(y(:));

% Michael Small's nonlinear prediction error code:
res = MS_nlpe(y,de,tau); % residuals

%% Get outputs
out.msqerr = mean(res.^2);
% out.rmserr = sqrt(mean(e.^2));
% out.mabserr = mean(abs(e));
% out.meanres = mean(e);

% Use MF_residanal on the residuals
% 1) Get statistics on residuals
residstats = MF_residanal(res);

% convert these to local outputs in quick loop
fields = fieldnames(residstats);
for k = 1:length(fields);
    eval(sprintf('out.%s = residstats.%s;',fields{k},fields{k}));
end


end