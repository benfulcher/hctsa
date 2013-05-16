function out = SC_bendfa(x,k,q,taustep)
%% THIS IS REDUNDANT!!!!
%% NOW SUBSUMED INTO SC_benfa -- use the 'dfa' option....
% fits kth order polynomials to y in a piecewise fashion
% (to reduce nonstationarity/low-frequency trends)
% q is the parameter in the fluctuation function (usually 2)
%  -- 2 gives RMS fluctuations
% ; then looks at fluctuations about these
% taustep is the step in the scale size, which is in the range
% suggeted by Peng (1995)
% Ben Fulcher September 2009

% Preliminaries
N = length(x);

% perform scaling over a range of tau, up to a fifth the length
% of the time series
taur = 5:taustep:floor(N/4); % Peng (1995) suggests 5:N/4
ntau = length(taur);

% 1) Compute integrated sequence

y = cumsum(x);

% 2) Compute the fluctuation function as follows
F = zeros(ntau,1);
% F is the fluctuation function
% each entry correponds to a given scale tau, and contains
% the fluctuation function at that scale

for i = 1:ntau
    % buffer the time series at the scale tau
    tau = taur(i); % the scale on which to compute fluctuations
    
    y_buff = buffer(y,tau);
    if size(y_buff,2)>floor(N/tau) % zero-padded, remove trailing set of points...
        y_buff = y_buff(:,1:end-1);
    end
    
    % analyzed length of time series (with trailing end-points removed)
    nn = size(y_buff,2)*tau;
    
    % Detrend in each (nonoverlapping) subsegment
    tt = (1:tau)'; % faux time range
    for j=1:size(y_buff,2);
        % fit a polynomial of order k in each subsegment
        p = polyfit(tt,y_buff(:,j),k);
        
        % remove the trend, store back in y_buff
        y_buff(:,j) = y_buff(:,j) - polyval(p,tt);
    end
    
    % reshape to a column vector, y_dt (detrended)
    y_dt = reshape(y_buff,nn,1);
    
%     plot(y_dt); input('heyheyhey');
    
    F(i) = (mean(y_dt.^q)).^(1/q);
    
end

[linfit stats] = robustfit(log(taur),log(F));

out.linfitint = linfit(1);
out.alpha = linfit(2);
out.stats_coeffcorr = abs(stats.coeffcorr(1,2)); % correlation of coefficient estimates
stats.se1 = stats.se(1); % standard error in intercept
stats.se2 = stats.se(2); % standard error in mean
stats.ssr = mean(stats.resid.^2); % mean squares residual


plot(log(taur),log(F));
% keyboard

end