function out = MF_arma_orders(y,pr,qr)
% Given ranges of AR order p and MA order q fits ARMA models and looks at
% the goodness of fit for all combinations. Returns statistics on
% appropriateness of the combinations of models
% Uses functions iddata, armax, and aic from Matlab's System Identification toolbox
% Ben Fulcher 17/2/2010

% ** Future improvements **
% (1) May want to quantify where a particular order starts to stand out, i.e.,
% may be quite sensitive to wanting p > 2, but may be quite indiscriminate
% when it comes to the MA order.
% (2) May want to do some prediction and get more statistics on quality of
% model rather than just in-sample FPE/AIC...

%% Check Inputs

% Convert y to time series object
y = iddata(y,[],1);

% ARMA(p,q): p range, pr
if nargin < 2 || isempty(pr)
   pr = 1:10;
end

% ARMA(p,q): q range, qr
if nargin < 3 || isempty(qr)
    qr = (1:5);
end


%% Preliminaries
fpes = zeros(length(pr),length(qr));
aics = zeros(length(pr),length(qr));

%% Fit the models
for i = 1:length(pr)
    p = pr(i);
    for j = 1:length(qr)
        q = qr(j);
        
        % Fit the ARMA(p,q) model
        m = armax(y,[p,q]);
        
        % Get statistics on it
        fpes(i,j) = m.EstimationInfo.FPE;
        aics(i,j) = aic(m);
    end 
end

% global minimum, aic
out.aic_min = min(aics(:));
[pi_opt qi_opt] = find(aics == min(aics(:)));
out.p_aic_opt = pr(pi_opt);
out.q_aic_opt = qr(qi_opt);

out.std_all_aics = std(aics(:)); % no idea why.
out.mean_all_aics = mean(aics(:)); % no idea why.

out.meanstd_aicsp = mean(std(aics));
out.meanstd_aicsq = mean(std(aics'));


end