% FC_LocalSimple
% 
% Does local forecasting using very simple predictors using the past $l$ values
% of the time series to predict its next value.
% 
% Three prediction methods are implemented:
% (i) 'mean': local mean prediction using the past ltrain time-series values,
% (ii) 'median': local median prediction using the past ltrain time-series values, and
% (iii) 'lfit': local linear prediction using the past ltrain time-series values.
% 
% INPUTS:
% y, the input time series
% fmeth, the forecasting method
% ltrain, the number of time-series values to use to forecast the next value
% 
% Outputs are the mean error, stationarity of residuals, Gaussianity of
% residuals, and their autocorrelation structure.

function out = FC_LocalSimple(y,fmeth,ltrain)
% Ben Fulcher, Nov 2009

% Check inputs
if nargin < 2 || isempty(fmeth)
    fmeth = 'mean';
end
if nargin < 2 || isempty(ltrain)
    ltrain = 3;
end

N = length(y); % time-series length

switch fmeth
    case 'mean'
        if strcmp(ltrain,'ac')
            lp = CO_fzcac(y); % make it tau
        else
            lp = ltrain; % the length of the subsegment preceeding to use to predict the subsequent value
        end
        evalr = lp+1:N; % range over which to evaluate the forecast
        res = zeros(length(evalr),1); % residuals
        for i = 1:length(evalr)
            res(i) = mean(y(evalr(i)-lp:evalr(i)-1)) - y(evalr(i)); % prediction-value
        end
        
    case 'median'
        if strcmp(ltrain,'ac')
            lp = CO_fzcac(y); % make it tau
        else
            lp = ltrain; % the length of the subsegment preceeding to use to predict the subsequent value
        end
        evalr = lp+1:N; % range over which to evaluate the forecast
        res = zeros(length(evalr),1); % residuals
        for i = 1:length(evalr)
            res(i) = median(y(evalr(i)-lp:evalr(i)-1)) - y(evalr(i)); % prediction-value
        end

%     case 'acf' % autocorrelation function
%         acl=ltrain; % autocorrelation 'length'
%         acc=zeros(acl,1); % autocorrelation coefficients
%         for i=1:acl, acc(i)=CO_autocorr(y,i); end
%         % normalize to a sum of 1 (so that operating on three mean values
%         % of the time series, also returns the mean value as output)
%         acc=acc/sum(acc);
%         
%         evalr=acl+1:N; % range over which to evaluate the forecast
%         res=zeros(length(evalr),1); % residuals
%         for i=1:length(evalr)
%             res(i)=sum(acc.*(y(evalr(i)-acl:evalr(i)-1))) - y(evalr(i)); % prediction-value
%         end

    case 'lfit'
        if strcmp(ltrain,'ac')
            lp = CO_fzcac(y); % make it tau
        else
            lp = ltrain; % the length of the subsegment preceeding to use to predict the subsequent value
        end
        evalr = lp+1:N; % range over which to evaluate the forecast
        res = zeros(length(evalr),1); % residuals
        for i = 1:length(evalr)
            % fit linear
            p = polyfit((1:lp)',y(evalr(i)-lp:evalr(i)-1),1);
            res(i) = polyval(p,lp+1) - y(evalr(i)); % prediction - value
        end
        
    otherwise
        error('Unknown forecasting method ''%s''',fmeth);
end

% out=res;
% plot(res);
% measures on the errors time series:

% mean error:
out.meanerr = mean(res);
out.rmserr = sqrt(mean(res.^2));
out.meanabserr = mean(abs(res));

% standard deviation of errors:
out.stderr = std(res);

% stationarity:
out.sws = SY_slidwin(res,'std','std',5,1);
out.swm = SY_slidwin(res,'mean','std',5,1);

% normality:
% out.chi2n=HT_disttests(res,'chi2gof','norm',10); % chi2
% out.ksn=HT_disttests(res,'ks','norm'); % Kolmogorov-Smirnov
tmp = DN_simplefit(y,'gauss1',0);
if ~isstruct(tmp) && isnan(tmp) % fitting failed
    out.gofnadjr2 = NaN;
else
    out.gofnadjr2 = tmp.adjr2; % degrees of freedom-adjusted rsqured
end

% autocorrelation structure:
out.ac1 = CO_autocorr(res,1);
out.ac2 = CO_autocorr(res,2);
out.taures = CO_fzcac(res);
out.tauresrat = CO_fzcac(res)/CO_fzcac(y);

end