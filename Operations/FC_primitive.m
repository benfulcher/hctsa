function out = FC_primitive(y,fmeth,fparam)
% Primitive forecasting
% **fmeth**
% 'mean': uses the mean of previous values to predict the next
%   <fparam> specifies the length of previous window with which to make
%   predictions; if 'ac' then uses tau
% 'median': uses median of previous values to predict the next
%   <fparam> specifies the length of the previous window with which to
%   make predictions
% 'acf': uses the autocorrelation coefficients at order given by <fparam>
% 'lfit': fit linear to a certain window of the time series; do linear
%           extrapolation
% Ben Fulcher, Nov 2009

switch fmeth
    case 'mean'
        if strcmp(fparam,'ac')
            lp = CO_fzcac(y); % make it tau
        else
            lp = fparam; % the length of the subsegment preceeding to use to predict the subsequent value
        end
        evalr = lp+1:length(y); % range over which to evaluate the forecast
        res = zeros(length(evalr),1); % residuals
        for i = 1:length(evalr)
            res(i) = mean(y(evalr(i)-lp:evalr(i)-1)) - y(evalr(i)); % prediction-value
        end
        
    case 'median'
        if strcmp(fparam,'ac')
            lp = CO_fzcac(y); % make it tau
        else
            lp = fparam; % the length of the subsegment preceeding to use to predict the subsequent value
        end
        evalr = lp+1:length(y); % range over which to evaluate the forecast
        res = zeros(length(evalr),1); % residuals
        for i = 1:length(evalr)
            res(i) = median(y(evalr(i)-lp:evalr(i)-1)) - y(evalr(i)); % prediction-value
        end

%     case 'acf' % autocorrelation function
%         acl=fparam; % autocorrelation 'length'
%         acc=zeros(acl,1); % autocorrelation coefficients
%         for i=1:acl, acc(i)=CO_autocorr(y,i); end
%         % normalize to a sum of 1 (so that operating on three mean values
%         % of the time series, also returns the mean value as output)
%         acc=acc/sum(acc);
%         
%         evalr=acl+1:length(y); % range over which to evaluate the forecast
%         res=zeros(length(evalr),1); % residuals
%         for i=1:length(evalr)
%             res(i)=sum(acc.*(y(evalr(i)-acl:evalr(i)-1))) - y(evalr(i)); % prediction-value
%         end

    case 'lfit'
        if strcmp(fparam,'ac')
            lp = CO_fzcac(y); % make it tau
        else
            lp = fparam; % the length of the subsegment preceeding to use to predict the subsequent value
        end
        evalr = lp+1:length(y); % range over which to evaluate the forecast
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
tmp = MF_M_mtlbfit(y,'gauss1',0);
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