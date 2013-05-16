function out = ST_cumul(y,wing,lorp,n)
% Compares the statistics from a given subset of the time series to that
% obtained from the full time series
% Ben Fulcher, September 2009

N = length(y);

% determine subset range to use: r
switch lorp
    case 'l'
        r=1:min(n,N); % takes first n points of time series
    case 'p'
    	r=1:round(N*n); % takes initial proportion n of time series
    case 'unicg'
        r=round(linspace(1,N,n)); % takes n uniformly distributed points in time series
    case 'randcg'
        r=randi(N,n,1); % takes n random points in time series; there could be repeats
        % This is quite unrobust, as it's taking just a single sample from
        % a test with a (possibly) large variance
end

% compare this subset to the full value
switch wing
    case 'mean'
        out = abs(mean(y(r))); %/mean(y); % Y SHOULD BE Z-SCORED;;
    case 'std'
        out = std(y(r)); %/std(y); % Y SHOULD BE Z-SCORED;;
    case 'median'
        out = median(y(r)); %/median(y); % if median is very small;; could be very noisy
    case 'iqr'
        out = abs(1-iqr(y(r))/iqr(y));
    case 'skewness'
        out = abs(1-skewness(y(r))/skewness(y)); % how far from true
    case 'kurtosis'
        out = abs(1-kurtosis(y(r))/kurtosis(y)); % how far from true
    case 'AC1'
        out = abs(1-CO_autocorr(y(r),1)/CO_autocorr(y,1)); % how far from true
    case 'SampEn1_01' % computationally expensive to calculate this full one each time...
        out = EN_sampenc(y(r),1,0.1)/EN_sampenc(y,1,0.1);
end


% N = length(y);
% cockburns = 10; % how far into the time series to look
% r=round(linspace(1,N,cockburns)); % the range of points in the time series
% cums = zeros(cockburns,1);
% switch wing
%     case 'mean'
%         for i=1:cockburns
%             cums(i)=mean(y(1:r(i)));
%         end
%         cums=cums/mean(y);
%     case 'std'
%         for i=1:cockburns
%             cums(i)=std(y(1:r(i)));
%         end
%         cums=cums/std(y);
%     case 'median'
%         for i=1:cockburns
%             cums(i)=median(y(1:r(i)));
%         end
%         cums=cums/median(y);
%     case 'iqr'
%         for i=1:cockburns
%             cums(i)=iqr(y(1:r(i)));
%         end
%         cums=cums/iqr(y);
%     case 'SampEn1_01'
%         for i=1:cockburns
%             cums(i)= EN_sampenc(y(1:r(i)),1,0.1); 
%         end
%         cums=cums/EN_sampenc(y,1,0.1);
% end
% 
% 
% plot(cums);



end