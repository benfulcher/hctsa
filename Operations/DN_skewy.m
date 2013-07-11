function out = DN_skewy(y,wh)
% Estimates custom skewness measures: either 'pearson' or 'bowley' for the input time series, y
% The 'bowley' method uses the quantile function from Matlab's Statistics Toolbox
% Ben Fulcher, 2009

switch wh
    case 'pearson'
        out = (3*mean(y)-median(y))./std(y);
    case 'bowley'
        qs = quantile(y,[0.25, 0.5, 0.75]);
        out = (qs(3)+qs(1) - 2*qs(2))./(qs(3)-qs(1));
        % Quartile skewness coefficient
    otherwise
        error('DN_skewy: unknown method %s',wh)
end

end