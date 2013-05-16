function out = DN_skewy(y,wh)

if strcmp(wh,'pearson')
    out = (3*mean(y)-median(y))./std(y);
elseif strcmp(wh,'bowley')
    qs = quantile(y,[0.25 0.5 0.75]);
    out = (qs(3)+qs(1) - 2*qs(2))./(qs(3)-qs(1)); % quartile skewness coefficient
else
    out = 0;
end

end