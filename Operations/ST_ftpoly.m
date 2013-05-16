function out = ST_ftpoly(y,n)
% Fits a polynomial to the time series and returns the mean RMS error of
% the fit.
% Ben Fulcher

x = 1:length(y);
if size(y,2)<size(y,1);
    y = y';
end
cf = polyfit(x,y,n);
f = polyval(cf,x);
out = sum((y-f).^2)/length(y); % mean RMS ERROR OF FIT

% % n=10;
% errs=zeros(n,1);
% x=1:length(y);
% for i=1:n
% cf=polyfit(x,y',i);
% f=polyval(cf,x);
% errs(i)=sum((y'-f).^2);
% % end
% % hold on;plot(f,'k')
% % plot(errs);

end