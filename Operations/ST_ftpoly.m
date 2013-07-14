function out = ST_ftpoly(y,k)
% Fits a polynomial of order k to the time series and returns the mean RMS error of
% the fit.
% Almost always kind of a stupid thing to do, but maybe it's somehow informative
% Ben Fulcher, 2009

N = length(y); % the length of the time series (number of samples)
t = (1:N)'; % Get a range for the time axis for time series y

% Supress the (valid!) warning from stupidly fitting a polynomial to a time series...
warning('off','MATLAB:polyfit:RepeatedPointsOrRescale');
cf = polyfit(t,y,k);
warning('on','MATLAB:polyfit:RepeatedPointsOrRescale');

f = polyval(cf,t);
out = sum((y-f).^2)/N; % mean RMS ERROR OF FIT

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