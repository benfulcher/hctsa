function out=MF_ARerr(y,p)

% fit AR model of order p

[a e]=arcov(y,p);
y_est=filter([0 -a(2:end)],1,y);
N=length(y);

switch n
	case 1 % rms errors
		err=sqrt(sum((y-y_est).^2)/length(y));
		out=err;
	case 2 % mean error
		err=y-y_est;
	    out=mean(err);
	case 3 % standard deviation of the residuals
		err=y-y_est;
	    out=std(err);
	case 4 % autocorrelation of residuals at lag 1
		err=y-y_est;
	    out=sum((err(1:N-1)-mean(err(1:N-1))).*(err(2:N)-mean(err(2:N))))/N/std(err(1:N-1))/std(err(2:N));
end