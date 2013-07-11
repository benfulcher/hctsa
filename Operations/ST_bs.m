function out = ST_bs(x,whatstat)
% Ben's Simple Stats (ST_bs)
N = length(x);

switch whatstat
    case 'zcross'
        % proportion of zero-crossings of the time series
        % (i.e., crosses its mean) should be higher for noise
        xch = x(1:end-1).*x(2:end);
        out = sum(xch<0)/N;
    case 'max'
        % proportion of local maxima in the time series
        % Changed by Ben Fulcher November 2009
%         dx=diff(x);
%         ddx=diff(dx);
%         sdsp=ddx(sgnchange(dx));
%         out=length(find(sdsp<0))/N;
        dx = diff(x);
        out = length(find(dx(1:end-1)>0 & dx(2:end)<0))/(N-1);
    case 'min'
        % proportion of local minima in the time series
        % Changed by Ben Fulcher November 2009
%         dx=diff(x);
%         ddx=diff(dx);
%         sdsp=ddx(sgnchange(dx));
%         out=length(find(sdsp>0))/length(x);
        dx = diff(x);
        out = length(find(dx(1:end-1)<0 & dx(2:end)>0))/(N-1);
    case 'pmcross'
        % ratio of times cross 1 to -1
        c1sig=length(sgnchange(x-1)); % num times cross 1
        c2sig=length(sgnchange(x+1)); % num times cross -1
        if c2sig == 0
            out=NaN;
        else
            out=c1sig/c2sig;
        end 
    case 'zsczcross'
        % ratio of zero crossings of raw to detrended time series
        % where the raw has zero mean
        x=zscore(x);
        xch=x(1:end-1).*x(2:end);
        h1=length(find(xch<0)); % num of zscross of raw series
        y=detrend(x);
        ych=y(1:end-1).*y(2:end);
        h2=length(find(ych<0)); % of detrended series
        if h1 == 0;
            out=NaN;
        else
            out=h2/h1;
        end
end

end 