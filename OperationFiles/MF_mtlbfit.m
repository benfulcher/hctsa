function out=MF_mtlbfit(x,whoa,fox,nbins)
% x: the (zscore) time series
% whoa: the model name
% fox: the output
% nbins: for distribution fits: number of bins in the histogram
%                                       (or, if nbins=0, uses ksdensity)
%% Fit the model
% Two cases: distribution fits and fits on the data
Distmods={'gauss1','gauss2','exp1','power1'}; % valid distribution models
TSmods={'sin1','sin2','sin3','fourier1','fourier2','fourier3'}; % valid time series models

if any(strcmp(Distmods,whoa)); % valid DISTRIBUTION model name
    if nargin<4; % haven't specified nbin
        disp('!! Bin count not specified, using histogram with 10 bins')
        nbins=10;
    end
    if nbins==0; % use ksdensity instead of a histogram
        [dny dnx]=ksdensity(x);
    else
        [dny dnx]=hist(x,nbins);
    end
    if size(dnx,2)>size(dnx,1);dnx=dnx';dny=dny';end % must be column vectors
    [cfun gof output]=fit(dnx,dny,whoa); % fit the model
elseif any(strcmp(TSmods,whoa)); % valid TIME SERIES model name
    if size(x,2)>size(x,1);x=x';end % x must be a column vector
    t=(1:length(x))'; % Assume equal sampling of the univariate time series
    [cfun gof output]=fit(t,x,whoa); % fit the model
else
    disp('Invalid model name'); return
end
%% Prepare the Output
switch fox
    case 'r2' % rsqured
        out=gof.rsquare;
    case 'adjr2' % degrees of freedom-adjusted rsqured
        out=gof.adjrsquare;
    case 'rmse' % root mean square error
        out=gof.rmse;
    case 'resAC1' % autocorrelation of residuals at lag 1
        out=CO_autocorr(output.residuals,1);
    case 'resAC2' % autocorrelation of residuals at lag 2
        out=CO_autocorr(output.residuals,2);
    case 'resruns' % runs test on residuals -- outputs p-value
        out=HT_hyptests(output.residuals,'runstest');
    case 'reslbq' % lbq test of randomness on residuals -- p-value
        % NEVER HAVE REQUIRED ECONOMETRICS TOOLBOX -- REMOVE THIS OPTION
        %         out=HT_hyptests(output.residuals,'lbq');
    otherwise
        disp('Invalid output setting');return
end
end