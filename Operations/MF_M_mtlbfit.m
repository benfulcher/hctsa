function out = MF_M_mtlbfit(x,dmodel,nbins)
% Fits a model simple time-series model or distribution to the data
% Uses the function fit from Matlab's Curve Fitting Toolbox
% An output uses the function runstest from Matlab's Statistics Toolbox
% x: the (z-scored) time series as a column vector
% dmodel: the model name
% nbins: for distribution fits: number of bins in the histogram
%                                       (or, if nbins=0, uses ksdensity)
% Ben Fulcher, 2009


%% Fit the model
% Two cases: distribution fits and fits on the data
Distmods = {'gauss1','gauss2','exp1','power1'}; % valid distribution models
TSmods = {'sin1','sin2','sin3','fourier1','fourier2','fourier3'}; % valid time series models

if any(strcmp(Distmods,dmodel)); % valid DISTRIBUTION model name
    if nargin < 3; % haven't specified nbin
        error('!! You must specify a bin count !!')
    end
    if nbins == 0; % use ksdensity instead of a histogram
        [dny, dnx] = ksdensity(x);
    else
        [dny, dnx] = hist(x,nbins);
    end
    if size(dnx,2) > size(dnx,1); dnx = dnx'; dny = dny'; end % must be column vectors
    
    try
        [cfun, gof, output] = fit(dnx,dny,dmodel); % fit the model
	catch emsg % this model can't even be fitted OR license problem...
        if strcmp(emsg.message,'NaN computed by model function.') ...
                || strcmp(emsg.message,'Inf computed by model function.') ...
                || strcmp(emsg.message,'Power functions cannot be fit to non-positive xdata.') ...
                || strcmp(emsg.identifier,'curvefit:fit:nanComputed')
            fprintf(1,'The model, %s, failed for this data -- returning NaNs for all fitting outputs\n',dmodel);
            out = NaN; return
        else
            error('MF_M_mtlbfit(x,%s,%u): Unexpected error fitting %s to the data distribution',dmodel,nbins,dmodel)
        end
	end
elseif any(strcmp(TSmods,dmodel)); % valid TIME SERIES model name
    if size(x,2)>size(x,1); x = x'; end % x must be a column vector
    t = (1:length(x))'; % Assume equal sampling of the univariate time series
    try
        [cfun, gof, output] = fit(t,x,dmodel); % fit the model
	catch emsg % this model can't even be fitted OR license problem
        if strcmp(emsg.message,'NaN computed by model function.') || strcmp(emsg.message,'Inf computed by model function.')
            disp(['The model ' dmodel ' failed for this data -- returning NaNs for all fitting outputs']);
            out = NaN;
            return
        else
            error(['MF_M_mtlbfit: Unexpected error fitting ' dmodel ' to the time series'])
        end
	end
else
    error('Invalid distribution or time-series model specified');
end

%% Prepare the Output
out.r2 = gof.rsquare; % rsquared
out.adjr2 = gof.adjrsquare; % degrees of freedom-adjusted rsqured
out.rmse = gof.rmse;  % root mean square error
out.resAC1 = CO_autocorr(output.residuals,1); % autocorrelation of residuals at lag 1
out.resAC2 = CO_autocorr(output.residuals,2); % autocorrelation of residuals at lag 2
out.resruns = HT_hyptests(output.residuals,'runstest'); % runs test on residuals -- outputs p-value

end