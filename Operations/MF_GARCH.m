function out = MF_GARCH(y,preproc,params)
% Tries to simulate a GARCH model fitting procedure:
% (1) Preprocessing: preprocesses the time series in an appropriate way to
%       remove strong trends (optional)
% (2) Pre-estimation: looks at correlation properties in the time series
%       (to motivate an appropriate GARCH model)
% (3) Fitting: fits a GARCH model, returning goodness of fit statistics and
%               parameters from the fitted model.
% (4) Post-estimation: returns statistics on residuals and standardized
%                       residuals

% Idea is that all of these stages can be pre-specified or skipped using
% arguments to the function.
% Also, I've minimal experience in actually fitting GARCH models...!
% This function should ideally be reformed by an expert.
% Uses functions from MATLAB's Econometrics Toolbox: archtest, lbqtest, autocorr, parcorr, garchset, garchfit, garchcount, aicbic
% Ben Fulcher 25/2/2010

%% Inputs

% Preprocessing settings: preproc:
%       (i) 'nothing' does nothing -- just takes in the raw time series
%       (ii) 'auto' [default] (returns the best)
%       (iii) 'diff' (differencing)
%       (iv) 'logreturns' (only for positive-only data)
%       (v) 'pwpoly' (fits piecewise polynomials to remove low-f trends)
if nargin < 2 || isempty(preproc)
    preproc = 'ar'; % do the preprocessing that maximizes stationarity/whitening
end

% % Preestimation settings: preest
% if nargin < 3 || isempty(preest)
%     preest = 'yep'; % do pre-estimation
% end

% Model fitting
if nargin < 3 || isempty(params)
    params = 'default'; % fit the default GARCH model
end


%% (1) Data preprocessing
y0 = y; % the original, unprocessed time series

switch preproc
    case 'nothing'
        % do nothing.
    case 'ar'
        % apply a number of standard preprocessings and return them in the
        % structure ypp. Also chooses the best preprocessing based on the worst fit
        % of an AR2 model to the processed time series.
        % has to beat doing nothing by 5% (error)
        % No spectral methods allowed...
        [ypp best] = benpp(y,'ar',2,0.05,0);
        eval(['y = ypp.' best ';']);
        disp(['Proprocessed according to AR(2) criterion using ' best]);
end

y = benzscore(y); % z-score the time series
N = length(y); % could be different to original (if choose a differencing, e.g.)

% Now have the preprocessed time series saved over y.
% The original, unprocessed time series is retained in y0.
% (Note that y=y0 is possible; when all preprocessings are found to be
%   worse at the given criterion).

%% (2) Data pre-estimation
% Aim is to return some statistics indicating the suitability of this class
% of modeling.
% Will use the statistics to motivate a GARCH model in the next
% section.
% Will use the statistics to compare to features of the residuals after
% modeling.

% (i) Engle's ARCH test
%       look at autoregressive lags 1:20
%       use the 10% significance level
[Engle_h_y, Engle_pValue_y, Engle_stat_y, Engle_cValue_y] = archtest(y,1:20,0.1);

% (ii) Ljung-Box Q-test
%       look at autocorrelation at lags 1:20
%       use the 10% significance level
%       departure from randomness hypothesis test
[lbq_h_y2, lbq_pValue_y2, lbq_stat_y2, lbq_cValue_y2] = lbqtest(y.^2,1:20,0.1);
% [lbq_h_y2, lbq_pValue_y2, lbq_stat_y2, lbq_cValue_y2] = lbqtest(y.^2,1:20,0.1,[]);


% (iii) Correlation in time series: autocorrelation
% autocorrs_y = CO_autocorr(y,1:20);
% autocorrs_var = CO_autocorr(y.^2,1:20);
[ACF_y,Lags_acf_y,bounds_acf_y] = autocorr(y,20,[],[]);
[ACF_var_y,Lags_acf_var_y,bounds_acf_var_y] = autocorr(y.^2,20,[],[]);

% (iv) Partial autocorrelation function: PACF
[PACF_y,Lags_pacf_y,bounds_pacf_y] = parcorr(y,20,[],[]);


%% (3) Create an appropriate GARCH model
switch params
    % Names of basic models
    case 'default'
        % (i) The default model:
        % a constant, C mean process
        % GARCH(1,1) conditionally Gaussian innovations
        
        % C: constant mean of mean-process
        % K: constant term in variance-process
        % GARCH(1): autoregressive parameter at lag 1 of variance
        % ARCH(1): regressive parameter at lag 1 of Gaussian noise process
        %           squared onto the variance process
        
        spec = garchset('P', 1, 'Q', 1); % equivalent to setting nothing in garchfit

    case 'auto'
        % automatically determines model from statistics above
        % This is not a great 'automatic' method, I think! But it's a
        % method...
        
        % AR component of model from ACF
        R = Lags_acf_y(find(abs(ACF_y) < bounds_acf_y(1),1,'first'))-1;
        % first time autocorrelation drops below significance level.
        if isempty(R), R = length(ACF_y)+1; end
        if R>4, R = 4; end
        % (**) R=0 implies that no AR component is required.
        
        
        % MA component of model from PACF
        M = Lags_pacf_y(find(abs(PACF_y) < bounds_pacf_y(1),1,'first'))-1;
        % first time partial autocorrelation drops below signifance level
        if isempty(M), M = Lags_pacf_y(end)+1; end
        if M>3, M=3; end % don't want to calculate massive mean models
        % (**) M=0 implies that no MA component is required.
        
        
        % I have no intutition as to how to add the GARCH component. I
        % guess if there's no evidence of heteroskedacity you wouldn't even
        % attempt a GARCH component, but for now let's fit one anyway. This
        % is a GARCH routine, after all...
        P = 1; % autoregression onto lagged conditional variance
        Q = 1; % regression of Gaussian noise process onto conditional variance
        
        
        % make an appropriate string
        garchpval = [R M P Q]; % 4-component vector of model orders
        garchpnam = {'R','M','P','Q'};
        s = '';
        for i=1:length(garchpval)
            if ~garchpval(i) == 0 % this should be a component to specify in 
                s = [s '''' garchpnam{i} ''', ' num2str(garchpval(i)) ', '];
            end
        end
        s = s(1:end-2); % remove the Oxford comma.
        % This string, s, should now specify an argument to garchset
        
        eval(['spec = garchset(' s ');']);
        
        
    otherwise
        % Specify the GARCH model as a string in the input
        % e.g., 'R, 2, M, 1, P, 1, Q, 1' will fit an ARMA(2,1) to mean
        % process and a GARCH(1,1) to the variance process
        try
            eval(['spec = garchset(' params ');']);
        catch emsg
           disp('Error in formatting input parameters specifying GARCH model. Exiting.')
           return
        end
    
end

spec = garchset(spec,'C', NaN); % fix C=0 -- zero-mean process
% In fact this gives quite different results, even when you C ends up being
% very close to zero...? Strangely not for the ARMA, but for the GARCH...

% Fit the model
try
    [coeff,errors,LLF,innovations,sigmas,summary] = garchfit(spec,y);
catch emsg
    disp('GARCH FIT FAILED BIG TIME');
    out = NaN; return
end

if all(isnan(innovations))
    disp('GARCH FIT FAILED');
    out = NaN; return
end

%% (4) Return statistics on fit

% (i) Return coefficients, and their standard errors as seperate statistics
% ___Mean_Process___
% --C--
% C=NaN --> C=0;

%   --AR--
if isfield(coeff,'AR')
    for i = 1:length(coeff.AR)
        eval(['out.AR_' num2str(i) ' = coeff.AR(' num2str(i) ');']);
        eval(['out.ARerr_' num2str(i) ' = errors.AR(' num2str(i) ');']);
    end 
end

%  --MA--
if isfield(coeff,'MA')
    for i = 1:length(coeff.MA)
        eval(['out.MA_' num2str(i) ' = coeff.MA(' num2str(i) ');']);
        eval(['out.MAerr_' num2str(i) ' = errors.MA(' num2str(i) ');']);
    end 
end

% __Variance_Process___
% -- K --
if isfield(coeff,'K')
    out.K = coeff.K;
end

% -- GARCH --
if isfield(coeff,'GARCH')
    for i = 1:length(coeff.GARCH)
        eval(['out.GARCH_' num2str(i) ' = coeff.GARCH(' num2str(i) ');']);
        eval(['out.GARCHerr_' num2str(i) ' = errors.GARCH(' num2str(i) ');']);
    end
end

% -- ARCH --
if isfield(coeff,'ARCH')
    for i = 1:length(coeff.ARCH)
        eval(['out.ARCH_' num2str(i) ' = coeff.ARCH(' num2str(i) ');']);
        eval(['out.ARCHerr_' num2str(i) ' = errors.ARCH(' num2str(i) ');']);
    end
end


% More statistics given from the fit

out.LLF = LLF; % log-likelihood function

out.summaryexitflag = summary.exitFlag; % whether the fit worked ok.
% This is just a record, really, since the numerical values are only
% symbolic.

nparams = garchcount(coeff); % number of parameters
out.nparams = nparams;

% use aicbic function
[AIC, BIC] = aicbic(LLF,nparams,N); % aic and bic of fit
out.aic = AIC;
out.bic = BIC;

%% Statistics on sigmas
% sigmas is the time series of conditional variance. For a time series with
% strong ARCH/GARCH effects, this will fluctuate; otherwise will be
% quite flat...
out.maxsigma = max(sigmas);
out.minsigma = min(sigmas);
out.rangesigma = max(sigmas)-min(sigmas);
out.stdsigma = std(sigmas);
out.meansigma = mean(sigmas);

%% Check residuals
e = innovations; % 'innovations' -- departures from mean process
stde = e./sigmas; % standardized residuals
stde2 = stde.^2;

% (i) Engle's ARCH test
%       look at autoregressive lags 1:20
%       use the 10% significance level
[Engle_h_stde, Engle_pValue_stde, Engle_stat_stde, Engle_cValue_stde] = archtest(stde,1:20,0.1);

% (ii) Ljung-Box Q-test
%       look at autocorrelation at lags 1:20
%       use the 10% significance level
%       departure from randomness hypothesis test
[lbq_h_stde2, lbq_pValue_stde2, lbq_stat_stde2, lbq_cValue_stde2] = lbqtest(stde2,1:20,0.1,[]);


% Ok, so now we've corrected for GARCH effects, how does this 'improve' the
% randomness/correlation in our signal. If the signal is much less
% correlated now, it is a signature that GARCH effects were significant in
% the original signal

% Mean improvement in Engle pValue
out.engle_mean_diff_p = mean(Engle_pValue_stde - Engle_pValue_y);
out.engle_max_diff_p = max(Engle_pValue_stde - Engle_pValue_y);

% Mean improvement in lbq pValue for squared time series
out.lbq_mean_diff_p = mean(lbq_pValue_stde2 - lbq_pValue_y2);
out.lbq_max_diff_p = max(lbq_pValue_stde2 - lbq_pValue_y2);

% Raw values:
out.engle_pval_stde_1 = Engle_pValue_stde(1);
out.engle_pval_stde_5 = Engle_pValue_stde(5);
out.engle_pval_stde_10 = Engle_pValue_stde(10);
out.minenglepval_stde = min(Engle_pValue_stde);
out.maxenglepval_stde = max(Engle_pValue_stde);

out.lbq_pval_stde_1 = lbq_pValue_stde2(1);
out.lbq_pval_stde_5 = lbq_pValue_stde2(5);
out.lbq_pval_stde_10 = lbq_pValue_stde2(10);
out.minlbqpval_stde2 = min(lbq_pValue_stde2);
out.maxlbqpval_stde2 = max(lbq_pValue_stde2);


% (iii) Correlation in time series: autocorrelation
% autocorrs_y = CO_autocorr(y,1:20);
% autocorrs_var = CO_autocorr(y.^2,1:20);
% [ACF_y,Lags_acf_y,bounds_acf_y] = autocorr(e,20,[],[]);
% [ACF_var,Lags_acf_var,bounds_acf_var] = autocorr(e.^2,20,[],[]);

% (iv) Partial autocorrelation function: PACF
% [PACF,Lags_pacf,bounds_pacf] = parcorr(e,20,[],[]);


% Use MF_residanal on the standardized innovations
% 1) Get statistics on standardized innovations
residout = MF_residanal(stde);

% convert these to local outputs in quick loop
fields = fieldnames(residout);
for k=1:length(fields);
    eval(['out.stde_' fields{k} ' = residout.' fields{k} ';']);
end


out.ac1_stde2 = CO_autocorr(stde2,1);
out.diff_ac1 = CO_autocorr(y.^2,1) - CO_autocorr(stde2,1);


%% (5) Comparison to other models
% e.g., does the additional heteroskedastic component improve the model fit
% over just the conditional mean component of the model.




end