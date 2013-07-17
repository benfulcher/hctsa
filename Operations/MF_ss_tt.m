function out = MF_ss_tt(y,model,order,howtosubset,samplep)
% Looks at robustness of training sets to fitted output parameters of a
% statespace model of given order (or 'best' order). Look at spread of
% parameters obtained (including in-sample goodness of fit statistics) --
% some indication of stationarity. Look at values of goodness of fit --
% some indication of model suitability.
% n.b., This code inherits strongly from this MF_ss_testset
% Ben Fulcher 12/2/2010


%% Preliminaries
N = length(y); % length of time series

%% Inputs

% (1) y: column vector time series
if nargin < 1 || isempty(y)
    disp('Give us a time series, ya mug'); return
end
% Convert y to time series object
y = iddata(y,[],1);

% (2) model, the type of model to fit
if nargin < 2 || isempty(model)
    model = 'ss'; % fit a state space model by default
end
    
% (3) order of model, order
if nargin < 3 || isempty(order)
    order = 2; % model of order 2 by default. Not very good defaults.
end

% (4) How to choose subsets from the time series, howtosubset
if nargin < 4 || isempty(howtosubset)
    howtosubset = 'rand'; % takes segments randomly from time series
end

% (5) Sampling parameters, samplep
if nargin < 5 || isempty(samplep)
    samplep = [20, 0.1]; % sample 20 times with 10%-length subsegments
end

% % (6) Predict some number of steps ahead in test sets, steps
% if nargin < 6 || isempty(steps)
%     steps = 2; % default: predict 2 steps ahead in test set
% end


%% Set the ranges beforehand
% Number of samples to take, npred
npred = samplep(1);
r = zeros(npred,2); % ranges

switch howtosubset
    case 'rand'
        if samplep(2) < 1 % specified a fraction of time series
            l = floor(N*samplep(2));
        else % specified an absolute interval
            l = samplep(2);
        end
        spts = randi(N-l+1,npred,1); % npred starting points
        r(:,1) = spts;
        r(:,2) = spts+l-1;
    case 'uniform'
        if length(samplep) == 1 % size will depend on number of unique subsegments
            spts = round(linspace(0,N,npred+1)); % npred+1 boundaries = npred portions
            r(:,1) = spts(1:npred)+1;
            r(:,2) = spts(2:end);
        else
            if samplep(2) < 1 % specified a fraction of time series
                l = floor(N*samplep(2));
            else % specified an absolute interval
                l = samplep(2);
            end
            spts = round(linspace(1,N-l+1,npred)); % npred+1 boundaries = npred portions
            r(:,1) = spts;
            r(:,2) = spts+l-1;
        end
    otherwise
        error('Unknown subset method ''%s''',howtosubset);
end

%% Fit the model to each training set
% model will be stored as a model object, m
% model is fitted using the entire dataset as the training set
% test sets will be smaller chunks of this.
% [could also fit multiple models using data not in multiple test sets, but
% this is messier]
switch model
    case 'arsbc'
        %% Fit AR models of 'best' order according to SBC, using arfit package
        % fit AR models of 'best' order, return statistics on how this best
        % order changes. The order input argument is not used for this
        % option.
        orders = zeros(npred,1);
        sbcs = zeros(npred,1);
        yy = y.y;
        for $1 $2 $3:npred
            % Use arfit software to retrieve the optimum AR(p) order by
            % Schwartz's Bayesian Criterion, SBC (or BIC); in the range
            % p = 1-10
            % Enforce zero mean level. This could be relaxed.
            try
                [west, Aest, Cest, SBC] = arfit(yy(r(i,1):r(i,2)), 1, 10, 'sbc', 'zero');
            catch emsg
                if strcmp(emsg.message,'Time series too short.')
                   fprintf(1,'Time Series is too short for ARFIT\n');
                   out = NaN; return
                end
            end
            orders(i) = length(Aest);
            sbcs(i) = min(SBC);
        end
        % Return statistics
        out.orders_mode = mode(orders);
        out.orders_mean = mean(orders);
        out.orders_std = std(orders);
        out.orders_max = max(orders);
        out.orders_min = min(orders);
        out.orders_range = range(orders);
        
        out.sbcs_mean = mean(sbcs);
        out.sbcs_std = std(sbcs);
        out.sbcs_range = range(sbcs);
        out.sbcs_min = min(sbcs);
        out.sbcs_max = max(sbcs);
        
    case 'ar'
        %% Fit AR model of specified order
        % Return statistics on parameters and goodness of fit
        fpes = zeros(npred,1);
        as = zeros(npred,order+1);
        for i = 1:npred
            % fit the ar model
            m = ar(y(r(i,1):r(i,2)),order);
            % get parameters and goodness of fit
            fpes(i) = m.EstimationInfo.FPE;
            as(i,:) = m.a;
        end
        
        % statistics on FPE
        out.fpe_std = std(fpes);
        out.fpe_mean = mean(fpes);
        out.fpe_max = max(fpes);
        out.fpe_min = min(fpes);
        out.fpe_range = range(fpes);
        
        % Statistics on fitted AR parameters
        for i = 1:order % first column will be ones
            eval(sprintf('out.a_%u_std = std(as(:,%u+1));',i,i));
            eval(sprintf('out.a_%u_mean = mean(as(:,%u+1));',i,i));
            eval(sprintf('out.a_%u_max = max(as(:,%u+1));',i,i));
            eval(sprintf('out.a_%u_min = min(as(:,%u+1));',i,i));
        end
    case 'ss'
        %% Fit state space models of specified order
        % Return statistics on goodness of fit
        % Could do parameters too, but I this would involve many outputs
        fpes = zeros(npred,1);
        for i = 1:npred
            try  m = n4sid(y(r(i,1):r(i,2)),order);
            catch
                % Some range of the time series is invalid for fitting the
                % model to.
                error('Couldn''t fit this state space model')
            end
            fpes(i) = m.EstimationInfo.FPE;
        end
        
        % statistics on FPE
        out.fpe_std = std(fpes);
        out.fpe_mean = mean(fpes);
        out.fpe_max = max(fpes);
        out.fpe_min = min(fpes);
        out.fpe_range = range(fpes);
        
    case 'arma'
        %% fit an ARMA model of specified order(s)
        % Note: order should be a two-component vector
        % Output parameters and goodness of fit
        fpes = zeros(npred,1);
        ps = zeros(npred,order(1)+1);
        qs = zeros(npred,order(2)+1);
        
        for i = 1:npred
            try
                m = armax(y(r(i,1):r(i,2)),order);
            catch emsg
                error('Couldn''t fit this ARMA model')
            end
            fpes(i) = m.EstimationInfo.FPE;
            ps(i,:) = m.a;
            qs(i,:) = m.c;
        end
        
        % statistics on FPE
        out.fpe_std = std(fpes);
        out.fpe_mean = mean(fpes);
        out.fpe_max = max(fpes);
        out.fpe_min = min(fpes);
        out.fpe_range = range(fpes);
        
        % Statistics on fitted AR parameters, p
        for i = 1:order % first column will be ones
            eval(sprintf('out.p_%u_std = std(ps(:,%u+1));',i,i));
            eval(sprintf('out.p_%u_mean = mean(ps(:,%u+1));',i,i));
            eval(sprintf('out.p_%u_max = max(ps(:,%u+1));',i,i));
            eval(sprintf('out.p_%u_min = min(ps(:,%u+1));',i,i));
        end
        
        % Statistics on fitted MA parameters, q
        for i = 1:order % first column will be ones
            eval(sprintf('out.q_%u_std = std(qs(:,%u+1));',i,i));
            eval(sprintf('out.q_%u_mean = mean(qs(:,%u+1));',i,i));
            eval(sprintf('out.q_%u_max = max(qs(:,%u+1));',i,i));
            eval(sprintf('out.q_%u_min = min(qs(:,%u+1));',i,i));
        end
    otherwise
        error('Unknown model ''%s''',model);
end


end