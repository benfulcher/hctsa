function out = GP_predict(y,covfunc,ntrain,ntest,npreds,pmode)
% Use Rasmussun GP code (from gaussianprocess.org) to optimize
% hyperparameters and predict the next data point
% Ben Fulcher 20/1/2010

doplot = 0;
%% Preliminaries
% length of the time series
N = length(y);

%% Check Inputs
if nargin < 2 || isempty(covfunc),
    fprintf(1,'Using a default covariance function: sum of squared exponential and noise\n')
    covfunc = {'covSum', {'covSEiso','covNoise'}};
end

if nargin < 3 || isempty(ntrain)
    ntrain = 20; % 20 previous data points to predict the next
end

if nargin < 4 || isempty(ntest)
    ntest = 5; % test on 5 data points into the future
end

if nargin < 5 || isempty(npreds) % number of predictions
    npreds = 10; % do it 10 times (equally-spaced) through the time series
end

if nargin < 6 || isempty(pmode)
    pmode = 'frombefore'; % predicts from previous ntrain datapoints
    % can also be 'randomgap' -- fills in random gaps in the middle of a string of
    % data
end

%% Set up loop
if ismember(pmode,{'frombefore','randomgap'})
    spns = floor(linspace(1,N-(ntest+ntrain),npreds)); % starting positions
elseif strcmp(pmode,'beforeafter')
    spns = floor(linspace(1,N-(ntest+ntrain*2),npreds)); % starting positions
end


mus = zeros(ntest,npreds); % predicted values
stderrs = zeros(ntest,npreds); % standard errors on predictions
yss = zeros(ntest,npreds); % test values
mlikelihoods = zeros(npreds,1); % marginal likelihoods of model

nhps = eval(feval(covfunc{:})); % number of hyperparameters
loghypers = zeros(nhps,npreds); % loghyperparameters

for i = 1:npreds
    %% (0) Set up test and training sets
    if strcmp(pmode,'frombefore')
        tt = (1:ntrain)'; % times (make from 1)
        rt = spns(i):spns(i)+ntrain-1; % training range
        yt = y(rt); % training data
        
        ts = (ntrain+1 : ntrain+1 + ntest-1)'; % times
        rs = spns(i)+ntrain : spns(i)+ntrain + ntest-1; % test range
        ys = y(rs); % test data
        
    elseif strcmp(pmode,'randomgap')
        t = (1:ntrain+ntest)';
        r = randperm(ntrain+ntest);
        yy = y(spns(i):spns(i)+ntrain+ntest-1);
        
        rt = sort(r(1:ntrain),'ascend');
        tt = t(rt);
        yt = yy(rt);
        
        rs = sort(r(ntrain+1:end),'ascend');
        ts = t(rs);
        ys = yy(rs);
        
    elseif strcmp(pmode,'beforeafter')
        t = (1:ntrain*2+ntest)';
        r = ntrain+1:ntrain+1 + ntest-1;
        yy = y(spns(i):spns(i)+2*ntrain+ntest-1);
        
        rt = [1:ntrain ntrain+ntest+1:ntrain*2+ntest];
        tt = t(rt);
        yt = yy(rt);
        
        rs = [ntrain+1 : ntrain+ntest];
        ts = t(rs);
        ys = yy(rs);
        
        
    end
    
    % Process to normalize scales
    ys = (ys-mean(yt))/std(yt); % same transformation as training set
    yt = (yt-mean(yt))/std(yt); % zscore training set
    
    %% (1) Learn hyperparameters from training set (t)
    
    loghyper = SUB_learnhyper(covfunc,-50,tt,yt);
    
    if isnan(loghyper)
        disp('Couldn''t learn hyperparameters');
        out = NaN; return
    end
    
    loghypers(:,i) = loghyper;
    
    % Get marginal likelihood for this model with hyperparameters optimized
    % over training data
    mlikelihoods(i) = - gpr(loghyper, covfunc, tt, yt);

    
    %% (2) Evaluate at test set (s)
    
    % evaluate at test points based on training time/data, predicting for
    % test times/data
    [mu S2] = gpr(loghyper, covfunc, tt, yt, ts);
    
    % Compare to actual test data --> store in row of errs
    mus(:,i) = mu; % ~predicted values for time series points
    stderrs(:,i) = 2*sqrt(S2); % ~errors on those predictions
    yss(:,i) = ys;
    
    % Plot
    if doplot
        if strcmp(pmode,'frombefore')
            plot(tt,yt,'.-k');
            hold on;
            plot(ts,ys,'.-b');
            errorbar(ts,mu,2*sqrt(S2),'m');
            hold off;
        else
            plot(tt,yt,'ok');
            hold on;
            plot(ts,ys,'ob');
            errorbar(ts,mu,2*sqrt(S2),'m');
            hold off;
        end
        keyboard
    end

%     for j=1:ntest
%         % set up structure output
%         err = abs(mu(j)-ys(j))/sqrt(S2(j)); % in units of std at this point
%         eval(['out.abserr' num2str(i) '_' num2str(j) ' = err;']);
%     end
    
end

% Ok, we're done.

%% Return statistics on how well it did
%% (1) PREDICTION ERROR MEASURES
allstderrs = abs(mus-yss)./stderrs; % differences between predictions and actual
                              % in units of standard error (95% confidence
                              % interval error bars)
allabserrs = abs(mus-yss);


% largest error:
out.maxstderr = max(allstderrs(:));
out.maxabserr = max(allabserrs(:));
% smallest error:
out.minstderr = min(allstderrs(:));
out.minabserr = min(allabserrs(:));
% mean error (across all):
out.meanstderr = mean(allstderrs(:));
out.meanabserr = mean(allabserrs(:));

% Summary of how it did on each run:
stderr_run = mean(abs(mus-yss)./stderrs); % per run
abserr_run = mean(abs(mus-yss));

% Mean error for a run
out.meanstderr_run = mean(stderr_run);
out.meanabserr_run = mean(abserr_run);

% Max error for a run
out.maxstderr_run = max(stderr_run);
out.maxabserr_run = max(abserr_run);

% Min error for a run
out.minstderr_run = min(stderr_run);
out.minabserr_run = min(abserr_run);

% Error bar stats
out.maxerrbar = max(stderrs(:)); % largest error bar
out.meanerrbar = mean(stderrs(:)); % mean error bar length
out.minerrbar = min(stderrs(:)); % minimum error bar length


%% (2) HYPERPARAMETER MEASURES
% mean and std for each hyperparameter
for i=1:nhps
    o1 = mean(loghypers(i,:));
    eval(['out.meanlogh' num2str(i) ' = o1;']);
    o2 = std(loghypers(i,:));
    eval(['out.stdlogh' num2str(i) ' = o2;']);
end

%% (3) Marginal likelihood measures
% Best marginal neg-log-likelihood attained
% Worst marginal neg-log-likelihood attained
% spread in marginal neg-log-likelihoods

out.maxmlik = max(mlikelihoods);
out.minmlik = min(mlikelihoods);
out.stdmlik = std(mlikelihoods);

%% Subfunctions

    function loghyper = SUB_learnhyper(covfunc,nfevals,t,y,init_loghyper)
        % nfevals--  negative: specifies maximum number of allowed
        % function evaluations
        % t: time
        % y: data
        
        if nargin < 5 || isempty(init_loghyper)
            % Use default starting values for parameters
            % How many hyperparameters
            s = feval(covfunc{:}); % string in form '2+1', ... tells how many
            % hyperparameters for each contribution to the
            % covariance function
            nhps = eval(s);
            init_loghyper = ones(nhps,1)*-1; % Initialize all log hyperparameters at -1
        end
%         init_loghyper(1) = log(mean(diff(t)));
        
        % Perform the optimization
        try
            loghyper = minimize(init_loghyper, 'gpr', nfevals, covfunc, t, y);
        catch emsg
            if strcmp(emsg.identifier,'MATLAB:posdef')
                disp('Error with lack of positive definite matrix for this function');
                loghyper = NaN; return
            elseif strcmp(emsg.identifier,'MATLAB:nomem')
                error('Out of memory fitting this Gaussian Process');
                % return as if a fatal error -- come back to this.
            else
                error(['Error fitting Gaussian Process: ' emsg.message])
            end
        end
        %        hyper = exp(loghyper);
        
    end

%     function loghyper = SUB_learnhyper(covfunc,nfevals,t,y)
%         % nfevals--  negative: specifies maximum number of allowed
%         % function evaluations
%         % t: time
%         % y: data
% 
%        % How many hyperparameters
%        s = feval(covfunc{:}); % string in form '2+1', ... tells how many
%        % hyperparameters for each contribution to the
%        % covariance function
%        nhps = eval(s);
%        init_loghyper = ones(nhps,1)*-1; % Initialize all log hyperparameters at -1
%        
%        % Perform the optimization
%        try
%            loghyper = minimize(init_loghyper, 'gpr', nfevals, covfunc, t, y);
%        catch emsg
%            if strcmp(emsg.identifier,'MATLAB:posdef')
%                disp('Error with lack of positive definite matrix for this function');
%                loghyper=NaN; return
%            elseif strcmp(emsg.identifier,'MATLAB:nomem')
%                disp('Out of memory. Leaving');
%                % return as if a fatal error -- come back to this.
%                return
%            else
%               keyboard 
%            end
%        end
% %        hyper = exp(loghyper); 
%         
%     end
    

end