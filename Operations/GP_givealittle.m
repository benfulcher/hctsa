function out = GP_givealittle(y,covfunc,npoints)
% Uses a Gaussian process to predict the time series in between
% equally-spaced points through the time series
% Ben Fulcher 22/1/2010

if nargin < 2 || isempty(covfunc)
    disp('GP_hps: We''re using default sum of SE and noise covariance function')
    covfunc = {'covSum', {'covSEiso','covNoise'}};
end

if nargin < 3 || isempty(npoints)
    npoints = 20;
end

doplot = 0; % set to 1 to visualize behavior

%% Get the points
N = length(y); % length of the time series
tt = floor(linspace(1,N,npoints))'; % time range (training)
yt = y(tt);

%% Optimize the GP parameters for the chosen covariance function

% Determine the number of hyperparameters, nhps
s = feval(covfunc{:}); % string in form '2+1', ... tells how many
% hyperparameters for each contribution to the
% covariance function
nhps = eval(s);


covfunc1 = covfunc{1};
covfunc2 = covfunc{2};
if strcmp(covfunc1,'covSum') && strcmp(covfunc2{1},'covSEiso') && strcmp(covfunc2{2},'covNoise')
    init_loghp = ones(3,1)*-1;
    % length parameter is in the ballpark of the difference between time
    % elements
    init_loghp(1) = log(mean(diff(tt)));
else
    init_loghp = ones(nhps,1)*-1;
end

loghyper = SUB_learnhyper(covfunc,-50,tt,yt,init_loghp);

if isnan(loghyper)
    out = NaN;
    return
end

%% Evaluate over the whole space now
% evaluate at test points based on training time/data, predicting for
% test times/data
if N <= 2000
    ts = (1:N)';
else % memory constraints force us to crudely resample
    ts = round(linspace(1,N,2000))';
end
[mu S2] = gpr(loghyper, covfunc, tt, yt, ts);


%% For Plotting
if doplot
    xstar = linspace(min(t),max(t),1000)';
    [mu S2] = gpr(loghyper, covfunc, t, y, ts);
    S2p = S2 - exp(2*loghyper(3)); % remove noise from predictions
    S2p = S2;
    figure('color','w');
    f = [mu+2*sqrt(S2p);flipdim(mu-2*sqrt(S2p),1)];
    fill([ts; flipdim(ts,1)], f, [6, 7, 7]/8, 'EdgeColor', [7, 7, 6]/8);
            % grayscale error bars
    hold on;
    plot(ts,mu,'k-','LineWidth',2); % mean function
    plot(ts,y(ts),'.-k'); % original data
end


%% Outputs
S = sqrt(S2); % standard deviation function, S
% rms error from mean function, mu
out.rmserr = mean(sqrt((y(ts)-mu).^2));
out.meanstderr = mean(abs(y(ts)-mu)./S);
out.stdmu = std(mu);
out.meanS = mean(S);
out.stdS = std(S);

% Marginal Likelihood
try
    out.mlikelihood = - gpr(loghyper, covfunc, ts, y(ts));
catch
    out.mlikelihood = NaN;
end

% Loghyperparameters
for i = 1:nhps
    % set up structure output
    eval(['out.logh' num2str(i) ' = loghyper(' num2str(i) ');']);
end

if strcmp(covfunc1,'covSum') && strcmp(covfunc2{1},'covSEiso') && strcmp(covfunc2{2},'covNoise')
   % Give extra output based on length parameter on length of time series
   out.h_lonN = exp(loghyper(1))/N;
end



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
                error('GP_givealittle: Out of memory');
                % return as if a fatal error -- come back to this.
            else
                error('Error fitting Gaussian Process to data')
            end
        end
        
    end


end