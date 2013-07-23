% GP_LearnHyperp
% 
% Function used by main Gaussian Process model fitting operations that learns
% Gaussian Process hyperparameters for the time series
% 
% INPUTS:
% covfunc, the covariance function, formated as gpml likes it
% nfevals, the number of function evaluations
% t, time
% y, data
% init_loghyper, inital values for hyperparameters
% 

function loghyper = GP_LearnHyperp(covfunc,nfevals,t,y,init_loghyper)
% Ben Fulcher
    
if nargin < 5 || isempty(init_loghyper)
    % Use default starting values for parameters
    % How many hyperparameters
    s = feval(covfunc{:}); % string in form '2+1', ... tells how many
    % hyperparameters for each contribution to the covariance function
    nhps = eval(s);
    init_loghyper = ones(nhps,1)*-1; % Initialize all log hyperparameters at -1
end
%         init_loghyper(1) = log(mean(diff(t)));
    
% Perform the optimization
try
    loghyper = minimize(init_loghyper, 'gpr', nfevals, covfunc, t, y);
catch emsg
    if strcmp(emsg.identifier,'MATLAB:posdef')
        fprintf(1,'Error: lack of positive definite matrix for this function');
        loghyper = NaN; return
    elseif strcmp(emsg.identifier,'MATLAB:nomem')
        error('Out of memory');
    else
        error('Error fitting Gaussian Process to data')
    end
end

end