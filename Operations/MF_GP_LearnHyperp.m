% MF_GP_LearnHyperp
% 
% Function used by main Gaussian Process model fitting operations that learns
% Gaussian Process hyperparameters for the time series.
% 
% References code 'minimize' from the GAUSSIAN PROCESS REGRESSION AND
% CLASSIFICATION Toolbox version 3.2, which is avilable at:
% http://gaussianprocess.org/gpml/code
% 
% INPUTS:
% 
% covfunc,       the covariance function, formated as gpml likes it
% nfevals,       the number of function evaluations
% t,             time
% y,             data
% init_loghyper, inital values for hyperparameters
% 
% ------------------------------------------------------------------------------
% Copyright (C) 2013,  Ben D. Fulcher <ben.d.fulcher@gmail.com>,
% <http://www.benfulcher.com>
%
% If you use this code for your research, please cite:
% B. D. Fulcher, M. A. Little, N. S. Jones., "Highly comparative time-series
% analysis: the empirical structure of time series and their methods",
% J. Roy. Soc. Interface 10(83) 20130048 (2010). DOI: 10.1098/rsif.2013.0048
%
% This function is free software: you can redistribute it and/or modify it under
% the terms of the GNU General Public License as published by the Free Software
% Foundation, either version 3 of the License, or (at your option) any later
% version.
% 
% This program is distributed in the hope that it will be useful, but WITHOUT
% ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
% FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
% details.
% 
% You should have received a copy of the GNU General Public License along with
% this program.  If not, see <http://www.gnu.org/licenses/>.
% ------------------------------------------------------------------------------

function loghyper = MF_GP_LearnHyperp(covfunc,nfevals,t,y,init_loghyper)
% Ben Fulcher, 2010
    
if nargin < 5 || isempty(init_loghyper)
    % Use default starting values for parameters
    % How many hyperparameters
    s = feval(covfunc{:}); % string in form '2+1', ... tells how many
    % hyperparameters for each contribution to the covariance function
    nhps = eval(s);
    init_loghyper = -1*ones(nhps,1); % Initialize all log hyperparameters at -1
end
%         init_loghyper(1) = log(mean(diff(t)));
    
% Perform the optimization
try
    loghyper = minimize(init_loghyper, 'gpr', nfevals, covfunc, t, y);
catch emsg
    switch emsg.identifier
    case 'MATLAB:posdef'
        fprintf(1,'Error: lack of positive definite matrix for this function');
        loghyper = NaN; return
    case 'MATLAB:nomem'
        error('Out of memory');
    otherwise
        error('Error fitting Gaussian Process to data')
    end
end

end