function out = IN_AutoMutualInfo(y,timeDelay,estMethod,extraParam)
% IN_AutoMutualInfo     Time-series automutual information
%
%---INPUTS:
%
% y: input time series
%
% timeDelay: time lag for automutual information calculation
%
% estMethod: the estimation method used to compute the mutual information:
%           (*) 'gaussian'
%           (*) 'kernel'
%           (*) 'kraskov1'
%           (*) 'kraskov2'
%
% cf. Kraskov, A., Stoegbauer, H., Grassberger, P., Estimating mutual
% information: http://dx.doi.org/10.1103/PhysRevE.69.066138

% ------------------------------------------------------------------------------
% Copyright (C) 2020, Ben D. Fulcher <ben.d.fulcher@gmail.com>,
% <http://www.benfulcher.com>
%
% If you use this code for your research, please cite the following two papers:
%
% (1) B.D. Fulcher and N.S. Jones, "hctsa: A Computational Framework for Automated
% Time-Series Phenotyping Using Massive Feature Extraction, Cell Systems 5: 527 (2017).
% DOI: 10.1016/j.cels.2017.10.001
%
% (2) B.D. Fulcher, M.A. Little, N.S. Jones, "Highly comparative time-series
% analysis: the empirical structure of time series and their methods",
% J. Roy. Soc. Interface 10(83) 20130048 (2013).
% DOI: 10.1098/rsif.2013.0048
%
% This function is free software: you can redistribute it and/or modify it under
% the terms of the GNU General Public License as published by the Free Software
% Foundation, either version 3 of the License, or (at your option) any later
% version.
%
% This program is distributed in the hope that it will be useful, but WITHOUT
% ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
% FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
% details.
%
% You should have received a copy of the GNU General Public License along with
% this program. If not, see <http://www.gnu.org/licenses/>.
% ------------------------------------------------------------------------------

% ------------------------------------------------------------------------------
% Check inputs:
% ------------------------------------------------------------------------------
if nargin < 2 || isempty(timeDelay)
    timeDelay = 1;
end
if ischar(timeDelay) && ismember(timeDelay,{'ac','tau'})
    timeDelay = CO_FirstCrossing(y,'ac',0,'discrete');
end

if nargin < 3 || isempty(estMethod)
    estMethod = 'kernel';
end

if nargin < 4
    extraParam = [];
end

N = length(y);
doPlot = false; % plot outputs to screen
minSamples = 5; % minimum 5 samples to compute a mutual information (could make higher?)

% ------------------------------------------------------------------------------
% Loop over time delays if a vector
numTimeDelays = length(timeDelay);
amis = nan(numTimeDelays,1);

if numTimeDelays > 1
    timeDelay = sort(timeDelay);
end

% Initialize miCalc object (needs to be reinitialized within the loop for kraskov):
if ~strcmp(estMethod,'gaussian')
    miCalc = IN_Initialize_MI(estMethod,extraParam,false); % NO ADDED NOISE!
end

for k = 1:numTimeDelays

    % Check enough samples to compute an automutual information
    if timeDelay(k) > N - minSamples
        % Time series is too short -- keep the remaining values as NaNs
        break
    end

    % Form the time-delay vectors y1 and y2
    y1 = y(1:end-timeDelay(k));
    y2 = y(1+timeDelay(k):end);

    if strcmp(estMethod,'gaussian')
        r = corr(y1,y2,'type','Pearson');
        amis(k) = -0.5*log(1-r^2);
    else
        % Reinitialize for Kraskov:
        miCalc.initialise(1,1);

        % Set observations to time-delayed versions of the time series:
        miCalc.setObservations(y1,y2);

        % Compute:
        amis(k) = miCalc.computeAverageLocalOfObservations();
    end

    % Plot:
    if doPlot
        plot(y1,y2,'.k')
        title(sprintf('ami = %.3f',amis(k)))
        pause(0.1)
    end
end

if any(isnan(amis))
    warning(['Time series (N=%u) is too short for automutual information calculations',...
                ' up to lags of %u'],N,max(timeDelay))
end

if doPlot
    plot(amis,'-k')
end

%-------------------------------------------------------------------------------
% Outputs:
%-------------------------------------------------------------------------------
if numTimeDelays == 1
    out = amis; % a scalar
else
    % A structure
    for k = 1:numTimeDelays
        out.(sprintf('ami%u',timeDelay(k))) = amis(k);
    end
end

end
