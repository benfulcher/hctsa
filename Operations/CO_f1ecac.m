function out = CO_f1ecac(y)
% CO_f1ecac     The 1/e correlation length
%
% Finds where autocorrelation function first crosses 1/e
%
%---INPUT:
% y, the input time series.

% ------------------------------------------------------------------------------
% Copyright (C) 2018, Ben D. Fulcher <ben.d.fulcher@gmail.com>,
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

N = length(y); % time-series length
thresh = 1/exp(1); % 1/e threshold

% With Fourier method you can compute the whole spectrum at once:
acf = CO_AutoCorr(y,[],'Fourier');

% First crossing of 1/e (always positive to negative, because a(1)=1)
firstCrossing = find((acf - thresh < 0),1,'first');

if ~isempty(firstCrossing)
    out = firstCrossing;
else
    out = N; % no crossing point anywhere across full ACF
end

%-------------------------------------------------------------------------------
% More general crossing events:
% crossedThreshold = ((acf(1:end-1)-thresh).*(acf(2:end)-thresh) < 0);
% if any(crossedThreshold)
    % out = find(crossedThreshold,1,'first');
% else
    % out = N; % no crossing point anywhere across full ACF
% end

%-------------------------------------------------------------------------------
% Plot to check:
% f = figure('color','w');
% hold on
% plot(0:length(acf)-1,acf,'o-k')
% plot([0,out],thresh*ones(2,1),'--r')

end
