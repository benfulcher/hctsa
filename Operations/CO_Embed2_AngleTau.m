function out = CO_Embed2_AngleTau(y,maxTau)
% CO_Embed2_AngleTau
%
% Investigates how the autocorrelation of angles between successive points in
% the two-dimensional time-series embedding change as tau varies from
% tau = 1, 2, ..., maxTau.
%
%---INPUTS:
% y, a column vector time series
% maxTau, the maximum time lag to consider

% ------------------------------------------------------------------------------
% Copyright (C) 2016, Ben D. Fulcher <ben.d.fulcher@gmail.com>,
% <http://www.benfulcher.com>
%
% If you use this code for your research, please cite:
% B. D. Fulcher, M. A. Little, N. S. Jones, "Highly comparative time-series
% analysis: the empirical structure of time series and their methods",
% J. Roy. Soc. Interface 10(83) 20130048 (2013). DOI: 10.1098/rsif.2013.0048
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

doPlot = 0;
tauRange = (1:1:maxTau);
numTau = length(tauRange);

% Ensure y is a column vector
if size(y,2) > size(y,1);
	y = y';
end

stats_store = zeros(3,numTau);

for i = 1:numTau
	tau = tauRange(i);

	m = [y(1:end-tau), y(1+tau:end)];

	theta = diff(m(:,2))./diff(m(:,1));
	theta = atan(theta); % measured as deviation from the horizontal

	stats_store(1,i) = CO_AutoCorr(theta,1,'Fourier');
	stats_store(2,i) = CO_AutoCorr(theta,2,'Fourier');
	stats_store(3,i) = CO_AutoCorr(theta,3,'Fourier');
end

if doPlot
    figure('color','w'); box('on');
    plot(stats_store');
end

% ------------------------------------------------------------------------------
% Compute lots of outputs statistics:
% ------------------------------------------------------------------------------
out.ac1_thetaac1 = CO_AutoCorr(stats_store(1,:),1,'Fourier');
out.ac1_thetaac2 = CO_AutoCorr(stats_store(2,:),1,'Fourier');
out.ac1_thetaac3 = CO_AutoCorr(stats_store(3,:),1,'Fourier');
out.mean_thetaac1 = mean(stats_store(1,:));
out.max_thetaac1 = max(stats_store(1,:));
out.min_thetaac1 = min(stats_store(1,:));
out.mean_thetaac2 = mean(stats_store(2,:));
out.max_thetaac2 = max(stats_store(2,:));
out.min_thetaac2 = min(stats_store(2,:));
out.mean_thetaac3 = mean(stats_store(3,:));
out.max_thetaac3 = max(stats_store(3,:));
out.min_thetaac3 = min(stats_store(3,:));
out.meanrat_thetaac12 = out.mean_thetaac1/out.mean_thetaac2;
out.diff_thetaac12 = sum(abs(stats_store(2,:)-stats_store(1,:)));

end
