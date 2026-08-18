function out = IN_MultivariateAMI(y, tauMethod, estMethod, extraParam)
% IN_MultivariateAMI   Multivariate automutual information, I(x_t; x_{t-tau}, x_{t-2tau})
%
% Measures how much information the past two points (spaced tau apart, with tau
% estimated from the series' own autocorrelation structure) jointly carry about
% the present -- and, via the synergy field, whether that joint information
% exceeds or falls short of the sum of the pairwise automutual informations at
% tau and 2*tau individually. Positive synergy means the two lagged points are
% jointly more informative than either alone would suggest (interaction
% structure); negative synergy means x_{t-2tau} is mostly redundant with
% x_{t-tau} once you already know it.
%
% ---INPUTS:
%
% y: input time series (column vector, expected z-scored)
%
% tauMethod: how to select the time delay, tau (as for BF_Embed):
%           (*) 'ac' (default): first zero-crossing of the autocorrelation function
%           (*) 'mi': first minimum of the (gaussian) automutual information
%           (*) a fixed positive integer
%
% estMethod: the estimation method used to compute the mutual information:
%           (*) 'gaussian' (default): closed-form, via the multiple correlation
%               coefficient -- fast, but blind to nonlinear/non-Gaussian structure
%           (*) 'kraskov2': nonparametric KSG estimator (Information Dynamics Toolkit)
%
% extraParam: number of nearest neighbors for the Kraskov estimator (default '4',
%             cf. IN_Initialize_MI.m)
%
% cf. Kraskov, A., Stoegbauer, H., Grassberger, P., Estimating mutual
% information: http://dx.doi.org/10.1103/PhysRevE.69.066138

% ------------------------------------------------------------------------------
% Copyright (C) 2013-2026, Ben D. Fulcher <ben.d.fulcher@gmail.com>,
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
if nargin < 2 || isempty(tauMethod)
	tauMethod = 'ac';
end
if nargin < 3 || isempty(estMethod)
	estMethod = 'gaussian';
end
if nargin < 4
	extraParam = [];
end

minSamples = 20; % need reasonably many samples for a stable bivariate MI estimate

% Ensure y is a column vector
if size(y, 2) > size(y, 1)
	y = y';
end

% ------------------------------------------------------------------------------
% Resolve tau and build the joint [x_{t-2tau}, x_{t-tau}, x_t] embedding in one
% step (reusing BF_Embed's own tau-resolution and windowing, exactly as used
% throughout hctsa's nonlinear/embedding-based operations):
% ------------------------------------------------------------------------------
y_embed = BF_Embed(y, tauMethod, 3);
if isnan(y_embed(1)) || size(y_embed, 1) < minSamples
	out = NaN; return
end

x_2tau = y_embed(:, 1); % x_{t-2tau}
x_tau = y_embed(:, 2);  % x_{t-tau}
x_now = y_embed(:, 3);  % x_t

if strcmp(tauMethod, 'ac')
	tau = CO_FirstCrossing(y, 'ac', 0, 'discrete');
elseif strcmp(tauMethod, 'mi')
	tau = CO_FirstMin(y, 'mi');
else
	tau = tauMethod;
end
out.tau = tau;

% ------------------------------------------------------------------------------
% Compute the joint and pairwise automutual informations, all on the identical
% aligned set of samples (for a fair synergy comparison):
% ------------------------------------------------------------------------------
switch estMethod
	case 'gaussian'
		% Multivariate Gaussian MI: I(Y;X) = -0.5*log(1 - R^2), where R^2 is the
		% squared multiple correlation coefficient of Y on X (exact generalization
		% of the univariate -0.5*log(1-r^2) formula used in IN_AutoMutualInfo.m):
		R = corr([x_now, x_2tau, x_tau]);
		corrVecXY = R(1, 2:3)';
		CorrXX = R(2:3, 2:3);
		Rsq = corrVecXY' * (CorrXX \ corrVecXY);
		Rsq = min(max(Rsq, 0), 1 - eps); % guard numerical over/undershoot
		out.multiAMI = -0.5 * log(1 - Rsq);

		out.ami_2tau = -0.5 * log(1 - R(1, 2)^2);
		out.ami_tau = -0.5 * log(1 - R(1, 3)^2);

	case {'kraskov1', 'kraskov2'}
		miCalc = IN_Initialize_MI(estMethod, extraParam, false, y); % tie-break-protected, no other added noise
		miCalc.initialise(2, 1);
		miCalc.setObservations([x_2tau, x_tau], x_now);
		out.multiAMI = miCalc.computeAverageLocalOfObservations();

		miCalc.initialise(1, 1);
		miCalc.setObservations(x_2tau, x_now);
		out.ami_2tau = miCalc.computeAverageLocalOfObservations();

		miCalc.initialise(1, 1);
		miCalc.setObservations(x_tau, x_now);
		out.ami_tau = miCalc.computeAverageLocalOfObservations();

	otherwise
		error('Unknown estimation method ''%s''', estMethod);
end

% The part of the joint information NOT explained by the two pairwise AMIs
% individually -- positive means synergistic (jointly more informative than
% the sum of parts), negative means redundant (x_{t-2tau} mostly duplicates
% what x_{t-tau} already tells you about x_t):
out.synergy = out.multiAMI - out.ami_tau - out.ami_2tau;

end
