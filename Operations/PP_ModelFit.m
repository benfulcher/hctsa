function out = PP_ModelFit(y, model, order, randomSeed)
% PP_ModelFit   Investigates whether AR model fit improves with different preprocessings.
%
% After performing the range of transformations to the time series, returns the
% in-sample root-mean-square (RMS) prediction errors for an AR model on each
% transformed time series as a ratio of the RMS prediction error of the original
% time series.
%
% PP_PreProcess.m is used to perform the preprocessings
%
% The AR model is fitted using the function ar and pe from Matlab's System
% Identification Toolbox
%
% Transformations performed include:
% (i) incremental differencing,
% (ii) filtering of the power spectral density function,
% (iii) removal of piece-wise polynomial trends, and
% (iv) rank mapping the values of the time series to a Gaussian distribution.
%
% Only one representative preprocessing per PP_PreProcess.m family is fit
% (d1, d2, p1_20, p2_5, rmgd): a redundancy check on 100 real series (Bonn
% EEG + Empirical1000) found the other candidates in PP_PreProcess.m's
% output each correlate at r>=0.95 with one of these, so fitting all of
% them just repeats the same handful of AR refits.
%
% ---INPUTS:
%
% y, the input time series
% model, the time-series model to fit to the transformed time series (currently
%           'ar' is the only option)
%
% order, the order of the AR model to fit to the data
%
% randomSeed, whether (and how) to reset the random seed, using BF_ResetSeed

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
%% Preliminaries
% ------------------------------------------------------------------------------
N = length(y); % length of the time series

% ------------------------------------------------------------------------------
%% Inputs
% ------------------------------------------------------------------------------
% Model: the model to fit preprocessed time series to
if nargin < 2 || isempty(model)
	model = 'ar';
end

% order: the order of model to fit
if nargin < 3 || isempty(order)
	order = 2;
end

% randomSeed: how to treat the randomization
if nargin < 4
	randomSeed = []; % default
end

% ------------------------------------------------------------------------------
%% Apply a range of preprocessings
% ------------------------------------------------------------------------------
yp = PP_PreProcess(y, '', [], [], [], randomSeed);
% Returns a structure, yp, with a range of time series in it, each a different
% transformation of the original, y.
%% ____________________FIT MODEL TO A CURATED SUBSET:_______________________ %%

% 'nothing' must stay first (it's the ratio denominator below); the rest are
% the 5 non-redundant representatives (see header note). 'peaks_08' was
% dropped 2026-08-15: r=0.985 with 'd2' on real EEG data, and 'd2' (simple
% 2nd-order differencing) is the more standard/interpretable of the pair.
fields = {'nothing', 'd1', 'd2', 'p1_20', 'p2_5', 'rmgd'};
numFields = length(fields);
% statstore = struct('fpes',{});

for i = 1:numFields
	% for each preprocessing, fit the model
	data = yp.(fields{i});
	% data is the current preprocessed data

	switch model % SO MANY OPTIONS! ;-)
		case 'ar'
			%% Check that a System Identification Toolbox license is available
			BF_CheckToolbox('identification_toolbox');

			data = zscore(data); % zscore the data from this preprocessing
			m = ar(data, order); % fit the model

			% Get statistics on fit
			%     () FPE
			statstore.fpe(i) = m.EstimationInfo.FPE;
			%     () in-sample prediction error
			e = pe(m, data);
			statstore.rmserr(i) = sqrt(mean(e.^2));
			statstore.mabserr(i) = mean(abs(e));
			statstore.ac1(i) = CO_AutoCorr(e, 1, 'Fourier');

		otherwise
			error('Unknown model ''%s''', model);
	end
end

% ------------------------------------------------------------------------------
%% Return statistics on statistics
% ------------------------------------------------------------------------------
% actually often as you make more stationary and remove trends it becomes
% harder to predict because these trends are very easy to predict, and
% making the series whiter will obviously decrease its predictability.

% (1) ratio of fpe of preprocessed to unprocessed time series
% I think just this is ok.
% for i=2:numFields
%     eval(['out.fperat_' fields{i} ' = ' num2str(statstore.fpe(i)/statstore.fpe(1)) ';']);
% end

% No, I'll just do in-sample rms error, for a single model no point fpeing
for i = 2:numFields
	out.(sprintf('stderat_%s', fields{i})) = statstore.rmserr(i) / statstore.rmserr(1);
end
% In fact, greater error in this case means a better detrending in some
% sense -- it's remobed more of the 'obvious' linear structure (assuming
% that's the aim).

% could also return statistics on other things like prediction error, but
% not alot of point, I think.

end
