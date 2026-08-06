function out = WL_fBM(y)
% WL_fBM   Parameters of fractional Gaussian noise/Brownian motion in a time series
%
% Uses the wfbmesti function from Matlab's Wavelet Toolbox
%
% ---INPUT:
% y, the time series to analyze.
%
% ---OUTPUTS: three Hurst-exponent estimates from wfbmesti's three internal
% estimators: a second-order-derivative estimate (H_deriv2), a wavelet-based
% version of the same (H_deriv2Wavelet, using a sym5 filter), and a
% wavelet-variance-vs-level regression estimate (H_varLevel).

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
%% Check that a Wavelet Toolbox license is available:
% ------------------------------------------------------------------------------
BF_CheckToolbox('wavelet_toolbox');

% Parameter estimation of fractional Brownian motion
hest = wfbmesti(y);
out.H_deriv2 = hest(1); % second-order discrete-derivative estimate
out.H_deriv2Wavelet = hest(2); % second-order discrete derivative, wavelet (sym5) version
out.H_varLevel = hest(3); % wavelet-variance-vs-level regression estimate

end
