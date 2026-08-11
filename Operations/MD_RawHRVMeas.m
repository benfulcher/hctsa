function out = MD_RawHRVMeas(x)
% MD_RawHRVMeas     Poincare plot measures used in HRV analysis.
%                   HRV: heart rate variability.
%
% Computes the triangular histogram index and Poincare plot measures to a time
% series assumed to measure sequences of consecutive RR intervals measured in
% milliseconds. Doesn't make much sense for other time series
%
% cf. "Do existing measures of Poincare plot geometry reflect nonlinear
%      features of heart rate variability?"
%      M. Brennan, et al., IEEE T. Bio.-Med. Eng. 48(11) 1342 (2001)
%
% Note that pNNx is not done here, but in MD_pNN.m
%
% This code is heavily derived from Max Little's hrv_classic.m code
% Max Little: http://www.maxlittle.net/

% ------------------------------------------------------------------------------
% Copyright (C) 2016, Max A. Little, <max.a.little@gmail.com>,
% <http://www.maxlittle.net/> and Ben D. Fulcher <ben.d.fulcher@gmail.com>,
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

N = length(x); % time-series length

% Triangular histogram index
out.tri10 = N / max(histcounts(x, 10));
out.tri20 = N / max(histcounts(x, 20));
out.trisqrt = N / max(histcounts(x, 'BinMethod', 'sqrt'));

% 'Poincare plot measures': see
% "Do Existing Measures ... ", Brennan et. al. (2001), IEEE Trans Biomed Eng 48(11)
% https://doi.org/10.1109/10.959330
diffx = diff(x);

% SD1: Standard deviation of first derivative (differences):
%       (variability in direction perpendicular to line of identity in (x_t,x_{t+1}))
out.SD1 = 1 / sqrt(2) * std(diffx) * 1000;

% SD2: (variability in direction of line of identity in (x_t,x_{t+1}))
% Mix of standard deviation of data and standard deviation of first derivative:
% High values for signal with high variance but low variance of differences
% Low values for signal with low variance but high variance of differences
out.SD2 = sqrt(2 * var(x) - 0.5 * std(diffx)^2) * 1000;

% NOTE: the SD1 measure is the same (up to a rescaling)
%       as the implementation in DN_Spread(diff(x),'std');

end
