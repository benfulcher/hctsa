function out = MD_rawHRVmeas(x)
% MD_rawHRVmeas     Heart rate variability (HRV) measures of a time series.
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

N = length(x); % time-series length

% Triangular histogram index
out.tri10 = N/max(histcounts(x,10));
out.tri20 = N/max(histcounts(x,20));
out.trisqrt = N/max(histcounts(x,'BinMethod','sqrt'));

% Poincare plot measures: see
% "Do Existing Measures ... ", Brennan et. al. (2001), IEEE Trans Biomed Eng 48(11)
diffx = diff(x);
out.SD1 = 1/sqrt(2) * std(diffx) * 1000;
out.SD2 = sqrt(2 * var(x) - (1/2) * std(diffx)^2) * 1000;

% NOTE: the SD1 measure is the same (up to a linear transformation)
%       as the implementation in DN_Spread(diff(x),'std');

end
