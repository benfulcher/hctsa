function out = SY_StdNthDer(y,n)
% SY_StdNthDer  Standard deviation of the nth derivative of the time series.
%
% Based on an idea by Vladimir Vassilevsky, a DSP and Mixed Signal Design
% Consultant in a Matlab forum, who stated that You can measure the standard
% deviation of the nth derivative, if you like".
%
% cf. http://www.mathworks.de/matlabcentral/newsreader/view_thread/136539
%
% The derivative is estimated very simply by simply taking successive increments
% of the time series; the process is repeated to obtain higher order
% derivatives.
%
% Note that this idea is popular in the heart-rate variability literature, cf.
% cf. "Do Existing Measures ... ", Brennan et. al. (2001), IEEE Trans Biomed Eng 48(11)
% (and function MD_hrv_classic)
%
%---INPUTS:
%
% y, time series to analyze
%
% n, the order of derivative to analyze

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

if nargin < 2 || isempty(n)
    n = 2;
end

yd = diff(y,n); % crude method of taking a derivative that could be improved
                % upon in future
out = std(yd);

end
