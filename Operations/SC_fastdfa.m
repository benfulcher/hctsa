function out = SC_fastdfa(y)
% SC_fastdfa   Matlab wrapper for Max Little's ML_fastdfa code
%
% Measures the scaling exponent of the time series using a fast implementation
% of detrended fluctuation analysis (DFA).
%
%---INPUT:
% y, the input time series, is fed straight into the fastdfa script.

% The original fastdfa code is by Max A. Little and publicly-available at
% http://www.maxlittle.net/software/index.php
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

if size(y,2) > size(y,1);
    y = y'; % Ensure input time series is a column vector
end

out = ML_fastdfa(y);

end
