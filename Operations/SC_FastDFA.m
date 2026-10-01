function out = SC_FastDFA(y)
% SC_FastDFA   Scaling exponent of the time series, from fast detrended fluctuation analysis.
%
% Detrended fluctuation analysis (DFA) measures how the size of the fluctuations
% of a series grows with the timescale. The series is integrated (cumulative
% sum) and cut into windows of s samples; a straight line is fitted to the
% integrated series in each window and subtracted; and the fluctuation F(s) is
% the root-mean-square of what is left, over the whole series. For a
% self-similar series F(s) ~ s^alpha, and the output is the exponent alpha: the
% slope of log F(s) against log s. White noise gives alpha of about 0.5,
% persistent (long-memory) series give alpha > 0.5 (about 1 for 1/f noise), and
% a random walk gives about 1.5.
%
% This is a wrapper for Max Little's ML_fastdfa (Toolboxes/Max_Little/fastdfa),
% which chooses the window sizes itself: the series length divided by 1, 2, 4, 8,
% ..., from the whole series down to windows of about 3 to 5 samples. The slope
% is a least-squares line through log10(F) against log10(s) over all of these
% scales, weighted equally.
%
% ---INPUTS:
% y, the input time series, is fed straight into ML_fastdfa as a column vector.
%
% ---OUTPUTS:
% a scalar: the DFA scaling exponent, alpha.
%
% ---NOTES:
% The original fastdfa code is by Max A. Little, publicly available at
% http://www.maxlittle.net/software/index.php (see the header of ML_fastdfa.m for
% how to cite it).

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

if size(y, 2) > size(y, 1);
	y = y'; % Ensure input time series is a column vector
end

out = ML_fastdfa(y);

end
