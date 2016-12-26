function out = EN_MS_LZcomplexity(y,n,preProc)
% EN_MS_LZcomplexity Lempel-Ziv complexity of a n-bit encoding of a time series
%
%---INPUTS:
% y, the input time series
% n, the (integer) number of bits to encode the data into
% preProc [opt], first apply a given preProcessing to the time series. For now,
%               just 'diff' is implemented, which zscores incremental
%               differences and then applies the complexity method.
%
%---OUTPUT: the normalized Lempel-Ziv complexity: i.e., the number of distinct
%           symbol sequences in the time series divided by the expected number
%           of distinct symbols for a noise sequence.

% Uses Michael Small's code: 'complexity' (renamed MS_complexity here).
%
% cf. M. Small, Applied Nonlinear Time Series Analysis: Applications in Physics,
% Physiology, and Finance (book) World Scientific, Nonlinear Science Series A,
% Vol. 52 (2005)
% Code is available at http://small.eie.polyu.edu.hk/matlab/
%
% The code is a wrapper for Michael Small's original code and uses the
% associated mex file compiled from complexitybs.c (renamed MS_complexitybs.c
% here).
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
    n = 2; % n-bit encoding
end
if nargin < 3
    preProc = []; % no preprocessing
end

% Apply some pre-processing to the time series before performing the analysis
if ischar(preProc)
    switch preProc
    case 'diff'
        y = zscore(diff(y));
    otherwise
        error('Unknown preprocessing setting ''%s''', preProc);
    end
end

% Run Michael Small's (mexed) code for calcaulting the Lempel-Ziv complexity:
out = MS_complexity(y,n);

end
