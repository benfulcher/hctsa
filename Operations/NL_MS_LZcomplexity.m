% NL_MS_LZcomplexity
% 
% Calculates the Lempel-Ziv complexity of a n-bit encoding of the time
% series using Michael Small's code: 'complexity' (renamed MS_complexity here).
% 
% cf. M. Small, Applied Nonlinear Time Series Analysis: Applications in Physics,
% Physiology, and Finance (book) World Scientific, Nonlinear Science Series A,
% Vol. 52 (2005)
% Code is available at http://small.eie.polyu.edu.hk/matlab/
% 
% The code is a wrapper for Michael Small's original code and uses the
% associated mex file compiled from complexitybs.c (renamed MS_complexitybs.c
% here).
% 
% INPUTS:
% y, the input time series
% n, the (integer) number of bits to encode the data into
% preproc [opt], first apply a given preprocessing to the time series. For now,
%               just 'diff' is implemented, which zscores incremental
%               differences and then applies the complexity method.
% 
% The function has a single output: the normalized Lempel-Ziv complexity: i.e.,
% the number of distinct symbol sequences in the time series divided by the
% expected number of distinct symbols for a noise sequence.
% 
% ------------------------------------------------------------------------------
% Copyright (C) 2013,  Ben D. Fulcher <ben.d.fulcher@gmail.com>,
% <http://www.benfulcher.com>
%
% If you use this code for your research, please cite:
% B. D. Fulcher, M. A. Little, N. S. Jones., "Highly comparative time-series
% analysis: the empirical structure of time series and their methods",
% J. Roy. Soc. Interface 10(83) 20130048 (2010). DOI: 10.1098/rsif.2013.0048
%
% This function is free software: you can redistribute it and/or modify it under
% the terms of the GNU General Public License as published by the Free Software
% Foundation, either version 3 of the License, or (at your option) any later
% version.
% 
% This program is distributed in the hope that it will be useful, but WITHOUT
% ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
% FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
% details.
% 
% You should have received a copy of the GNU General Public License along with
% this program.  If not, see <http://www.gnu.org/licenses/>.
% ------------------------------------------------------------------------------

function out = NL_MS_LZcomplexity(y,n,preproc)
% Ben Fulcher, 19/2/2010

if nargin < 2 || isempty(n)
    n = 2; % n-bit encoding
end
if nargin < 3
    preproc = []; % no preprocessing
end

% Apply some pre-processing to the time series before performing the analysis
if ischar(preproc)
    switch preproc
    case 'diff'
        y = BF_zscore(diff(y));
    otherwise
        error('Unknown preprocessing setting ''%s''', preproc);
    end
end

% Run Michael Small's (mexed) code for calcaulting the Lempel-Ziv complexity:
out = MS_complexity(y,n);

end