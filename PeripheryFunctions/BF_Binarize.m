function yBin = BF_Binarize(y,binarizeHow)
% BF_Binarize    Converts an input vector into a binarized version

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

%-------------------------------------------------------------------------------
% Check inputs, set defaults:
%-------------------------------------------------------------------------------
if nargin < 2 || isempty(binarizeHow)
    binarizeHow = 'diff';
end

%-------------------------------------------------------------------------------
% Do the binary transformation:
%-------------------------------------------------------------------------------

switch binarizeHow
    case 'diff'
        % Binary signal: 1 for stepwise increases, 0 for stepwise decreases
        yBin = ((sign(diff(y)))+1)/2;

    case 'mean'
        % Binary signal: 1 for above mean, 0 for below mean
        yBin = (sign(y - mean(y))+1)/2;

    case 'median'
        % Binary signal: 1 for above median, 0 for below median
        yBin = (sign(y - median(y))+1)/2;

    case 'iqr'
        % Binary signal: 1 if inside interquartile range, 0 otherwise
        iqr = quantile(y,[0.25, 0.75]);
        iniqr = (y > iqr(1) & y <= iqr(2));
        yBin = zeros(length(y),1);
        yBin(iniqr) = 1;

    otherwise
        error('Unknown binary transformation setting ''%s''',binarizeHow)
end

end
