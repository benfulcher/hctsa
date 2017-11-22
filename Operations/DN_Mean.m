function out = DN_Mean(y,meanType)
% DN_Mean   A given measure of location of a data vector.
%
%---INPUTS:
%
% y, the input data vector
%
% meanType, (i) 'norm' or 'arithmetic', arithmetic mean
%           (ii) 'median', median
%           (iii) 'geom', geometric mean
%           (iv) 'harm', harmonic mean
%           (v) 'rms', root-mean-square
%           (vi) 'iqm', interquartile mean
%           (vii) 'midhinge', midhinge

% ------------------------------------------------------------------------------
% Copyright (C) 2017, Ben D. Fulcher <ben.d.fulcher@gmail.com>,
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
% Check Inputs
% ------------------------------------------------------------------------------
if nargin < 2 || isempty(meanType)
    meanType = 'arithmetic'; % normal mean
end

N = length(y); % time-series length

switch meanType
	case {'norm','arithmetic'} % mean
		out = mean(y);

    case 'median' % median
        out = median(y);

	case 'geom' % geometric mean
		out = geomean(y); %(prod(y))^(1/N);

	case 'harm' % harmonic mean
		out = harmmean(y); %N/sum(y.^(-1));

	case 'rms' % rms
		out = sqrt(sum(y.^2)/N);

    case 'iqm' % interquartile mean, cf. DN_TrimmedMean
        p = prctile(y, [25; 75]);
        out = mean(y(y >= p(1) & y <= p(2)));

    case 'midhinge' % average of 1st and third quartiles
        p = prctile(y, [25; 75]);
        out = mean(p);

    otherwise
        error('Unknown mean type ''%s''', meanType);
end

end
