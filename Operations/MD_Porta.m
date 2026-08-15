function out = MD_Porta(x,numLevels)
% MD_Porta      Porta's symbolic-dynamics word-type indices.
%
% Quantizes the time series into a small number of levels and classifies
% consecutive length-3 "words" of symbols by their pattern of variation:
%   0V   -- no variation (all three symbols equal)
%   1V   -- one variation (exactly one of the two transitions is flat)
%   2LV  -- two like variations (both transitions move the same direction)
%   2UV  -- two unlike variations (transitions move in opposite directions)
%
% Originally developed for heart-rate-variability analysis, quantifying the
% complexity/regularity of the symbolic dynamics of RR interval sequences.
%
% cf. A. Porta et al., "Quantifying the strength of the linear and
% nonlinear relationships between heart period and arterial pressure",
% IEEE Trans. Biomed. Eng. 45(8) 1017 (1998)
%
% ---INPUTS:
% x, the input time series
% numLevels, the number of quantization levels (default: 6, as in the
%            original papers)
%
% ---OUTPUTS:
% out.pV0, out.pV1, out.pV2LV, out.pV2UV: percentage of length-3 words of
% each type

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

if nargin < 2 || isempty(numLevels)
    numLevels = 6; % as in the original papers
end

if std(x) == 0
    % Constant series: quantization is undefined
    out.pV0 = NaN;
    out.pV1 = NaN;
    out.pV2LV = NaN;
    out.pV2UV = NaN;
    return
end

sym = discretize(x, numLevels); % quantize into 1:numLevels

d = diff(sym); % transitions between consecutive symbols
d1 = d(1:end-1); % first transition of each length-3 word
d2 = d(2:end);   % second transition of each length-3 word
numWords = length(d1);

is0V = (d1 == 0) & (d2 == 0);
is1V = xor(d1 == 0, d2 == 0);
is2V = (d1 ~= 0) & (d2 ~= 0);
is2LV = is2V & (sign(d1) == sign(d2));
is2UV = is2V & (sign(d1) ~= sign(d2));

out.pV0 = 100 * sum(is0V) / numWords;
out.pV1 = 100 * sum(is1V) / numWords;
out.pV2LV = 100 * sum(is2LV) / numWords;
out.pV2UV = 100 * sum(is2UV) / numWords;

end
