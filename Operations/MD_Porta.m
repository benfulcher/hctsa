function out = MD_Porta(x,numLevels)
% MD_Porta   Porta's symbolic-dynamics word-type indices.
%
% Quantizes the time series into a small number of levels (numLevels equal-width
% bins spanning its range, with explicit edges from BF_HistEdges) and classifies consecutive length-3
% "words" of symbols by their pattern of variation:
%   0V   -- no variation (all three symbols equal)
%   1V   -- one variation (exactly one of the two transitions is flat)
%   2LV  -- two like variations (both transitions move the same direction)
%   2UV  -- two unlike variations (transitions move in opposite directions)
%
% Originally developed for heart-rate-variability analysis, quantifying the
% complexity/regularity of the symbolic dynamics of RR interval sequences.
%
% ---INPUTS:
% x, the input time series
% numLevels, the number of quantization levels (default: 6, as in the original
%    papers)
%
% ---OUTPUTS:
% A structure with fields:
% pV0, the percentage of length-3 words with no variation
% pV1, the percentage of length-3 words with one variation
% pV2LV, the percentage of length-3 words with two like variations
% pV2UV, the percentage of length-3 words with two unlike variations
% All four are NaN for a constant series.
%
% ---REFERENCES:
% A. Porta, S. Guzzetti, N. Montano, R. Furlan, M. Pagani, A. Malliani and S. Cerutti,
% "Entropy, entropy rate, and pattern classification as tools to typify complexity
% in short heart period variability series", IEEE Trans. Biomed. Eng. 48(11),
% 1282-1291 (2001). DOI: 10.1109/10.959324
% (defines the 0V, 1V, 2LV and 2UV classes of three-beat patterns.)
%
% A. Porta, E. Tobaldini, S. Guzzetti, R. Furlan, N. Montano and T. Gnecchi-Ruscone,
% "Assessment of cardiac autonomic modulation during graded head-up tilt by
% symbolic analysis of heart rate variability", Am. J. Physiol. Heart Circ.
% Physiol. 293(1), H702-H708 (2007). DOI: 10.1152/ajpheart.00006.2007
%
% ---NOTES:
% Earlier versions of this docstring cited "A. Porta et al., Quantifying the
% strength of the linear and nonlinear relationships between heart period and
% arterial pressure, IEEE Trans. Biomed. Eng. 45(8) 1017 (1998)". No paper with
% that title was found; IEEE Trans. Biomed. Eng. 45(8), 1017-1023 (1998) is a
% different paper (Cammarota and Onaral, DOI: 10.1109/10.704870). The references
% above are the papers that define and use the word classes implemented here. The
% default of six quantization levels is the value reported for this method in
% secondary sources; it was not checked in the two papers above.

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

sym = discretize(x, BF_HistEdges(x, numLevels)); % quantize into 1:numLevels equal-width levels spanning the data

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
