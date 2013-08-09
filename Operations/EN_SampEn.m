% EN_sampen
% 
% Estimates the Sample Entropy of the time series, SampEn(m,r), by referencing
% code from PhysioNet.
% 
% The publicly-available PhysioNet code, sampenc (renamed here to RN_sampenc) is
% available from:
% http://www.physionet.org/physiotools/sampen/matlab/1.1/sampenc.m
%
% cf. "Physiological time-series analysis using approximate entropy and sample
% entropy", J. S. Richman and J. R. Moorman, Am. J. Physiol. Heart Circ.
% Physiol., 278(6) H2039 (2000)
% 
% This function can also calculate the SampEn of successive increments of time
% series, i.e., we using an incremental differencing pre-processing, as
% used in the so-called Control Entropy quantity:
% 
% "Control Entropy: A complexity measure for nonstationary signals"
% E. M. Bollt and J. Skufca, Math. Biosci. Eng., 6(1) 1 (2009)
% 
% INPUTS:
% y, the input time series
% M, the embedding dimension
% r, the threshold
% preprocess [opt], (i) 'diff1', incremental difference preprocessing.
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

function out = EN_SampEn(y,M,r,preprocess)
% Ben Fulcher, November 2009

if nargin < 4
    preprocess = ''; % don't apply preprocessing
end

if ~isempty(preprocess)
    switch preprocess
    case 'diff1'
        % First do an incremental differencing of the time series
        % thus yielding the 'Control Entropy'
        y = BF_zscore(diff(y));
    otherwise
        error('Unknown preprocessing setting ''%s''',preprocess);
    end
end

% Use the physionet code to calculate the Sample Entropy using these parameters:
[e, p, ~, ~] = PN_sampenc(y,M,r);

%% Give outputs
for i = 1:M
    eval(sprintf('out.p%u = p(%u);',i,i));
    eval(sprintf('out.sampen%u = e(%u);',i,i));
end

out.meanchsampen = mean(diff(e));
out.meanchp = mean(diff(p));

end