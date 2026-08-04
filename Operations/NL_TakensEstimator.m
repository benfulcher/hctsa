function out = NL_TakensEstimator(y, Nref, rad, past, embedParams, randomSeed)
% NL_TakensEstimator   Taken's estimator for correlation dimension.
%
% cf. "Detecting strange attractors in turbulence", F. Takens.
% Lect. Notes Math. 898 p366 (1981)
%
% Uses the TISEAN routines d2 and c2t (Takens' estimator from correlation sum
% data) rather than TSTOOL's takens_estimator, which this operation used
% previously. d2 estimates the correlation sum across embedding dimensions
% 1:m (only the dimension-m block is used here); c2t then reads off Takens'
% maximum-likelihood dimension estimate as a function of length scale from
% that correlation-sum data, and the value at an upper length scale of rad
% standard deviations of y is taken as the output (matching Kantz &
% Schreiber's recommendation of using half a standard deviation, already
% used the same way in NL_TISEAN_d2.m's takens05 -- TSTOOL's own rad
% parameter was instead defined as a proportion of "attractor size", so this
% is an equivalent-in-spirit but not numerically identical length-scale
% convention).
%
% ---INPUTS:
% y, the input time series
% Nref, the number of reference points (can be -1 to use all points)
% rad, the upper length scale to read off the dimension estimate, in standard
%       deviations of y (cf. TSTOOL's rad, a proportion of attractor size)
% past, the Theiler window
% embedParams, the embedding parameters for BF_Embed, in the form {tau,m}
% randomSeed, whether (and how) to reset the random seed, using BF_ResetSeed
%               (relevant if an embedding-dimension method requiring
%               randomization is used)
%
% ---OUTPUT: the Taken's estimator of the correlation dimension, d2.

% ------------------------------------------------------------------------------
% Copyright (C) 2020, Ben D. Fulcher <ben.d.fulcher@gmail.com>,
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
%% Check inputs
% ------------------------------------------------------------------------------
N = length(y); % time-series length

% 1) Nref
if nargin < 2 || isempty(Nref)
	Nref = -1; % use all points
end

% 2) Upper length scale (standard deviations of y) at which to read off the
% dimension estimate from the correlation-sum data:
if nargin < 3 || isempty(rad)
	rad = 0.05;
end

% 3) Theiler window
if nargin < 4 || isempty(past)
	past = 1; % just exclude current point
end
if (past > 0) && (past < 1)
	past = floor(N * past); % specify a fraction of the time series length...
end

% 4) Embedding parameters
if nargin < 5 || isempty(embedParams)
	embedParams = {'ac', 'fnn'};
	fprintf(1, 'Using default time-delay embedding using autocorrelation and fnn-mar\n');
else
	if length(embedParams) ~= 2
		error('Embedding parameters are incorrectly formatted, we need {tau,m}')
	end
end

% 5) randomSeed: how to treat the randomization
if nargin < 6
	randomSeed = []; % default
end

% ------------------------------------------------------------------------------
%% Resolve the embedding parameters (tau, m)
% ------------------------------------------------------------------------------
tm = BF_Embed(y, embedParams{1}, embedParams{2}, true, randomSeed);
tau = tm(1);
m = tm(2);

% ------------------------------------------------------------------------------
%% Write the file for TISEAN to work with
% ------------------------------------------------------------------------------
filePath = BF_WriteTempFile(y);

% Map TSTOOL's Nref convention (-1 = use all points) onto TISEAN's d2
% (-N 0 = use all pairs):
if Nref == -1
	NrefTISEAN = 0;
else
	NrefTISEAN = Nref;
end

% ------------------------------------------------------------------------------
%% Run the TISEAN codes, d2 then c2t
% ------------------------------------------------------------------------------
% d2 over embedding dimensions 1:m (only the m-th is used below). Note:
% "-M<m>,<m>" (i.e. asking for a single, fixed embedding dimension) triggers
% a bug in this TISEAN build where it reports "0 lines read" and produces no
% output at all; "-M1,<m>" (a genuine range, as NL_TISEAN_d2.m already uses)
% works correctly, so that's used here too and the dimension-m block is
% picked out afterwards. The swept range of length scales is left at d2's
% default (spans the full data interval, so comfortably covers the
% rad-standard-deviations cutoff used below):
[~, res] = system(sprintf('d2 -d%u -M1,%u -t%u -N%u %s', tau, m, past, NrefTISEAN, filePath));
if exist([filePath '.stat'], 'file'), delete([filePath '.stat']); end
if exist([filePath '.d2'], 'file'), delete([filePath '.d2']); end
if exist([filePath '.h2'], 'file'), delete([filePath '.h2']); end

if isempty(res) || ~isempty(regexp(res, 'command not found', 'once'))
	if exist([filePath '.c2'], 'file'), delete([filePath '.c2']); end
	error('Call to TISEAN function ''d2'' failed.');
end

[~, res] = system(sprintf('c2t %s.c2', filePath));
delete([filePath '.c2']);

if isempty(res) || ~isempty(regexp(res, 'command not found', 'once'))
	error('Call to TISEAN function ''c2t'' failed.');
end

% ------------------------------------------------------------------------------
%% Parse the c2t output and read off the estimate at rad standard deviations
% ------------------------------------------------------------------------------
s = textscan(res, '%[^\n]'); s = s{1};
wi = strmatch('writing to stdout', s);
if isempty(wi)
	error('TISEAN routine ''c2t'' returned unexpected output.');
end
s = s(wi + 1:end);

% There should be one '#m=' block per embedding dimension 1:m; take the
% last one (dimension m, the one actually requested):
w = strmatch('#m=', s);
if length(w) ~= m
	error('TISEAN routine ''c2t'' returned an unexpected number of data blocks.');
end
w(end + 1) = length(s) + 1;
ss = s(w(m) + 1:w(m + 1) - 1);

rc = zeros(length(ss), 2); % [length scale r, Takens estimate]
nn = 0;
for jj = 1:length(ss)
	tmp = textscan(ss{jj}, '%f%f');
	if all(cellfun(@isempty, tmp))
		break % a trailing comment line
	end
	nn = nn + 1;
	rc(nn, :) = horzcat(tmp{:});
end
rc = rc(1:nn, :);

eup = rad * std(y); % upper length scale, in standard deviations of y
theIndex = find(rc(:, 1) > eup, 1, 'first');
if isempty(theIndex)
	out = NaN;
else
	out = rc(theIndex, 2);
end

end
