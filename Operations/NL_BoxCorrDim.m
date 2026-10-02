function out = NL_BoxCorrDim(y, numBins, embedParams)
% NL_BoxCorrDim   How the box-counting entropy of a delay embedding grows with embedding dimension.
%
% Time-delay embeds the series in d = 1,...,m dimensions and partitions the
% space into boxes of side epsilon, using TISEAN's 'boxcount' (this operation
% previously used TSTOOL's 'corrdim'). With p_i the fraction of embedded points
% in box i, boxcount gives the order-2 Renyi (collision) entropy
% H(epsilon,d) = -log(sum_i p_i^2) for a sweep of numBins box sizes, from the
% full range of the series downward, and the increment over the
% (d-1)-dimensional embedding, I(epsilon,d) = H(epsilon,d) - H(epsilon,d-1)
% (defined for d = 2,...,m; at d = 1, boxcount reports H itself, which is not an
% increment, so d = 1 is excluded from all summaries). The matrix I (length scales
% by embedding dimensions 2,...,m) is summarized by its mean, median and minimum
% across length scales for each dimension, across dimensions for each length
% scale, and overall.
%
% ---INPUTS:
% y, column vector of time series data
% numBins, number of length-scale (epsilon) values in the box-counting sweep
%          (default: 100). TSTOOL's "maximum number of partitions per axis" has no
%          exact TISEAN equivalent; this is the closest analogue.
% embedParams [opt], embedding parameters as {tau,m} in a 2-entry cell, a
%          time delay, tau, and embedding dimension, m, as inputs to BF_Embed
%          (default: {'ac','fnn'})
%
% ---OUTPUTS: a structure of summaries of the matrix I(epsilon,d), with d the
% embedding dimension and r the index of the length scale (r = 1 is the full
% range of the series, larger r are finer scales):
% meand<d>, mediand<d>, mind<d>: mean, median and minimum of I over length
%          scales, at embedding dimension d = 2,...,m
% meanr<r>, medianr<r>, minr<r>: mean, median and minimum of I over embedding
%          dimensions 2,...,m, at length scale r = 2,...,numBins
% meanchr<r>: mean change of I from one embedding dimension to the next
%          (d = 2,...,m), at length scale r = 2,...,numBins
% stdmean, stdmedian: standard deviation, across embedding dimensions 2,...,m, of
%          the mean (or median) of I over length scales
% medianstretch, minstretch, iqrstretch: median, minimum and interquartile
%          range of I over all length scales and embedding dimensions 2,...,m
%
% ---NOTES:
% The increment I approaches the entropy rate of the process (the K2 entropy,
% per delay step) rather than a slope against log(epsilon), so these features
% are entropy-rate-like even though the function is named for the correlation
% dimension.

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

doPlot = false; % plot outputs to a figure

% ------------------------------------------------------------------------------
%% Check inputs, preliminaries
% ------------------------------------------------------------------------------
% (1) Maxmum number of partitions per axis, numBins
if nargin < 2 || isempty(numBins)
	numBins = 100; % default number of bins per axis is 100
end

% (2) Set embedding parameters to defaults
if nargin < 3 || isempty(embedParams)
	embedParams = {'ac', 'fnn'};
else
	if length(embedParams) ~= 2
		error('Embedding parameters should be formatted like {tau,m}')
	end
end

% ------------------------------------------------------------------------------
%% Resolve the embedding parameters (tau, m)
% ------------------------------------------------------------------------------
tm = BF_Embed(y, embedParams{1}, embedParams{2}, true);
tau = tm(1);
if isnan(tau)
	warning('Could not determine embedding parameters for this time series');
	out = NaN; return
end
mMax = tm(2);

% ------------------------------------------------------------------------------
%% Run the TISEAN code, boxcount
% ------------------------------------------------------------------------------
filePath = BF_WriteTempFile(y);
outFilePath = [filePath '.box'];

[~, res] = BF_TiseanSystem(sprintf('boxcount -M1,%u -d%u -Q2.0 -#%u -o %s %s', ...
						  mMax, tau, numBins, outFilePath, filePath));

if isempty(res) || ~isempty(regexp(res, 'command not found', 'once'))
	if exist(outFilePath, 'file'), delete(outFilePath); end
	error('Call to TISEAN function ''boxcount'' failed.');
end

if ~exist(outFilePath, 'file')
	error('TISEAN function ''boxcount'' did not produce a .box output file.');
end

fid = fopen(outFilePath);
fileLines = textscan(fid, '%[^\n]');
fclose(fid);
delete(outFilePath);
fileLines = fileLines{1};

% One '#component = 1 embedding = <d>' block per embedding dimension d = 1:mMax:
w = strmatch('#component', fileLines);
if length(w) ~= mMax
	% Data-dependent: TISEAN couldn't produce output at every requested embedding
	% dimension for this series (e.g. too short relative to mMax/tau).
	warning('TISEAN function ''boxcount'' returned an unexpected number of data blocks.');
	out = NaN; return
end
w(end + 1) = length(fileLines) + 1;

rs = zeros(numBins, mMax); % local dimension estimate at each (length scale, embedding dim)
for d = 1:mMax
	ss = fileLines(w(d) + 1:w(d + 1) - 1);
	nn = 0;
	for jj = 1:length(ss)
		tmp = textscan(ss{jj}, '%f%f%f');
		if all(cellfun(@isempty, tmp))
			break % a trailing blank/comment line
		end
		nn = nn + 1;
		rs(nn, d) = tmp{3}; % the increment over the (d-1)-dim embedding (see header comment)
	end
	if nn ~= numBins
		error('TISEAN function ''boxcount'' returned an unexpected number of length scales.');
	end
end

% Contains ldr as rows for embedding dimensions 1:m as columns;
if doPlot
	figure('color', 'w'); box('on');
	plot(rs, 'k');
end

% ------------------------------------------------------------------------------
%% Output Statistics
% ------------------------------------------------------------------------------
% These statistics are just from intuition

m = size(rs, 2); % number of embedding dimensions (= mMax)
ldr = size(rs, 1); % number of length scales (= numBins)

if m < 2
	% The increment I is only defined from d = 2 (d = 1 holds H itself)
	warning('Embedding dimension m = %u is too low for a box-counting entropy increment', m);
	out = NaN; return
end

for i = 2:m
	out.(sprintf('meand%u', i)) = mean(rs(:, i));
	out.(sprintf('mediand%u', i)) = median(rs(:, i));
	out.(sprintf('mind%u', i)) = min(rs(:, i));
end

for i = 2:ldr
	out.(sprintf('meanr%u', i)) = mean(rs(i, 2:end));
	out.(sprintf('medianr%u', i)) = median(rs(i, 2:end));
	out.(sprintf('minr%u', i)) = min(rs(i, 2:end));
	out.(sprintf('meanchr%u', i)) = mean(diff(rs(i, 2:end)));
end

out.stdmean = std(mean(rs(:, 2:end)));
out.stdmedian = std(median(rs(:, 2:end)));

rsstretch = rs(:, 2:end);
rsstretch = rsstretch(:);
out.medianstretch = median(rsstretch);
out.minstretch = min(rsstretch); % same as at maximum embedding dimension, m, or usually at maximum ldr (18)
out.iqrstretch = iqr(rsstretch);

end
