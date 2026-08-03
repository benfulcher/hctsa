function out = NL_BoxCorrDim(y, numBins, embedParams)
% NL_BoxCorrDim  Correlation dimension of a time series.
%
% Estimates the correlation dimension of a time-delay embedded time series
% using a box-counting approach, via TISEAN's 'boxcount' (this operation
% previously used TSTOOL's 'corrdim').
%
% 'boxcount' estimates the Renyi entropy of order Q (Q = 2.0 here, giving
% the box-counting correlation entropy/dimension) via partitioning, for
% every embedding dimension 1:m and a sweep of length scales, writing (for
% each embedding dimension d) both the raw entropy H_Q(epsilon,d) and its
% increment over the (d-1)-dimensional embedding, H_Q(epsilon,d) -
% H_Q(epsilon,d-1) -- it is this increment (which approaches the
% correlation dimension itself as d grows) that plays the role of TSTOOL's
% corrdim matrix, whose columns/rows this operation's output statistics
% below already summarise as per-embedding-dimension and per-length-scale
% local dimension estimates.
%
% ---INPUTS:
% y, column vector of time series data
% numBins, number of length-scale (epsilon) values in the box-counting sweep
%          (TSTOOL's own "maximum number of partitions per axis" doesn't
%          have an exact TISEAN equivalent; this is the closest analogue --
%          it controls the resolution of the length-scale sweep the same
%          way numBins previously did).
% embedParams [opt], embedding parameters as {tau,m} in 2-entry cell for a
%                   time-delay, tau, and embedding dimension, m. As inputs to BF_Embed.
%
% ---OUTPUTS: Simple summaries of the outputs from boxcount.

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
	embedParams = {'ac', 'fnnmar'};
else
	if length(embedParams) ~= 2
		error('Embedding parameters should be formatted like {tau,m}')
	end
end

% ------------------------------------------------------------------------------
%% Resolve the embedding parameters (tau, m)
% ------------------------------------------------------------------------------
tm = BF_Embed(y, embedParams{1}, embedParams{2}, 2);
tau = tm(1);
mMax = tm(2);

% ------------------------------------------------------------------------------
%% Run the TISEAN code, boxcount
% ------------------------------------------------------------------------------
filePath = BF_WriteTempFile(y);
outFilePath = [filePath '.box'];

[~, res] = system(sprintf('boxcount -M1,%u -d%u -Q2.0 -#%u -o %s %s', ...
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
	error('TISEAN function ''boxcount'' returned an unexpected number of data blocks.');
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

for i = 2:m
	out.(sprintf('meand%u', i)) = mean(rs(:, i));
	out.(sprintf('mediand%u', i)) = median(rs(:, i));
	out.(sprintf('mind%u', i)) = min(rs(:, i));
end

for i = 2:ldr
	out.(sprintf('meanr%u', i)) = mean(rs(i, :));
	out.(sprintf('medianr%u', i)) = median(rs(i, :));
	out.(sprintf('minr%u', i)) = min(rs(i, :));
	out.(sprintf('meanchr%u', i)) = mean(diff(rs(i, :)));
end

out.stdmean = std(mean(rs));
out.stdmedian = std(median(rs));

rsstretch = rs(:);
out.medianstretch = median(rsstretch);
out.minstretch = min(rsstretch); % same as at maximum embedding dimension, m, or usually at maximum ldr (18)
out.iqrstretch = iqr(rsstretch);

end
