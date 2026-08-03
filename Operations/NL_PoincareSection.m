function out = NL_PoincareSection(y, ref, embedParams)
% NL_PoincareSection   Poincare section analysis of a time series.
%
% Time-delay embeds the time series and computes a Poincare section using
% TISEAN's 'poincare', which cuts the trajectory on a fixed embedding
% coordinate (the last, by convention) held at its own mean, in a single
% crossing direction.
%
% ---INPUTS:
% y, the input time series
%
% ref: 'max' or 'min', selecting the crossing direction. This operation
%      previously used TSTOOL's 'poincare', which cut a hyperplane
%      orthogonal to the local tangent vector at a chosen reference point in
%      the time series -- a construction TISEAN has no equivalent for.
%      'ref' is repurposed here to instead pick which of the two possible
%      crossing directions to use: 'max' selects crossings heading toward a
%      local maximum (ascending through the mean, TISEAN's "from below",
%      -C0); 'min' selects crossings heading toward a local minimum
%      (descending through the mean, "from above", -C1).
%
% embedParams: the usual thing to give BF_Embed for the time-delay embedding, as
%               {tau,m}. m is forced to 3 -- i.e., embed in a 3
%               dimensional space so that the Poincare section is 2-dimensional.
%
% ---OUTPUTS: include statistics on the x- and y- components of these vectors on the
% Poincare surface, on distances between adjacent points, distances from the
% mean position, and the entropy of the vector cloud.
%
% Uses TISEAN's 'poincare' (this operation previously used TSTOOL's
% 'poincare', a different construction -- see 'ref' above).
% TISEAN: https://www.pks.mpg.de/tisean/

% Another thing that could be cool to do is to analyze variation in the plots as
% ref changes... (not done here)
%
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

if nargin < 2 || isempty(ref)
	ref = 'max';
end

if nargin < 3 || isempty(embedParams)
	embedParams = {'mi', 3};
	fprintf(1, 'Using default embedding settings: minimum of the automutual information for tau and m = 3\n');
end

if embedParams{2} ~= 3
	embedParams{2} = 3;
	fprintf(1, 'Three-dimensional embedding\n');
end

% TISEAN's poincare has no reference-point concept (see header comment);
% 'ref' is repurposed to select the crossing direction:
switch ref
	case 'max'
		direction = 0; % crossing from below (heading toward a local maximum)
	case 'min'
		direction = 1; % crossing from above (heading toward a local minimum)
	otherwise
		error(['NL_PoincareSection: ref must be ''max'' or ''min'' ' ...
			   '(TISEAN''s poincare has no reference-point concept, only ' ...
			   'a choice of crossing direction).']);
end

doPlot = 0; % plot outputs to a figure

% ------------------------------------------------------------------------------
%% Run TISEAN's poincare
% ------------------------------------------------------------------------------
N = length(y); % length of the time series
tm = BF_Embed(y, embedParams{1}, embedParams{2}, 2);
tau = tm(1);
m = tm(2);

filePath = BF_WriteTempFile(y);
outFilePath = [filePath '.poin'];

% Cut on the last (m-th) embedding coordinate, at TISEAN's own default
% threshold (that coordinate's mean):
[~, res] = system(sprintf('poincare -d%u -m%u -q%u -C%u -o %s %s', ...
						  tau, m, m, direction, outFilePath, filePath));
delete(filePath); % remove the temporary time-series data file

if isempty(res) || ~isempty(regexp(res, 'command not found', 'once'))
	if exist(outFilePath, 'file'), delete(outFilePath); end
	error('Call to TISEAN function ''poincare'' failed.');
end

if ~exist(outFilePath, 'file')
	error('TISEAN function ''poincare'' did not produce a .poin output file.');
end

% An empty output file means no crossings were found in this direction
% (dlmread errors on an empty file, so check first):
fileInfo = dir(outFilePath);
if fileInfo.bytes == 0
	delete(outFilePath);
	fprintf(1, 'No section points found to run NL_PoincareSection\n');
	out = NaN; return
end

% Columns are the two uncut embedding coordinates, followed by the
% (interpolated) crossing time -- only the first two are point coordinates:
v = dlmread(outFilePath);
delete(outFilePath);
v = v(:, 1:2);
NN = size(v, 1);

if NN < 2
	fprintf(1, 'No section points found to run NL_PoincareSection\n');
	out = NaN; return
end

% Labeling poincare surface plane x-y
x = v(:, 1);
y = v(:, 2);

if doPlot
	figure('color', 'w'); box('on');
	plot(x, y, '.k'); axis equal
end

% ------------------------------------------------------------------------------
%% Get statistics out
% ------------------------------------------------------------------------------

% Basic statistics:
out.pcross = NN / N; % proportion of time series that crosses poincare surface

out.maxx = max(x);
out.minx = min(x);
out.stdx = std(x);
out.iqrx = iqr(x);
out.meanx = mean(x);
out.ac1x = CO_AutoCorr(x, 1, 'Fourier');
out.ac2x = CO_AutoCorr(x, 2, 'Fourier');
out.tauacx = CO_FirstCrossing(x, 'ac', 0, 'continuous');

out.maxy = max(y);
out.miny = min(y);
out.stdy = std(y);
out.iqry = iqr(y);
out.meany = mean(y);
out.ac1y = CO_AutoCorr(y, 1, 'Fourier');
out.ac2y = CO_AutoCorr(y, 2, 'Fourier');
out.tauacy = CO_FirstCrossing(y, 'ac', 0, 'continuous');

out.boxarea = range(x) * range(y);

% Statistics on distance between adjacent points, ds
vdiff = v(2:end, :) - v(1:end - 1, :);
ds = sqrt(vdiff(:, 1).^2 + vdiff(:, 2).^2);

% Probability that next point in series is within radius r of current point in
% the poincare section:
out.pwithinr01 = sum(ds < 0.1) / (NN - 1);
out.pwithin02 = sum(ds < 0.2) / (NN - 1);
out.pwithin03 = sum(ds < 0.3) / (NN - 1);
out.pwithin05 = sum(ds < 0.5) / (NN - 1);
out.pwithin1 = sum(ds < 1) / (NN - 1);
out.pwithin2 = sum(ds < 2) / (NN - 1);
out.meands = mean(ds);
out.maxds = max(ds);
out.minds = min(ds);
out.iqrds = iqr(ds);

% Now normalize both axes and look for structure in the cloud of points
% don't normalize for standard deviation -- this probably reveals some
% structure...? But location is already noted.
x = x - mean(x);
y = y - mean(y);

% Statistics on distance on Poincare surface from (mean,mean)
D = sqrt(x.^2 + y.^2); % distance from (mean,mean)
out.maxD = max(D);
out.minD = min(D);
out.stdD = std(D);
out.iqrD = iqr(D);
out.meanD = mean(D);
out.ac1D = CO_AutoCorr(D, 1, 'Fourier');
out.ac2D = CO_AutoCorr(D, 2, 'Fourier');
out.tauacD = CO_FirstCrossing(D, 'ac', 0, 'continuous');

% ------------------------------------------------------------------------------
%% Statistics of the boxed distribution:
% ------------------------------------------------------------------------------

numPartitions = [5, 10];
% (i) 5 partitions per axis
% (ii) 10 partitions per axis

for i = 1:length(numPartitions)
	boxcounts = subcountboxes(x, y, numPartitions(i));
	pbox = boxcounts / NN;

	out.(sprintf('maxpbox%u', numPartitions(i))) = max(pbox(:));
	out.(sprintf('minpbox%u', numPartitions(i))) = min(pbox(:));
	out.(sprintf('zerospbox%u', numPartitions(i))) = sum(pbox(:) == 0);
	out.(sprintf('meanpbox%u', numPartitions(i))) = mean(pbox(:));
	out.(sprintf('rangepbox%u', numPartitions(i))) = range(pbox(:));
	% This probably needs to be normalized:
	out.(sprintf('hboxcounts%u', numPartitions(i))) = -sum(pbox(pbox > 0) .* log(pbox(pbox > 0)));
	out.(sprintf('tracepbox%u', numPartitions(i))) = sum(diag(pbox)); % trace
end

% NOTE: max and range provide almost the same information on real data.

% can imagine doing many more things; like seeing different slabs of the
% space, or finding the line whos vicinity includes many points, etc. but I
% think this is enough for now.

% ------------------------------------------------------------------------------
function boxcounts = subcountboxes(x, y, nbox)
	boxcounts = zeros(nbox);
	% boxes are quantiles along each axis
	xbox = quantile(x, linspace(0, 1, nbox + 1));
	ybox = quantile(y, linspace(0, 1, nbox + 1));
	xbox(end) = xbox(end) + 1;
	ybox(end) = ybox(end) + 1;

	for ii = 1:nbox % x
		rx = (x >= xbox(ii) & x < xbox(ii + 1)); % these x are in range
		for jj = 1:nbox % y
			% only need to look at those ys for which the xs are in range
			boxcounts(ii, jj) = sum(y(rx) >= ybox(jj) & y(rx) < ybox(jj + 1));
		end
	end
end
% ------------------------------------------------------------------------------

end
