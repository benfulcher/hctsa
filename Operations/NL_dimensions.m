function out = NL_dimensions(y, numBins, embedParams)
% NL_dimensions Box counting, information, and correlation dimension of a time series.
%
% Computes the box counting (D0) and correlation (D2) dimension of a
% time-delay embedded time series across a range of embedding dimensions,
% using TISEAN's 'boxcount' (Renyi entropy of order Q=0.0, giving raw
% ln(N(epsilon)), the log box-count -- the direct analogue of D0) and
% 'd2' (the classic Grassberger-Procaccia pair-counting correlation sum,
% giving ln(C(epsilon)) -- the direct analogue of D2 and, unlike
% 'boxcount' at Q=2.0, the same pair-counting construction TSTOOL's own
% correlation-dimension estimate used). This operation previously used
% TSTOOL's 'dimensions'. This function contains extensive code for
% estimating the best scaling range to estimate the dimension using a
% penalized regression procedure -- all of that downstream analysis code
% is embedding-agnostic (it only consumes (logr, logN) / (logr, logC)
% matrices) and is unchanged here.
%
% The information dimension (D1, Q=1) was already disabled in the
% previous TSTOOL-based version of this operation ("there's not extra
% information in it by these estimates") and, since that dead code
% depended entirely on TSTOOL's now-removed 'dimensions' output, has been
% removed rather than carried forward as an unreachable branch.
%
% ---INPUTS:
% y, column vector of time series data
% numBins, maximum number of partitions per axis
% embedParams, embedding parameters to feed BF_Embed() for embedding the
%              signal in the form {tau,m}
%
% ---OUTPUTS:
% A range of statistics are returned about how each dimension estimate changes
% with m, the scaling range in r, and the embedding dimension at which the best
% fit is obtained.

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
%% Preliminaries, check inputs
% ------------------------------------------------------------------------------
doPlot = false; % plot outputs to screen

% (1) Maximum number of bins, numBins
if nargin < 2 || isempty(numBins)
	numBins = 50; % 50 points
	fprintf(1, 'Using a default of 50 bins per axis.\n');
end

% (2) Set embedding parameters to defaults
if nargin < 3 || isempty(embedParams)
	embedParams = {'ac', 'fnnmar'};
	fprintf(1, 'Using default time-delay embedding parameters: autocorrelation and cao');
else
	if length(embedParams) ~= 2
		error('Embedding parameters are incorrectly formatted -- need {tau,m}')
	end
end

% ------------------------------------------------------------------------------
%% Resolve the embedding parameters (tau, m)
% ------------------------------------------------------------------------------
tm = BF_Embed(y, embedParams{1}, embedParams{2}, true);
tau = tm(1);
mopt = tm(2); % the resolved embedding dimension, possibly < 3

% dimension curves are always computed for m = 1:M (M >= 3, matching the
% original re-embed-at-3 behavior, since bc_logN(:,3) etc. are used below):
M = max(mopt, 3);

filePath = BF_WriteTempFile(y);

% ------------------------------------------------------------------------------
%% Run the TISEAN code, boxcount (Q=0.0), for the box-counting dimension (D0)
% ------------------------------------------------------------------------------
outFilePath = [filePath '.box'];
[~, res] = system(sprintf('boxcount -M1,%u -d%u -Q0.0 -#%u -o %s %s', ...
						  M, tau, numBins, outFilePath, filePath));
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

w = strmatch('#component', fileLines);
if length(w) ~= M
	error('TISEAN function ''boxcount'' returned an unexpected number of data blocks.');
end
w(end + 1) = length(fileLines) + 1;

bc_logr = [];
bc_logN = zeros(numBins, M); % ln(N(r)), the raw box count (Q=0 Renyi entropy) at each (length scale, embedding dim)
for d = 1:M
	ss = fileLines(w(d) + 1:w(d + 1) - 1);
	r = zeros(numBins, 1); logN = zeros(numBins, 1);
	nn = 0;
	for jj = 1:length(ss)
		tmp = textscan(ss{jj}, '%f%f%f');
		if all(cellfun(@isempty, tmp))
			break
		end
		nn = nn + 1;
		r(nn) = tmp{1};
		logN(nn) = tmp{2}; % raw H_0(epsilon) = ln(N(epsilon)) already
	end
	if nn ~= numBins
		error('TISEAN function ''boxcount'' returned an unexpected number of length scales.');
	end
	if d == 1
		bc_logr = log(r); % raw (linearly-spaced-ish) epsilon -> ln(epsilon)
	end
	bc_logN(:, d) = logN;
end

% ------------------------------------------------------------------------------
%% Run the TISEAN code, d2, for the correlation dimension (D2)
% ------------------------------------------------------------------------------
% "-M1,M" (a genuine range) rather than "-M<M>,<M>" works around a real bug
% in this TISEAN build where a fixed single embedding dimension reports "0
% lines read" and produces no output (see NL_GPCorrSum.m):
% Length scale, in standard deviations of y scaled by sqrt(M) (the same
% embedding-dimension correction validated in NL_GPCorrSum.m -- a
% pairwise distance in an M-dimensional embedding scales as std(y)*sqrt(M),
% not std(y) alone). Unlike NL_GPCorrSum (which only ever needs one
% target dimension and can tolerate an occasional honest NaN there), every
% embedding dimension 1:M needs usable data here, and the highest dimension
% M is by far the most likely to have zero pairs at TISEAN's own default
% minimum epsilon (data interval/1000) -- confirmed empirically: for
% structureless noise embedded at a representative M=7, TISEAN's own
% default left 25/50 bins with zero pairs (i.e. ln(C(r)) = -Inf) at that
% dimension, and even minEps = maxEps/100 or /30 still left 25 or 16/50
% zero; minEps = maxEps/10 was the first tested ratio to leave zero:
maxEps = std(y) * sqrt(M);
minEps = maxEps / 10;
[~, res] = system(sprintf('d2 -d%u -M1,%u -t0 -N0 -r%g -R%g -#%u %s', ...
						  tau, M, minEps, maxEps, numBins, filePath));
if isempty(res) || ~isempty(regexp(res, 'command not found', 'once'))
	if exist([filePath '.c2'], 'file'), delete([filePath '.c2']); end
	error('Call to TISEAN function ''d2'' failed.');
end
if exist([filePath '.stat'], 'file'), delete([filePath '.stat']); end
if exist([filePath '.d2'], 'file'), delete([filePath '.d2']); end
if exist([filePath '.h2'], 'file'), delete([filePath '.h2']); end
if ~exist([filePath '.c2'], 'file')
	error('TISEAN function ''d2'' did not produce a .c2 output file.');
end

fid = fopen([filePath '.c2']);
fileLines = textscan(fid, '%[^\n]');
fclose(fid);
delete([filePath '.c2']);
fileLines = fileLines{1};

w = strmatch('#dim=', fileLines);
if length(w) ~= M
	error('TISEAN function ''d2'' returned an unexpected number of data blocks.');
end
w(end + 1) = length(fileLines) + 1;

co_logr = [];
co_logC = zeros(numBins, M); % ln(C(r)), the raw correlation sum at each (length scale, embedding dim)
for d = 1:M
	ss = fileLines(w(d) + 1:w(d + 1) - 1);
	r = zeros(numBins, 1); Cr = zeros(numBins, 1);
	nn = 0;
	for jj = 1:length(ss)
		tmp = textscan(ss{jj}, '%f%f');
		if all(cellfun(@isempty, tmp))
			break
		end
		nn = nn + 1;
		r(nn) = tmp{1};
		Cr(nn) = tmp{2};
	end
	if nn ~= numBins
		error('TISEAN function ''d2'' returned an unexpected number of length scales.');
	end
	if d == 1
		co_logr = log(r); % .c2 stores raw (geometrically-spaced) r and C(r), both need log()
	end
	co_logC(:, d) = log(Cr);
end

% A zero-pair bin (Cr = 0, so log(Cr) = -Inf) would silently poison every
% downstream mean/min/range/polyfit statistic that touches it. The minEps
% floor above was tuned to avoid this in practice (see comment there), but
% as a safety net for configurations it doesn't fully cover, treat a
% residual non-finite value as an honest "not enough data to estimate a
% correlation-dimension curve at every requested embedding dimension" and
% bail out, rather than silently propagating -Inf/NaN into the fits below:
if any(~isfinite(bc_logN(:))) || any(~isfinite(co_logC(:)))
	fprintf(1, 'No good outputs obtained from the box-counting/correlation dimension curves.\n');
	out = NaN; return
end

if doPlot
	plot(bc_logr, bc_logN, 'o-')
	input('BC')
	plot(co_logr, co_logC, 'o-')
	input('CO')
end

% *** We now have to look for scaling regimes in each of these dimensions

% ------------------------------------------------------------------------------
%% Basic statistics on curves
% ------------------------------------------------------------------------------
out = struct;

% -------------------------------------------------------------------------------
%% How do curves change with m?
% -------------------------------------------------------------------------------
% Use SUB_mch

% Box counting dimension:
out = SUB_mch(bc_logr, bc_logN, 'bc', out);

% Correlation dimension:
out = SUB_mch(co_logr, co_logC, 'co', out);

% ------------------------------------------------------------------------------
%% What is the scaling range in r?
% ------------------------------------------------------------------------------
% ... and how good is the fit over this range?
% Use SUB_ScalingRange

% Box counting dimension, m = 1
out = SUB_ScalingRange(bc_logr, bc_logN(:, 1), 'scr_bc_m1', out);

% Box counting dimension m = 2
out = SUB_ScalingRange(bc_logr, bc_logN(:, 2), 'scr_bc_m2', out);

% Box counting dimension m = 3
out = SUB_ScalingRange(bc_logr, bc_logN(:, 3), 'scr_bc_m3', out);

% Box counting dimension m = chosen/given
out = SUB_ScalingRange(bc_logr, bc_logN(:, mopt), 'scr_bc_mopt', out);

% Correlation dimension, m = 1
out = SUB_ScalingRange(co_logr, co_logC(:, 1), 'scr_co_m1', out);

% Correlation dimension m = 2
out = SUB_ScalingRange(co_logr, co_logC(:, 2), 'scr_co_m2', out);

% Correlation dimension m = 3
out = SUB_ScalingRange(co_logr, co_logC(:, 3), 'scr_co_m3', out);

% Correlation dimension m = chosen/given
out = SUB_ScalingRange(co_logr, co_logC(:, mopt), 'scr_co_mopt', out);

% ------------------------------------------------------------------------------
%% What m gives best fit?
% ------------------------------------------------------------------------------
% Use SUB_bestm

% Box counting dimension
out = SUB_bestm(bc_logr, bc_logN, 'bc', out);

% Correlation dimension
out = SUB_bestm(co_logr, co_logC, 'co', out);

% ------------------------------------------------------------------------------
% ------------------------------------------------------------------------------
function out = SUB_mch(logr, logN, prefix, out)
	% looks at how changes with m. Since m will in general be different for each
	% different time series (i.e., if choosing an automatic method for
	% determining the embedding parameters), we have that m is at least
	% 3 here so that we can do statistics on at least these ones...

	% (i) on average the raw means at each m up to m = 3
	out.([prefix, '_meanm1']) = mean(logN(:, 1));
	out.([prefix, '_meanm2']) = mean(logN(:, 2));
	out.([prefix, '_meanm3']) = mean(logN(:, 3));
	out.([prefix, '_meanmmax']) = mean(logN(:, end));

	% (ii) raw minimum at each m up to m = 3
	out.([prefix, '_minm1']) = min(logN(:, 1));
	out.([prefix, '_minm2']) = min(logN(:, 2));
	out.([prefix, '_minm3']) = min(logN(:, 3));
	out.([prefix, '_minmmax']) = min(logN(:, end));

	% (iii) range at each m up to m = 3
	out.([prefix, '_range1']) = range(logN(:, 1));
	out.([prefix, '_range2']) = range(logN(:, 2));
	out.([prefix, '_range3']) = range(logN(:, 3));
	out.([prefix, '_rangemmax']) = range(logN(:, end));

	% (iv) increments with m
	out.([prefix, '_mindiff']) = mean([min(logN(:, 2)) - min(logN(:, 1)), min(logN(:, 3)) - min(logN(:, 2))]);
	out.([prefix, '_meandiff']) = mean([mean(logN(:, 2)) - mean(logN(:, 1)), mean(logN(:, 3)) - mean(logN(:, 2))]);

	% (v) slopes and goodness of fit across whole r range
	% logr/logN are both columns here (this used to rely on TSTOOL's
	% spacing() returning logr as a row, transposing logN(:,k) to match; with
	% logr now a column too, that stray transpose inside subsublinfit turned
	% "y - pfit" into an N-by-N broadcast instead of an N-by-1 residual, so
	% it's dropped):
	[out.([prefix, '_lfitm1']), out.([prefix, '_lfitb1']), out.([prefix, '_lfitmeansqdev1'])] = subsublinfit(logr, logN(:, 1));
	[out.([prefix, '_lfitm2']), out.([prefix, '_lfitb2']), out.([prefix, '_lfitmeansqdev2'])] = subsublinfit(logr, logN(:, 2));
	[out.([prefix, '_lfitm3']), out.([prefix, '_lfitb3']), out.([prefix, '_lfitmeansqdev3'])] = subsublinfit(logr, logN(:, 3));
	[out.([prefix, '_lfitmmax']), out.([prefix, '_lfitbmax']), out.([prefix, '_lfitmeansqdevmax'])] = subsublinfit(logr, logN(:, end));

	function [m, b, meansqdev] = subsublinfit(x, y)
		p1 = polyfit(x, y, 1);
		pfit = p1(1) * x + p1(2);
		res = y - pfit;
		m = p1(1); % gradient
		b = p1(2); % intercept
		meansqdev = mean(res.^2);
	end
end

% ------------------------------------------------------------------------------
function out = SUB_ScalingRange(logr, logN, prefix, out)
	% determines the scaling range in r for some m
	% we remove points from either extreme in r until minimize some
	% error measure
	% two dimensional optimization: over starting point and ending
	% point.

	l = length(logr);
	stptr = 1:floor(l / 2) - 1; % must be in the first half (not necessarily, but for here)
	endptr = ceil(l / 2) + 1:l; % must be in second half (not necessarily, but for here)
	mybad = zeros(length(stptr), length(endptr));
	for i = 1:length(stptr)
		for j = 1:length(endptr)
			% logr/logN both columns here -- see SUB_mch's comment on the
			% same stray-transpose issue:
			mybad(i, j) = lfitbadness(logr(stptr(i):endptr(j)), logN(stptr(i):endptr(j)));
		end
	end
	[a, b] = find(mybad == min(min(mybad))); % this defines the 'best' scaling range
	%         plot(logr,logN,'o-b'); hold on; plot(logr(stptr(a):endptr(b)),logN(stptr(a):endptr(b)),'o-r');
	%         hold off
	%         disp(['keep from ' num2str(stptr(a)) ' to ' num2str(endptr(b))])

	out.([prefix, '_logrmin']) = logr(stptr(a)); % minimum of scaling range
	out.([prefix, '_logrmax']) = logr(endptr(b)); % maximum of scaling range
	out.([prefix, '_logrrange']) = logr(endptr(b)) - logr(stptr(a)); % range of scaling... range
	out.([prefix, '_pgone']) = (stptr(a) - 1 + l - endptr(b)) / length(logr); % number of points removed in process
	% of choosing the optimum scaling range

	% Do the optimum fit again
	x = logr(stptr(a):endptr(b));
	y = logN(stptr(a):endptr(b));
	p = polyfit(x, y, 1);
	pfit = p(1) * x + p(2);
	res = pfit - y;
	out.([prefix, '_meanabsres']) = mean(abs(res));
	out.([prefix, '_meansqres']) = mean(res.^2);
	out.([prefix, '_scaling_exp']) = p(1);
	out.([prefix, '_scaling_int']) = p(2);
	out.([prefix, '_minbad']) = min(min(mybad));

	function badness = lfitbadness(x, y)
		gamma = 0.02; % reguralization parameter gamma selected empirically, could be tweaked in future work
		p = polyfit(x, y, 1);
		pfit = p(1) * x + p(2);
		res = pfit - y;
		badness = mean(abs(res)) - gamma * length(x); % want to still maximize length(x)
	end
end

% ------------------------------------------------------------------------------
function out = SUB_bestm(logr, logNN, prefix, out)
	% logNN is a matrix... logN is a vector for a given m
	% determines the scaling range in r for some m
	% we remove points from either extreme in r until minimize some
	% error measure
	% two dimensional optimization: over starting point and ending
	% point.

	store_scalingexps = zeros(size(logNN, 2), 1);
	store_meansqres = zeros(size(logNN, 2), 1);
	for k = 1:size(logNN, 2);
		logN = logNN(:, k); % take this element

		l = length(logr);
		stptr = 1:floor(l / 2) - 1; % must be in the first half (not necessarily, but for here)
		endptr = ceil(l / 2) + 1:l; % must be in second half (not necessarily, but for here)
		mybad = zeros(length(stptr), length(endptr));
		for i = 1:length(stptr)
			for j = 1:length(endptr)
				% logr/logN both columns here -- see SUB_mch's comment on the
				% same stray-transpose issue:
				mybad(i, j) = lfitbadness(logr(stptr(i):endptr(j)), logN(stptr(i):endptr(j)));
			end
		end
		[a, b] = find(mybad == min(min(mybad))); % this defines the 'best' scaling range

		% Do the optimum fit again
		x = logr(stptr(a):endptr(b));
		y = logN(stptr(a):endptr(b));
		p = polyfit(x, y, 1);
		pfit = p(1) * x + p(2);
		res = pfit - y;
		% subout.meanabsres = mean(abs(res));
		% subout.meansqres = mean(res.^2);
		% subout.scaling_exp = p(1);
		% subout.scaling_int = p(2);
		% subout.minbad = min(min(mybad));

		store_scalingexps(k) = p(1);
		store_meansqres(k) = mean(res.^2);
	end

	out.([prefix, '_minscalingexp']) = min(store_scalingexps);
	out.([prefix, '_meanscalingexp']) = mean(store_scalingexps);
	out.([prefix, '_maxscalingexp']) = max(store_scalingexps);
	out.([prefix, '_mbestfit']) = find(store_meansqres == min(store_meansqres), 1, 'first');

	function badness = lfitbadness(x, y)
		gamma = 0.02; % reguralization parameter gamma selected empirically, could be tweaked in future work
		p = polyfit(x, y, 1);
		pfit = p(1) * x + p(2);
		res = pfit - y;
		badness = mean(abs(res)) - gamma * length(x); % want to still maximize length(x)
	end
end

end
