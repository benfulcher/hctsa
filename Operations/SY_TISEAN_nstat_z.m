function out = SY_TISEAN_nstat_z(y,numSeg,embedParams)
% SY_TISEAN_nstat_z     Cross-forecast errors of zeroth-order time-series models
%
% Uses the nstat_z routine from the TISEAN package for nonlinear time-series
% analysis to calculate cross-forecast errors of zeroth-order models for the
% time-delay embedded time series.
%
% The program looks for nonstationarity in a time series by dividing it
% into a number of segments and calculating the cross-forecast errors
% between the different segments. The model used for the forecast is
% zeroth order model as proposed by Schreiber.
%
% cf. "Practical implementation of nonlinear time series methods: The TISEAN
% package", R. Hegger, H. Kantz, and T. Schreiber, Chaos 9(2) 413 (1999)
%
% Available here:
% http://www.mpipks-dresden.mpg.de/~tisean/Tisean_3.0.1/index.html
%
% The TISEAN routines are performed in the command line using 'system' commands
% in Matlab, and require that TISEAN is installed and compiled, and able to be
% executed in the command line.
%
%---INPUTS:
%
% y, the input time series
%
% numSeg, the number of equally-spaced segments to divide the time series into,
%       and used to predict the other time series segments
%
% embedParams, in the form {tau,m}, as usual for BF_embed. That is, for an
%               embedding dimension, tau, and embedding dimension, m. E.g.,
%               {1,3} has a time-delay of 1 and embedding dimension of 3.
%
%
%---OUTPUTS: include the trace of the cross-prediction error matrix, the mean,
% minimum, and maximum cross-prediction error, the minimum off-diagonal
% cross-prediction error, and eigenvalues of the cross-prediction error matrix.

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

% ------------------------------------------------------------------------------
%% Check Inputs / Set defaults
% ------------------------------------------------------------------------------

if nargin < 1
    error('Input a time series')
end

if nargin < 2 || isempty(numSeg)
    numSeg = 5; % divide the data into 5 segments by default
end

if nargin < 3
    embedParams = {1,3};
    fprintf(1,'Using default embedding using tau = 1 and m = 3\n');
end

N = length(y); % length of the time series

% ------------------------------------------------------------------------------
%% Write the file to disk for TISEAN to work with
% ------------------------------------------------------------------------------

% Write a temporary file in the system temp directory:
filePath = BF_WriteTempFile(y);
% fprintf(1,'Wrote temporary data file ''%s'' for TISEAN.\n',filePath)

% Get embedding parameters:
tm = BF_embed(y,embedParams{1},embedParams{2},2);
tau = tm(1); % time delay
m = tm(2); % embedding dimension

% ------------------------------------------------------------------------------
% Do some preliminary checks:
% ------------------------------------------------------------------------------
clength = (N-(m-1)*tau)/numSeg;
incStep = 1; % step increment
minNeighbors = 30; % minimum number of neighbors for fit

% Don't understand the c++ source code to work out exactly where the
% find_neighbors routine is crashing, but it's definitely to do with not being
% able to find enough neighbors and getting stuck in a while loop...
% Try this heuristic (multiplying the actual minimum number by 1.5):
if (clength-(m-1)*tau-incStep) <= minNeighbors*1.5
    delete(filePath); % remove the temporary file
    warning('Not enough neighbors to reliably estimate prediction errors with these settings');
    out = NaN; return
end

% ------------------------------------------------------------------------------
%% Do the calculation in the commandline
% ------------------------------------------------------------------------------

[~, res] = system(sprintf('nstat_z -# %u -d%u -m%u %s',numSeg,tau,m,filePath));
delete(filePath) % remove the temporary file filePath
if isempty(res), error('Call to TISEAN function ''nstat_z'' failed.'), end

% ------------------------------------------------------------------------------
%% Read the output from TISEAN
% ------------------------------------------------------------------------------
s = textscan(res,'%[^\n]'); s = s{1};
wi = strmatch('Writing to stdout',s);
if isempty(wi)
    error('TISEAN routine ''nstat_z'' returned unexpected output...');
end
s = s(wi+1:end);

xperr = zeros(numSeg); % cross prediction error from using segment i to forecast segment j

for i = 1:numSeg
    for j = 1:numSeg
        tmp = textscan(s{(i-1)*numSeg+j},'%n%n%n');
        xperr(i,j) = tmp{3};
    end
end

% pcolor(xperr)

% ------------------------------------------------------------------------------
%% Output statistics
% ------------------------------------------------------------------------------
% diagonal elements are using a segment to predict itself -- ought to be
% pretty good:
out.trace = sum(diag(xperr)); % trace
out.mean = mean(xperr(:));
out.median = median(xperr(:));

% minimum prediction error: the best you can do
out.min = min(xperr(:));
% maximum prediction error: the worst you can do
out.max = max(xperr(:));
% measures of spread of prediction error: stationarity
out.iqr = iqr(xperr(:));
out.std = std(xperr(:));
out.range = range(xperr(:));

% minimum prediction error not on diagonal
lowertri = tril(xperr,-1); lowertri = lowertri(lowertri>0);
uppertri = triu(xperr,1); uppertri = uppertri(uppertri>0);
offdiag = [lowertri; uppertri];
if isempty(lowertri)
    out.minlower = NaN;
else
    out.minlower = min(lowertri);
end
if isempty(uppertri)
    out.minupper = NaN;
else
    out.minupper = min(uppertri);
end
if isempty(offdiag)
    out.minoffdiag = NaN;
    out.iqroffdiag = NaN;
    out.stdoffdiag = NaN;
    out.rangeoffdiag = NaN;
else
    out.minoffdiag = min(offdiag);
    % measures of spread: non-stationarity
    out.iqroffdiag = iqr(offdiag);
    out.stdoffdiag = std(offdiag);
    out.rangeoffdiag = range(offdiag);
end

% Comparing columns/rows
out.stdmean = std(mean(xperr));
out.rangemean = range(mean(xperr));
out.stdmedian = std(median(xperr));
out.rangemedian = range(median(xperr));
out.rangerange = range(range(xperr));
out.stdrange = std(range(xperr));
out.rangestd = range(std(xperr));
out.stdstd = std(std(xperr));

% Eigenvalues
eigs = eig(xperr);
imagEigs = imag(eigs);
out.maximageig = max(imagEigs);
out.minimageig = min(imagEigs); % covaries (negatively) with maximageig

realEigs = real(eigs);
out.rangeeig = range(realEigs); % range of real parts of eigenvalues
out.stdeig = std(realEigs);
out.mineig = min(realEigs);
out.maxeig = max(realEigs);


end
