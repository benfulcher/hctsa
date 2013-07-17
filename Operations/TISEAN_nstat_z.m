function out = TISEAN_nstat_z(y,nseg,embedparams)

% Exploits the TISEAN routine nstat_z:
% http://www.mpipks-dresden.mpg.de/~tisean/Tisean_3.0.1/index.html
% This program looks for nonstationarity in a time series by dividing it
% into a number of segments and calculating the cross-forecast errors
% between the different segments. The model used for the forecast is
% zeroth order model as proposed by Schreiber.

% Inputs:
% y: the column of time series data
% nseg: number of subsegments to divide the data into
% embedparams: embedding parameters in usual format for BF_embed

% Ben Fulcher 17/11/2009

if nargin < 1
    error('Input a time series')
end

if nargin < 2
    nseg = 5; % divide the data into 5 segments by default
end

%% Check Inputs / Set defaults
if nargin < 3
    embedparams = {1,3};
    fprintf(1,'Using default embedding using tau = 1 and m = 3\n')
end
tm = BF_embed(y,embedparams{1},embedparams{2},2);

%% Write the file
tnow = datestr(now,'yyyymmdd_HHMMSS_FFF');
% to the millisecond (only get double-write error for same function called in same millisecond
fn = sprintf('tisean_temp_nstat_z_%s.dat',tnow);
dlmwrite(fn,y);
fprintf(1,'Just written temporary file %s for TISEAN\n',fn)

%% Do the calculation
N = length(y); % length of the time series
if N/tm(1) < nseg*8 % heuristic
    % it may be more tm(1) itself rather than compared to nseg...?
    fprintf(1,'Time delay tau = %u too large for time series length, N = %u\n with %u segments',tm(1),N,nseg);
    delete(fn) % remove the temporary file fn
    out = NaN; return
end

% disp(['H:\bin\','stp -d' num2str(tm(1)) ' -m' num2str(tm(2)) ...
%                   ' -%' num2str(flevel) ' -t' num2str(tsteps) ' ' fn]);
[~, res] = system(sprintf('nstat_z -#%u -d%u -m%u %s',nseg,tm(1),tm(2),fn));
% [~, res] = system(sprintf('nstat_z -#' num2str(nseg) ' -d' num2str(tm(1)) ' -m' num2str(tm(2)) ' ' fn]);
delete(fn) % remove the temporary file fn
if isempty(res), error('Call to TISEAN function ''nstat_z'' failed.'), end


%% Read the input
s = textscan(res,'%[^\n]'); s = s{1};
wi = strmatch('Writing to stdout',s);
if isempty(wi)
    error('TISEAN routine ''nstat_z'' didn''t return what I expected...');
end
s = s(wi+1:end);

xperr = zeros(nseg); % cross prediction error from using segment i to forecast segment j

for i = 1:nseg
    for j = 1:nseg
        tmp = textscan(s{(i-1)*nseg+j},'%n%n%n');
        xperr(i,j) = tmp{3};
    end
end

% pcolor(xperr)


%% Output statistics
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


% comparing columns/rows
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
out.maximageig = max(imag(eigs));
out.minimageig = min(imag(eigs));

eigs = real(eigs);
out.rangeeig = range(eigs); % range of real parts of eigenvalues
out.stdeig = std(eigs);
out.mineig = min(eigs);
out.maxeig = max(eigs);


end