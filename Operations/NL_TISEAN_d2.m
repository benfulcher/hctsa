function out = NL_TISEAN_d2(y, tau, maxm, theilerWin)
% NL_TISEAN_d2  d2 routine from the TISEAN package.
%
% The function estimates the correlation sum, the correlation dimension and
% the correlation entropy of a given time series, y. Our code uses the outputs
% from this algorithm to return a set of informative features about the results.
%
%---INPUTS:
%
% y, input time series
%
% tau, time-delay (can be 'ac' or 'mi' for first zero-crossing of
%       autocorrelation function, or first minimum of the automutual
%       information)
%
% maxm, the maximum embedding dimension
%
% theilerWin, the Theiler window

% cf. "Practical implementation of nonlinear time series methods: The TISEAN
% package", R. Hegger, H. Kantz, and T. Schreiber, Chaos 9(2) 413 (1999)
%
% The TISEAN package is available here:
% http://www.mpipks-dresden.mpg.de/~tisean/Tisean_3.0.1/index.html
%
% The TISEAN routines are performed in the command line using 'system' commands
% in Matlab, and require that TISEAN is installed and compiled, and able to be
% executed in the command line.
%
% cf. "Spurious dimension from correlation algorithms applied to limited
% time-series data", J. Theiler, Phys. Rev. A, 34(3) 2427 (1986)
%
% cf. "Nonlinear Time Series Analysis", Cambridge University Press, H. Kantz
% and T. Schreiber (2004)
%
% Taken's estimator is computed for the correlation dimension, as well as related
% statistics, including other dimension estimates by finding appropriate scaling
% ranges, and searching for a flat region in the output of TISEAN's h2
% algorithm, which indicates determinism/deterministic chaos.
%
% To find a suitable scaling range, a penalized regression procedure is used to
% determine an optimal scaling range that simultaneously spans the greatest
% range of scales and shows the best fit to the data, and return the range, a
% goodness of fit statistic, and a dimension estimate.

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

N = length(y); % data length (number of samples)

% ------------------------------------------------------------------------------
%% Check inputs
% ------------------------------------------------------------------------------
% time delay, tau
if nargin < 2 || isempty(tau)
    tau = 1;
end
if strcmp(tau,'ac')
    tau = CO_FirstZero(y,'ac');
elseif strcmp(tau,'mi')
    tau = CO_FirstMin(y,'mi');
end

% Maximum embedding dimension
if nargin < 3 || isempty(maxm)
    maxm = 10;
end

% Theiler window
if nargin < 4 || isempty(theilerWin)
   theilerWin = 0.01; % Set a Theiler window of 1%% of the data length
end
if (theilerWin > 0) && (theilerWin < 1) % specify proportion of time-series length
    theilerWin = round(theilerWin*N);
end

% ------------------------------------------------------------------------------
%% Write the file
% ------------------------------------------------------------------------------
filePath = BF_WriteTempFile(y);
fprintf(1,'Wrote the input time series (N = %u) to the temporary file ''%s'' for TISEAN\n',length(y),filePath);

% ------------------------------------------------------------------------------
%% Run the TISEAN code, d2
% ------------------------------------------------------------------------------
[~, res] = system(sprintf('d2 -d%u -M1,%u -t%u %s',tau,maxm,theilerWin,filePath));
delete(filePath) % remove the temporary time-series data file
%  * extension .stat: This file shows the current status of the estimate.
if exist([filePath '.stat'],'file')
    delete([filePath '.stat']); % perhaps this file has something useful in it, but it's probably not for us...
end

if isempty(res) || ~isempty(regexp(res,'command not found', 'once')) % nothing came out??
    if exist([filePath '.c2'],'file'), delete([filePath '.c2']); end
    if exist([filePath '.d2'],'file'), delete([filePath '.d2']); end
    if exist([filePath '.h2'],'file'), delete([filePath '.h2']); end
    if isempty(res)
        error('Call to TISEAN function ''d2'' failed.')
    else
        error('Call to TISEAN function ''d2'' failed: %s',res)
    end
elseif strfind(res,'dyld: Library not loaded')
    error('DYLD library not found -- try recompiling TISEAN:\n%s',res);
end

% Check that all required files were generated (could not be due to problems with path or filename?)
if ~exist([filePath '.c2'],'file') || ~exist([filePath '.d2'],'file') || ~exist([filePath '.h2'],'file')
    % Delete all temporary files that exist:
    if exist([filePath '.c2'],'file'), delete([filePath '.c2']); end
    if exist([filePath '.d2'],'file'), delete([filePath '.d2']); end
    if exist([filePath '.h2'],'file'), delete([filePath '.h2']); end
    % Then throw an error:
    error([filePath,'.c2/.d2/.h2 not generated?']);
end


% this creates files in the local directory:
%  * extension .c2: This file contains the correlation sums for all treated length scales and embedding dimensions.
%  * extension .d2: This file contains the local slopes of the logarithm of the correlation sum, the correlation dimension.
%  * extension .h2: This file contains the correlation entropies.

%% Retreve output from file

% % (1) --------- C2 -----------
% fid_c2 = fopen([filePath '.c2']);
% s = textscan(fid_c2,'%[^\n]');
% if isempty(s)
%     disp(['Error reading TISEAN output file ' filePath '.c2'])
%     fclose(fid_c2); % close the file
%     delete([filePath '.c2']) % delete this file
%     delete([filePath '.d2']) % delete this file
%     delete([filePath '.h2']) % delete this file
%     return;
% end
% s = s{1};
% w = zeros(maxm+1,1);
% c2dat = SUB_readTISEANout(s,maxm,'#dim=',2);
% fclose(fid_c2); % close the file
%
% % c2dat now contains the correlation sums (second column)
% % as a function of the length scale epsilon (first column) for
% % the different embedding dimension (cell components 1--10)

% ------- GAUSSIAN KERNEL CORRELATION INTEGRAL -----------
% Now implement Gaussian Kernel Correlation integral
[~, res] = system(sprintf('c2g %s.c2',filePath));
% output is in res -- check it
s = textscan(res,'%[^\n]'); s = s{1};
wi = strmatch('writing to stdout',s);
s = s(wi+1:end);
if isempty(s) % TISEAN did produce valid output
    delete([filePath,'.c2']); delete([filePath,'.d2']); delete([filePath,'.h2']) % just in case these files were generated...
    error('TISEAN d2 produced invalid output (perhaps due to long tau = %u, N = %u)\n: %s',tau,N,res)
end
try
    c2gdat = SUB_readTISEANout(s,maxm,'#m=',3);
catch
    delete([filePath '.c2']); delete([filePath '.d2']); delete([filePath '.h2'])
    error('There are probably some Inf and NAN values in the TISEAN output files...?')
end

% c2gdat contains r (1), the Gaussian kernel correlation integral (2), and its
% logarithmic derivative with respect to r (3)

% ----- TAKENS MAXIMUM LIKELIHOOD ESTIMATOR FROM CORRELATION SUMS ----
% The integral is computed from the discrete values of C(r) by assuming an
% exact power law between the available points.
[~, res] = system(sprintf('c2t %s.c2',filePath));
% output is in res
s = textscan(res,'%[^\n]'); s = s{1};
wi = strmatch('writing to stdout',s);
s = s(wi+1:end);
c2tdat = SUB_readTISEANout(s,maxm,'#m=',2);
delete([filePath '.c2']);

% c2tdat contains upper length scale r (1), and the takens estimator (2)
% The integral is computed from the discrete values of C(r) by assuming an
% exact power law between the available points


% (2) --------- D2 ------------
fid_d2 = fopen([filePath '.d2']);
s = textscan(fid_d2,'%[^\n]');
s = s{1};
% FEED THIS INTO SUBROUTINE
d2dat = SUB_readTISEANout(s,maxm,'#dim=',2);
fclose(fid_d2); % close the file
delete([filePath '.d2']); % delete the file

% d2dat now contains the local slopes of the logarithm of correlation sums
% (second column) as a function of the length scale epsilon (first column) for
% the different embedding dimension (cell components 1--10)


% (3) -------------- H2 -------------

fid_h2 = fopen([filePath '.h2']);
s = textscan(fid_d2,'%[^\n]');
s = s{1};
% FEED THIS INTO SUBROUTINE
h2dat = SUB_readTISEANout(s,maxm,'#dim=',2);
fclose(fid_h2); % close the file
delete([filePath '.h2']); % delete the file

% h2dat now contains the correlation entropies (second column)

%% Time to obtain something useful from all this data

% ------------------------------------------------------------------------------
%% (1) TAKENS ESTIMATOR
% ------------------------------------------------------------------------------
% correlation dimension at upper length scale of eup
% (for z-scored time series, std = 1...; in units of this)
% Kantz & Shreiber recommend taking at half the std of the signal
takens05 = SUB_takens(c2tdat,0.5);
out.takens05_mean = mean(takens05);
out.takens05_median = median(takens05);
out.takens05_max = max(takens05);
out.takens05_min = min(takens05);
out.takens05_std = std(takens05);
out.takens05_iqr = iqr(takens05);

% Find outliers as means of inferring m_min
% looking for approaching a constant for m>m_min
mmintakens05 = SUB_findmmin(takens05);
out.takens05mmin_ri = mmintakens05.ri1; % minimum dimension to observe a scaling range
out.takens05mmin_goodness = mmintakens05.goodness;
out.takens05mmin_stabled = mmintakens05.stabled;
out.takens05mmin_linrmserr = mmintakens05.linrmserr;


% ------------------------------------------------------------------------------
%% (2) D2, local slopes of correlation integral
% ------------------------------------------------------------------------------
% semilogx(d2dat{1}(:,1),d2dat{1}(:,2))

% (2i) Estimate dimensions using Ben's method
% convert cell to matrix, taking second column in each case:
if all(cellfun(@isempty,d2dat))
    error('No data...')
end
[d2dat_v, d2dat_M] = SUB_celltomat(d2dat,2);

try
    benfindd2 = findscalingr_ind(d2dat_M);
catch
    error('Error finding scaling range')
end

% rows: increasing embedding m
% columns: stpt, endpt, goodness, dim
out.bend2_mindim = min(benfindd2(:,4));
out.bend2_maxdim = max(benfindd2(:,4));
out.bend2_meandim = mean(benfindd2(:,4));
out.bend2_meangoodness = mean(benfindd2(:,3));

mminfulcherd2 = SUB_findmmin(benfindd2(:,4));
if isempty(mminfulcherd2.ri1)
    out.benmmind2_logminl = NaN;
else
    out.benmmind2_logminl = log(d2dat_v(mminfulcherd2.ri1)); % minimum scale to observe a scaling range
end
out.benmmind2_goodness = mminfulcherd2.goodness;
out.benmmind2_stabledim = mminfulcherd2.stabled;
out.benmmind2_linrmserr = mminfulcherd2.linrmserr;

% Fulcher reshaped (frs) -- only for large enough m (as determined by
% criteria above):
d2dat_M_frs = d2dat_M(mminfulcherd2.ri1:end,:);
% find scaling region across m for a saturated range m.
scd2 = findscalingr(d2dat_M_frs);

if isempty(scd2.ri1);
    out.d2_logminscr = NaN;
else
    out.d2_logminscr = log(d2dat_v(scd2.ri1));
end
if isempty(scd2.ri2);
    out.d2_logmaxscr = NaN;
else
    out.d2_logmaxscr = log(d2dat_v(scd2.ri2));
end
if isnan(out.d2_logmaxscr) || isnan(out.d2_logminscr)
    out.d2_logscr = NaN;
else
    out.d2_logscr = out.d2_logmaxscr-out.d2_logminscr;
end

out.d2_goodness = scd2.goodness;
out.d2_dimest = scd2.dimest;
out.d2_dimstd = scd2.dimstd;


% hold off;
% semilogx(d2dat_v,d2dat_M,'o-m')
% hold on;
% semilogx(d2dat_v,d2dat_M_frs,'o-k')
% hold off;

% semilogx(c2tdat{1}(:,1),c2tdat{1}(:,2));

% ------------------------------------------------------------------------------
%% (3) Use Gaussian-smoothed estimates
% ------------------------------------------------------------------------------
% c2gdat
% we have the local slopes (d2) in the third column.
% Do all the same stuff as d2.

% (2i) Estimate dimensions using Ben's method
% convert cell to matrix, taking second column in each case:
[d2gdat_v, d2gdat_M] = SUB_celltomat(c2gdat,3);

try
    benfindd2g = findscalingr_ind(d2gdat_M);
catch
    error('Error finding scaling range')
    % out = NaN;
    % return
end

% rows: increasing embedding m
% columns: stpt, endpt, goodness, dim
out.bend2g_mindim = min(benfindd2g(:,4));
out.bend2g_maxdim = max(benfindd2g(:,4));
out.bend2g_meandim = mean(benfindd2g(:,4));
out.bend2g_meangoodness = mean(benfindd2g(:,3));

mminfulcherd2g = SUB_findmmin(benfindd2g(:,4));
if isempty(mminfulcherd2g.ri1)
    out.benmmind2g_logminl = NaN;
else
    out.benmmind2g_logminl = log(d2gdat_v(mminfulcherd2g.ri1)); % minimum scale to observe a scaling range
end
out.benmmind2g_goodness = mminfulcherd2g.goodness;
out.benmmind2g_stabledim = mminfulcherd2g.stabled;
out.benmmind2g_linrmserr = mminfulcherd2g.linrmserr;

% Fulcher reshaped (frs) -- only for large enough m (as determined by
% criteria above):
d2gdat_M_frs = d2gdat_M(mminfulcherd2g.ri1:end,:);
% find scaling region across m for a saturated range m.
scd2g = findscalingr(d2gdat_M_frs);

if isempty(scd2g.ri1)
    out.d2g_logminscr = NaN;
else
    out.d2g_logminscr = log(d2gdat_v(scd2g.ri1));
end
if isempty(scd2g.ri2)
    out.d2g_logmaxscr = NaN;
else
    out.d2g_logmaxscr = log(d2gdat_v(scd2g.ri2));
end
if any(isempty([scd2g.ri1, scd2g.ri2]))
    out.d2g_logscr = NaN;
else
    out.d2g_logscr = out.d2g_logmaxscr-out.d2g_logminscr;
end
out.d2g_goodness = scd2g.goodness;
out.d2g_dimest = scd2g.dimest;
out.d2g_dimstd = scd2g.dimstd;


% ------------------------------------------------------------------------------
%% (4) H2
% ------------------------------------------------------------------------------
% h2dat
% A flat region in this indicates determinism/deterministic chaos
[h2dat_v, h2dat_M] = SUB_celltomat(h2dat,2);
% semilogx(h2dat_v,h2dat_M,'ok')
% keyboard
% semilogx(h2dat_v(1:end-1),diff(h2dat_M(1,:)),'ok')
% keyboard
try h2results = SUB_getslopes(h2dat_v,h2dat_M);
catch emsg
    if strcmp(emsg.identifier,'MATLAB:subsassigndimmismatch')
        out = NaN;
        return
    end
end
% gets slopes for each dimension
slopesh2 = h2results(:,4);

% What are the (robust, mid-range) slopes like?
% plot(slopesh2); % dominant mid-range slope
findch_h2 = SUB_findmmin(slopesh2);
if isempty(findch_h2.ri1)
    out.slopesh2_ri1 = NaN;
else
    out.slopesh2_ri1 = findch_h2.ri1;
end
out.slopesh2_goodness = findch_h2.goodness;
out.slopesh2_stabled = findch_h2.stabled;
out.slopesh2_linrmserr = findch_h2.linrmserr;

% Are the any intermediate flat regions (signature of deterministic chaos)?
flattens = SUB_doesflatten(h2dat_v, h2dat_M);
out.h2meangoodness = mean(flattens(:,1)); % how close to having intermediate 'flat' regions
out.h2bestgoodness = min(flattens(:,1)); % best you can do
out.h2besth2 = flattens(find(flattens(:,1) == min(flattens(:,1)),1,'first'),2);
out.meanh2 = mean(flattens(:,2));
out.medianh2 = median(flattens(:,2));
flatsh2min = SUB_findmmin(flattens(:,2));
if isempty(flatsh2min.ri1)
    out.flatsh2min_ri1 = NaN;
else
    out.flatsh2min_ri1 = flatsh2min.ri1;
end
out.flatsh2min_goodness = flatsh2min.goodness;
out.flatsh2min_stabled = flatsh2min.stabled;
out.flatsh2min_linrmserr = flatsh2min.linrmserr;


% can look for local slopes using
% c2d:
% [pop res] = system(['H:\bin\c2d ' filePath '.c2']);


% ------------------------------------------------------------------------------
    function dimdat = SUB_readTISEANout(s,maxm,blocker,nc)
        % blocker the string distinguishing sections of output
        % nc number of columns in string

        w = strmatch(blocker,s);
        if length(w)~=maxm
            error('error reading TISEAN output');
        end
        w(end+1) = length(s)+1; % as if there were another marker at the entry after the last data row

        dimdat = cell(maxm,1); % stores data for each embedding dimension
        for ii = 1:maxm
            ss = s(w(ii)+1:w(ii+1)-1);
            nn = zeros(length(ss),nc);
            for jj = 1:length(ss)
                if nc == 2
                    tmp = textscan(ss{jj},'%f%f');
                elseif nc == 3
                    tmp = textscan(ss{jj},'%f%f%f');
                end
                if all(cellfun(@isempty,tmp))
                    % Ben Fulcher, 2015-03-06
                    % Sometimes a comment at the bottom of the output file
                    nn = nn(1:jj-1,:);
                    break
                else
                    nn(jj,:) = horzcat(tmp{:});
                end
            end
            dimdat{ii} = nn;
        end

    end
% ------------------------------------------------------------------------------

% ------------------------------------------------------------------------------
    function takensp = SUB_takens(dat,eup)
        % dat is takens estimator data, cell with a component corresponding to each
        % embedding dimension
        % eup is the cutoff length scale
        % returns a vector containing the dimension estimate at eup with
        % each element corresponding to an embedding dimension m up to the
        % maximum
        mmax = length(dat);
        takensp = zeros(mmax,1);
        for ii = 1:mmax
            theindex = find(dat{ii}(:,1)>eup,1,'first');
            if ~isempty(theindex)
                takensp(ii) = dat{ii}(theindex,2);
            else
                takensp(ii) = NaN;
            end
        end

    end
% ------------------------------------------------------------------------------

% ------------------------------------------------------------------------------
    function out = findscalingr(x)
        % finds constant regions in matrix x
        % if x a matrix, finds scaling regions requiring all columns to
        % match up. (i.e., to exhibit scaling at the same time)
        % starting point must be in first half of data
        % end point must be in last half of data

        l = size(x,2);
        gamma = 0.002; % regularization parameter selected empirically
        % for a consistent regularizer, need to data in [0,1] so weights
        % are consistent...
%         x=(x-min(x(:)))./(max(x(:))-min(x(:)));
        stptr = 1:floor(l/2)-1; % must be in the first half
        endptr = ceil(l/2)+1:l; % must be in second half
        mybad = zeros(length(stptr),length(endptr));
        for i = 1:length(stptr)
            for j = 1:length(endptr)
                % mean in this range -- mean (across range) of mean of
                % points (at each point)
                mm = mean(mean(x(:,stptr(i):endptr(j)))); % middle value for this range: exponent estimate
                nock = endptr(j)-stptr(i)+1; % extent of scaling range
                spreads = zeros(nock,1); % mean square deviation from mm at each point
                for k = 1:nock
                    spreads(k) = mean((x(:,stptr(i)+k-1)-mm).^2);
                end

                mybad(i,j) = mean(spreads)-gamma*nock; % want to still maximize length(x)
            end
        end
        [a, b] = find(mybad == min(min(mybad)),1,'first'); % this defines the 'best' scaling range
        ri1 = stptr(a); % minimum index of scaling range
        ri2 = endptr(b); % maximum index of scaling range
        out.ri1 = ri1;
        out.ri2 = ri2;
        out.goodness = min(mybad(:));
        out.dimest = mean(mean(x(:,ri1:ri2)));
        out.dimstd = std(mean(x(:,ri1:ri2)));


%         hold off;
%         plot(1:l,x,'o-k');
%         hold on;
%         plot(ri1:ri2,mean(mean(x(:,ri1:ri2)))*ones(ri2-ri1+1),'--r');
%         hold off;
%         keyboard
    end
% ------------------------------------------------------------------------------

% ------------------------------------------------------------------------------
    function results = findscalingr_ind(x)
        % AS ABOVE EXCEPT LOOKS FOR SCALING RANGES FOR INDIVIDUAL DIMENSIONS
        % finds constant regions in matrix x
        % if x a matrix, finds scaling regions requiring all columns to
        % match up. (i.e., to exhibit scaling at the same time)
        % starting point must be in first half of data
        % end point must be in last half of data

        l = size(x,2);      % number of distance/scaling points per dimension
        ndim = size(x,1);   % number of dimensions
        gamma = 1E-3;       % regularization parameter selected 'empirically'

        stptr = 1:floor(l/2)-1; % must be in the first half
        endptr = ceil(l/2)+1:l; % must be in second half
        results = zeros(ndim,4); %stpt, endpt, goodness, dim

        for c = 1:ndim
            mybad = zeros(length(stptr),length(endptr));
            v = x(c,:); % the vector of data for length scales
            vnorm = (v-min(v))./(max(v)-min(v)); % normalize regardless of range
            for i = 1:length(stptr)
                for j = 1:length(endptr)
                    mybad(i,j) = std(vnorm(stptr(i):endptr(j)))-gamma*(endptr(j)-stptr(i)+1);
                end
            end
            [a, b] = find(mybad == min(mybad(:)),1,'first'); % this defines the 'best' scaling range
            results(c,1) = stptr(a);
            results(c,2) = endptr(b);
            results(c,3) = min(mybad(:));
            results(c,4) = mean(v(stptr(a):endptr(b)));

%             hold off;
%             plot(1:l,v,'o-k');
%             hold on;
%             plot(stptr(a):endptr(b),mean(v(:,stptr(a):endptr(b)))*ones(endptr(b)-stptr(a)+1),'--r');
%             hold off;
%             keyboard
        end

    end
% ------------------------------------------------------------------------------

% ------------------------------------------------------------------------------
    function [thevector, thematrix] = SUB_celltomat(thecell,thecolumn)
        % converts cell to matrix, where each (specified) column in cell
        % becomes a column in the new matrix
%         thecelltest = thecell;
        % But higher dimensions may not reach low enough length scales
        % rescale range to greatest common span
        nn = length(thecell);
        mini = min(thecell{1}(:,1));
        maxi = max(thecell{1}(:,1));
        for ii = 2:nn
            mini = max([mini min(thecell{ii}(:,1))]);
            maxi = min([maxi max(thecell{ii}(:,1))]);
        end
        for ii = 1:nn % rescales each dimension so all share common scale
            thecell{ii} = thecell{ii}(thecell{ii}(:,1) >= mini & thecell{ii}(:,1) <= maxi,:);
        end

        thevector = thecell{1}(:,1);
        ee = length(thevector);

        goodones = cellfun(@(x)length(x),thecell) == ee;
        if ~all(goodones)
            % there's a bug in TISEAN where sometimes there are repeated 'x'
            % values -- check for this
            theproblems = find(goodones == 0);
            for ii = 1:length(theproblems)
                [~, m] = unique(thecell{theproblems(ii)}(:,1));
                thecell{theproblems(ii)} = thecell{theproblems(ii)}(m,:);
            end
        end

        thematrix = zeros(nn,ee); % across the rows for dimensions; across columns for lengths/epsilons
        for ii = 1:nn
            try thematrix(ii,:) = thecell{ii}(:,thecolumn);
            catch
                return
            end
        end

    end
% ------------------------------------------------------------------------------

% ------------------------------------------------------------------------------
    function out = SUB_findmmin(ds)
        % estimated dimensions for d = 1, ... , m
        % estimates when they stabilize to a limiting value

        % algorithm leaves out starting at the start, progressively,
        % to minimize variance

        l = length(ds);
        gamma = 0.1; % regularizer: CHOSEN AD HOC!! (maybe it's nicer to say 'empirically'...)
        dsraw = ds; % before normalization
        ds = (ds-min(ds))./(max(ds)-min(ds)); % rescale data to [0,1] so weights are consistent
        stptr = 1:l-1; % need at least two points
        mybad = zeros(length(stptr),1);
        for ii = 1:length(stptr)
                % mean in this range -- mean (across range) of mean of
                % points (at each point)
                mybad(ii) = std(ds(stptr(ii):end))-gamma*(l-stptr(ii)+1);
                % (*gamma): want to still maximize extent of constant region
        end
        a = find(mybad == min(mybad),1,'first'); % this defines the 'best' scaling range
        out.ri1 = stptr(a); % minimum index of scaling range
        out.goodness = min(mybad);
        out.stabled = mean(dsraw(a:end));

%         hold off;
%         plot(1:l,ds,'o-k');
%         hold on;
%         plot(stptr(a):l,mean(ds(a:end))*ones(l-stptr(a)+1),'--r');
%         hold off;
%         keyboard
        % how linear is it?
        p = polyfit((1:l)',ds,1);
        pfit = polyval(p,1:l);
        out.linrmserr = sqrt(mean((ds'-pfit).^2));

    end
% ------------------------------------------------------------------------------

% ------------------------------------------------------------------------------
    function results = SUB_getslopes(x,Y)
        dx = log10(x(2))-log10(x(1));
%         dx = x(2) - x(1);
        ndim = size(Y,1); % number of embeding dimensions
        gamma = 2E-3; % regularizer, chosen 'empirically' (i.e., ad hoc)
        l = size(Y,2)-1; % number of distance/scaling points per dimension
        stptr = 1:floor(l/2)-1; % must be in the first half
        endptr = ceil(l/2)+1:l; % must be in second half
        results = zeros(ndim,4); %stpt, endpt, goodness, dim

        for c = 1:ndim
            % () find best scaling region in which to estimate gradient

            mybad = zeros(length(stptr),length(endptr));
            v = diff(Y(c,:)).*dx; % make transformation to vector of local gradients
            vnorm = (v-min(v))./(max(v)-min(v)); % normalize regardless of range
            for i = 1:length(stptr)
                for j = 1:length(endptr)
                    mybad(i,j) = std(vnorm(stptr(i):endptr(j)))-gamma*(endptr(j)-stptr(i)+1);
                end
            end
            [a, b] = find(mybad == min(mybad(:)),1,'first'); % this defines the 'best' scaling range
            results(c,1) = stptr(a);
            results(c,2) = endptr(b);
            results(c,3) = min(mybad(:));
            results(c,4) = mean(v(stptr(a):endptr(b)));

%             hold off;
%             plot(v,'.k')
%             hold on;
%             plot(stptr(a):endptr(b),mean(v(stptr(a):endptr(b)))*ones(endptr(b)-stptr(a)+1),'--r');
%             hold off;
%             keyboard
        end
    end
% ------------------------------------------------------------------------------

% ------------------------------------------------------------------------------
    function results = SUB_doesflatten(x,Y)
        % look for region of zero gradient amidst regions of negative
        % gradient -- e.g., by two moving boundaries and a t-test between
        % them... (this would be a better way, perhaps)
        % for each embedding dimension (size(Y,1)), returns the goodness of
        % flatness in first column, and the mean of the quantity in Y
        % across this best range

        dx = log10(x(2))-log10(x(1));
        ndim = size(Y,1); % number of embeding dimensions
        l = size(Y,2)-1; % number of distance/scaling points per dimension
        stptr = 5:floor(l/2)-1; % must be in the first half
        endptr = ceil(l/2)+1:l-5; % must be in second half
        results = zeros(ndim,2); % h2, goodness

        for c = 1:ndim
            % regions that deviate least from zero
            mybad = zeros(length(stptr),length(endptr));
            v = diff(Y(c,:)).*dx; % make transformation to vector of local gradients
            vnorm = abs(v)./max(abs(v));
            for i = 1:length(stptr)
                for j = 1:length(endptr)
                    mybad(i,j) = abs(mean(vnorm(stptr(i):endptr(j)))) ... % deviations from zero in middle region
                                 - abs(mean(vnorm(1:stptr(i)))) ... % minus deviations from outside regions
                                 - abs(mean((vnorm(endptr(j):end))));
%                                         -gamma*(endptr(j)-stptr(i)+1); %
%                                         bonus for longer intermediate
%                                         regions
                end
            end
            [a, b] = find(mybad == min(mybad(:)),1,'first'); % this defines the 'best' scaling range

            % best range for being near-flat: ri1--ri2
            ri1 = stptr(a);
            ri2 = endptr(b);

            % goodness function: is this range close to zero compared to rest?
            results(c,1) = min(mybad(:)); % goodness
            results(c,2) = mean(Y(c,ri1:ri2)); % only relevant if goodness is small enough

        end
    end
% ------------------------------------------------------------------------------

end
