function out = TSTL_corrsum(y,Nref,r,thwin,nbins,embedparams,dotwo)
% Uses TSTOOL code corrsum (or corrsum2)
% uses Grassberger-Proccacia algorithm using fast nearest neighbour search
% y: column vector of time series data
% Nref: number of (randomly-chosen) reference points (-1: use all points,
%       if a decimal, then use this fraction of the time series length)
% r: 0<r<1, maximum search radius relative to attractor size
% thwin: number of samples to exclude before and after each reference index
% (~ Theiler window)
% nbins: number of partitioned bins
% embedparams: embedding parameters to feed benembed.m for embedding the
% signal in the form {tau,m}
% dotwo: if this is set to 2, will use corrsum2 instead of corrsum, in which
% case n specifies the number of pairs per bin. Default is 1, i.e., to use
% corrsum.
% Ben Fulcher November 2009

%% Preliminaries
N = length(y); % length of time series

% (1) Number of reference points, Nref
if nargin < 2 || isempty(Nref)
    Nref = 500; % 500 points
end
if Nref < 1 && Nref > 0
    Nref = round(N*Nref); % specify a proportion of time series length
end
if Nref >= N
    Nref = -1; % Number of reference points capped at time series length
end

% (2) Maximum relative search radius, r
if nargin < 3 || isempty(r)
    r = 0.05; % 5% of attractor radius
end

% (3) Remove spurious correlations of adjacent points, thwin
if nargin < 4 || isempty(thwin)
    thwin = 10; % default window length
end

% (4) Number of bins, nbins
if nargin < 5 || isempty(nbins)
    nbins = 20; % defulat number of bins
end

% (5) Set embedding parameters to defaults
if nargin < 6 || isempty(embedparams)
    embedparams = {'ac','cao'};
else
    if length(embedparams)~=2
        disp('given embedding parameters incorrectly formatted -- need {tau,m}')
    end
end

if nargin < 7 || isempty(dotwo)
    dotwo = 1; % use corrsum rather than corrsum2
end

if Nref == -1 && dotwo == 2
    % we need a *number* of pairs for corrsum2, round down from 50% of time series
    % length
    Nref = floor(N*0.5);
end

%% Embed the signal
% convert to embedded signal object for TSTOOL
s = benembed(y,embedparams{1},embedparams{2},1);

if ~strcmp(class(s),'signal') && isnan(s); % embedding failed
    error('Embedding failed')
elseif length(data(s)) < thwin
    fprintf(1,'Embedded time series (N = %u, m = %u, tau = %u) too short to do a correlation sum\n',N,embedparams{1},embedparams{2})
    out = NaN; return
end

%% Run
me = []; % error catch
if dotwo == 1 % use corrsum
    try
        rs = corrsum(s,Nref,r,thwin,nbins);
    catch me % DEAL WITH ERROR MESSAGE BELOW
    end
elseif dotwo == 2 % use corrsum2
    try
        rs = corrsum2(s,Nref,r,thwin,nbins);
    catch me
    end
end

if ~isempty(me) && strcmp(me.message,'Maximal search radius must be greater than starting radius')
    disp('Max search radius less than starting radius. Returning NaNs.')
    out = NaN;
    return
elseif ~isempty(me) && strcmp(me.message,'Cannot find an interpoint distance greater zero, maybe ill-conditioned data set given')
    disp('Cannot find an interpoint distance greater than zero. Shit. Returning NaNs.')
    out = NaN;
    return
elseif ~isempty(me) && strcmp(me.message,'Reference indices out of range')
    disp('Reference indicies out of range. Returning NaNs.')
    out = NaN;
    return
end

lnr = spacing(rs);
lnCr = data(rs);
% plot(lnr,lnCr);
% keyboard
% Contains ln(r) in rows and values are ln(C(r));
% keyboard

%% remove any Infs in lnCr
rgood = find(isfinite(lnCr));
if isempty(rgood)
    out = NaN; return
end
lnCr = lnCr(rgood);
lnr = lnr(rgood);

%% Output Statistics
% basic
out.minlnr = min(lnr);
out.maxlnr = max(lnr);
out.minlnCr = min(lnCr);
out.maxlnCr = max(lnCr);
out.rangelnCr = range(lnCr);
out.meanlnCr = mean(lnCr);


% fit linear to log-log plot (full range)
enoughpoints=1;
try
    [a, stats] = robustfit(lnr,lnCr);
catch me
    if strcmp(me.message,'Not enough points to perform robust estimation.')
        enoughpoints = 0;
    end
end

if enoughpoints
    out.robfit_a1 = a(1);
    out.robfit_a2 = a(2);
    out.robfit_sigrat = stats.ols_s/stats.robust_s;
    out.robfit_s = stats.s;
    out.robfit_sea1 = stats.se(1);
    out.robfit_sea2 = stats.se(2);


    fit_lnCr = a(2)*lnr+a(1);
    % hold on;plot(lnr,fit_lnCr,'r');hold off
    res = lnCr-fit_lnCr';
    out.robfitresmeanabs = mean(abs(res));
    out.robfitresmeansq = mean(res.^2);
    out.robfitresac1 = CO_autocorr(res,1);
else
    out.robfit_a1 = NaN;
    out.robfit_a2 = NaN;
    out.robfit_sigrat = NaN;
    out.robfit_s = NaN;
    out.robfit_sea1 = NaN;
    out.robfit_sea2 = NaN;
    
    out.robfitresmeanabs = NaN;
    out.robfitresmeansq = NaN;
    out.robfitresac1 = NaN;
end
    

% now non-robust linear fit
% [p S] = polyfit(lnr',lnCr,1);


%     function out = SUB_allNaNs
%         % return real NaNs
%         out.robfit_a1 = NaN;
%         out.robfit_a2 = NaN;
%         out.robfit_sigrat = NaN;
%         out.robfit_s = NaN;
%         out.robfit_sea1 = NaN;
%         out.robfit_sea2 = NaN;
%         out.robfitresmeanabs = NaN;
%         out.robfitresmeansq = NaN;
%         out.robfitresac1 = NaN;
%         out.minlnr = NaN;
%         out.maxlnr = NaN;
%         out.minlnCr = NaN;
%         out.maxlnCr = NaN;
%         out.rangelnCr = NaN;
%         out.meanlnCr = NaN; 
%     end

end