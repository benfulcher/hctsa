function out = TSTL_fracdims(y,kmin,kmax,Nref,gstart,gend,past,steps,embedparams)
% Uses TSTOOL code fracdims
% Computes fractal dimension spectrum D(q) using moments of neighbour
% distances for time-delay reconstructed time series y.
% y: column vector of time series data
% kmin: minimum number of neighbours for each reference point
% kmax: maximum number of neighbours for each reference point
% Nref: number of randomly-chosen reference points (-1: use all points)
% gstart: starting value for moments
% gend: end value for moments
% past [opt]: number of samples to exclude before an after each reference
% index (default=0)
% steps [opt]: number of moments to calculate (default=32);
% embedparams: how to embed the time series using a time-delay
% reconstruction

% Ben Fulcher November 2009

%% Preliminaries
N = length(y); % length of time series

% (1) Minimum number of neighbours, kmin
if nargin < 2 || isempty(kmin)
    kmin = 3; % default
    fprintf(1,'Using default, minimum number of neighbours, kmin = %u\n',kmin)
end

% (2) Maximum number of neighbours, kmax
if nargin < 3 || isempty(kmax)
    kmax = 10; % default
    fprintf(1,'Using default maximum number of neighbours, kmax = %u\n',kmax)
end

% (3) Number of randomly-chosen reference points, Nref
if nargin < 4 || isempty(Nref)
    Nref = 0.2; % default:  20% of the time series length
    fprintf(1,'Using default number of reference points: Nref = %f\n',Nref))
end
if Nref < 1 && Nref > 0
    Nref = round(N*Nref); % specify a proportion of time series length
end

% (4) moment starting value, gstart
if nargin < 5 || isempty(gstart)
    gstart = 1; % default
    fprintf(1,'Using default moment starting value, gstart = %u\n', gstart)
end

% (5) moment ending value, gend
if nargin < 6 || isempty(gend)
    gend = 10; % default
    fprintf(1,'Using default moment ending value, gend = %u\n', gend)
end

% (6) past
if nargin < 7 || isempty(past)
    past = 10; % default
    fprintf(1,'Using default past correlation exclusion window value, past = %u\n')
end

% (7) steps
if nargin < 8 || isempty(steps)
    steps = 32;
end

% (8) Embedding parameters
if nargin < 9 || isempty(embedparams)
    embedparams = {'ac','cao'};
    fprintf(1,'Using default embedding parameters of autocorrelation for tau and cao method for m')
end


%% Embed the signal
% convert scalar time series y to embedded signal object s for TSTOOL
s = BF_embed(y,embedparams{1},embedparams{2},1);

if ~strcmp(class(s),'signal') && isnan(s); % embedding failed
    error('Embedding failed')
end


%% Run
try
    rs = fracdims(s,kmin,kmax,Nref,gstart,gend,past,steps);
catch me
    if strcmp(me.message,'Fast nearest neighbour searcher : To many neighbors for each query point are requested')
        out = NaN;
        return
    end
end
% view(rs);

Dq = data(rs);
q = spacing(rs);
% plot(q,dq);

%% Get output stats
out.rangeDq = range(Dq);
out.minDq = min(Dq);
out.maxDq = max(Dq);
out.meanDq = mean(Dq);

out.minq = min(q);
out.maxq = max(q);
out.rangeq = range(q);
out.meanq = mean(q);

% fit linear
p = polyfit(q,Dq',1);
p_fit = q*p(1) + p(2);
res = p_fit - Dq';
out.linfit_a = p(1);
out.linfit_b = p(2);
out.linfit_rmsqres = sqrt(mean(res.^2));

% fit exponential
s = fitoptions('Method','NonlinearLeastSquares','StartPoint',[range(Dq) -0.5 min(Dq)]);
f = fittype('a*exp(b*x)+c','options',s);
[c, gof] = fit(q',Dq,f);
out.expfit_a = c.a;
out.expfit_b = c.b;
out.expfit_c = c.c;
out.expfit_r2 = gof.rsquare; % this is more important!
out.expfit_adjr2 = gof.adjrsquare;
out.expfit_rmse = gof.rmse;


end