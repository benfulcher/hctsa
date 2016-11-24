function out = NL_TSTL_FractalDimensions(y,kmin,kmax,Nref,gstart,gend,past,steps,embedParams)
% NL_TSTL_FractalDimensions    Fractal dimension spectrum, D(q), of a time series.
%
%---INPUTS:
% y, column vector of time series data
% kmin, minimum number of neighbours for each reference point
% kmax, maximum number of neighbours for each reference point
% Nref, number of randomly-chosen reference points (-1: use all points)
% gstart, starting value for moments
% gend, end value for moments
% past [opt], number of samples to exclude before an after each reference
%             index (default=0)
% steps [opt], number of moments to calculate (default=32);
% embedParams, how to embed the time series using a time-delay reconstruction
%
%---OUTPUTS: include basic statistics of D(q) and q, statistics from a linear fit,
% and an exponential fit of the form D(q) = Aexp(Bq) + C.

% Computes the fractal dimension spectrum, D(q), using moments of neighbor
% distances for time-delay embedded time series by referencing the code,
% fracdims, from the TSTOOL package.
%
% TSTOOL: http://www.physik3.gwdg.de/tstool/
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
% Check a curve-fitting toolbox license is available:
% ------------------------------------------------------------------------------
BF_CheckToolbox('curve_fitting_toolbox');

% ------------------------------------------------------------------------------
%% Preliminaries
% ------------------------------------------------------------------------------
N = length(y); % Length of time series
doPlot = 0; % Don't plot results by default

% ------------------------------------------------------------------------------
%% Check inputs
% ------------------------------------------------------------------------------
% (1) Minimum number of neighbours, kmin
if nargin < 2 || isempty(kmin)
    kmin = 3; % default
    fprintf(1,'Using default, minimum number of neighbours, kmin = %u\n',kmin);
end

% (2) Maximum number of neighbours, kmax
if nargin < 3 || isempty(kmax)
    kmax = 10; % default
    fprintf(1,'Using default maximum number of neighbours, kmax = %u\n',kmax);
end

% (3) Number of randomly-chosen reference points, Nref
if nargin < 4 || isempty(Nref)
    Nref = 0.2; % default:  20% of the time series length
    fprintf(1,'Using default number of reference points: Nref = %f\n',Nref);
end
if (Nref > 0) && (Nref < 1)
    Nref = round(N*Nref); % specify a proportion of time series length
end

% (4) moment starting value, gstart
if nargin < 5 || isempty(gstart)
    gstart = 1; % default
    fprintf(1,'Using default moment starting value, gstart = %u\n',gstart);
end

% (5) moment ending value, gend
if nargin < 6 || isempty(gend)
    gend = 10; % default
    fprintf(1,'Using default moment ending value, gend = %u\n',gend);
end

% (6) past
if nargin < 7 || isempty(past)
    past = 10; % default
    fprintf(1,'Using default past correlation exclusion window value, past = %u\n',past);
end

% (7) steps
if nargin < 8 || isempty(steps)
    steps = 32;
end

% (8) Embedding parameters
if nargin < 9 || isempty(embedParams)
    embedParams = {'ac','fnnmar'};
    fprintf(1,'Using default embedding parameters of autocorrelation for tau and cao method for m\n');
end


% ------------------------------------------------------------------------------
%% Embed the signal
% ------------------------------------------------------------------------------
% Convert the scalar time series, y, to embedded signal object s for TSTOOL
s = BF_embed(y,embedParams{1},embedParams{2},1);

if ~isa(s,'signal') && isnan(s); % embedding failed
    error('Embedding of the %u-sample time series failed',N)
end

% ------------------------------------------------------------------------------
%% Run the TSTOOL code, fracdims:
% ------------------------------------------------------------------------------
% Checks that tstoolbox/@signal/fracdims exists
if ~exist(fullfile('tstoolbox','@signal','fracdims'),'file')
    error(['Cannot find the code ''fracdims'' from the TSTOOL package. ' ...
            'Is it installed and in the Matlab path?']);
end
try
    rs = fracdims(s,kmin,kmax,Nref,gstart,gend,past,steps);
catch me
    if strcmp(me.message,['Fast nearest neighbour searcher : ' ...
            'To many neighbors for each query point are requested'])
        out = NaN; return
    else
        error('Error occurred calling fracdims: %s',me.message);
    end
end

Dq = data(rs);
q = spacing(rs);

% ------------------------------------------------------------------------------
% Plot the results in a figure:
% ------------------------------------------------------------------------------
if doPlot
    figure('color','w'); box('on');
    subplot(2,1,1); view(rs);
    subplot(2,1,2); plot(q,dq);
end

% ------------------------------------------------------------------------------
%% Get output statistics:
% ------------------------------------------------------------------------------
out.rangeDq = range(Dq);
out.maxDq = max(Dq);
out.meanDq = mean(Dq);

% out.minq = min(q);
out.maxq = max(q);
out.rangeq = range(q);
out.meanq = mean(q);

%---Fit linear
p = polyfit(q,Dq',1);
p_fit = q*p(1) + p(2);
res = p_fit - Dq';
out.linfit_a = p(1);
out.linfit_b = p(2);
out.linfit_rmsqres = sqrt(mean(res.^2));

%---Fit exponential
% s = fitoptions('Method','NonlinearLeastSquares','StartPoint',[range(Dq) -0.5 min(Dq)]);
% f = fittype('a*exp(b*x)+c','options',s);
% [c, gof] = fit(q',Dq,f);
% out.expfit_a = c.a;
% out.expfit_b = c.b;
% out.expfit_c = c.c;
% out.expfit_r2 = gof.rsquare; % I reckon this one is the most important!
% out.expfit_adjr2 = gof.adjrsquare;
% out.expfit_rmse = gof.rmse;


end
