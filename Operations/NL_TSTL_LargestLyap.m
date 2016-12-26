function out = NL_TSTL_LargestLyap(y,Nref,maxtstep,past,NNR,embedParams)
% NL_TSTL_LargestLyap   Largest Lyapunov exponent of a time series.
%
%---INPUTS:
% y, the time series to analyze
% Nref, number of randomly-chosen reference points (-1 == all)
% maxtstep, maximum prediction length (samples)
% past, exclude -- Theiler window idea
% NNR, number of nearest neighbours
% embedParams, input to BF_embed, how to time-delay-embed the time series, in
%               the form {tau,m}, where string specifiers can indicate standard
%               methods of determining tau or m.
%
%---OUTPUTS: a range of statistics on the outputs from this function, including
% a penalized linear regression to the scaling range in an attempt to fit to as
% much of the range of scales as possible while simultaneously achieving the
% best possible linear fit.

% Uses the TSTOOL code 'largelyap'.
%
% TSTOOL: http://www.physik3.gwdg.de/tstool/
%
% The algorithm used (using formula (1.5) in Parlitz Nonlinear Time Series
% Analysis book) is very similar to the Wolf algorithm:
% "Determining Lyapunov exponents from a time series", A. Wolf et al., Physica D
% 16(3) 285 (1985)
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
%% Check a curve-fitting toolbox license is available:
% ------------------------------------------------------------------------------
BF_CheckToolbox('curve_fitting_toolbox');

doPlot = 0; % whether to plot outputs to a figure

% ------------------------------------------------------------------------------
%% Check inputs, preliminaries:
% ------------------------------------------------------------------------------
N = length(y); % length of time series

% (1) Nref: number of randomly-chosen reference points
if nargin < 2 || isempty(Nref)
    Nref = 0.5; % use half the length of the time series
end
if Nref < 1 && Nref > 0
    Nref = round(N*Nref); % specify a proportion of time series length
end
if Nref > -1
    fprintf('This algorithm relies on stochastic sampling of %u data points.\n',Nref)
end

% (2) maxtstep: maximum prediction length
if nargin < 3 || isempty(maxtstep)
    maxtstep = 0.1; % 10% length of time series
end
if maxtstep < 1 && maxtstep > 0
    maxtstep = round(N*maxtstep); % specify a proportion of time series length
end
if maxtstep < 10;
    maxtstep = 10; % minimum prediction length; for output stats purposes...
end
if maxtstep > 0.5*N;
    maxtstep = 0.5*N; % can't look further than half the time series length, methinks
end

% (3) past/theiler window
if nargin < 4 || isempty(past)
    past = 40;
end
if past < 1 && past > 0
    past = floor(past*N);
    if past == 0, past = 1; end
end

% (4) Number of neighest neighbours
if nargin < 5 || isempty(NNR)
    NNR = 3;
end

% (5) Embedding parameters, embedParams
if nargin < 6 || isempty(embedParams)
    embedParams = {'ac','fnnmar'};
    disp('using default embedding using autocorrelation and cao')
else
    if length(embedParams) ~= 2
        error('Embedding parameters formatted incorrectly -- should be {tau,m}')
    end
end

% ------------------------------------------------------------------------------
%% Embed the signal
% ------------------------------------------------------------------------------
% convert to embedded signal object for TSTOOL
s = BF_embed(y,embedParams{1},embedParams{2},1);

if ~isa(s,'signal') && isnan(s); % embedding failed
    error('Embedding failed');
end

% ------------------------------------------------------------------------------
%% Run the TSTOOL code, largelyap (which is stochastic):
% ------------------------------------------------------------------------------
try
    rs = largelyap(s,Nref,maxtstep,past,NNR);
catch
   disp('Error evaluating the TSTOOL method ''largelyap''')
   out = NaN;
   return
end

p = data(rs);
t = spacing(rs);

if doPlot
    figure('color','w'); box('on');
    plot(t,p,'.-k')
end

% we have the prediction error p as a function of the prediction length...?
% * function file says: output - vector of length taumax+1, x(tau) = 1/Nref *
%                                sum(log2(dist(reference point + tau, nearest neighbor +
%                                tau)/dist(reference point, nearest neighbor)))

% ------------------------------------------------------------------------------
%% Get output stats
% ------------------------------------------------------------------------------

if all(p == 0)
    out = NaN; return
end

% p at lags up to 5
% (note that p1 = 0, so not so useful to record)
for i = 1:5
    % evaluate p(1), p(2), ..., p(5) for the output structure
    out.(sprintf('p%u',i)) = p(i);
end
out.maxp = max(p);

% Number/proportion of crossings at 80% and 90% of maximum
ncrossx = @(x) sum((p(1:end-1)-x*max(p)).*(p(2:end)-x*max(p)) < 0);

out.ncross09maxold = sum((p(1:end-1)-0.9*max(p)).*(p(2:end)-0.9*max(p)) < 0);

out.ncross08max = ncrossx(0.8);
out.pcross08max = ncrossx(0.8)/(length(p)-1);

out.ncross09max = ncrossx(0.9);
out.pcross09max = ncrossx(0.9)/(length(p)-1);
% out.pcross08max = sum((p(1:end-1)-0.8*max(p)).*(p(2:end)-0.8*max(p)) < 0)/(length(p)-1);
% out.pcross09max = sum((p(1:end-1)-0.9*max(p)).*(p(2:end)-0.9*max(p)) < 0)/(length(p)-1);

% Time taken to get to n% maximum
ttomaxx = @(x) find(p > x*max(p),1,'first')-1;
out.to095max = ttomaxx(0.95);
% out.to095max = find(p > 0.95*max(p),1,'first')-1;
if isempty(out.to095max), out.to095max = NaN; end
out.to09max = ttomaxx(0.9);
% out.to09max = find(p > 0.9*max(p),1,'first')-1;
if isempty(out.to09max), out.to09max = NaN; end
out.to08max = ttomaxx(0.8);
% out.to08max = find(p > 0.8*max(p),1,'first')-1;
if isempty(out.to08max), out.to08max = NaN; end
out.to07max = ttomaxx(0.7);
% out.to07max = find(p > 0.7*max(p),1,'first')-1;
if isempty(out.to07max), out.to07max = NaN; end
out.to05max = ttomaxx(0.5);
% out.to05max = find(p > 0.5*max(p),1,'first')-1;
if isempty(out.to05max), out.to05max = NaN; end


% ------------------------------------------------------------------------------
%% Find scaling region:
% ------------------------------------------------------------------------------
% fit from zero to 95% of maximum...
imax = find(p > 0.95*max(p),1,'first');

if imax <= 3
    % not a suitable range for finding scaling
    % return NaNs for these
    out.vse_meanabsres = NaN;
    out.vse_rmsres = NaN;
    out.vse_gradient = NaN;
    out.vse_intercept = NaN;
    out.vse_minbad = NaN;

    out.ve_meanabsres = NaN;
    out.ve_rmsres = NaN;
    out.ve_gradient = NaN;
    out.ve_intercept = NaN;
    out.ve_minbad = NaN;
else
    t_scal = t(1:imax);
    p_scal = p(1:imax);
%     pp = polyfit(t_scal,p_scal',1); pfit = pp(1)*t_scal+pp(2);
    % hold on; plot(t_scal,p_scal,'.-r'); hold off
    % hold on; plot(t_scal,pfit,'-r'); hold off;
    % keyboard

    % ------------------------------------------------------------------------------
    %% Adjust start and end times for best scaling
    % ------------------------------------------------------------------------------

    l = imax; % = length(t_scal)
    stptr = 1:floor(l/2)-1; % start point must be in the first half (not necessarily, but for here)
    endptr = ceil(l/2)+1:l; % end point must be in second half (not necessarily, but for here)
    mybad = zeros(length(stptr),length(endptr));
    for i = 1:length(stptr)
        for j = 1:length(endptr)
            mybad(i,j) = lfitbadness(t_scal(stptr(i):endptr(j)),p_scal(stptr(i):endptr(j))');
        end
    end
    [a, b] = find(mybad == min(mybad(:))); % this defines the 'best' scaling range

    % Do the optimum fit again
    t_opt = t_scal(stptr(a):endptr(b));
    p_opt = p_scal(stptr(a):endptr(b))';
    pp = polyfit(t_opt,p_opt,1);
    pfit = pp(1)*t_opt+pp(2);
    res = pfit - p_opt;

    % hold on; plot(t_opt,p_opt,'og'); hold off;
    % hold on; plot(t_opt,pfit,'-g'); hold off;
    % vse == vary start and end times
    out.vse_meanabsres = mean(abs(res));
    out.vse_rmsres = sqrt(mean(res.^2));
    out.vse_gradient = pp(1);
    out.vse_intercept = pp(2);
    out.vse_minbad = min(mybad(:));
    if isempty(out.vse_minbad), out.vse_minbad = NaN; end

    %% Adjust just end time for best scaling
    imin = find(p > 0.50*max(p),1,'first');

    endptr = imin:imax; % end point is at least at 50% mark of maximum
    mybad = zeros(length(endptr),1);
    for i = 1:length(endptr)
        mybad(i) = lfitbadness(t_scal(1:endptr(i)),p_scal(1:endptr(i))');
    end
    b = find(mybad == min(mybad(:))); % this defines the 'best' scaling range

    % Do the optimum fit again
    t_opt = t_scal(1:endptr(b));
    p_opt = p_scal(1:endptr(b))';
    pp = polyfit(t_opt,p_opt,1);
    pfit = pp(1)*t_opt+pp(2);
    res = pfit-p_opt;

    % hold on; plot(t_opt,p_opt,'om'); hold off;
    % hold on; plot(t_opt,pfit,'-m'); hold off;
    out.ve_meanabsres = mean(abs(res));
    out.ve_rmsres = sqrt(mean(res.^2));
    out.ve_gradient = pp(1);
    out.ve_intercept = pp(2);
    out.ve_minbad = min(mybad(:));
    if isempty(out.ve_minbad), out.ve_minbad = NaN; end

end

% Fit exponential
s = fitoptions('Method','NonlinearLeastSquares','StartPoint',[max(p) -0.5]);
f = fittype('a*(1-exp(b*x))','options',s);
fitWorked = 1;
try
    [c, gof] = fit(t',p,f);
catch
    fitWorked = 0;
end
if fitWorked
    out.expfit_a = c.a;
    out.expfit_b = c.b;
    out.expfit_r2 = gof.rsquare;
    out.expfit_adjr2 = gof.adjrsquare;
    out.expfit_rmse = gof.rmse;
else
    out.expfit_a = NaN;
    out.expfit_b = NaN;
    out.expfit_r2 = NaN;
    out.expfit_adjr2 = NaN;
    out.expfit_rmse = NaN;
end

if doPlot
    hold on
    plot(t,c.a*(1-exp(c.b*t)),':r');
    hold off
end

    % ------------------------------------------------------------------------------
    function badness = lfitbadness(x,y,gamma)
        if nargin < 3,
            gamma = 0.006; % regularization parameter, gamma, chosen empirically
        end
        pp = polyfit(x,y,1);
        pfit = pp(1)*x + pp(2);
        res = pfit - y;
        badness = mean(abs(res)) - gamma*length(x); % want to still maximize length(x)
    end
    % ------------------------------------------------------------------------------

end
