function out = NL_TSTL_PoincareSection(y,ref,embedParams)
% NL_TSTL_PoincareSection   Poincare sectino analysis of a time series.
%
% Obtains a Poincare section of the time-delay embedded time series, producing a
% set of vector points projected orthogonal to the tangential vector at the
% specified index using TSTOOL code 'poincare'.
%
%---INPUTS:
% y, the input time series
%
% ref: the reference point. Can be an absolute number (2 takes the second point
%      in the (embedded) time series) or a string like 'max' or 'min' that takes
%      the first maximum, ... in the (scalar) time series, ...
%
% embedParams: the usual thing to give BF_embed for the time-delay embedding, as
%               {tau,m}. A common choice for m is 3 -- i.e., embed in a 3
%               dimensional space so that the Poincare section is 2-dimensional.
%
%---OUTPUTS: include statistics on the x- and y- components of these vectors on the
% Poincare surface, on distances between adjacent points, distances from the
% mean position, and the entropy of the vector cloud.
%
% TSTOOL: http://www.physik3.gwdg.de/tstool/

% Another thing that could be cool to do is to analyze variation in the plots as
% ref changes... (not done here)
%
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
%% Check inputs
% ------------------------------------------------------------------------------

if nargin < 2 || isempty(ref)
    ref = 'max'; % reference point is the first maximum of the time series
end

if nargin < 3 || isempty(embedParams)
    embedParams = {'mi', 3};
    fprintf(1,'Using default embedding settings: minimum of the automutual information for tau and m = 3\n');
end

if embedParams{2} ~= 3
    embedParams{2} = 3;
    fprintf(1,'Three-dimensional embedding\n');
end

% set ref
if ischar(ref)
    switch ref
        case 'max'
            % first local maximum
            dydt = diff(y);
            ref = find(dydt(1:end-1)>=0 & dydt(2:end)<0,1,'first')+1;

        case 'min'
            % first local minimum
            dydt = diff(y);
            ref = find(dydt(1:end-1)<=0 & dydt(2:end)>0,1,'first')+1;

        otherwise
            error('Unknown reference setting ''%s''',ref);
    end
end

if isempty(ref), out = NaN; return; end % ridiculous
if ref < 2, ref = 2; end % gives an error, uses previous value in algorithm

doPlot = 0; % plot outputs to a figure

% ------------------------------------------------------------------------------
%% Do your magic, TSTOOL!:
% ------------------------------------------------------------------------------
% time-delay embed the signal:
N = length(y); % length of the time series
s = BF_embed(y,embedParams{1},3,1);

if ~isa(s,'signal') && isnan(s); % embedding failed
    error('Embedding failed');
end

% Run external TSTOOL code, poincare:
try
    rs = poincare(s,ref);
catch me
    if strcmp(me.message,'No section points found') ...
            || strcmp(me.identifier,'MATLAB:badsubscript')
        fprintf(1,'No section points found to run NL_TSTL_PoincareSection\n');
        out = NaN; return
    else
        error(me.message)
    end
end

% Convert back to Matlab forms
v = data(rs); % vectors on poincare surface
NN = length(v);
% Labeling poincare surface plane x-y
x = (v(:,1));
y = (v(:,2));

if doPlot
    figure('color','w'); box('on');
    plot(x,y,'.k'); axis equal
end

% ------------------------------------------------------------------------------
%% Get statistics out
% ------------------------------------------------------------------------------

% Basic statistics:
out.pcross = NN/N; % proportion of time series that crosses poincare surface

out.maxx = max(x);
out.minx = min(x);
out.stdx = std(x);
out.iqrx = iqr(x);
out.meanx = mean(x);
out.ac1x = CO_AutoCorr(x,1,'Fourier');
out.ac2x = CO_AutoCorr(x,2,'Fourier');
out.tauacx = CO_FirstZero(x,'ac');

out.maxy = max(y);
out.miny = min(y);
out.stdy = std(y);
out.iqry = iqr(y);
out.meany = mean(y);
out.ac1y = CO_AutoCorr(y,1,'Fourier');
out.ac2y = CO_AutoCorr(y,2,'Fourier');
out.tauacy = CO_FirstZero(y,'ac');

out.boxarea = range(x)*range(y);

% Statistics on distance between adjacent points, ds
vdiff = v(2:end,:)-v(1:end-1,:);
ds = sqrt(vdiff(:,1).^2 + vdiff(:,2).^2);

% Probability that next point in series is within radius r of current point in
% the poincare section:
out.pwithinr01 = sum(ds<0.1)/(NN-1);
out.pwithin02 = sum(ds<0.2)/(NN-1);
out.pwithin03 = sum(ds<0.3)/(NN-1);
out.pwithin05 = sum(ds<0.5)/(NN-1);
out.pwithin1 = sum(ds<1)/(NN-1);
out.pwithin2 = sum(ds<2)/(NN-1);
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
D = sqrt(x.^2+y.^2); % distance from (mean,mean)
out.maxD = max(D);
out.minD = min(D);
out.stdD = std(D);
out.iqrD = iqr(D);
out.meanD = mean(D);
out.ac1D = CO_AutoCorr(D,1,'Fourier');
out.ac2D = CO_AutoCorr(D,2,'Fourier');
out.tauacD = CO_FirstZero(D,'ac');

% ------------------------------------------------------------------------------
%% Statistics of the boxed distribution:
% ------------------------------------------------------------------------------

numPartitions = [5,10];
% (i) 5 partitions per axis
% (ii) 10 partitions per axis

for i = 1:length(numPartitions)
    boxcounts = subcountboxes(x,y,numPartitions(i));
    pbox = boxcounts/NN;

    out.(sprintf('maxpbox%u',numPartitions(i))) = max(pbox(:));
    out.(sprintf('minpbox%u',numPartitions(i))) = min(pbox(:));
    out.(sprintf('zerospbox%u',numPartitions(i))) = sum(pbox(:) == 0);
    out.(sprintf('meanpbox%u',numPartitions(i))) = mean(pbox(:));
    out.(sprintf('rangepbox%u',numPartitions(i))) = range(pbox(:));
    % This probably needs to be normalized:
    out.(sprintf('hboxcounts%u',numPartitions(i))) = -sum(pbox(pbox > 0).*log(pbox(pbox > 0)));
    out.(sprintf('tracepbox%u',numPartitions(i))) = sum(diag(pbox)); % trace
end

% NOTE: max and range provide almost the same information on real data.

% can imagine doing many more things; like seeing different slabs of the
% space, or finding the line whos vicinity includes many points, etc. but I
% think this is enough for now.

% ------------------------------------------------------------------------------
function boxcounts = subcountboxes(x,y,nbox)
    boxcounts = zeros(nbox);
    % boxes are quantiles along each axis
    xbox = quantile(x,linspace(0,1,nbox+1));
    ybox = quantile(y,linspace(0,1,nbox+1));
    xbox(end) = xbox(end)+1;
    ybox(end) = ybox(end)+1;

    for ii = 1:nbox % x
        rx = (x >= xbox(ii) & x < xbox(ii+1)); % these x are in range
        for jj = 1:nbox % y
            % only need to look at those ys for which the xs are in range
            boxcounts(ii,jj) = sum(y(rx) >= ybox(jj) & y(rx) < ybox(jj+1));
        end
    end
end
% ------------------------------------------------------------------------------

end
