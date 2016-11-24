function out = NL_TSTL_dimensions(y,nbins,embedParams)
% NL_TSTL_dimensions box counting, information, and correlation dimension of a time series.
%
% Computes the box counting, information, and correlation dimension of a
% time-delay embedded time series using the TSTOOL code 'dimensions'.
% This function contains extensive code for estimating the best scaling range to
% estimate the dimension using a penalized regression procedure.
%
%---INPUTS:
% y, column vector of time series data
% nbins, maximum number of partitions per axis
% embedParams, embedding parameters to feed BF_embed.m for embedding the
%              signal in the form {tau,m}
%
%---OUTPUTS:
% A range of statistics are returned about how each dimension estimate changes
% with m, the scaling range in r, and the embedding dimension at which the best
% fit is obtained.

% cf. TSTOOL, http://www.physik3.gwdg.de/tstool/
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
%% Preliminaries, check inputs
% ------------------------------------------------------------------------------
doPlot = 0; % plot outputs to screen

% (1) Maximum number of bins, nbins
if nargin < 2 || isempty(nbins)
    nbins = 50; % 50 points
    fprintf(1,'Using a default of 50 bins per axis\n');
end

% (2) Set embedding parameters to defaults
if nargin < 3 || isempty(embedParams)
    embedParams = {'ac','fnnmar'};
    fprintf(1,'Using default time-delay embedding parameters: autocorrelation and cao');
else
    if length(embedParams) ~= 2
        error('Embedding parameters are incorrectly formatted -- need {tau,m}')
    end
end

% ------------------------------------------------------------------------------
%% Embed the signal
% ------------------------------------------------------------------------------
% Convert to embedded signal object for TSTOOL
s = BF_embed(y,embedParams{1},embedParams{2},1);

if ~isa(s,'signal') && isnan(s); % embedding failed
    error('Time-delay embedding for TSTOOL failed')
end

if size(data(s),2) < 3 % embedded with dimension < 3
    % note the 'true' predicted embedding dimension
    mopt = size(data(s),2);
    % embed with dimension m = 3
    s = BF_embed(y,embedParams{1},3,1);
    fprintf(1,'Re-embedded with embedding dimension 3\n');
else
	mopt = size(data(s),2);
end

% ------------------------------------------------------------------------------
%% Run the TSTOOL function:
% ------------------------------------------------------------------------------
% This looks for the dimensions file in the tstoolbox/@signal/dimensions directory
if ~exist(fullfile('tstoolbox','@signal','dimensions'),'file')
    error('Cannot find the code ''dimensions'' from the TSTOOL package. Is it installed and in the Matlab path?');
end
try
    [bc, ~, co] = dimensions(s,nbins);
catch me
    error('Error running TSTOOL code dimensions: %s',me.message);
end

% We now have the scaling of the boxcounting dimension, D0, the information
% dimension D1, and the correlation dimension D2.

% It seems like there's not extra information in the in dimension, D1, by these
% estimates -- focus on the bc, D0 and correlation, D2.
% Can switch this here:
compute_in = 0;

% ------------------------------------------------------------------------------
%% Convert output to vectors
% ------------------------------------------------------------------------------
% calculations for each dimension up to the maximum:
% Seems to be in units of log_2 -- log2, so actually doesn't span a very wide
% range of length scales... -- although I think the maximum length is at 1,
% so I think it might be referring to fractions of the 'attractor size'...

% (1) Boxcounting dimension (BC)
bc_logN = data(bc); % this is log(N(r)) -- number within radius (look for this to scale linearly)
bc_logr = spacing(bc); % this is log(r) -- length scale
bc_logNlogr = bc_logN./(ones(size(bc_logN,2),1)*bc_logr)'; % look for this to be constant

% (2) Information dimension (IN)
if compute_in
    in_logl = data(in); % I think this is log(l(r)), look for this to scale linearly
    in_logr = spacing(in);
    in_logllogr = in_logl./(ones(size(in_logl,2),1)*in_logr)'; % look for this to be constant
end

% (3) Correlation dimension (CO)
co_logC = data(co); % look for this to scale linearly
co_logr = spacing(co);
co_logClogr = co_logC./(ones(size(co_logC,2),1)*co_logr)'; % look for this to be constant


if doPlot
    plot(bc_logr,bc_logNlogr,'o-')
    plot(bc_logr,bc_logN,'o-')
    input('BC')
    if compute_in
        plot(in_logr,in_logllogr,'o-')
        plot(in_logr,in_logl,'o-')
        input('IN')
    end
    plot(co_logr,co_logClogr,'o-')
    plot(co_logr,co_logC,'o-')
    input('CO')
end

% *** We now have to look for scaling regimes in each of these dimensions


% ------------------------------------------------------------------------------
%% Basic statistics on curves
% ------------------------------------------------------------------------------
out = struct;

%-------------------------------------------------------------------------------
%% How do curves change with m?
%-------------------------------------------------------------------------------
% Use SUB_mch

% Box counting dimension:
out = SUB_mch(bc_logr,bc_logN,'bc',out);

% Information dimension:
if compute_in
    out = SUB_mch(bc_logr,bc_logN,'in',out);
end

% Correlation dimension:
out = SUB_mch(co_logr,co_logC,'co',out);

% ------------------------------------------------------------------------------
%% What is the scaling range in r?
% ------------------------------------------------------------------------------
% ... and how good is the fit over this range?
% Use SUB_scr

% Box counting dimension, m = 1
out = SUB_scr(bc_logr,bc_logN(:,1),'scr_bc_m1',out);

% Box counting dimension m = 2
out = SUB_scr(bc_logr,bc_logN(:,2),'scr_bc_m2',out);

% Box counting dimension m = 3
out = SUB_scr(bc_logr,bc_logN(:,3),'scr_bc_m3',out);

% Box counting dimension m = chosen/given
out = SUB_scr(bc_logr,bc_logN(:,mopt),'scr_bc_mopt',out);

if compute_in
    % Information dimension, m = 1
    out = SUB_scr(in_logr,in_logl(:,1),'scr_in_m1',out);

    % Information dimension m = 2
    out = SUB_scr(in_logr,in_logl(:,2),'scr_in_m2',out);

    % Information dimension m = 3
    out = SUB_scr(in_logr,in_logl(:,3),'scr_in_m3',out);

    % Information dimension m = chosen/given
    out = SUB_scr(in_logr,in_logl(:,mopt),'scr_in_mopt',out);
end

% Correlation dimension, m = 1
out = SUB_scr(co_logr,co_logC(:,1),'scr_co_m1',out);

% Correlation dimension m = 2
out = SUB_scr(co_logr,co_logC(:,2),'scr_co_m2',out);

% Correlation dimension m = 3
out = SUB_scr(co_logr,co_logC(:,3),'scr_co_m3',out);

% Correlation dimension m = chosen/given
out = SUB_scr(co_logr,co_logC(:,mopt),'scr_co_mopt',out);

% ------------------------------------------------------------------------------
%% What m gives best fit?
% ------------------------------------------------------------------------------
% Use SUB_bestm

% Box counting dimension
out = SUB_bestm(bc_logr,bc_logN,'bc',out);

if compute_in
    % Information dimension
    out = SUB_bestm(in_logr,in_logl,'in',out);
end

% Correlation dimension
out = SUB_bestm(co_logr,co_logC,'co',out);

% ------------------------------------------------------------------------------
% ------------------------------------------------------------------------------
function out = SUB_mch(logr,logN,prefix,out)
    % looks at how changes with m. Since m will in general be different for each
    % different time series (i.e., if choosing an automatic method for
    % determining the embedding parameters), we have that m is at least
    % 3 here so that we can do statistics on at least these ones...

    % (i) on average the raw means at each m up to m = 3
    out.([prefix,'_meanm1']) = mean(logN(:,1));
    out.([prefix,'_meanm2']) = mean(logN(:,2));
    out.([prefix,'_meanm3']) = mean(logN(:,3));
    out.([prefix,'_meanmmax']) = mean(logN(:,end));

    % (ii) raw minimum at each m up to m = 3
    out.([prefix,'_minm1']) = min(logN(:,1));
    out.([prefix,'_minm2']) = min(logN(:,2));
    out.([prefix,'_minm3']) = min(logN(:,3));
    out.([prefix,'_minmmax']) = min(logN(:,end));

    % (iii) range at each m up to m = 3
    out.([prefix,'_range1']) = range(logN(:,1));
    out.([prefix,'_range2']) = range(logN(:,2));
    out.([prefix,'_range3']) = range(logN(:,3));
    out.([prefix,'_rangemmax']) = range(logN(:,end));

    % (iv) increments with m
    out.([prefix,'_mindiff']) = mean([min(logN(:,2))-min(logN(:,1)),min(logN(:,3))-min(logN(:,2))]);
    out.([prefix,'_meandiff']) = mean([mean(logN(:,2))-mean(logN(:,1)),mean(logN(:,3))-mean(logN(:,2))]);

    % (v) slopes and goodness of fit across whole r range
    [out.([prefix,'_lfitm1']), out.([prefix,'_lfitb1']), out.([prefix,'_lfitmeansqdev1'])] = subsublinfit(logr,logN(:,1)');
    [out.([prefix,'_lfitm2']), out.([prefix,'_lfitb2']), out.([prefix,'_lfitmeansqdev2'])] = subsublinfit(logr,logN(:,2)');
    [out.([prefix,'_lfitm3']), out.([prefix,'_lfitb3']), out.([prefix,'_lfitmeansqdev3'])] = subsublinfit(logr,logN(:,3)');
    [out.([prefix,'_lfitmmax']), out.([prefix,'_lfitbmax']), out.([prefix,'_lfitmeansqdevmax'])] = subsublinfit(logr,logN(:,end)');

    function [m, b, meansqdev] = subsublinfit(x,y)
        p1 = polyfit(x,y,1);
        pfit = p1(1)*x + p1(2);
        res = y-pfit;
        m = p1(1); % gradient
        b = p1(2); % intercept
        meansqdev = mean(res.^2);
    end
end

% ------------------------------------------------------------------------------
function out = SUB_scr(logr,logN,prefix,out)
    % determines the scaling range in r for some m
    % we remove points from either extreme in r until minimize some
    % error measure
    % two dimensional optimization: over starting point and ending
    % point.

    l = length(logr);
    stptr = 1:floor(l/2)-1; % must be in the first half (not necessarily, but for here)
    endptr = ceil(l/2)+1:l; % must be in second half (not necessarily, but for here)
    mybad = zeros(length(stptr),length(endptr));
    for i = 1:length(stptr)
        for j = 1:length(endptr)
            mybad(i,j) = lfitbadness(logr(stptr(i):endptr(j)),logN(stptr(i):endptr(j))');
        end
    end
    [a, b] = find(mybad == min(min(mybad))); % this defines the 'best' scaling range
%         plot(logr,logN,'o-b'); hold on; plot(logr(stptr(a):endptr(b)),logN(stptr(a):endptr(b)),'o-r');
%         hold off
%         disp(['keep from ' num2str(stptr(a)) ' to ' num2str(endptr(b))])


    out.([prefix,'_logrmin']) = logr(stptr(a)); % minimum of scaling range
    out.([prefix,'_logrmax']) = logr(endptr(b)); % maximum of scaling range
    out.([prefix,'_logrrange']) = logr(endptr(b))-logr(stptr(a)); % range of scaling... range
    out.([prefix,'_pgone']) = (stptr(a)-1 + l - endptr(b))/length(logr); % number of points removed in process
                                               				  % of choosing the optimum scaling range

	% Do the optimum fit again
	x = logr(stptr(a):endptr(b));
	y = logN(stptr(a):endptr(b))';
    p = polyfit(x,y,1);
    pfit = p(1)*x+p(2);
    res = pfit-y;
    out.([prefix,'_meanabsres']) = mean(abs(res));
	out.([prefix,'_meansqres']) = mean(res.^2);
    out.([prefix,'_scaling_exp']) = p(1);
	out.([prefix,'_scaling_int']) = p(2);
	out.([prefix,'_minbad']) = min(min(mybad));

    function badness = lfitbadness(x,y)
        gamma = 0.02; % reguralization parameter gamma selected empirically, could be tweaked in future work
        p = polyfit(x,y,1);
        pfit = p(1)*x+p(2);
        res = pfit-y;
        badness = mean(abs(res)) - gamma*length(x); % want to still maximize length(x)
    end
end

% ------------------------------------------------------------------------------
function out = SUB_bestm(logr,logNN,prefix,out)
    % logNN is a matrix... logN is a vector for a given m
	% determines the scaling range in r for some m
    % we remove points from either extreme in r until minimize some
    % error measure
    % two dimensional optimization: over starting point and ending
    % point.

	store_scalingexps = zeros(size(logNN,2),1);
	store_meansqres = zeros(size(logNN,2),1);
	for k = 1:size(logNN,2);
		logN = logNN(:,k); % take this element

        l = length(logr);
        stptr = 1:floor(l/2)-1; % must be in the first half (not necessarily, but for here)
        endptr = ceil(l/2)+1:l; % must be in second half (not necessarily, but for here)
        mybad = zeros(length(stptr),length(endptr));
        for i = 1:length(stptr)
            for j = 1:length(endptr)
                mybad(i,j) = lfitbadness(logr(stptr(i):endptr(j)),logN(stptr(i):endptr(j))');
            end
        end
        [a, b] = find(mybad == min(min(mybad))); % this defines the 'best' scaling range

		% Do the optimum fit again
		x = logr(stptr(a):endptr(b));
		y = logN(stptr(a):endptr(b))';
        p = polyfit(x,y,1);
        pfit = p(1)*x+p(2);
        res = pfit-y;
		% subout.meanabsres = mean(abs(res));
		% subout.meansqres = mean(res.^2);
		% subout.scaling_exp = p(1);
		% subout.scaling_int = p(2);
		% subout.minbad = min(min(mybad));

		store_scalingexps(k) = p(1);
		store_meansqres(k) = mean(res.^2);
	end

	out.([prefix,'_minscalingexp']) = min(store_scalingexps);
	out.([prefix,'_meanscalingexp']) = mean(store_scalingexps);
	out.([prefix,'_maxscalingexp']) = max(store_scalingexps);
	out.([prefix,'_mbestfit']) = find(store_meansqres == min(store_meansqres),1,'first');

    function badness = lfitbadness(x,y)
        gamma = 0.02; % reguralization parameter gamma selected empirically, could be tweaked in future work
        p = polyfit(x,y,1);
        pfit = p(1)*x + p(2);
        res = pfit - y;
        badness = mean(abs(res)) - gamma*length(x); % want to still maximize length(x)
    end
end

end
