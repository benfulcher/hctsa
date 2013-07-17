function out = TSTL_dimensions(y,nbins,embedparams)
% Uses TSTOOL code dimensions
% Computes boxcounting, information, and correlation dimensions of the time
% series embedded using embedparams.
% y: column vector of time series data
% nbins: maximum number of partitions per axis.
% embedparams: embedding parameters to feed BF_embed.m for embedding the
% signal in the form {tau,m}

% Ben Fulcher November 2009

%% Preliminaries
N = length(y); % length of time series

% (1) Maximum number of bins, nbins
if nargin < 2 || isempty(nbins)
    nbins = 50; % 50 points
    fprintf(1,'Using a default of 50 bins per axis\n');
end

% (2) Set embedding parameters to defaults
if nargin < 3 || isempty(embedparams)
    embedparams = {'ac','cao'};
    fprintf(1,'Using default time-delay embedding parameters: autocorrelation and cao')
else
    if length(embedparams) ~= 2
        error('Embedding parameters are incorrectly formatted -- need {tau,m}')
    end
end

%% Embed the signal
% convert to embedded signal object for TSTOOL
s = BF_embed(y,embedparams{1},embedparams{2},1);

if ~strcmp(class(s),'signal') && isnan(s); % embedding failed
    error('Time-delay embedding for TSTOOL failed')
end

if size(data(s),2) < 3 % embedded with dimension < 3
    % note the 'true' predicted embedding dimension
    mopt = size(data(s),2);
    % embed with dimension m = 3
    s = BF_embed(y,embedparams{1},3,1);
    fprintf(1,'Re-embedded with embedding dimension 3\n')
else
	mopt = size(data(s),2);
end

%% Run

[bc, in, co] = dimensions(s,nbins);

% we now have the scaling of the boxcounting dimension, D0, the information
% dimension D1, and the correlation dimension D2.

%% Convert output to vectors
% calculations for each dimension up to the maximum:
% Seems to be in units of log_2 -- log2, so actually doesn't span a very wide
% range of length scales... -- although I think the maximum length is at 1,
% so I think it might be referring to fractions of the 'attractor size'...

% (1) Boxcounting dimension (BC)
bc_logN = data(bc); % this is log(N(r)) -- number within radius (look for this to scale linearly)
bc_logr = spacing(bc); % this is log(r) -- length scale
bc_logNlogr = bc_logN./(ones(size(bc_logN,2),1)*bc_logr)'; % look for this to be constant

% (2) Information dimension (IN)
in_logl = data(in); % I think this is log(l(r)), look for this to scale linearly
in_logr = spacing(in);
in_logllogr = in_logl./(ones(size(in_logl,2),1)*in_logr)'; % look for this to be constant

% (3) Correlation dimension (CO)
co_logC = data(co); % look for this to scale linearly
co_logr = spacing(co);
co_logClogr = co_logC./(ones(size(co_logC,2),1)*co_logr)'; % look for this to be constant


% plot(bc_logr,bc_logNlogr,'o-')
% plot(bc_logr,bc_logN,'o-')
% input('BC')
% plot(in_logr,in_logllogr,'o-')
% plot(in_logr,in_logl,'o-')
% input('IN')
% plot(co_logr,co_logClogr,'o-')
% plot(co_logr,co_logC,'o-')
% input('CO')


% we now have to look for scaling regimes in each of these dimensions
%% Basic statistics on curves

%% How do curves change with m?
% Box counting dimension
outmch_bc = SUB_mch(bc_logr,bc_logN);
out.bc_meanm1 = outmch_bc.meanm1;
out.bc_meanm2 = outmch_bc.meanm2;
out.bc_meanm3 = outmch_bc.meanm3;
out.bc_meanmmax = outmch_bc.meanmmax;
out.bc_minm1 = outmch_bc.minm1;
out.bc_minm2 = outmch_bc.minm2;
out.bc_minm3 = outmch_bc.minm3;
out.bc_minmmax = outmch_bc.minmmax;
out.bc_range1 = outmch_bc.range1;
out.bc_range2 = outmch_bc.range2;
out.bc_range3 = outmch_bc.range3;
out.bc_rangemmax = outmch_bc.rangemmax;
out.bc_mindiff = outmch_bc.mindiff;
out.bc_meandiff = outmch_bc.meandiff;
out.bc_lfitm1 = outmch_bc.lfitm1;
out.bc_lfitm2 = outmch_bc.lfitm2;
out.bc_lfitm3 = outmch_bc.lfitm3;
out.bc_lfitmmax = outmch_bc.lfitmmax;
out.bc_lfitb1 = outmch_bc.lfitb1;
out.bc_lfitb2 = outmch_bc.lfitb2;
out.bc_lfitb3 = outmch_bc.lfitb3;
out.bc_lfitbmax = outmch_bc.lfitbmax;
out.bc_lfitmeansqdev1 = outmch_bc.lfitmeansqdev1;
out.bc_lfitmeansqdev2 = outmch_bc.lfitmeansqdev2;
out.bc_lfitmeansqdev3 = outmch_bc.lfitmeansqdev3;
out.bc_lfitmeansqdevmax = outmch_bc.lfitmeansqdevmax;

% Information dimension
outmch_in = SUB_mch(in_logr,in_logl);
out.in_meanm1 = outmch_in.meanm1;
out.in_meanm2 = outmch_in.meanm2;
out.in_meanm3 = outmch_in.meanm3;
out.in_meanmmax = outmch_in.meanmmax;
out.in_minm1 = outmch_in.minm1;
out.in_minm2 = outmch_in.minm2;
out.in_minm3 = outmch_in.minm3;
out.in_minmmax = outmch_in.minmmax;
out.in_range1 = outmch_in.range1;
out.in_range2 = outmch_in.range2;
out.in_range3 = outmch_in.range3;
out.in_rangemmax = outmch_in.rangemmax;
out.in_mindiff = outmch_in.mindiff;
out.in_meandiff = outmch_in.meandiff;
out.in_lfitm1 = outmch_in.lfitm1;
out.in_lfitm2 = outmch_in.lfitm2;
out.in_lfitm3 = outmch_in.lfitm3;
out.in_lfitmmax = outmch_in.lfitmmax;
out.in_lfitb1 = outmch_in.lfitb1;
out.in_lfitb2 = outmch_in.lfitb2;
out.in_lfitb3 = outmch_in.lfitb3;
out.in_lfitbmax = outmch_in.lfitbmax;
out.in_lfitmeansqdev1 = outmch_in.lfitmeansqdev1;
out.in_lfitmeansqdev2 = outmch_in.lfitmeansqdev2;
out.in_lfitmeansqdev3 = outmch_in.lfitmeansqdev3;
out.in_lfitmeansqdevmax = outmch_in.lfitmeansqdevmax;

% Correlation dimension
outmch_co = SUB_mch(co_logr,co_logC);
out.co_meanm1 = outmch_co.meanm1;
out.co_meanm2 = outmch_co.meanm2;
out.co_meanm3 = outmch_co.meanm3;
out.co_meanmmax = outmch_co.meanmmax;
out.co_minm1 = outmch_co.minm1;
out.co_minm2 = outmch_co.minm2;
out.co_minm3 = outmch_co.minm3;
out.co_minmmax = outmch_co.minmmax;
out.co_range1 = outmch_co.range1;
out.co_range2 = outmch_co.range2;
out.co_range3 = outmch_co.range3;
out.co_rangemmax = outmch_co.rangemmax;
out.co_mindiff = outmch_co.mindiff;
out.co_meandiff = outmch_co.meandiff;
out.co_lfitm1 = outmch_co.lfitm1;
out.co_lfitm2 = outmch_co.lfitm2;
out.co_lfitm3 = outmch_co.lfitm3;
out.co_lfitmmax = outmch_co.lfitmmax;
out.co_lfitb1 = outmch_co.lfitb1;
out.co_lfitb2 = outmch_co.lfitb2;
out.co_lfitb3 = outmch_co.lfitb3;
out.co_lfitbmax = outmch_co.lfitbmax;
out.co_lfitmeansqdev1 = outmch_co.lfitmeansqdev1;
out.co_lfitmeansqdev2 = outmch_co.lfitmeansqdev2;
out.co_lfitmeansqdev3 = outmch_co.lfitmeansqdev3;
out.co_lfitmeansqdevmax = outmch_co.lfitmeansqdevmax;


%% What is the scaling range in r?
% ... and how good is the fit over this range?


% Box counting dimension, m = 1
outscr_bc_m1 = SUB_scr(bc_logr,bc_logN(:,1));
out.scr_bc_m1_minbad = outscr_bc_m1.minbad;
out.scr_bc_m1_logrmin = outscr_bc_m1.logrmin;
out.scr_bc_m1_logrmax = outscr_bc_m1.logrmax;
out.scr_bc_m1_logrrange = outscr_bc_m1.logrrange;
out.scr_bc_m1_pgone = outscr_bc_m1.pgone;
out.scr_bc_m1_meanabsres = outscr_bc_m1.meanabsres;
out.scr_bc_m1_meansqres = outscr_bc_m1.meansqres;
out.scr_bc_m1_scaling_exp = outscr_bc_m1.scaling_exp;
out.scr_bc_m1_scaling_int = outscr_bc_m1.scaling_int;
out.scr_bc_m1_minbad = outscr_bc_m1.minbad;
% Box counting dimension m = 2
outscr_bc_m2 = SUB_scr(bc_logr,bc_logN(:,2));
out.scr_bc_m2_minbad = outscr_bc_m2.minbad;
out.scr_bc_m2_logrmin = outscr_bc_m2.logrmin;
out.scr_bc_m2_logrmax = outscr_bc_m2.logrmax;
out.scr_bc_m2_logrrange = outscr_bc_m2.logrrange;
out.scr_bc_m2_pgone = outscr_bc_m2.pgone;
out.scr_bc_m2_meanabsres = outscr_bc_m2.meanabsres;
out.scr_bc_m2_meansqres = outscr_bc_m2.meansqres;
out.scr_bc_m2_scaling_exp = outscr_bc_m2.scaling_exp;
out.scr_bc_m2_scaling_int = outscr_bc_m2.scaling_int;
out.scr_bc_m2_minbad = outscr_bc_m2.minbad;
% Box counting dimension m = 3
outscr_bc_m3 = SUB_scr(bc_logr,bc_logN(:,3));
out.scr_bc_m3_minbad = outscr_bc_m3.minbad;
out.scr_bc_m3_logrmin = outscr_bc_m3.logrmin;
out.scr_bc_m3_logrmax = outscr_bc_m3.logrmax;
out.scr_bc_m3_logrrange = outscr_bc_m3.logrrange;
out.scr_bc_m3_pgone = outscr_bc_m3.pgone;
out.scr_bc_m3_meanabsres = outscr_bc_m3.meanabsres;
out.scr_bc_m3_meansqres = outscr_bc_m3.meansqres;
out.scr_bc_m3_scaling_exp = outscr_bc_m3.scaling_exp;
out.scr_bc_m3_scaling_int = outscr_bc_m3.scaling_int;
out.scr_bc_m3_minbad = outscr_bc_m3.minbad;
% Box counting dimension m = chosen/given
outscr_bc_mopt = SUB_scr(bc_logr,bc_logN(:,mopt));
out.scr_bc_mopt_minbad = outscr_bc_mopt.minbad;
out.scr_bc_mopt_logrmin = outscr_bc_mopt.logrmin;
out.scr_bc_mopt_logrmax = outscr_bc_mopt.logrmax;
out.scr_bc_mopt_logrrange = outscr_bc_mopt.logrrange;
out.scr_bc_mopt_pgone = outscr_bc_mopt.pgone;
out.scr_bc_mopt_meanabsres = outscr_bc_mopt.meanabsres;
out.scr_bc_mopt_meansqres = outscr_bc_mopt.meansqres;
out.scr_bc_mopt_scaling_exp = outscr_bc_mopt.scaling_exp;
out.scr_bc_mopt_scaling_int = outscr_bc_mopt.scaling_int;
out.scr_bc_mopt_minbad = outscr_bc_mopt.minbad;


% Information dimension, m = 1
outscr_in_m1 = SUB_scr(in_logr,in_logl(:,1));
out.scr_in_m1_minbad = outscr_in_m1.minbad;
out.scr_in_m1_logrmin = outscr_in_m1.logrmin;
out.scr_in_m1_logrmax = outscr_in_m1.logrmax;
out.scr_in_m1_logrrange = outscr_in_m1.logrrange;
out.scr_in_m1_pgone = outscr_in_m1.pgone;
out.scr_in_m1_meanabsres = outscr_in_m1.meanabsres;
out.scr_in_m1_meansqres = outscr_in_m1.meansqres;
out.scr_in_m1_scaling_exp = outscr_in_m1.scaling_exp;
out.scr_in_m1_scaling_int = outscr_in_m1.scaling_int;
out.scr_in_m1_minbad = outscr_in_m1.minbad;
% Information dimension m = 2
outscr_in_m2 = SUB_scr(in_logr,in_logl(:,2));
out.scr_in_m2_minbad = outscr_in_m2.minbad;
out.scr_in_m2_logrmin = outscr_in_m2.logrmin;
out.scr_in_m2_logrmax = outscr_in_m2.logrmax;
out.scr_in_m2_logrrange = outscr_in_m2.logrrange;
out.scr_in_m2_pgone = outscr_in_m2.pgone;
out.scr_in_m2_meanabsres = outscr_in_m2.meanabsres;
out.scr_in_m2_meansqres = outscr_in_m2.meansqres;
out.scr_in_m2_scaling_exp = outscr_in_m2.scaling_exp;
out.scr_in_m2_scaling_int = outscr_in_m2.scaling_int;
out.scr_in_m2_minbad = outscr_in_m2.minbad;
% Information dimension m = 3
outscr_in_m3 = SUB_scr(in_logr,in_logl(:,3));
out.scr_in_m3_minbad = outscr_in_m3.minbad;
out.scr_in_m3_logrmin = outscr_in_m3.logrmin;
out.scr_in_m3_logrmax = outscr_in_m3.logrmax;
out.scr_in_m3_logrrange = outscr_in_m3.logrrange;
out.scr_in_m3_pgone = outscr_in_m3.pgone;
out.scr_in_m3_meanabsres = outscr_in_m3.meanabsres;
out.scr_in_m3_meansqres = outscr_in_m3.meansqres;
out.scr_in_m3_scaling_exp = outscr_in_m3.scaling_exp;
out.scr_in_m3_scaling_int = outscr_in_m3.scaling_int;
out.scr_in_m3_minbad = outscr_in_m3.minbad;
% Information dimension m = chosen/given
outscr_in_mopt = SUB_scr(in_logr,in_logl(:,mopt));
out.scr_in_mopt_minbad = outscr_in_mopt.minbad;
out.scr_in_mopt_logrmin = outscr_in_mopt.logrmin;
out.scr_in_mopt_logrmax = outscr_in_mopt.logrmax;
out.scr_in_mopt_logrrange = outscr_in_mopt.logrrange;
out.scr_in_mopt_pgone = outscr_in_mopt.pgone;
out.scr_in_mopt_meanabsres = outscr_in_mopt.meanabsres;
out.scr_in_mopt_meansqres = outscr_in_mopt.meansqres;
out.scr_in_mopt_scaling_exp = outscr_in_mopt.scaling_exp;
out.scr_in_mopt_scaling_int = outscr_in_mopt.scaling_int;
out.scr_in_mopt_minbad = outscr_in_mopt.minbad;



% coformation dimension, m = 1
outscr_co_m1 = SUB_scr(co_logr,co_logC(:,1));
out.scr_co_m1_minbad = outscr_co_m1.minbad;
out.scr_co_m1_logrmin = outscr_co_m1.logrmin;
out.scr_co_m1_logrmax = outscr_co_m1.logrmax;
out.scr_co_m1_logrrange = outscr_co_m1.logrrange;
out.scr_co_m1_pgone = outscr_co_m1.pgone;
out.scr_co_m1_meanabsres = outscr_co_m1.meanabsres;
out.scr_co_m1_meansqres = outscr_co_m1.meansqres;
out.scr_co_m1_scaling_exp = outscr_co_m1.scaling_exp;
out.scr_co_m1_scaling_int = outscr_co_m1.scaling_int;
out.scr_co_m1_minbad = outscr_co_m1.minbad;
% coformation dimension m = 2
outscr_co_m2 = SUB_scr(co_logr,co_logC(:,2));
out.scr_co_m2_minbad = outscr_co_m2.minbad;
out.scr_co_m2_logrmin = outscr_co_m2.logrmin;
out.scr_co_m2_logrmax = outscr_co_m2.logrmax;
out.scr_co_m2_logrrange = outscr_co_m2.logrrange;
out.scr_co_m2_pgone = outscr_co_m2.pgone;
out.scr_co_m2_meanabsres = outscr_co_m2.meanabsres;
out.scr_co_m2_meansqres = outscr_co_m2.meansqres;
out.scr_co_m2_scaling_exp = outscr_co_m2.scaling_exp;
out.scr_co_m2_scaling_int = outscr_co_m2.scaling_int;
out.scr_co_m2_minbad = outscr_co_m2.minbad;
% coformation dimension m = 3
outscr_co_m3 = SUB_scr(co_logr,co_logC(:,3));
out.scr_co_m3_minbad = outscr_co_m3.minbad;
out.scr_co_m3_logrmin = outscr_co_m3.logrmin;
out.scr_co_m3_logrmax = outscr_co_m3.logrmax;
out.scr_co_m3_logrrange = outscr_co_m3.logrrange;
out.scr_co_m3_pgone = outscr_co_m3.pgone;
out.scr_co_m3_meanabsres = outscr_co_m3.meanabsres;
out.scr_co_m3_meansqres = outscr_co_m3.meansqres;
out.scr_co_m3_scaling_exp = outscr_co_m3.scaling_exp;
out.scr_co_m3_scaling_int = outscr_co_m3.scaling_int;
out.scr_co_m3_minbad = outscr_co_m3.minbad;
% coformation dimension m = chosen/given
outscr_co_mopt = SUB_scr(co_logr,co_logC(:,mopt));
out.scr_co_mopt_minbad = outscr_co_mopt.minbad;
out.scr_co_mopt_logrmin = outscr_co_mopt.logrmin;
out.scr_co_mopt_logrmax = outscr_co_mopt.logrmax;
out.scr_co_mopt_logrrange = outscr_co_mopt.logrrange;
out.scr_co_mopt_pgone = outscr_co_mopt.pgone;
out.scr_co_mopt_meanabsres = outscr_co_mopt.meanabsres;
out.scr_co_mopt_meansqres = outscr_co_mopt.meansqres;
out.scr_co_mopt_scaling_exp = outscr_co_mopt.scaling_exp;
out.scr_co_mopt_scaling_int = outscr_co_mopt.scaling_int;
out.scr_co_mopt_minbad = outscr_co_mopt.minbad;


%% What m gives best fit?
% box counting dimension
bestm_bc = SUB_bestm(bc_logr,bc_logN);
out.bc_minscalingexp = bestm_bc.minscalingexp;
out.bc_maxscalingexp = bestm_bc.maxscalingexp;
out.bc_meanscalingexp = bestm_bc.meanscalingexp;
out.bc_mbestfit = bestm_bc.mbestfit;

% information dimension
bestm_in = SUB_bestm(in_logr,in_logl);
out.in_minscalingexp = bestm_in.minscalingexp;
out.in_maxscalingexp = bestm_in.maxscalingexp;
out.in_meanscalingexp = bestm_in.meanscalingexp;
out.in_mbestfit = bestm_in.mbestfit;

% correlation dimension
bestm_co = SUB_bestm(co_logr,co_logC);
out.co_minscalingexp = bestm_co.minscalingexp;
out.co_maxscalingexp = bestm_co.maxscalingexp;
out.co_meanscalingexp = bestm_co.meanscalingexp;
out.co_mbestfit = bestm_co.mbestfit;


    function subout = SUB_mch(logr,logN)
        % looks at how changes with m. Since m will in general be different for each
        % different time series (i.e., if choosing an automatic method for
        % determining the embedding parameters), we have that m is at least
        % 3 here so that we can do statistics on at least these ones...
        
        % (i) on average the raw means at each m up to m = 3
        subout.meanm1 = mean(logN(:,1));
        subout.meanm2 = mean(logN(:,2));
        subout.meanm3 = mean(logN(:,3));
        subout.meanmmax = mean(logN(:,end));
        
        % (ii) raw minimum at each m up to m = 3
        subout.minm1 = min(logN(:,1));
        subout.minm2 = min(logN(:,2));
        subout.minm3 = min(logN(:,3));
        subout.minmmax = min(logN(:,end));
        
        % (iii) range at each m up to m = 3
        subout.range1 = range(logN(:,1));
        subout.range2 = range(logN(:,2));
        subout.range3 = range(logN(:,3));
        subout.rangemmax = range(logN(:,end));
        
        % (iv) increments with m
        subout.mindiff = mean([min(logN(:,2))-min(logN(:,1)),min(logN(:,3))-min(logN(:,2))]);
        subout.meandiff = mean([mean(logN(:,2))-mean(logN(:,1)),mean(logN(:,3))-mean(logN(:,2))]);
        
        % (v) slopes and goodness of fit across whole r range
        [subout.lfitm1 subout.lfitb1 subout.lfitmeansqdev1] = subsublinfit(logr,logN(:,1)');
        [subout.lfitm2 subout.lfitb2 subout.lfitmeansqdev2] = subsublinfit(logr,logN(:,2)');
        [subout.lfitm3 subout.lfitb3 subout.lfitmeansqdev3] = subsublinfit(logr,logN(:,3)');
        [subout.lfitmmax subout.lfitbmax subout.lfitmeansqdevmax] = subsublinfit(logr,logN(:,end)');
        
        
        function [m, b, meansqdev] = subsublinfit(x,y)
            p1 = polyfit(x,y,1);
            pfit = p1(1)*x + p1(2);
            res = y-pfit;
            m = p1(1); % gradient
            b = p1(2); % intercept
            meansqdev = mean(res.^2);
        end
        
        
    end

    function subout = SUB_scr(logr,logN)
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


        subout.logrmin = logr(stptr(a)); % minimum of scaling range
        subout.logrmax = logr(endptr(b)); % maximum of scaling range
        subout.logrrange = logr(endptr(b))-logr(stptr(a)); % range of scaling... range
        subout.pgone = (stptr(a)-1 + l - endptr(b))/length(logr); % number of points removed in process 
                                                   				  % of choosing the optimum scaling range

		% Do the optimum fit again
		x = logr(stptr(a):endptr(b));
		y = logN(stptr(a):endptr(b))';
        p = polyfit(x,y,1);
        pfit = p(1)*x+p(2);
        res = pfit-y;
        subout.meanabsres = mean(abs(res));
		subout.meansqres = mean(res.^2);
        subout.scaling_exp = p(1);
		subout.scaling_int = p(2);
		subout.minbad = min(min(mybad));
        
        function badness = lfitbadness(x,y)
            gamma = 0.02; % reguralization parameter gamma selected empirically, could be tweaked in future work
            p = polyfit(x,y,1);
            pfit = p(1)*x+p(2);
            res = pfit-y;
            badness = mean(abs(res)) - gamma*length(x); % want to still maximize length(x)
            
        end
    end

    function subout = SUB_bestm(logr,logNN)
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

		subout.minscalingexp = min(store_scalingexps);
		subout.meanscalingexp = mean(store_scalingexps);
		subout.maxscalingexp = max(store_scalingexps);
		subout.mbestfit = find(store_meansqres == min(store_meansqres),1,'first');
		
        function badness = lfitbadness(x,y)
            gamma = 0.02; % reguralization parameter gamma selected empirically, could be tweaked in future work
            p = polyfit(x,y,1);
            pfit = p(1)*x + p(2);
            res = pfit - y;
            badness = mean(abs(res)) - gamma*length(x); % want to still maximize length(x)

        end

    end

end