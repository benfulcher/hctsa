% PP_Iterate
% 
% Iteratively applies a transformation to the time series and measures how
% various properties of the time series change as the transformation is
% iteratively applied to it.
% 
% INPUTS:
% 
% y, the input time series
% 
% detrndmeth, the detrending method to apply:
%           (i) 'spline' removes a spine fit,
%           (ii) 'diff' takes incremental differences,
%           (iii) 'medianf' applies a median filter,
%           (iv) 'rav' applies a running mean filter,
%           (v) 'resampleup' progressively upsamples the time series,
%           (vi) 'resampledown' progressively downsamples the time series.
% 
% ------------------------------------------------------------------------------
% Copyright (C) 2013,  Ben D. Fulcher <ben.d.fulcher@gmail.com>,
% <http://www.benfulcher.com>
%
% If you use this code for your research, please cite:
% B. D. Fulcher, M. A. Little, N. S. Jones., "Highly comparative time-series
% analysis: the empirical structure of time series and their methods",
% J. Roy. Soc. Interface 10(83) 20130048 (2010). DOI: 10.1098/rsif.2013.0048
%
% This function is free software: you can redistribute it and/or modify it under
% the terms of the GNU General Public License as published by the Free Software
% Foundation, either version 3 of the License, or (at your option) any later
% version.
% 
% This program is distributed in the hope that it will be useful, but WITHOUT
% ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
% FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
% details.
% 
% You should have received a copy of the GNU General Public License along with
% this program.  If not, see <http://www.gnu.org/licenses/>.
% ------------------------------------------------------------------------------

function out = PP_Iterate(y,detrndmeth)
% Ben Fulcher, 10/7/09

%% FOREPLAY
N = length(y); % length of time series

%% Determine the number of times the processing will be progressively performed
% this information is stored in nr
switch detrndmeth
    case 'spline', nr = (1:20);
    case 'diff', nr = (1:5);
    case 'medianf', nr = round(linspace(1,N/25,25));
    case 'rav', nr = round(linspace(1,N/25,25));
    case 'resampleup', nr = (1:20);
    case 'resampledown', nr = (1:20);
    otherwise
        error('Unknown detrending method ''%s''',detrndmeth);
end

%% Do the progessive processing with running statistical evaluation
outmat = zeros(length(nr),10);
for q = 1:length(nr)
    n = nr(q);
    switch detrndmeth
        case 'spline' % Spline detrend
            %% Check that a Curve-Fitting Toolbox license is available:
            BF_CheckToolbox('curve_fitting_toolbox')
            
            nknots = n; % progressively make more knots
            intp = 4; % cubic interpolants
            spline = spap2(nknots,intp,[1:N]',y); % just a single middle knot with cubic interpolants
            y_spl = fnval(spline,1:N); % evaluate at the 1:N time intervals
            y_d = y - y_spl';
            
        case 'diff' % Differencing
            ndiffs = n; % progressively difference
            y_d = diff(y,ndiffs);
            
        case 'medianf' % Median Filter; n is the order of filtering
            y_d = medfilt1(y,n);
            
        case 'rav' % Running Average; n is the window size
            y_d = filter(ones(1,n)/n,1,y);
            
        case 'resampleup' % upsample
            y_d = resample(y,n,1);
            
        case 'resampledown' % downsample
            y_d = resample(y,1,n);
    end
    outmat(q,:) = doyourcalcthing(y,y_d);
end
% out1 = outmat;


%% DETERMINE THE STATISTICS ON EACH TEST
% four statistics on each test

stats = zeros(10,4);
for t = 1:10;
    if any(~isfinite(stats(t,:))),
        fprintf(1,'This is a bad statistic\n')
    end
    stats(t,:) = doyourtestthing(outmat(:,t)); 
end

%% WRITE OUTPUT
% 1) StatAv_5
out.statav5_trend = stats(1,1);
out.statav5_jump = stats(1,2);
out.statav5_lin = stats(1,3);
out.statav5_exp = stats(1,4);

% 2) Sliding Window Mean
out.swms5_2_trend = stats(2,1);
out.swms5_2_jump = stats(2,2);
out.swms5_2_lin = stats(2,3);
out.swms5_2_exp = stats(2,4);

% 3) Sliding Window std
out.swss5_2_trend = stats(3,1);
out.swss5_2_jump = stats(3,2);
out.swss5_2_lin = stats(3,3);
out.swss5_2_exp = stats(3,4);

% 4) Gaussian kernel density fit, rmse
out.gauss1_kd_trend = stats(4,1);
out.gauss1_kd_jump = stats(4,2);
out.gauss1_kd_lin = stats(4,3);
out.gauss1_kd_exp = stats(4,4);

% 5) Gaussian 10-bin histogram fit, rmse
out.gauss1_h10_trend = stats(5,1);
out.gauss1_h10_jump = stats(5,2);
out.gauss1_h10_lin = stats(5,3);
out.gauss1_h10_exp = stats(5,4);

% 6) Compare normal fit
out.norm_kscomp_trend = stats(6,1);
out.norm_kscomp_jump = stats(6,2);
out.norm_kscomp_lin = stats(6,3);
out.norm_kscomp_exp = stats(6,4);

% 7) Outliers: ben test
out.ol_trend = stats(7,1);
out.ol_jump = stats(7,2);
out.ol_lin = stats(7,3);
out.ol_exp = stats(7,4);

% 8) Cross correlation to original signal (-1)
out.xcn1_trend = stats(8,1);
out.xcn1_jump = stats(8,2);
out.xcn1_lin = stats(8,3);
out.xcn1_exp = stats(8,4);

% 9) Cross correlation to original signal (+1)
out.xc1_trend = stats(9,1);
out.xc1_jump = stats(9,2);
out.xc1_lin = stats(9,3);
out.xc1_exp = stats(9,4);

% 10) Norm of differences to original and processed signals
out.normdiff_trend = stats(10,1);
out.normdiff_jump = stats(10,2);
out.normdiff_lin = stats(10,3);
out.normdiff_exp = stats(10,4);


    %% TESTS
    function f = doyourcalcthing(y,y_d)
        y = BF_zscore(y); y_d = BF_zscore(y_d);
        f = zeros(10,1); % vector of features to output
        % 1) Stationarity
        % (a) StatAv
        f(1) = SY_StatAv(y_d,'seg',5);
        % (b) Sliding window mean
        f(2) = SY_SlidingWindow(y_d,'mean','std',5,2);
        % (c) Sliding window std
        f(3) = SY_SlidingWindow(y_d,'std','std',5,2)/SY_SlidingWindow(y,'std','std',5,2);
        
        % 2) Gaussianity
        %   (a) kernel density estimation method
        me1 = DN_SimpleFit(y_d,'gauss1',0); % kernel density fit to 1-peak gaussian
        
        if ~isstruct(me1) && isnan(me1)
            f(4) = NaN;
        else
            f(4) = me1.rmse;
        end
        
        %   (b) histogram 10 bins
        me1 = DN_SimpleFit(y_d,'gauss1',10); % 10-bin histogram fit to 1-peak gaussian
        if ~isstruct(me1) && isnan(me1)
            f(5) = NaN;
        else
            f(5) = me1.rmse;
        end
        
        %   (c) compare distribution to fitted normal distribution
        me1 = DN_CompareKSFit(y_d,'norm');
        if ~isstruct(me1) && isnan(me1)
            f(6) = NaN;
        else
            f(6) = me1.adiff;
        end
        
        % 3) Outliers
        f(7) = DN_OutlierTest(y_d,5,'mean');
        
        % Cross Correlation to original signal
        if length(y) == length(y_d)
            xc = xcorr(y,y_d,1,'coeff');
            f(8) = xc(1);
            f(9) = xc(3);

            % Norm of differences between original and randomized signals
            f(10) = norm(y-y_d)/length(y);
        else
           f(8:10) = NaN; % like for differencing where lose some points
        end
    end

    function g = doyourtestthing(f)
        if ~all(isfinite(f))
            g = NaN*ones(4,1); return
        else
            g = zeros(4,1);
        end
        f = zscore(f);
        % for each return a set of simple little tests:
        % (1) Is it increasing/decreasing?
        %       do this by the sum of differences; could take sign?
        g(1) = sum(diff(f));
        
        % (2) Is there a jump anywhere? If so, where?
        % this is done crudely -- needs to be a sudden jump in a single
        % step
        % find running differnce between mean before and mean after
        mbfatd = zeros(length(f),4);
        for j = 1:length(f)
            mbfatd(j,1) = mean(f(1:j));
            mbfatd(j,2) = mean(f(j:end));
            mbfatd(j,3) = std(f(1:j))/sqrt(length(1:j));
            mbfatd(j,4) = std(f(j:end))/sqrt(length(j:length(f)));
        end

        % Pick the maximum difference
        [c, ind] = max(abs(mbfatd(:,1)-mbfatd(:,2)));

        % t-statistic at this point
        tstat = abs((mbfatd(ind,1)-mbfatd(ind,2))/sqrt(mbfatd(ind,3)^2+mbfatd(ind,4)^2));
        if tstat > 15; %2 && ind < length(f)-2 && mbfatd(ind-1) < mbfatd(ind) && mbfatd(ind+1) < mbfatd(ind)
            g(2) = ind;
        else g(2) = NaN;
        end
        
        % (3) is it linear?
        try
            [cfun, gof] = fit((1:length(f))',f,'poly1');
        catch emsg
            if ~(strcmp(emsg.message,'Inf computed by model function.') || strcmp(emsg.message,'NaN computed by model function.'))
                return
            end
        end
        
        if gof.rsquare > 0.9
            g(3) = cfun.p1; % returns the gradient of the line
        else
            g(3) = NaN;
        end
        
        % (4) is it exponential?
        
        if g(1) > 0 % increasing, fit saturating exponential
            stpt = [-f(end), -0.1, f(end)];
        else % decreasing, fit decaying exponential
            stpt = [f(end), -0.1, f(end)];
        end
        
        fopt = fitoptions('Method','NonlinearLeastSquares','StartPoint',stpt);
        ftyp = fittype('a*exp(b*x)+c','options',fopt);
        try 
            [c, gof] = fit((1:length(f))',f,ftyp);
            if gof.rsquare > 0.9
                g(4) = c.b; % return the decay coefficient
            else
                g(4) = NaN;
            end
        catch emsg
            if (strcmp(emsg.message,'Inf computed by model function.') || strcmp(emsg.message,'NaN computed by model function.'))
                g(4) = NaN; % error fitting -- return a NaN
            else % unable to run 'fit' command -- fatal error
                return
            end
        end

        
    end
end
     