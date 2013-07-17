function out = CO_acfshape(y)
% This function analyzes the shape of the autocorrelation function
% Ben Fulcher 2009

% Only look up to when two consecutive values are under the threshold for
% significance:

N = length(y); % length of the time series
th = 2/sqrt(N);
acf = zeros(N,1);

% Calculate the autocorrelation function, up to a maximum lag, length of
% time series (hopefully it's cropped by then)
for i = 1:N
    acf(i) = CO_autocorr(y,i-1); % *** NOTE THIS! *** acf vector indicies are not lags
    if i > 2 && abs(acf(i)) < th && abs(acf(i-1)) < th
       acf = acf(1:i-2);
       break
    end
end

Nac = length(acf);

out.Nac = length(acf); % the distance the acf lasts until significance is 'drowned out' (by my definition)

% Count peaks
dacf = diff(acf);
ddacf = diff(dacf);
extrr = BF_sgnchange(dacf,1);
sdsp = ddacf(extrr);
maxr = extrr(sdsp < 0);
minr = extrr(sdsp > 0);
nmaxr = length(maxr);
nminr = length(minr);

% Number of local minima
out.nminima = sum(sdsp > 0);
out.meanminima = mean(sdsp(sdsp > 0));

% Proportion of local maxima
out.nmaxima = sum(sdsp < 0);
out.meanmaxima = abs(mean(sdsp(sdsp < 0))); % must be negative: make it positive

% Proportion of extrema
out.nextrema = length(sdsp);
out.pextrema = length(sdsp)/Nac;

% Mean of the ACF
out.meanacf = mean(acf);
out.meanabsacf = mean(abs(acf));

% Correlations between extrema
if nmaxr > 4 % need at least 5 points to do this
    out.maximaspread = std(diff(maxr)); % spread of inter-maxima intervals
    out.ac1maxima = CO_autocorr(acf(maxr),1);
else % less than 5 points, return NaNs:
    out.maximaspread = NaN;
    out.ac1maxima = NaN;
end
if nminr > 4 % need at least 5 points to do this
    out.minimaspread = std(diff(minr)); % spread of inter-minima intervals
    out.ac1minima = CO_autocorr(acf(minr),1);
else % less than 5 points, return NaNs:
    out.minimaspread = NaN;
    out.ac1minima = NaN;
end

% Autocorrelation of the ACF
out.ac1 = CO_autocorr(acf,1);
out.ac2 = CO_autocorr(acf,2);
out.ac3 = CO_autocorr(acf,3);
out.actau = CO_autocorr(acf,CO_fzcac(acf));


if Nac > 3 % Need at least four points to fit exponential
    
    %% Fit exponential decay to absolute ACF:
    s = fitoptions('Method','NonlinearLeastSquares','StartPoint',[1, -0.5]);
    f = fittype('a*exp(b*x)','options',s);
    b = 1;
    try
        [c, gof] = fit((1:Nac)',abs(acf),f);
    catch
        b = 0;
    end
    if b == 1
        out.fexpabsacf_a = c.a;
        out.fexpabsacf_b = c.b; % this is important
        out.fexpabsacf_r2 = gof.rsquare; % this is more important!
        out.fexpabsacf_adjr2 = gof.adjrsquare;
        out.fexpabsacf_rmse = gof.rmse;
    
        expfit = c.a*exp(c.b*[1:Nac]');
        res = abs(acf)-expfit;
        out.fexpabsacf_varres = var(res);
    else % fit failed -- return NaNs
        out.fexpabsacf_a = NaN;
        out.fexpabsacf_b = NaN;
        out.fexpabsacf_r2 = NaN;
        out.fexpabsacf_adjr2 = NaN;
        out.fexpabsacf_rmse = NaN;
        out.fexpabsacf_varres = NaN;
    end
    
    %% fit linear to local maxima
    s = fitoptions('Method','NonlinearLeastSquares','StartPoint',[-0.1 1]);
    f = fittype('a*x+b','options',s);
%     plot(maxr,acf(maxr),'ok');
    b = 1;
    try [c, gof] = fit(maxr,acf(maxr),f);
    catch
        b = 0;
    end
    if b == 1; % Fit was successful
        out.flinlmxacf_a = c.a;
        out.flinlmxacf_b = c.b;
        out.flinlmxacf_r2 = gof.rsquare;
        out.flinlmxacf_adjr2 = gof.adjrsquare;
        out.flinlmxacf_rmse = gof.rmse;
    else % Fit failed -- return NaNs
        out.flinlmxacf_a = NaN;
        out.flinlmxacf_b = NaN;
        out.flinlmxacf_r2 = NaN;
        out.flinlmxacf_adjr2 = NaN;
        out.flinlmxacf_rmse = NaN;
    end
else
    out.fexpabsacf_a = NaN;
    out.fexpabsacf_b = NaN;
    out.fexpabsacf_r2 = NaN;
    out.fexpabsacf_adjr2 = NaN;
    out.fexpabsacf_rmse = NaN;
    out.fexpabsacf_varres = NaN;
    out.flinlmxacf_a = NaN;
    out.flinlmxacf_b = NaN; % this is an important statistic
    out.flinlmxacf_r2 = NaN; % this is a more important statistic!
    out.flinlmxacf_adjr2 = NaN;
    out.flinlmxacf_rmse = NaN;
end


end