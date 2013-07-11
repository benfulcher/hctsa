function out = PD_wang07(y)
% Implements the periodicity extraction measure proposed in Wang (2007).
% IEEE International Conference on Data Mining DOI 10.1109/ICDM.2007.103
% Code by Ben Fulcher 9/7/09
% I think that the threshold being 0.01 is questionable
% So I've implemented for 0,0.01,0.1,0.2 and [1,5,10]/sqrt(N)
% Thus output is for each of the following thresholds:
% (0 0.01 0.1 0.2 1/sqrt(N) 5/sqrt(N) 10/sqrt(N))

% Inputs:
% y is the (univariate) time series vector
% Check the time series is zscored
if ~BF_iszscored(y)
    warning('The input time series should be z-scored for EN_progranz')
end

%% Foreplay
N = length(y); % length of the time series
ths = [0, 0.01, 0.1, 0.2, 1/sqrt(N), 5/sqrt(N), 10/sqrt(N)]; % the thresholds with which to count a peak
nths = length(ths); % the number of thresholds

%% 1: Detrend using a regression spline with 3 knots
% I'm not quite sure how to do this, but I'm doing it like this:
% y_or=y; % the original series
% r=linspace(1,N,3);% range for spline (3 knots)
% y_sp=spline(1:N,y,r); % fit the spline on the data, y
% respline=spline(r,y_sp,1:N);
% y=y-respline'; % the detrended series

spline = spap2(2,4,1:N,y); % just a single middle knot with cubic interpolants
y_spl = fnval(spline,1:N); % evaluated at the 1:N time intervals
y = y - y_spl';
% plot(y_or,'k'); hold on; plot(y,'r'); hold off
% input('is this ok, mdear?')

%% 2. Compute autocorrelations up to 1/3 the length of the time series.
acmax = ceil(N/3); % compute autocorrelations up to this lag
acf = zeros(acmax,1); % the autocorrelation function
for i = 1:acmax % i is the \tau, the AC lag
    acf(i) = mean(y(1:N-i).*y(i+1:N));
end
% plot(acf)

%% 3. Frequency is the first peak satisfying the following conditions:
% (a) a trough before it
% (b) difference between peak and trough is at least 0.01
% (c) peak corresponds to positive correlation

% (i) find peaks and troughs in ACF
diffac = diff(acf); % differenced time series
sgndiffac = sign(diffac); % sign of differenced time series
bath = diff(sgndiffac); % differenced, signed, differenced time series
troughs = find(bath == 2)+1; % finds troughs
peaks = find(bath == -2)+1; % finds peaks
npeaks = length(peaks);

thefreqs = zeros(nths,1);
for k = 1:nths
    thefreqs(k) = 1;
    for i = 1:npeaks % search through all peaks for one that meets the condition
        ipeak = peaks(i); % acf lag at which there is a peak
        thepeak = acf(ipeak); % acf at the peak
        ftrough = find(troughs < ipeak,1,'last');
        if isempty(ftrough); continue; end
        itrough = troughs(ftrough); % acf lag at which there is a trough (the first one preceeding the peak)
        thetrough = acf(itrough); % acf at the trough
        
        % (a) a trough before it: should be implicit in the ftrough bit above
        %     if troughs(1)>ipeak % the first trough is after it
        %         continue
        %     end
        
        % (b) difference between peak and trough is at least 0.01
        if thepeak - thetrough < ths(k)
            continue
        end
        
        % (c) peak corresponds to positive correlation
        if thepeak < 0
            continue
        end
        
        % we made it! Use this frequency!
        thefreqs(k) = ipeak; break
    end
end

%% Convert vector into a structure for output
% output entries are out.th1, out.th2, ..., out.th7:
for i = 1:nths
    eval(sprintf('out.th%u = thefreqs(%u);',i,i))
end

end