function out = SY_slidwin(y,windowstat,acrosswindowstat,nseg,nmov)
% Calculate windowstat in each window, and compute acrosswindowstat for the 
% set of statistics calculated in each window.
% Ben Fulcher

if nargin < 2 || isempty(windowstat)
    windowstat = 'mean'; % do a sliding window mean
end

wlen = floor(length(y)/nseg); % size of window
inc = floor(wlen/nmov); % increment to move at each step
if inc == 0; inc = 1; end % increment rounded down to zero, prop it up

nsteps = (floor((length(y)-wlen)/inc)+1);
qs = zeros(nsteps,1);

switch windowstat
    case 'mean' % n=1: sliding window mean
        for i = 1:nsteps
            qs(i) = mean(y((i-1)*inc + 1:(i-1)*inc + wlen));
        end
    case 'std' % n=2: sliding window std
        for i = 1:nsteps
            qs(i) = std(y((i-1)*inc + 1:(i-1)*inc + wlen));
        end
    case 'ent' %n=3: sliding window distributional entropy
        for i = 1:nsteps
            ksstats = DN_kssimp(y((i-1)*inc + 1:(i-1)*inc + wlen),'entropy');
            qs(i) = ksstats.entropy;
        end
    case 'apen' % sliding window ApEn
        for i = 1:nsteps
            qs(i) = EN_ApEn(y((i-1)*inc + 1:(i-1)*inc + wlen),1,0.2);
        end
    case 'mom3'
        for i = 1:nsteps
            qs(i) = DN_moments(y((i-1)*inc + 1:(i-1)*inc + wlen),3);
        end
    case 'mom4'
        for i = 1:nsteps
            qs(i) = DN_moments(y((i-1)*inc + 1:(i-1)*inc + wlen),4);
        end
    case 'mom5' % don't worry about this... too much
        for i = 1:nsteps
            qs(i) = DN_moments(y((i-1)*inc + 1:(i-1)*inc + wlen),5);
        end
    case 'lillie'
        for i = 1:nsteps
            qs(i) = HT_disttests(y((i-1)*inc + 1:(i-1)*inc + wlen),'lillie','norm');
        end
    case 'AC1'
        for i = 1:nsteps
            qs(i) = CO_autocorr(y((i-1)*inc + 1:(i-1)*inc + wlen),1);
        end
    otherwise
        error('Unknown statistic ''%s''',windowstat)
end

%         plot(round(wlen/2):inc:(nsteps-1)*inc+round(wlen/2),qs,'r');
switch acrosswindowstat
    case 'std'
        out = std(qs)/std(y);
    case 'apen'
        out = EN_ApEn(qs,1,0.2); % ApEn of the sliding window measures
    case 'ent'
        kssimpouts = DN_kssimp(qs); % get a load of statistics from kernel-smoothed distribution
        out = kssimpouts.entropy; % distributional entropy
    otherwise
        error('Unknown statistic: ''%s''',acrosswindowstat)
end

end