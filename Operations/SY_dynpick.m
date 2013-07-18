function out = SY_dynpick(y,l)
% Picks a number (100) of time series segments of a given length;
% returns statistics on measures calculated on these subsegments
% Kind of a bootstrap sort of idea
% Ben Fulcher August 2009

doplot = 0; % set to 1 to plot outputs to figure
if nargin < 2 || isempty(l)
    l = 100; % by default use 100 samples
end


if strcmp(l,'ac2')
    taug = CO_fzcac(y); % tau (global)
    l = 2*taug;
elseif strcmp(l,'ac5')
    taug = CO_fzcac(y); % tau (global)
    l = 5*taug;
end

N = length(y); % the length of the time series

if l > 0.9*N % not suitable -- too short
	fprintf(1,'This time series (N = %u) is too short to use l = %u\n',N,l)
    out = NaN; % NaN means not suitable
    return
end

% nseg segments, each of length segl data points

nfeat = 9; % number of features
nseg = 100; % number of segments
fs = zeros(nseg,nfeat);
qs = zeros(nseg,nfeat);

for j = 1:nseg
    % pick a range
    % in this implementation, ranges CAN overlap
    ist = randint(1,1,[1, N-1-l]); % random start point (not exceeding the endpoint)
    ifh = ist+l-1; % finish index
    rs = ist:ifh; % sample range (from starting to finishing index)
    ysub = y(rs); % subsection

    taul = CO_fzcac(ysub);
    
    qs(j,1) = mean(ysub); % mean
    qs(j,2) = std(ysub); % standard deviation
    qs(j,3) = skewness(ysub); % skewness
    qs(j,4) = kurtosis(ysub); % kurtosis
    qs(j,5) = EN_ApEn(ysub,1,0.2); % ApEn_1
    qs(j,6) = EN_sampenc(ysub,1,0.2); % SampEn_1
    qs(j,7) = CO_autocorr(ysub,1); % AC1
    qs(j,8) = CO_autocorr(ysub,2); % AC2
    qs(j,9) = taul;
end
    
if doplot
    figure('color','w');
    subplot(2,1,1); hold on; plot(y,'k'); plot(ists,y(ists),'.r'); title('time series')
    subplot(2,1,2); plot(qs(:,1),'b'); title('local means')
    input('what do you think?')
end

% Can think of this as a big bootstrapped distribution of the timeseries at
% a scale given by the length l
fs(1:nfeat) = mean(qs); % the mean value of the feature across subsegments of the time series
fs(nfeat+1:nfeat*2) = std(qs); % the spread of the feature across subsegments of the time series
    
%     fs(i,nfeat+1:2*nfeat)=std(qs);

% switch meattray
%     case 'std'
%         out=std(qs)/std(y);
%     case 'apen'
%         out=EN_ApEn(qs,1,0.2); % ApEn of the sliding window measures
%     case 'ent'
%         out=DN_kssimp(qs,'entropy'); % distributional entropy
%     case 'lbq' % lbq test for randomness
%         [h, p] = lbqtest(y);
%         out=p;
% end
% end

% plot(fs)

out.meanmean = fs(1);
out.meanstd = fs(2);
out.meanskew = fs(3);
out.meankurt = fs(4);
out.meanapen1_02 = fs(5);
out.meansampen1_02 = fs(6);
out.meanac1 = fs(7);
out.meanac2 = fs(8);
out.meantaul = fs(9);

out.stdmean = fs(10);
out.stdstd = fs(11);
out.stdskew = fs(12);
out.stdkurt = fs(13);
out.stdapen1_02 = fs(14);
out.stdsampen1_02 = fs(15);
out.stdac1 = fs(16);
out.stdac2 = fs(17);
out.stdtaul = fs(18);

end