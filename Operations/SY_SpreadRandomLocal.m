% SY_SpreadRandomLocal
% 
% Implements a bootstrap-based stationarity measure: nseg time-series segments
% of length l are selected at random from the time series and in each
% segment a local quantity is calculated: mean, standard deviation, skewness,
% kurtosis, ApEn(1,0.2), SampEn(1,0.2), AC(1), AC(2), and the first
% zero-crossing of the autocorrelation function.
% 
% INPUTS:
% y, the input time series
% 
% l, the length of local time-series segments to analyze as a positive integer.
%    Can also be a specified character string:
%       (i) 'ac2': twice the first zero-crossing of the autocorrelation function
%       (ii) 'ac5': five times the first zero-crossing of the autocorrelation function
% 
% nseg, the number of randomly-selected local segments to analyze
% 
% Outputs are the mean and also the standard deviation of this set of 100 local
% estimates.
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

function out = SY_SpreadRandomLocal(y,l,nseg)
% Ben Fulcher, August 2009

doplot = 0; % set to 1 to plot outputs to figure

if nargin < 2 || isempty(l)
    l = 100; % by default use 100 samples
end

if ischar(l)
    switch l
    case 'ac2'
        taug = CO_FirstZero(y,'ac'); % tau (global)
        l = 2*taug;
    case 'ac5'
        taug = CO_FirstZero(y,'ac'); % tau (global)
        l = 5*taug;
    otherwise
        error('Unknown specifier ''%s''',l);
    end
end

if nargin < 3 || isempty(nseg)
    nseg = 100; % 100 segments by default
end

N = length(y); % the length of the time series

if l > 0.9*N % not suitable -- too short
	fprintf(1,'This time series (N = %u) is too short to use l = %u\n',N,l)
    out = NaN; return % NaN means not suitable
end

% nseg segments, each of length segl data points

nfeat = 9; % number of features
fs = zeros(nseg,nfeat);
qs = zeros(nseg,nfeat);

for j = 1:nseg
    % pick a range
    % in this implementation, ranges CAN overlap
    % ist = randint(1,1,[1, N-1-l]); % random start point (not exceeding the endpoint)
    ist = randi(N-1-l,1); % random start point (not exceeding the endpoint)
    ifh = ist+l-1; % finish index
    rs = ist:ifh; % sample range (from starting to finishing index)
    ysub = y(rs); % subsection

    taul = CO_FirstZero(ysub,'ac');
    
    qs(j,1) = mean(ysub); % mean
    qs(j,2) = std(ysub); % standard deviation
    qs(j,3) = skewness(ysub); % skewness
    qs(j,4) = kurtosis(ysub); % kurtosis
    qs(j,5) = EN_ApEn(ysub,1,0.2); % ApEn_1
    qs(j,6) = PN_sampenc(ysub,1,0.2,1); % SampEn_1
    qs(j,7) = CO_AutoCorr(ysub,1); % AC1
    qs(j,8) = CO_AutoCorr(ysub,2); % AC2
    qs(j,9) = taul;
end
    
if doplot
    figure('color','w');
    subplot(2,1,1); hold on; plot(y,'k'); plot(ists,y(ists),'.r'); title('time series')
    subplot(2,1,2); plot(qs(:,1),'b'); title('local means')
    input('What do you think?')
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
%         out=DN_FitKernelSmooth(qs,'entropy'); % distributional entropy
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