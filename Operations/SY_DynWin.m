% SY_DynWin
% 
% Analyzes how stationarity estimates depend on the number of segments used to
% segment up the time series.
% 
% Specifically, variation in a range of local measures are implemented: mean,
% standard deviation, skewness, kurtosis, ApEn(1,0.2), SampEn(1,0.2), AC(1),
% AC(2), and the first zero-crossing of the autocorrelation function.
% 
% The standard deviation of local estimates of these quantities across the time
% series are calculated as an estimate of the stationarity in this quantity as a
% function of the number of splits, n_{seg}, of the time series.
% 
% INPUTS:
% 
% y, the input time series
% 
% maxnseg, the maximum number of segments to consider. Will sweep from 2
%           segments to maxnseg.
% 
% 
% Outputs of the operation are the standard deviation of this set of
% 'stationarity' estimates across these window sizes.
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

function out = SY_DynWin(y,maxnseg)
% Ben Fulcher, August 2009

if nargin < 2 || isempty(maxnseg)
    maxnseg = 10;
end

nsegr = (2:1:maxnseg); % range of nseg to sweep across
nmov = 1; % controls window overlap

nfeat = 11; % number of features
fs = zeros(length(nsegr),nfeat);
taug = CO_FirstZero(y,'ac'); % global tau

for i = 1:length(nsegr)
    nseg = nsegr(i);
    wlen = floor(length(y)/nseg); % window length
    inc = floor(wlen/nmov); % increment to move at each step
    if inc == 0; inc = 1; end % increment rounded down to zero, prop it up

    nsteps = (floor((length(y)-wlen)/inc)+1);
    % qs = struct;
    qs = zeros(nsteps,nfeat);    
    
    for j = 1:nsteps
        ysub = y((j-1)*inc+1:(j-1)*inc+wlen);
        taul = CO_FirstZero(ysub,'ac');
        
        % qs.mean(j) = mean(ysub); % mean
        % qs.std(j) = std(ysub); % standard deviation
        % qs.skew(j) = skewness(ysub); % skewness
        % qs.kurt(j) = kurtosis(ysub); % kurtosis
        % qs.apen(j) = EN_ApEn(ysub,1,0.2); % ApEn_1
        % qs.sampen(j) = EN_sampenc(ysub,1,0.2); % SampEn_1
        % qs.ac1(j) = CO_AutoCorr(ysub,1); % AC1
        % qs.ac2(j) = CO_AutoCorr(ysub,2); % AC2
        % qs.tauglob(j) = CO_AutoCorr(ysub,taug); % AC_glob_tau
        % qs.tauloc(j) = CO_AutoCorr(ysub,taul); % AC_loc_tau
        % qs.taul(j) = taul;

        qs(j,1) = mean(ysub); % mean
        qs(j,2) = std(ysub); % standard deviation
        qs(j,3) = skewness(ysub); % skewness
        qs(j,4) = kurtosis(ysub); % kurtosis
        qs(j,5) = EN_ApEn(ysub,1,0.2); % ApEn_1
        qs(j,6) = PN_sampenc(ysub,1,0.2,1); % SampEn_1
        qs(j,7) = CO_AutoCorr(ysub,1); % AC1
        qs(j,8) = CO_AutoCorr(ysub,2); % AC2
        qs(j,9) = CO_AutoCorr(ysub,taug); % AC_glob_tau
        qs(j,10) = CO_AutoCorr(ysub,taul); % AC_loc_tau
        qs(j,11) = taul;
    end
    % plot(qs,'o-');
    % input('what do you think?')
    
    % fs(i,1:nfeat) = structfun(@(x)std(x),qs,'UniformOutput',1);
    % fs(i) = structfun(@(x)std(x),qs,'UniformOutput',0);
    fs(i,1:nfeat) = std(qs);
end

% fs contains std of quantities at all different 'scales' (segment lengths)

fs = std(fs); % how much does the 'std stationarity' vary over different scales?

% plot(fs)
out.stdmean = fs(1);
out.stdstd = fs(2);
out.stdskew = fs(3);
out.stdkurt = fs(4);
out.stdapen1_02 = fs(5);
out.stdsampen1_02 = fs(6);
out.stdac1 = fs(7);
out.stdac2 = fs(8);
out.stdactaug = fs(9);
out.stdactaul = fs(10);
out.stdtaul = fs(11);


end