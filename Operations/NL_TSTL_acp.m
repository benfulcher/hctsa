function out = NL_TSTL_acp(y,tau,past,maxDelay,maxDim,Nref,randomSeed)
% NL_TSTL_acp   acp function in TSTOOL
%
% The documentation isn't crystal clear, but this function seems to be related
% to cross-prediction.
%
%---INPUTS:
% y, time series
% tau, delay time
% past, number of samples to exclude before and after each index (to avoid
%               correlation effects ~ Theiler window)
% maxDelay, maximal delay (<< length(y))
% maxDim, maximal dimension to use
% Nref, number of reference points
% randomSeed, whether (and how) to reset the random seed, using BF_ResetSeed
%
%---OUTPUTS: statistics summarizing the output of the routine.

% TSTOOL: http://www.physik3.gwdg.de/tstool/
%
% May in future want to also make outputs normalized by first value; so get
% metrics on both absolute values at each dimension but also some
% indication of the shape
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

try
    s = signal(y);
catch
    error('Error running ''signal'' on the input time series -- has TSTOOL been installed?')
end
if ~isa(s,'signal')
    error('Error making a signal class of the input time series')
end
N = length(y); % length of the time series

% ------------------------------------------------------------------------------
%% Check inputs
% ------------------------------------------------------------------------------
% (*) tau
if nargin < 2
    tau = 'ac'; % use first zero-crossing of autocorrelation function as default
end
if strcmp(tau,'mi')
    tau = CO_FirstMin(y,'mi');
elseif strcmp(tau,'ac')
    tau = CO_FirstZero(y,'ac');
end
% time delay can't be more than 1/20th of time series length
if tau > N/20
    tau = floor(N/20);
end

% (*) past
if nargin < 3 || isempty(past)
    past = 1;
end
if (past > 0) && (past < 1)
	past = floor(past*N); % specify a proportion of the time series length
end

% (*) maxDelay
if nargin < 4 || isempty(maxDelay)
    maxDelay = ceil(N/4);
end
if (maxDelay > 0) && (maxDelay < 1)
	maxDelay = ceil(maxDelay*N); % specify a proportion of the time series length
end

% (*) maxDim
if nargin < 5 || isempty(maxDim)
    maxDim = 10;
end

% (*) Nref
if nargin < 6 || isempty(Nref)
    Nref = ceil(N/10); % 1/10 of the time series length (I think should be greater -- 50%)
end
if (Nref > 0) && (Nref < 1)
	Nref = ceil(Nref*N); % specify a proportion of the time series length
end

% (*) randomSeed: how to treat the randomization
if nargin < 7
    randomSeed = []; % default
end

% ------------------------------------------------------------------------------
% Run and get data output from TSTOOL function acp:
% ------------------------------------------------------------------------------

% Control the random seed (for reproducibility):
BF_ResetSeed(randomSeed);

% Run TSTOOL's acp function
acpf = data(acp(s,tau,past,maxDelay,maxDim,Nref));

% ------------------------------------------------------------------------------
%% Get outputs
% ------------------------------------------------------------------------------
% Now, the documentation is pretty vague, in fact does not mention at all
% what form the output is in... which is a bit shit. But I gather that each
% column corresponds to an embedding dimension (size(acpf,2)=maxDim), and
% each row corresponds to an increasing time delay up to maxDelay
% (size(acpf,1) = maxDelay+1)...
% I propose a returning some measures here

macpf = mean(acpf); % mean vector of length maxDim
sacpf = std(acpf); % std vector of length maxDim
iqracpf = iqr(acpf); % iqr vector of length maxDim

% How does the mean crossprediction function decay with m?
dmacpf = diff(macpf);
out.mmacpfdiff = mean(abs(dmacpf));
out.stdmacpfdiff = std(abs(dmacpf));
out.propdecmacpf = sum(dmacpf < 0)/length(dmacpf);

for i = 1:maxDim-1
    % Give proportion drop at each increase in m
    out.(sprintf('macpfdrop_%u',i)) = abs(macpf(i)/macpf(i+1)-1);
end

% output statistics on the acp at each dimension
for i = 1:maxDim
    % mean acp at each dimension
    out.(sprintf('macpf_%u',i)) = macpf(i);
    % std of acp at each dimension
    out.(sprintf('sacpf_%u',i)) = sacpf(i);
    % iqr of acp at each dimension
    out.(sprintf('iqracpf_%u',i)) = iqracpf(i);
    % abs(AC1) of acp at each dimension
    out.(sprintf('ac1_acpf_%u',i)) = abs(CO_AutoCorr(acpf(:,i),1,'Fourier'));
end

% plot(macpf)

end
