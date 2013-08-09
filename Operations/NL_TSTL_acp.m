% NL_TSTL_acp
% 
% Implements the TSTOOL routine acp using a time lag, tau, a Theiler window,
% past, maximum delay, maxdelay, maximum embedding dimension, maxdim, and number
% of reference points, Nref.
% 
% The documentation isn't crystal clear, but I think this function has to do
% with cross-prediction.
% 
% TSTOOL: http://www.physik3.gwdg.de/tstool/
%
% INPUTS:
% 
% y, time series
% 
% tau, delay time
% 
% past, number of samples to exclude before and after each index (to avoid
%               correlation effects ~ Theiler window)
% 
% maxdelay, maximal delay (<< length(y))
% 
% maxdim, maximal dimension to use
% 
% Nref, number of reference points
% 
% Outputs are statistics summarizing the output of the routine.
% 
% May in future want to also make outputs normalized by first value; so get
% metrics on both absolute values at each dimension but also some
% indication of the shape
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

function out = NL_TSTL_acp(y,tau,past,maxdelay,maxdim,Nref)
% Ben Fulcher, October 2009

try
    s = signal(y);
catch
    error('Error running ''signal'' on the input time series -- has TSTOOL been installed?')
end
if ~strcmp(class(s),'signal')
    error('Error making a signal class of the input time series')
end
N = length(y); % length of the time series

%% Check inputs
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


% (*) maxdelay
if nargin < 4 || isempty(maxdelay)
    maxdelay = ceil(N/4);
end
if (maxdelay > 0) && (maxdelay < 1)
	maxdelay = ceil(maxdelay*N); % specify a proportion of the time series length
end

% (*) maxdim
if nargin < 5 || isempty(maxdim)
    maxdim = 10;
end

% (*) Nref
if nargin < 6 || isempty(Nref)
    Nref = ceil(N/10); % 1/10 of the time series length (I think should be greater -- 50%)
end
if (Nref > 0) && (Nref < 1)
	Nref = ceil(Nref*N); % specify a proportion of the time series length
end


% Run and get data output from TSTOOL function acp:
acpf = data(acp(s,tau,past,maxdelay,maxdim,Nref));


%% Get outputs
% Now, the documentation is pretty vague, in fact does not mention at all
% what form the output is in... which is a bit shit. But I gather that each
% column corresponds to an embedding dimension (size(acpf,2)=maxdim), and
% each row corresponds to an increasing time delay up to maxdelay
% (size(acpf,1) = maxdelay+1)...
% I propose a returning some measures here

macpf = mean(acpf); % mean vector of length maxdim
sacpf = std(acpf); % std vector of length maxdim
iqracpf = iqr(acpf); % iqr vector of length maxdim

% How does the mean crossprediction function decay with m?
dmacpf = diff(macpf);
out.mmacpfdiff = mean(abs(dmacpf));
out.stdmacpfdiff = std(abs(dmacpf));
out.propdecmacpf = sum(dmacpf < 0)/length(dmacpf);

for i = 1:maxdim-1
    % give proportion drop at each increase in m
    drophere = abs(macpf(i)/macpf(i+1)-1);
    eval(sprintf('out.macpfdrop_%u = drophere;',i));
end

% output statistics on the acp at each dimension
for i = 1:maxdim
    % mean acp at each dimension
    eval(sprintf('out.macpf_%u = macpf(%u);',i,i));
    % std of acp at each dimension
    eval(sprintf('out.sacpf_%u = sacpf(%u);',i,i));
    % iqr of acp at each dimension
    eval(sprintf('out.iqracpf_%u = iqracpf(%u);',i,i));
    % AC1 of acp at each dimension
    ac1 = abs(CO_AutoCorr(acpf(:,i),1));
    eval(sprintf('out.ac1_acpf_%u = ac1;',i));
end

% plot(macpf)

end