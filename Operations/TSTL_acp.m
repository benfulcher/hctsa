function out = TSTL_acp(y,tau,past,maxdelay,maxdim,nref)
% Uses the TSTOOL package: acp
% y: column vector of data
% tau: delay time
% past: number of samples to exclude before and after each index
%       (to avoid correlation effects)
% maxdelay: maximal delay (<< length(y))
% maxdim: maximal dimension to use
% nref: number of reference points
% may in future want to also make outputs normalized by first value; so get
% metrics on both absolute values at each dimension but also some
% indication of the shape
% Ben Fulcher October 2009

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
    tau = CO_firstmin(y,'mi');
elseif strcmp(tau,'ac')
    tau = CO_fzcac(y);
end
% time delay can't be more than 1/20th of time series length
if tau > N/20
    tau = floor(N/20);
end


% (*) past
if nargin < 3 || isempty(past)
    past = 1;
end
if past > 0 && past < 1
	past = floor(past*N); % specify a proportion of the time series length
end


% (*) maxdelay
if nargin < 4 || isempty(maxdelay)
    maxdelay = floor(N/4);
end
if maxdelay > 0 && maxdelay < 1
	maxdelay = floor(maxdelay*N); % specify a proportion of the time series length
end

% (*) maxdim
if nargin < 5 || isempty(maxdim)
    maxdim = 10;
end

% (*) nref
if nargin < 6 || isempty(nref)
    nref = floor(N/10); % 1/10 of the time series length (I think should be greater -- 50%)
end
if nref > 0 && nref < 1
	nref = floor(nref*N); % specify a proportion of the time series length
end

acpf = data(acp(s,tau,past,maxdelay,maxdim,nref));

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
    ac1 = abs(CO_autocorr(acpf(:,i),1));
    eval(sprintf('out.ac1_acpf_%u = ac1;',i));
end

% plot(macpf)

end