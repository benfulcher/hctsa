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

s = signal(y);
N = length(y);

%% Check inputs
% (*) tau
if strcmp(tau,'mi')
    tau = CO_fmmi(y);
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

try
    acpf = data(acp(s,tau,past,maxdelay,maxdim,nref));
catch emsg
    if strcmp(emsg.message,'Average interpoint-distance in data set B seems to be zero');
        out = NaN;
        return
    else
		disp(['Unexpected error in TSTL_acp: ' emsg])
	end
end

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
out.propdecmacpf = length(find(dmacpf<0))/length(dmacpf);

for i = 1:maxdim-1
    % give proportion drop at each increase in m
    drophere = abs(macpf(i)/macpf(i+1)-1);
    eval(['out.macpfdrop_' num2str(i) ' = drophere;']);
end

% output statistics on the acp at each dimension
for i = 1:maxdim
    % mean acp at each dimension
    eval(['out.macpf_' num2str(i) ' = macpf(' num2str(i) ');']);
    % std of acp at each dimension
    eval(['out.sacpf_' num2str(i) ' = sacpf(' num2str(i) ');']);
    % iqr of acp at each dimension
    eval(['out.iqracpf_' num2str(i) ' = iqracpf(' num2str(i) ');']);
    % AC1 of acp at each dimension
    ac1 = abs(CO_autocorr(acpf(:,i),1));
    eval(['out.ac1_acpf_' num2str(i) ' = ac1;']);
end

% plot(macpf)

end