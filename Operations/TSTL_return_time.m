function out = TSTL_return_time(y,NNR,maxT,past,Nref,embedparams)
% Uses TSTOOL code return_time to calculate a histogram of return times.
% INPUTS:
% y: scalar time series
% NNR: number of nearest neighbours
% maxT: maximum return time to consider
% past: Theiler window
% Nref: number of reference indicies
% embedparams: to feed into BF_embed
% Ben Fulcher 12/11/2009

N = length(y); % length of the input time series

%% Check Inputs
% Number of nearest neighbours, NNR
if nargin < 2 || isempty(NNR)
    NNR = 5;
end
if NNR > 0 && NNR < 1 % specify a proportion of time series length
    NNR = floor(NNR*N); if NNR == 0, NNR = 1; end
end

% Maximum return time, maxT
if nargin < 3 || isempty(maxT)
    maxT = 0.1;
end
if maxT > 0 && maxT <= 1 % specify a proportion
    maxT = floor(N*maxT);
    if maxT = 0, maxT = 1; end
end

% Theiler window, past
if nargin < 4 || isempty(past) 
    past = 10;
end
if past > 0 && past < 1 % specify a proportion
    past = floor(N*past);
    if past == 0, past = 1; end
end

% Number of reference points
if nargin < 5 || isempty(Nref)
    Nref = -1; % use all available points
end

% embed parameters
if nargin < 6 || isempty(embedparams)
    embedparams = {'ac','cao'};
    fprintf(1,'Using default embedding using autocorrelation and cao\n')
end


%% Embed the signal
s = BF_embed(y,embedparams{1},embedparams{2},1);

if ~strcmp(class(s),'signal') && isnan(s); % embedding failed
    out = NaN; return
end

%% Run the code
try
    rs = return_time(s, NNR, maxT, past, Nref);
catch emsg
    if strcmp(emsg.message,'Index exceeds matrix dimensions.')
        fprintf(1,'Error evaluating return_time\n');
        out = NaN; return
    end
end

Trett = data(rs);
NN = length(Trett);
% plot(Trett),keyboard

%% Quantify structure in output
out.max = max(Trett);
out.std = std(Trett);
out.pzeros = sum(Trett == 0)/NN;
out.pg05 = sum(Trett>max(Trett)*0.5)/NN;
out.iqr = iqr(Trett);

% recurrent peaks:
icross05 = find((Trett(1:end-1)-0.5*max(Trett)).*(Trett(2:end)-0.5*max(Trett)) < 0);
if ~isempty(icross05) && length(icross05) > 2
    difficross05 = diff(icross05);
    difficross05 = difficross05(difficross05 > 0.4*max(difficross05)); % remove small entries, crossing peaks
    
    out.meanpeaksep = mean(difficross05)/NN;
    out.maxpeaksep = max(difficross05)/NN;
    out.minpeaksep = min(difficross05)/NN;
    out.rangepeaksep = range(difficross05)/NN;
    out.stdpeaksep = std(difficross05)/sqrt(NN);
else
    out.meanpeaksep = NaN;
    out.maxpeaksep = NaN;
    out.minpeaksep = NaN;
    out.rangepeaksep = NaN;
    out.stdpeaksep = NaN;
end

out.statrtys = std(Trett(1:floor(end/2)))/std(Trett(floor(end/2)+1:end));
out.statrtym = mean(Trett(1:floor(end/2)))/mean(Trett(floor(end/2)+1:end));

out.hhist = -sum(Trett(Trett>0).*log(Trett(Trett>0)));


%% course-grain a bit, to 20 bins
nbins = 20;
cglav = zeros(nbins,1);
inds = round(linspace(0,NN,21));
for i = 1:nbins
    cglav(i) = sum(Trett(inds(i)+1:inds(i+1)));
end
% plot(cglav,'k'), keyboard
out.hcgdist = -sum(cglav(cglav>0).*log(cglav(cglav>0)));
out.rangecgdist = range(cglav);
out.pzeroscgdist = sum(cglav == 0)/nbins;


%% Get distribution of distribution of return times
[nhist, xout] = hist(Trett,30);
nhist = nhist/sum(nhist);
% plot(xout,nhist,'o-k'), keyboard
out.maxhisthist = max(nhist);
out.phisthistmin = nhist(1);
out.hhisthist = -sum(nhist(nhist > 0).*log(nhist(nhist > 0)));

end