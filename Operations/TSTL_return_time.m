function out = TSTL_return_time(y,NNR,maxT,past,Nref,embedparams)
% Uses TSTOOL code return_time to calculate a histogram of return times.
% INPUTS:
% y: scalar time series
% NNR: number of nearest neighbours
% maxT: maximum return time to consider
% past: Theiler window
% Nref: number of reference indicies
% embedparams: to feed into benembed
% Ben Fulcher 12/11/2009

N = length(y);

%% Check Inputs
% Number of nearest neighbours, NNR
if nargin < 2 || isempty(NNR)
    NNR = 5;
end
if NNR>0 && NNR<1 % specify a proportion of time series length
    NNR = floor(NNR*N);
end

% Maximum return time, maxT
if nargin < 3 || isempty(maxT)
    maxT = 0.1;
end
if maxT>0 && maxT<=1 % specify a proportion
    maxT = floor(N*maxT);
end

% Theiler window, past
if nargin < 4 || isempty(past) 
    past = 10;
end
if past>0 && past<1 % specify a proportion
    past = floor(N*past);
end

% Number of reference points
if nargin < 5 || isempty(Nref)
    Nref = -1; % use all available points
end

% embed parameters
if nargin < 6 || isempty(embedparams)
    embedparams={'ac','cao'};
    disp('using default embedding using autocorrelation and cao')
end


%% Embed the signal
s = benembed(y,embedparams{1},embedparams{2},1);

if ~strcmp(class(s),'signal') && isnan(s); % embedding failed
    out = NaN;
    return
end
% s = signal(y);

%% Run the code
try
    rs = return_time(s, NNR, maxT, past, Nref);
catch emsg
    if strcmp(emsg.message,'Index exceeds matrix dimensions.')
        disp('Error evaluating return_time');
        out = NaN;
        return
    end
end
% view(rs)

lavery = data(rs);
NN = length(lavery);
% plot(lavery),keyboard

%% Quantify structure in output
out.max = max(lavery);
out.std = std(lavery);
out.pzeros = sum(lavery == 0)/NN;
out.pg05 = sum(lavery>max(lavery)*0.5)/NN;
out.iqr = iqr(lavery);

% recurrent peaks:
icross05 = find((lavery(1:end-1)-0.5*max(lavery)).*(lavery(2:end)-0.5*max(lavery))<0);
if ~isempty(icross05) && length(icross05)>2
    difficross05 = diff(icross05);
    difficross05 = difficross05(difficross05>0.4*max(difficross05)); % remove small entries, crossing peaks
%     plot(difficross05); keyboard
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

out.statrtys = std(lavery(1:floor(end/2)))/std(lavery(floor(end/2)+1:end));
out.statrtym = mean(lavery(1:floor(end/2)))/mean(lavery(floor(end/2)+1:end));

out.hhist = -sum(lavery(lavery>0).*log(lavery(lavery>0)));


%% course-grain a bit, to 20 bins
nbins = 20;
cglav = zeros(nbins,1);
inds = round(linspace(0,NN,21));
for i=1:nbins
    cglav(i) = sum(lavery(inds(i)+1:inds(i+1)));
end
% plot(cglav,'k'), keyboard
out.hcgdist = -sum(cglav(cglav>0).*log(cglav(cglav>0)));
out.rangecgdist = range(cglav);
out.pzeroscgdist = sum(cglav == 0)/nbins;


%% Get distribution of distribution of return times
[maria xout] = hist(lavery,30);
maria = maria/sum(maria);
% plot(xout,maria,'o-k'), keyboard
out.maxhisthist = max(maria);
out.phisthistmin = maria(1);
out.hhisthist = -sum(maria(maria>0).*log(maria(maria>0)));

end