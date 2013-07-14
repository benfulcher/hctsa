function out = TSTL_amutual(y,maxtau,nbins)
% Uses amutual from the TSTOOL package
% maxtau: the maximum lag for which to calculate the auto mutual information
% nbins: the number of bins for histogram calculation
% Ben Fulcher October 2009

%% Preliminaries
N = length(y); % length of time series
s = signal(y); % convert to signal object for TSTOOL

%% Check Inputs
if nargin < 2 || isempty(maxtau)
    maxtau = ceil(N/4);
end

if nargin < 3 || isempty(nbins)
    nbins = round(sqrt(N/10)); % this is an arbitrary choice (!!) ;-)
end

%% Run
ami = data(amutual(s,maxtau,nbins));
lami = length(ami);
% plot(ami)

% change ami vector to a structure for output
for i = 1:maxtau+1
    eval(sprintf('out.ami%u = ami(%u);',i,i));
end

% mean mutual information over this lag range
out.mami = mean(ami);
out.stdami = std(ami);

% first miniimum of mutual information across range
dami = diff(ami);
extremai = find(dami(1:end-1).*dami(2:end)<0);
out.pextrema = length(extremai)/(lami-1);
if isempty(extremai)
   out.fmmi = lami; % actually represents lag, because indexes don't but diff delays by 1
else
    out.fmmi = extremai(1);
end

% Look for periodicities in local maxima
maximai = find(dami(1:end-1) > 0 & dami(2:end) < 0) + 1;
dmaximai = diff(maximai);
% is there a big peak in dmaxima?
 % (no need to normalize since a given method inputs its range; but do it anyway... ;-))
out.pmaxima = length(dmaximai)/floor(lami/2);
out.modeperiodmax = mode(dmaximai);
out.pmodeperiodmax = sum(dmaximai == mode(dmaximai))/length(dmaximai);


% hold on; plot(maximai,ami(maximai),'or'); hold off


end