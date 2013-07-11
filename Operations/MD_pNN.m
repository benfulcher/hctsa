function out = MD_pNN(x)
% Packages up classic HRV operations -- adapted from code emailed from Max 
% Little on 26/1/2009. Implemented on 12/4/2010. Just the usual efficient
% turn-around time.
% Ben Fulcher 12/4/2010.
% Ben Fulcher 24/2/2011 -- repackaged old MD_hrv_classic.m becuase it did not make medical sense to do PNN on a z-scored time series. But now PSD doesn't make too much sense. Just do the pNN measures then.

diffx = diff(x); % successive increments of the time series
N = length(x); % length of the time series

%% Calculate pNNx percentage
% pNNx: recommendation as per Mietus et. al. 2002, "The pNNx files: ...", Heart
% strange to do this for a z-scored time series...
% pnntime = 20;

Dx = abs(diffx) * 1000; % assume milliseconds as for RR intervals
pnns = [5, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100];
for i = 1:length(pnns)
    thispnn = sum(Dx > pnns(i)) / (N-1);
    eval(sprintf('out.pnn%u =  thispnn;', pnns(i)))
end

end