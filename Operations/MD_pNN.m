% MD_pNN
% 
% Applies pNNx measures to time series assumed to represent sequences of
% consecutive RR intervals measured in milliseconds.
% 
% cf. "The pNNx files: re-examining a widely used heart rate variability
%           measure", J.E. Mietus et al., Heart 88(4) 378 (2002)
% 
% INPUTS,
% x, the input time series
% 
% This code is derived from MD_hrv_classic.m becuase it doesn't make medical
% sense to do PNN on a z-scored time series.
% 
% But now PSD doesn't make too much sense, so we just evaluate the pNN measures.
% 
% Code is heavily derived from that provided by Max A. Little:
% http://www.maxlittle.net/
%

function out = MD_pNN(x)
% Ben Fulcher 24/2/2011

diffx = diff(x); % successive increments of the time series
N = length(x); % length of the time series

%% Calculate pNNx percentage
% strange to do this for a z-scored time series...
% pnntime = 20;

Dx = abs(diffx) * 1000; % assume milliseconds as for RR intervals
pnns = [5, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100];
for i = 1:length(pnns)
    thispnn = sum(Dx > pnns(i)) / (N-1);
    eval(sprintf('out.pnn%u =  thispnn;', pnns(i)))
end

end