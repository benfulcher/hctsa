function out = MD_pNN(x)
% Packages up classic HRV operations -- adapted from code emailed from Max 
% Little on 26/1/2009. Implemented on 12/4/2010. Just the usual efficient
% turn-around time.
% Ben Fulcher 12/4/2010.
% Ben Fulcher 24/2/2011 -- repackaged old MD_hrv_classic.m becuase was
% medically silly to do pnn on zscored time series. But now psd doesn't
% make too much sense. Just do the pNN measures then.

% function [ tri, SD1, SD2] = hrv_classic(series)
% Calculates classical HRV analysis measures from a given NN/RR time series
% in units of seconds.

diffx = diff(x);
N = length(x);

%% Calculate pNNx percentage
% pNNx: recommendation as per Mietus et. al. 2002, "The pNNx files: ...", Heart
% strange to do this for a z-scored time series...
% pnntime = 20;

Dx = abs(diffx) * 1000; % assume milliseconds as for RR intervals
pnns = [5,10,20,30,40,50,60,70,80,90,100];
for i = 1:length(pnns)
    gaga = sum(Dx > pnns(i))/(N-1);
    eval(['out.pnn' num2str(pnns(i)) ' =  gaga;'])
end



end