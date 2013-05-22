function out = MD_rawHRVmeas(x)
% Ben Fulcher 24/2/2011 -- want to do this on raw RR intervals; probably will not make much sense for other time series
% Note that pNNx is done in MD_pNN.m
% This code is very much derived from Max Little's hrv_classic.m code

N = length(x);

% Triangular histogram index
out.tri10 = N/max(hist(x,10));
out.tri20 = N/max(hist(x,20));
out.trisqrt = N/max(hist(x,sqrt(N)));

% Poincare plot measures: see
% "Do Existing Measures ... ", Brennan et. al. (2001), IEEE Trans Biomed Eng 48(11)
diffx = diff(x);
out.SD1 = 1/sqrt(2) * std(diffx) * 1000;
out.SD2 = sqrt(2 * var(x) - (1/2) * std(diffx)^2) * 1000;

end
