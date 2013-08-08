% MD_rawHRVmeas
% 
% Evaluates the triangular histogram index and Poincare plot measures to a time
% series assumed to measure sequences of consecutive RR intervals measured in
% milliseconds. Doesn't make much sense for other time series
% 
% cf. "Do existing measures of Poincare plot geometry reflect nonlinear
%      features of heart rate variability?"
%      M. Brennan, et al., IEEE T. Bio.-Med. Eng. 48(11) 1342 (2001)
% 
% Note that pNNx is not done here, but in MD_pNN.m
% 
% This code is heavily derived from Max Little's hrv_classic.m code
% Max Little: http://www.maxlittle.net/
% 

function out = MD_rawHRVmeas(x)
% Ben Fulcher 24/2/2011 -- want to do this on raw RR intervals

N = length(x); % time-series length

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
