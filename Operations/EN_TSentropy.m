% TSentropy
% 
% Estimates the Tsallis entropy of a signal using a parameter q that
% measures the non-extensivity of the system; q = 1 recovers the Shannon
% entropy.
% 
% INPUTS:
% x, the time series
% q, the non-extensivity parameter
% 
% Uses code written by D. Tolstonogov and obtained from
% http://download.tsresearchgroup.com/all/tsmatlablink/TSentropy.m.
%

function out = EN_TSentropy(x, q)
% Wrapper for TS_entropy
% Ben Fulcher 2009

if nargin < 2
    q = 1; % Shannon entropy by default
end

out = TS_entropy(x,q);

end