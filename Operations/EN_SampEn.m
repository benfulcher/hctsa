% EN_sampen
% 
% Estimates the Sample Entropy of the time series, SampEn(m,r), using Physionet code
% 
% "Physiological time-series analysis using approximate entropy and sample entropy"
% J. S. Richman and J. R. Moorman, Am. J. Physiol. Heart Circ. Physiol., 278(6)
% H2039 (2000)
% 
% The function can also calculate the SampEn of successive increments of time
% series, i.e., we using an incremental differencing pre-processing, as
% used in the so-called Control Entropy quantity:
% "Control Entropy: A complexity measure for nonstationary signals"
% E. M. Bollt and J. Skufca, Math. Biosci. Eng., 6(1) 1 (2009)
% 
% Our implementation is based on publicly-available physionet code from
% http://www.physionet.org/physiotools/sampen/matlab/1.1/sampenc.m
% 
% The function has been renamed RN_sampenc here
% 
% INPUTS:
% y, the input time series
% M, the embedding dimension
% r, the threshold
% preprocess [opt], incremental difference preprocessing, 'diff1'
% 

function out = EN_SampEn(y,M,r,preprocess)
% Ben Fulcher, November 2009

if nargin < 4
    preprocess = ''; % don't apply preprocessing
end

if ~isempty(preprocess)
    switch preprocess
    case 'diff1'
        % First do an incremental differencing of the time series
        % thus yielding the 'Control Entropy'
        y = BF_zscore(diff(y));
    otherwise
        error('Unknown preprocessing setting ''%s''',preprocess);
    end
end

% Use the physionet code to calculate the Sample Entropy using these parameters:
[e, p, ~, ~] = PN_sampenc(y,M,r);


%% Give outputs
for i = 1:M
    eval(sprintf('out.p%u = p(%u);',i,i));
    eval(sprintf('out.sampen%u = e(%u);',i,i));
end

out.meanchsampen = mean(diff(e));
out.meanchp = mean(diff(p));

end