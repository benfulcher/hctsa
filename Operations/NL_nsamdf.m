function out = NL_nsamdf(x,fs,winLenRel,shiftLenRel,lagRel,degree,doPlot)
% NL_nsamdf computes the nonlinearity measure L through nsAMDF
% (nonlinear average magnitude difference function), developed by:
%
% Ozkurt et al. (2020), "Identification of nonlinear features in cortical and
% subcortical signals of Parkinson's Disease patients via a novel efficient
% measure", NeuroImage.
%
% Please refer to and cite this paper if you use this function in your work.
%
%---INPUTS:
%
% data = One dimensional input time-series (it can be raw, but it is a good idea
%           to low-pass filter it to get rid of high frequency nuisance)
%           For the neural data (LFP, MEG) in the aferomentioned paper, we
%           low-passed the raw data for 40 Hz.
%           The data should be long enough for proper estimation of nonlinearity.
%
% winLenRel = window length (a long enough segment is important to estimate the nonlinearity)
%
% shiftLenRel = This amounts to window length - overlap length btw windows
%
% lagRel = TMaximum lag for nsAMDF, we chose it as 1.
%
% degree = The chosen degree p should ideally be large enough to capture the
%           highest order of nonlinearity within the data.
%           We chose p=7 for in our case of Parkinsonian data in the paper.
%
% doPlot = true to plot nsAMDF sequences, otherwise just assign it false.
%
%---OUTPUTS:
%
% L: nonlinearity measure
%
% s2: normalized nsAMDF for the degree 2
%
% sd: normalized nsAMDF for the chosen degree greater than 2
%
%
% Required subfunctions are NormedSingleCurveLengthWindowed.m & NormedSingleCurveLength.m
%
%
%   Authored by Tolga Esat Ozkurt, 2020. (tolgaozkurt@gmail.com)
%   Edits by Ben Fulcher for incorporating into hctsa.

%-------------------------------------------------------------------------------
% Set defaults:
if nargin < 2
    fs = 1;
end
if nargin < 3
    winLenRel = 14;
end
windowLength = winLenRel*fs;
if nargin < 4
    shiftLenRel = 0.5;
end
shiftLength = shiftLenRel*windowLength;
if nargin < 5
    lagRel = 1;
end
lag = fs*lagRel;
if nargin < 6
    degree = 7;
end
if nargin < 7
    doPlot = false
end

%-------------------------------------------------------------------------------
% nsAMDF for p = 2:
s2 = NormedSingleCurveLengthWindowed(x,windowLength,shiftLength,lag,fs,2);
out.s2 = s2 ./ max(s2); % normalized

% nsAMDF for p = degree:
sd = NormedSingleCurveLengthWindowed(x,windowLength,shiftLength,lag,fs,degree);
out.sd = sd ./ max(sd); % normalized

% If you like, you can bandpass filter s2 and sd for the specific frequency band
% of nonlinear effect both to compute L and plot them as such
out.L = norm(s2 - sd);

%-------------------------------------------------------------------------------
if doPlot
    figure
    plot(s2,'b')
    hold on
    plot(sd,'g')
end

end
