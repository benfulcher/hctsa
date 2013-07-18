% Performs fast detrended fluctuation analysis on a nonstationary input signal to
% obtain an estimate for the scaling exponent.
%
% Useage:
% [alpha, intervals, flucts] = fastdfa(x)
% [alpha, intervals, flucts] = fastdfa(x, intervals)
% Inputs
%    x          - input signal: must be a COLUMN! vector
% Optional inputs
%    intervals  - List of sample interval widths at each scale
%                 (If not specified, then a binary subdivision is constructed)
%
% Outputs:
%    alpha      - Estimated scaling exponent
%    intervals  - List of sample interval widths at each scale
%    flucts     - List of fluctuation amplitudes at each scale
%
% (c) 2006 Max Little. If you use this code, please cite:
% M. Little, P. McSharry, I. Moroz, S. Roberts (2006),
% Nonlinear, biophysically-informed speech pathology detection
% in Proceedings of ICASSP 2006, IEEE Publishers: Toulouse, France.

function out = SC_fastdfa(y)
% Matlab wrapper for Max Little's fastdfa code

if size(y,1) < size(y,2);
    y = y'; % ensure input time series is a column vector
end

out = fastdfa(y);

end