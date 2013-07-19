% Measures the scaling exponent of the time series using a fast implementation
% of detrended fluctuation analysis (DFA).
% Coded by Max A. Little and
% publicly-available at
% http://www.physics.ox.ac.uk/users/littlem/software/index.html.
%
% INPUTS,
% y, the input time series, is fed straight into the fastdfa script.

function out = SC_fastdfa(y)
% Matlab wrapper for Max Little's ML_fastdfa code

if size(y,1) < size(y,2);
    y = y'; % ensure input time series is a column vector
end

out = ML_fastdfa(y);

end