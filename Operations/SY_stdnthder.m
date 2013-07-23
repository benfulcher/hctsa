% SY_StdNthDer
% 
% Estimates the standard deviation of the nth derivative of the time series.
% 
% Based on an idea by Vladimir Vassilevsky, a DSP and Mixed Signal Design
% Consultant in a Matlab forum, who stated that You can measure the standard
% deviation of the nth derivative, if you like".
% 
% cf. http://www.mathworks.de/matlabcentral/newsreader/view_thread/136539
% 
% The derivative is estimated very simply by simply taking successive increments
% of the time series; the process is repeated to obtain higher order
% derivatives.
% 
% INPUTS:
% 
% y, time series to analyze
% 
% n, the order of derivative to analyze
% 

function out = SY_StdNthDer(y,n)
% Ben Fulcher, 2010

if nargin < 2 || isempty(n)
    n = 2;
end

yd = diff(y,n); % crude method of taking a derivative that could be improved
                % upon in future
out = std(yd);

end