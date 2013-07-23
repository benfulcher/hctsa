% RN_Gaussian
% 
% Returns a random number drawn from a normal distribution.
% 
% This operation was sometimes used as a control, as a kind of 'null operation'
% with which to compare to the performance of other operations that use
% informative properties of the time-series data, but is otherwise useless for
% any real analysis.

function out = RN_Gaussian(y)
% Ben Fulcher, 2009

out = randn(1,1);

end