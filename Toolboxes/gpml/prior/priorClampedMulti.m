function [lp,dlp] = priorClampedMulti(x)

% Dummy hyperparameter prior distribution to fix the value of a hyperparameter.
% The function is not intended to be evaluated but exists merely to make
% the user aware of the possibility to use it. The function is equivalent
% to priorDeltaMulti.
%
% For more help on design of priors, try "help priorDistributions".
%
% Copyright (c) by Roman Garnett and Hannes Nickisch, 2014-12-06.
%
% See also priorDistributions.m.

error('The function is not intended to be called directly.')