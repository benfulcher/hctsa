function [lp,dlp] = priorSameMulti(x)

% Dummy hyperparameter prior distribution to have a group of hyperparameter
% share the same value.
% The function is not intended to be evaluated but exists merely to make
% the user aware of the possibility to use it. The function is equivalent
% to priorEqualMulti.
%
% For more help on design of priors, try "help priorDistributions".
%
% Copyright (c) by Hannes Nickisch, 2016-10-26.
%
% See also priorDistributions.m.

error('The function is not intended to be called directly.')