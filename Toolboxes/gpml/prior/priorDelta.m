function [lp,dlp] = priorDelta(x)
% PRIORDELTA fix the value of a hyperparameter.
%
% This is a summy hyperparameter prior distribution.
% The function is not intended to be evaluated but exists merely to make
% the user aware of the possibility to use it. The function is equivalent
% to priorClamped.
%
% For more help on design of priors, try "help priorDistributions".
%
% See also PRIORDISTRIBUTIONS

% Copyright (c) by Roman Garnett and Hannes Nickisch, 2014-12-06.

error('The function is not intended to be called directly.')
