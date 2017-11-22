function BF_ResetSeed(resetHow)
% BF_ResetSeed      Reset the random seed generator in Matlab
%
% Allows functions using random numbers to produce repeatable results with a
% consistent syntax.

% ------------------------------------------------------------------------------
% Copyright (C) 2017, Ben D. Fulcher <ben.d.fulcher@gmail.com>,
% <http://www.benfulcher.com>
%
% If you use this code for your research, please cite the following two papers:
%
% (1) B.D. Fulcher and N.S. Jones, "hctsa: A Computational Framework for Automated
% Time-Series Phenotyping Using Massive Feature Extraction, Cell Systems 5: 527 (2017).
% DOI: 10.1016/j.cels.2017.10.001
%
% (2) B.D. Fulcher, M.A. Little, N.S. Jones, "Highly comparative time-series
% analysis: the empirical structure of time series and their methods",
% J. Roy. Soc. Interface 10(83) 20130048 (2013).
% DOI: 10.1098/rsif.2013.0048
%
% This work is licensed under the Creative Commons
% Attribution-NonCommercial-ShareAlike 4.0 International License. To view a copy of
% this license, visit http://creativecommons.org/licenses/by-nc-sa/4.0/ or send
% a letter to Creative Commons, 444 Castro Street, Suite 900, Mountain View,
% California, 94041, USA.
% ------------------------------------------------------------------------------

% Check defaults
if nargin < 1 || isempty(resetHow)
    resetHow = 'default';
end

% ------------------------------------------------------------------------------
% Reset the seed
switch resetHow
case 'default'
    % Reset to the default (the Mersenne Twister with seed 0)
    rng(0,'twister');
case 'none' % don't change the seed
    return
otherwise
    error('Now sure how to reset using ''%s''',resetHow);
end

end
