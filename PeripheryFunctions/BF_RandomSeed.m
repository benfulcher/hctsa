function seed = BF_RandomSeed(randomSeed)
% BF_RandomSeed     The integer seed (for BF_Random) that a randomSeed input stands for
%
% Functions that draw random numbers from BF_Random (the portable generator that gives
% the same numbers in every language) accept a randomSeed input, with the same
% options as BF_ResetSeed.
%
%---INPUTS:
% randomSeed, one of:
%           'default' or empty -- the fixed seed 0
%           a numeric scalar -- that seed (rounded, made non-negative)
%           'none' -- a seed drawn from MATLAB's global random stream (so repeated
%                     calls differ and depend on MATLAB's generator)
%
%---OUTPUTS:
% seed, a non-negative integer below 4e9

% ------------------------------------------------------------------------------
% Copyright (C) 2013-2026, Ben D. Fulcher <ben.d.fulcher@gmail.com>,
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

if nargin < 1 || isempty(randomSeed)
    randomSeed = 'default';
end

if isnumeric(randomSeed) && isscalar(randomSeed)
    seed = mod(round(abs(randomSeed)), 4e9);
    return
end

switch randomSeed
case 'default'
    seed = 0;
case 'none'
    seed = floor(4e9 * rand);
otherwise
    error('Not sure how to interpret the random seed ''%s''', randomSeed);
end

end
