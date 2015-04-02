% ------------------------------------------------------------------------------
% BF_ResetSeed
% ------------------------------------------------------------------------------
% 
% Resets the random seed generator in Matlab, so that operations using random
% numbers produce repeatable results
% 
%---HISTORY:
% Ben Fulcher, 2015-03-15
% 
% ------------------------------------------------------------------------------
% Copyright (C) 2013,  Ben D. Fulcher <ben.d.fulcher@gmail.com>,
% <http://www.benfulcher.com>
% 
% If you use this code for your research, please cite:
% B. D. Fulcher, M. A. Little, N. S. Jones, "Highly comparative time-series
% analysis: the empirical structure of time series and their methods",
% J. Roy. Soc. Interface 10(83) 20130048 (2010). DOI: 10.1098/rsif.2013.0048
% 
% This work is licensed under the Creative Commons
% Attribution-NonCommercial-ShareAlike 3.0 Unported License. To view a copy of
% this license, visit http://creativecommons.org/licenses/by-nc-sa/3.0/ or send
% a letter to Creative Commons, 444 Castro Street, Suite 900, Mountain View,
% California, 94041, USA.
% ------------------------------------------------------------------------------

function BF_ResetSeed(resetHow)

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