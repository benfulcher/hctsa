function [INP_mops,INP_ops] = TS_GiveMeInputFiles(whatFeatureSet)
% TS_GiveMeInputFiles  Gives a set of input files corresponding to a named feature set
%
%---INPUT:
% whatFeatureSet: The name of a feature set.
%
%---OUTPUTS:
% INP_mops: the filename for a master operations input file.
% INP_ops: the filename for an operations input file.

%-------------------------------------------------------------------------------
% Copyright (C) 2022, Ben D. Fulcher <ben.d.fulcher@gmail.com>,
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
%-------------------------------------------------------------------------------

if nargin < 1
    whatFeatureSet = 'hctsa';
end

switch whatFeatureSet
case 'hctsa'
    INP_mops = 'INP_mops_hctsa.txt';
    INP_ops = 'INP_ops_hctsa.txt';
case 'catch22'
    INP_mops = 'INP_mops_catch22.txt';
    INP_ops = 'INP_ops_catch22.txt';
otherwise
    error('Unknown feature set: ''%s''',whatFeatureSet)
end

end
