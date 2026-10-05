function perm = BF_SpreadPerm(N)
% BF_SpreadPerm     A fixed ordering of 1..N whose first k elements are evenly spread, for every k
%
% A deterministic substitute for a random permutation when the first k elements are
% to be used as a subsample (reference points, pairs, segment starts, ...): the
% indices are ordered by the fractional parts of j*phi, j = 1..N (phi = the golden
% ratio conjugate), the low-discrepancy Weyl sequence. The first k elements of the
% result therefore cover 1..N almost uniformly for any k, without the clumps and gaps
% of a random subsample, and without aliasing with periodicities of the series (as
% a regular stride has). The result does not depend on any random seed.
%
%---INPUTS:
% N, the number of elements
%
%---OUTPUTS:
% perm, a column vector containing each of the integers 1..N once: perm(j) is the
%       rank of frac(j*phi) among frac(1*phi), ..., frac(N*phi)

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

perm = zeros(N, 1);
[~, order] = sort(mod((1:N)' * 0.6180339887498949, 1));
perm(order) = 1:N;

end
