function data = DVV_dvv(x, m, Nsub, nd, Ntv)
% Delay Vector Variance method for real and complex signals
%
%
% USAGE: C = dvv (X, m, Nsub, nd, Ntv)
%	X       original real-valued or complex time series
%	m       delay embedding dimension
%	Ntv     number of points on horizontal axes
%	Nsub	number of reference DVs to consider
%	nd      Span over which to perform DVV
%
%
%   A Delay Vector Variance (DVV) toolbox for MATLAB
%   (c) Copyright Danilo P. Mandic 2008
%   http://www.commsp.ee.ic.ac.uk/~mandic/dvv.htm
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   This program is free software; you can redistribute it and/or modify
%   it under the terms of the GNU General Public License as published by
%   the Free Software Foundation; either version 2 of the License, or
%   (at your option) any later version.
%
%   This program is distributed in the hope that it will be useful,
%   but WITHOUT ANY WARRANTY; without even the implied warranty of
%   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%   GNU General Public License for more details.
%
%   You can obtain a copy of the GNU General Public License from
%   http://www.gnu.org/copyleft/gpl.html or by writing to
%   Free Software Foundation, Inc.,675 Mass Ave, Cambridge, MA 02139, USA.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% ------------------------------------------------------------------------------
% Check parameters
% ------------------------------------------------------------------------------
if nargin < 1
	error('Not enough input arguments');
end
if nargin < 2 || isempty(m)
	m = 3;
end
if nargin < 3 || isempty(Nsub)
	Nsub = 200;
end
if nargin < 4 || isempty(nd)
	nd = 2.0;
end
if nargin < 5 || isempty(Ntv)
	Ntv = 25*nd;
end
if nargin < 6 || isempty(numSurr)
    numSurr = 10;
end

% ------------------------------------------------------------------------------
% Initial Conditions
% ------------------------------------------------------------------------------
N = length(x);              % Length of input vector
tau = 1;                    % Time delay parameter
d = zeros(N-m*tau, Nsub);
y = zeros(Ntv,1);

% Make input vector x a column vector
if size(x,2) > size(x,1)
    x = x';
end

% Generate Nsub subset from existing DV's, randomly
ref = randsample(N - m*tau,Nsub) + m*tau;

% Compute pairwise distances between reference DVs and all DVs
% (hctsa modification: the original evaluated norm() once per (reference,
% point) pair -- Nsub x (N - m*tau) scalar calls, and this function is run
% once per surrogate too -- which dominated NL_DVV's cost. Build the delay
% vectors once and use pdist2; d, avg and variance are exactly as before.)
Nd = N - m*tau;
Y = zeros(Nd, m); % delay vectors: row r holds x(r : tau : r + (m-1)*tau), i.e., the DV of point j = r + m*tau
for k = 1:m
    Y(:,k) = x((1:Nd) + (k-1)*tau);
end
d = pdist2(Y, Y(ref - m*tau, :)); % Nd x Nsub: d(j - m*tau, i) = ||DV_ref(i) - DV_j||

% Mean and std variation calculation of input data, over all (reference,
% point) pairs except each reference's pairing with itself (distance 0):
count = Nsub*Nd - Nsub;
avg = sum(d(:))/count;
variance = sqrt((sum((d(:)-avg).^2) - Nsub*avg^2)/(count-1));

% Calculates the range vector consisting of Ntv equally spaced regions
n = (1:Ntv)-1;
rd = avg-nd*variance + (2*nd*variance*n)/(Ntv-1);

% Creates sets of DV's, for each ref element of subset and value rd, which have norms closer than distance rd to ref
% (hctsa modification: the original looped over every (rd, reference) pair
% with a find() and a var() over the selected targets. Equivalent here:
% for each reference, sort its distances once and read the variance of
% the targets within each rd off cumulative sums. The original excluded the
% target with raw index k (the loop counter, not ref(k)) from each set; that
% is reproduced exactly so results are unchanged.)
tot = zeros(1, Ntv);
count = zeros(1, Ntv);
rdRow = rd(:)';
for k = 1:Nsub
    [ds, o] = sort(d(:,k));
    xs = x(o + m*tau);
    S1 = cumsum(xs);
    S2 = cumsum(xs.^2);
    c = sum(ds <= rdRow, 1); % number of targets within each rd (1 x Ntv)
    s1 = zeros(1, Ntv); s2 = zeros(1, Ntv);
    s1(c > 0) = S1(c(c > 0)); s2(c > 0) = S2(c(c > 0));
    % Exclude target index k, as the original did (only possible when k is a
    % valid target index, i.e., k > m*tau, and it lies within rd):
    if k > m*tau
        isIn = d(k - m*tau, k) <= rdRow;
        c(isIn) = c(isIn) - 1;
        s1(isIn) = s1(isIn) - x(k);
        s2(isIn) = s2(isIn) - x(k)^2;
    end
    % Only those variance values are considered for which the corresponding
    % sets have atleast 30 DVs
    ok = (c >= 30) & (rdRow > 0);
    v = (s2(ok) - s1(ok).^2 ./ c(ok)) ./ (c(ok) - 1); % sample variance of x over the set
    tot(ok) = tot(ok) + v;
    count(ok) = count(ok) + 1;
end
y = nan(Ntv, 1);
y(count > 0) = (tot(count > 0) ./ (count(count > 0) * var(x)))';

% Horizontal axis
T = (rd'-avg)/variance;

% DVV Output
data = [T,y];

end
