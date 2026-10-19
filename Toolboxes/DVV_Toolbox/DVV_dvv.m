function data = DVV_dvv(x, m, Nsub, nd, Ntv, tau, refDV)
% Delay Vector Variance method for real and complex signals
%
%
% USAGE: C = dvv (X, m, Nsub, nd, Ntv, tau)
%	X       original real-valued or complex time series
%	m       delay embedding dimension
%	Ntv     number of points on horizontal axes
%	Nsub	number of reference DVs to consider
%	nd      Span over which to perform DVV
%	tau     time delay between delay-vector elements, in samples (default 1).
%	        hctsa modification: the original fixed tau = 1.
%	refDV   (optional) indices (1 to N - m*tau) of the Nsub reference delay vectors.
%	        hctsa modification: by default they are drawn at random (randsample),
%	        as in the original.
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
if nargin < 6 || isempty(tau)
    tau = 1; % time delay (samples) between the elements of each delay vector
end

% ------------------------------------------------------------------------------
% Initial Conditions
% ------------------------------------------------------------------------------
N = length(x);              % Length of input vector
d = zeros(N-m*tau, Nsub);
y = zeros(Ntv,1);

% Make input vector x a column vector
if size(x,2) > size(x,1)
    x = x';
end

% Generate Nsub subset from existing DV's, randomly
if nargin < 7 || isempty(refDV)
    ref = randsample(N - m*tau,Nsub) + m*tau;
else
    ref = refDV(:) + m*tau; % (the target index of each reference delay vector)
end

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
% for each reference, sort its distances once; the set within each rd is then
% a leading run of the sorted targets, whatever the order among equal
% distances. The set around a reference excludes the reference vector itself;
% the original instead excluded the target with raw index k, the loop counter
% rather than ref(k).)
% (hctsa modification: the variance of each set is computed in two passes,
% the mean and then the mean squared deviation from it (Chan, Golub and
% LeVeque, Am. Stat. 37, 242 (1983)), which cannot be negative and is exactly
% 0 for equal targets. A variance below tol^2, with tol = 1e-12 of the range of
% x, is that of targets equal up to rounding, and is 0.)
tol = 1e-12 * (max(x) - min(x));
tot = zeros(1, Ntv);
count = zeros(1, Ntv);
rdRow = rd(:)';
% (hctsa modification: the variance of each leading run of sorted targets comes
% from prefix sums of the mean-shifted targets, and the counts within each rd
% from one histogram of the sorted distances. Prefix sums lose accuracy when a
% set's variance is small relative to the spread of x, so below 1e-4 of the
% mean square of the shifted x (error ~1e-12 relative) the two-pass var() is
% used, as before; this also keeps the tol clamp exact.)
xc = x - mean(x);
vSmall = 1e-4 * mean(xc.^2);
rdUp = rdRow + eps(rdRow); % next double above rd (rd > 0 where used): ds <= rd  <=>  ds < rdUp
for k = 1:Nsub
    [ds, o] = sort(d(:,k));
    notSelf = (o ~= ref(k) - m*tau); % (the reference vector itself)
    ds = ds(notSelf);
    xs = x(o(notSelf) + m*tau);
    % number of targets within each rd (1 x Ntv), counting ds <= rd
    c = cumsum(histcounts(ds, [-Inf, rdUp, Inf]));
    c = c(1:Ntv);
    % Only those variance values are considered for which the corresponding
    % sets have atleast 30 DVs
    ok = find((c >= 30) & (rdRow > 0));
    if isempty(ok), continue; end
    xk = xc(o(notSelf) + m*tau);
    cn = c(ok)';
    cs1 = cumsum(xk(1:max(cn)));
    cs2 = cumsum(xk(1:max(cn)).^2);
    vv = ((cs2(cn) - cs1(cn).^2 ./ cn) ./ (cn - 1))';
    for jj = find(~(vv >= vSmall)) % small or NaN: two-pass
        vv(jj) = var(xs(1:cn(jj)));
    end
    vv(vv <= tol^2) = 0;
    tot(ok) = tot(ok) + vv;
    count(ok) = count(ok) + 1;
end
y = nan(Ntv, 1);
y(count > 0) = (tot(count > 0) ./ (count(count > 0) * var(x)))';

% Horizontal axis
T = (rd'-avg)/variance;

% DVV Output
data = [T,y];

end
