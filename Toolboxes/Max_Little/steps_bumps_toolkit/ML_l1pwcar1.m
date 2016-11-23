% Performs discrete correlated total variation denoising (CTVD) using a
% primal-dual interior-point solver. It minimizes the following discrete
% functional:
%
%  E=(1/2)||y_0-ay_1-x||_2^2+lambda*||Dx||_1,
%
% over the variable x, given the input signal y, according to each
% value of the regularization parameter lambda > 0. a is the
% autocorrelation of y at time lag 1. D is the first difference matrix.
% Uses hot-restarts from each value of lambda to speed up convergence for
% subsequent values: best use of this feature is made by ensuring that the
% chosen lambda values are close to each other.
%
% Usage:
% [x, E, s] = ML_l1pwcar1(y, lambda, a, display, stoptol, maxiter)
%
% Input arguments:
% - y          Original signal to denoise, size N x 1.
% - lambda     A vector of positive regularization parameters, size L x 1.
%              TVD will be applied to each value in the vector.
% - a          Autocorrelation parameter.
% - display    (Optional) Set to 0 to turn off progress display, 1 to turn
%              on. If not specifed, defaults to progress display on.
% - stoptol    (Optional) Precision as determined by duality gap tolerance,
%              if not specified, defaults to 1e-3.
% - maxiter    (Optional) Maximum interior-point iterations, if not
%              specified defaults to 60.
%
% Output arguments:
% - x          Denoised output signal for each value of lambda, size N x L.
% - E          Objective functional at minimum for each lambda, size L x 1.
% - s          Optimization result, 1 = solved, 0 = maximum iterations
%              exceeded before reaching duality gap tolerance, size L x 1.
%
% (c) Max Little, 2010. Based around code originally written by 
% S.J. Kim, K. Koh, S. Boyd and D. Gorinevsky.
% If you use this code for your research, please cite:
% M.A. Little, B.C. Steel, F. Bai, Y. Sowa, T. Bilyard, D.M. Mueller,
% R.M. Berry, N.S. Jones (2011)
% Steps and bumps: precision extraction of discrete states of molecular machines
% Biophysical Journal, 101(2):477-485
%
% This code is released under the terms of GNU General Public License as
% published by the Free Software Foundation; version 2 or later.
% 

function [x, E, s] = ML_l1pwcar1(y, lambda, a, display, stoptol, maxiter)

narginchk(3,6);

if (nargin < 4)
    display = 1;
end
if (nargin < 5)
    stoptol = 1e-3;
end
if (nargin < 6)
    maxiter = 60;
end

y = y(:);

N = length(y);
y1 = y(2:N);
y0 = y(1:N-1);
yd = y1-a*y0;

[mu, E, s] = ML_l1pwc(yd, lambda, display, stoptol, maxiter);
x = mu/(1-a);
