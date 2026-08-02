function K = apxSparse(cov, xu, hyp, x, z)

% apxSparse - Covariance function for sparse approximations.
%
% This is a special covariance function which handles some non-standard
% requirements of sparse covariances
% 1) it holds the inducing inputs in xu, and
% 2) it returns covariances for prediction purposes. Note, that sparse
% inference methods don't actually call this function.
%
% Any sparse approximation to the posterior Gaussian process is equivalent to
% using the approximate covariance function:
%   Kt = Q + s*diag(g); g = diag(K-Q); Q = Ku'*inv(Kuu+diag(snu2))*Ku;
% where Ku and Kuu are covariances w.r.t. to inducing inputs xu and
% snu2 is a vector with the noise variance of the inducing inputs.
% As a default, we fix the standard deviation of the inducing inputs
% snu = sfu/1e3 * ones(nu,1) to be a one per mil of the signal standard
% deviation sfu^2 = trace(Kuu)/nu of the inducing inputs. Alternatively, the
% noise of the inducing inputs can be treated as a hyperparameter by means
% of the field hyp.snu that can contain a scalar value log(snu) or a vector
% of inducing point specific log noise standard deviations.
%
% The parameter s from [0,1] (see [1]) interpolates between the two limiting
% cases s=1: FITC [2] and s=0 VFE [3].
%
% The function is designed to be used with infGaussLik and infLaplace.
%
% The implementation exploits the Woodbury matrix identity
%   inv(Kt+inv(W)) = D - D*Ku'*inv(Kuu+snu2*eye(nu)+Ku*D*Ku')*Ku*D, where
%                    D = diag(d), and d = W./(s*g.*W+1)
% and the Cholesky decomposition in order to be applicable to large datasets.
% The computational complexity is O(n nu^2) where n is the number of data
% points x and nu the number of inducing inputs in xu.
%
% The inducing points can be specified through
%  1) the 2nd apxSparse parameter or by
%  2) providing a hyp.xu hyperparameters to the inference function.
% Note that 2) has priority over 1).
% In case 2) is provided and derivatives dnlZ are requested, there will also be
% a dnlZ.xu field allowing to optimise w.r.t. to the inducing points xu.
%
% [1] Bui, Yan & Turner, A Unifying Framework for Sparse GP Approximation
%     using Power EP, 2016, https://arxiv.org/abs/1605.07066.
% [2] Snelson & Ghahramani, Sparse GPs using pseudo-inputs, NIPS, 2006.
% [3] Titsias, Var. Learning of Inducing Variables in Sparse GPs, AISTATS, 2009
%
% Copyright (c) by Carl Edward Rasmussen and Hannes Nickisch, 2016-12-13.
%
% See also covFunctions.m, cov/apx.m, inf/infLaplace.m, inf/infGaussLik.m.

if nargin<4, K = feval(cov{:}); return, end

if size(xu,2) ~= size(x,2)
  error('Dimensionality of inducing inputs must match training inputs')
end

if strcmp(z, 'diag')
  K = feval(cov{:}, hyp, x,  z);                % normal evaluation for diagonal
else
  K = feval(cov{:}, hyp, xu, z);               % substitue inducing inputs for x
end