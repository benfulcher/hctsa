% Covariance functions to be use by Gaussian process functions. There are two
% different kinds of covariance functions: simple and composite:
%
% 1) Elementary and standalone covariance functions:
%   covZero       - zero covariance function
%   covEye        - unit covariance function
%   covOne        - unit constant covariance function
%   covDiscrete   - precomputed covariance for discrete data
%
% 2) Composite covariance functions:
%   covScale      - scaled version of a covariance function
%   covSum        - sums of covariance functions
%   covProd       - products of covariance functions
%   covMask       - mask some dimensions of the data
%   covPref       - difference covariance for preference learning
%   covPER        - make stationary covariance periodic
%   covADD        - additive covariance function
%   covWarp       - warp input to a covariance function
%
% 3) Mahalanobis distance based covariances and their modes
%   covMaha       - generic "mother" covariance
%    * eye        - unit length scale
%    * iso        - isotropic length scale
%    * ard        - automatic relevance determination
%    * prot       - (low-rank) projection in input space
%    * fact       - factor analysis covariance
%    * vlen       - spatially varying length scale
%   covGE         - Gamma exponential covariance
%   covMatern     - Matern covariance function with nu=1/2, 3/2 or 5/2
%   covPP         - piecewise polynomial covariance function (compact support)
%   covRQ         - rational quadratic covariance function
%   covSE         - squared exponential covariance function
%
% 4) Dot product based covariances and their modes
%   covDot        - generic "mother" covariance
%    * eye        - unit length scale
%    * iso        - isotropic length scale
%    * ard        - automatic relevance determination
%    * pro        - (low-rank) projection in input space
%    * fac        - factor analysis covariance
%   covLIN        - linear covariance function
%   covPoly       - polynomial covariance function
%
% 5) Time series covariance functions on the positive real line
%   covFBM        - fractional Brownian motion covariance
%   covULL        - underdamped linear Langevin process covariance
%   covW          - i-times integrated Wiener process covariance
%   covOU         - i-times integrated Ornstein-Uhlenbeck process covariance
%
% 6) Standalone covariances
%   covNNone      - neural network covariance function
%   covLINone     - linear covariance function with bias
%   covPeriodic   - smooth periodic covariance function (1d)
%   covPeriodicNoDC - as above but with zero DC component and properly scaled
%   covCos        - sine periodic covariance function (1d) with unit period
%   covGabor      - Gabor covariance function
%
% 7) Shortcut covariances assembled from library
%   covConst      - covariance for constant functions
%   covNoise      - independent covariance function (i.e. white noise)
%   covPERiso     - make isotropic stationary covariance periodic
%   covPERard     - make ARD stationary covariance periodic
%   covMaterniso  - Matern covariance function with nu=1/2, 3/2 or 5/2
%   covMaternard  - Matern covariance function with nu=1/2, 3/2 or 5/2 with ARD
%   covPPiso      - piecewise polynomial covariance function (compact support)
%   covPPard      - piecewise polynomial covariance function (compact support)
%   covRQiso      - isotropic rational quadratic covariance function
%   covRQard      - rational quadratic covariance function with ARD
%   covSEiso      - isotropic squared exponential covariance function
%   covSEisoU     - same as above but without latent scale
%   covSEard      - squared exponential covariance function with ARD
%   covSEvlen     - spatially varying lengthscale squared exponential
%   covSEproj     - projection squared exponential covariance function
%   covLINiso     - linear covariance function
%   covLINard     - linear covariance function with ARD
%   covGaborard   - Gabor covariance function with ARD
%   covGaborsio   - isotropic Gabor covariance function
%   covSM         - spectral mixture covariance function
%
% 8) Special purpose (approximation) covariance functions
%   apxSparse     - sparse approximation: to be used for large scale inference
%                   problems with inducing points aka FITC
%   apxGrid       - grid interpolation:   to be used for large scale inference
%                   problems with Kronecker/Toeplitz/BTTB covariance matrix
%
% The covariance functions are written according to a special convention where
% the exact behaviour depends on the number of input and output arguments
% passed to the function. If you want to add new covariance functions, you
% should follow this convention if you want them to work with the function gp.
% There are four different ways of calling the covariance functions:
%
% 1) With no (or one) input argument(s):
%
%    s = cov
%
% The covariance function returns a string s telling how many hyperparameters it
% expects, using the convention that "D" is the dimension of the input space.
% For example, calling "covRQard" returns the string '(D+2)'.
%
% 2) With two input arguments and one output argument:
%
%    K = cov(hyp, x) equivalent to K = cov(hyp, x, [])
%
% The function computes and returns the covariance matrix where hyp are
% the hyperparameters and x is an n by D matrix of cases, where
% D is the dimension of the input space. The returned covariance matrix is of
% size n by n.
%
% 3) With three input arguments and one output argument:
%
%    Kz = cov(hyp, x, z)
%    kx = cov(hyp, x, 'diag')
%
% The function computes test set covariances; kx is a vector of self covariances
% for the test cases in x (of length n) and Kz is an (n by nz) matrix of cross
% covariances between training cases x and test cases z.
%
% 4) With two output arguments:
%
%     [K,dK] = cov(hyp, x) equivalent to [K,dK] = cov(hyp, x, [])
%     [K,dK] = cov(hyp, x, z)
%     [K,dK] = cov(hyp, x, 'diag')
%
% The function computes and returns the covariances K as in 3) above.
% In addition to that, the (linear) directional derivative **function** dK is
% returned. The two possible calls dhyp = dK(Q) and [dhyp,dx] = dK(Q) for a
% direction Q of the same size as K are possible. The return arguments dhyp
% and dx are the directional derivatives dhyp = d trace(Q'*K) / d hyp and
% dx = d trace(Q'*K) / d x are of the same size as the hyperparameter
% vector hyp and the input data x, respectively. The components of dhyp and
% dx are defined as follows: dhyp(i) = trace(Q'*( d K / d hyp(i) ))
% and dx(i,j) = trace(Q'*( d K / d x(i,j) )).
%
% Covariance functions can be specified in two ways: either as a string
% containing the name of the covariance function or using a cell array. For
% example:
%
%   cov = 'covRQard';
%   cov = {'covRQard'};
%   cov = {@covRQard};
%
% are supported. Only the second and third form using the cell array can be used
% for specifying composite covariance functions, made up of several
% contributions. For example:
%
%        cov = {'covScale', {'covRQiso'}};
%        cov = {'covSum', {'covRQiso','covSEard','covNoise'}};
%        cov = {'covProd',{'covRQiso','covSEard','covNoise'}};
%        cov = {'covMask',mask,'covSEiso'}
%   q=1; cov = {'covPPiso',q};
%   d=3; cov = {'covPoly',d};
%        cov = {'covADD',{[1,2],'covSEiso'}};
%        cov = {@apxSparse, {@covSEiso}, u};  where u are the inducing inputs
%
% specifies a covariance function which is the sum of three contributions. To
% find out how many hyperparameters this covariance function requires, we do:
%
%   feval(cov{:})
% 
% which returns the string '3+(D+1)+1' (i.e. the 'covRQiso' contribution uses
% 3 parameters, the 'covSEard' uses D+1 and 'covNoise' a single parameter).
%
% See also USAGECOV

% Copyright (c) by Carl Edward Rasmussen and Hannes Nickisch, 2022-04-04.
%                                      File automatically generated using noweb.

help covFunction
