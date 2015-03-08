% covariance functions to be use by Gaussian process functions. There are two
% different kinds of covariance functions: simple and composite:
%
% simple covariance functions:
%   covConst      - covariance for constant functions
%   covCos        - sine periodic covariance function (1d) with unit period
%   covLIN        - linear covariance function without parameters
%   covLINard     - linear covariance function with ARD
%   covLINiso     - linear covariance function
%   covLINone     - linear covariance function with bias
%   covMaternard  - Matern covariance function with nu=1/2, 3/2 or 5/2 with ARD
%   covMaterniso  - Matern covariance function with nu=1/2, 3/2 or 5/2
%   covNNone      - neural network covariance function
%   covNoise      - independent covariance function (i.e. white noise)
%   covPeriodic   - smooth periodic covariance function (1d)
%   covPeriodicNoDC - as above but with zero DC component and properly scaled
%   covPoly       - polynomial covariance function
%   covPPard      - piecewise polynomial covariance function (compact support)
%   covPPiso      - piecewise polynomial covariance function (compact support)
%   covRQard      - rational quadratic covariance function with ARD
%   covRQiso      - isotropic rational quadratic covariance function
%   covSEard      - squared exponential covariance function with ARD
%   covSEiso      - isotropic squared exponential covariance function
%   covSEisoU     - same as above but without latent scale
%   covSEvlen     - spatially varying lengthscale squared exponential
%   covSEfact     - factor analysis squared exponential covariance function
%   covSM         - spectral mixture covariance function
%   covGaborard   - Gabor covariance function with ARD
%   covGaborsio   - isotropic Gabor covariance function
%   covDiscrete   - precomputed covariance for discrete data
%
% composite (meta) covariance functions (see explanation at the bottom):
%   covScale      - scaled version of a covariance function
%   covProd       - products of covariance functions
%   covSum        - sums of covariance functions
%   covADD        - additive covariance function
%   covMask       - mask some dimensions of the data
%   covPERard     - make ARD stationary covariance periodic
%   covPERiso     - make isotropic stationary covariance periodic
%   covPref       - difference covariance for preference learning
%
% special purpose (wrapper) covariance functions
%   covFITC       - to be used in conjunction with infFITC* for large scale 
%                   inference problems; any covariance can be wrapped by
%                   covFITC such that the FITC approximation is applicable
%   covGrid       - to be used in conjunction with infGrid* for large scale 
%                   inference problems on grids resulting Kronecker structure
%
% Naming convention: all covariance functions are named "cov/cov*.m". A trailing
% "iso" means isotropic, "ard" means Automatic Relevance Determination, and
% "one" means that the distance measure is parameterized by a single parameter.
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
% 2) With two input arguments:
%
%    K = cov(hyp, x) equivalent to K = cov(hyp, x, [])
%
% The function computes and returns the covariance matrix where hyp are
% the hyperparameters and x is an n by D matrix of cases, where
% D is the dimension of the input space. The returned covariance matrix is of
% size n by n.
%
% 3) With three input arguments:
%
%    Ks  = cov(hyp, x, xs)
%    kss = cov(hyp, xs, 'diag')
%
% The function computes test set covariances; kss is a vector of self covariances
% for the test cases in xs (of length ns) and Ks is an (n by ns) matrix of cross
% covariances between training cases x and test cases xs.
%
% 4) With four input arguments:
%
%     dKi   = cov(hyp, x, [], i)
%     dKsi  = cov(hyp, x, xs, i)
%     dkssi = cov(hyp, xs, 'diag', i)
%
% The function computes and returns the partial derivatives of the
% covariance matrices with respect to hyp(i), i.e. with
% respect to the hyperparameter number i.
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
%        cov = {'covMask',{mask,'covSEiso'}}
%   q=1; cov = {'covPPiso',q};
%   d=3; cov = {'covPoly',d};
%        cov = {'covADD',{[1,2],'covSEiso'}};
%        cov = {@covFITC, {@covSEiso}, u};  where u are the inducing inputs
%
% specifies a covariance function which is the sum of three contributions. To 
% find out how many hyperparameters this covariance function requires, we do:
%
%   feval(cov{:})
% 
% which returns the string '3+(D+1)+1' (i.e. the 'covRQiso' contribution uses
% 3 parameters, the 'covSEard' uses D+1 and 'covNoise' a single parameter).
%
% See also doc/usageCov.m.
%
% Copyright (c) by Carl Edward Rasmussen and Hannes Nickisch, 2014-12-08.
%                                      File automatically generated using noweb.
