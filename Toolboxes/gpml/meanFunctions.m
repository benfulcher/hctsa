% mean functions to be use by Gaussian process functions. There are two
% different kinds of mean functions: simple and composite:
%
% simple mean functions:
%
%   meanZero      - zero mean function
%   meanOne       - one mean function
%   meanConst     - constant mean function
%   meanLinear    - linear mean function
%   meanPoly      - polynomial mean function
%   meanDiscrete  - precomputed mean for discrete data
%   meanGP        - predictive mean of another GP
%   meanGPexact   - predictive mean of a regression GP
%   meanNN        - nearest neighbor mean function
%
% composite covariance functions (see explanation at the bottom):
%
%   meanScale     - scaled version of a mean function
%   meanPow       - power of a mean function
%   meanProd      - products of mean functions
%   meanSum       - sums of mean functions
%   meanMask      - mask some dimensions of the data
%   meanPref      - difference mean for preference learning
%
% Naming convention: all mean functions are named "mean/mean*.m".
%
%
% 1) With no or only a single input argument:
%
%    s = meanNAME  or  s = meanNAME(hyp)
%
% The mean function returns a string s telling how many hyperparameters hyp it
% expects, using the convention that "D" is the dimension of the input space.
% For example, calling "meanLinear" returns the string 'D'.
%
% 2) With two input arguments:
%
%    m = meanNAME(hyp, x) 
%
% The function computes and returns the mean vector where hyp are the 
% hyperparameters and x is an n by D matrix of cases, where D is the dimension
% of the input space. The returned mean vector is of size n by 1.
%
% 3) With three input arguments:
%
%    dm = meanNAME(hyp, x, i)
%
% The function computes and returns the n by 1 vector of partial derivatives
% of the mean vector w.r.t. hyp(i) i.e. hyperparameter number i.
%
% See also doc/usageMean.m.
%
% Copyright (c) by Carl Edward Rasmussen and Hannes Nickisch, 2014-12-08.
%                                      File automatically generated using noweb.
