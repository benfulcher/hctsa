% Mean functions to be use by Gaussian process functions. There are two
% different kinds of mean functions: simple and composite:
%
% Simple mean functions:
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
%   meanWSPC      - weighted sum of projected cosines
%
% Composite mean functions (see explanation at the bottom):
%
%   meanScale     - scaled version of a mean function
%   meanSum       - sum of mean functions
%   meanProd      - product of mean functions
%   meanPow       - power of a mean function
%   meanMask      - mask some dimensions of the data
%   meanPref      - difference mean for preference learning
%   meanWarp      - warped mean function
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
% 2) With two input arguments and one output argument:
%
%    m = meanNAME(hyp, x)
%
% The function computes and returns the mean vector m with components
% m(i) = m(x(i,:)) where hyp are the hyperparameters and x is an n by D matrix
% of data, where D is the dimension of the input space. The returned mean
% vector m is of size n by 1.
%
% 3) With two input arguments and two output arguments:
%
%    [m, dm] = meanNAME(hyp, x)
%
% The function computes and returns the mean vector m as in 2) above.
% In addition to that, the (linear) directional derivative function dm is
% returned. The call dhyp = dm(q) for a direction vector q of size n by 1
% returns a vector of directional derivatives dhyp = d (q'*m(x)) / d hyp of
% the same size as the hyperparameter vector hyp. The components of dhyp are
% defined as follows: dhyp(i) = q'*( d m(x) / d hyp(i) ).
%
% See also USAGEMEAN

% Copyright (c) by Carl Edward Rasmussen and Hannes Nickisch, 2022-04-04.
%                                      File automatically generated using noweb.

help meanFunctions
