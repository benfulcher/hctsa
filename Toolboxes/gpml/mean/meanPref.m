function A = meanPref(mean, hyp, x, varargin)

% meanPref - mean function for preference learning.
%
% m(x) = m_0(x1)-m_0(x2), where x1=x(:,1:D), x2=x(:,D+1:2*D), D = size(x,2)/2.
%
% The hyperparameters are:
%
% hyp = [ hyp_m0 ]
%
% This function doesn't actually compute very much on its own, it merely does
% some bookkeeping, and calls another mean function to do the actual work.
%
% See Collaborative Gaussian Processes for Preference Learning, NIPS 2014.
%
% Copyright (c) by Hannes Nickisch and Roman Garnett, 2014-11-01.
%
% See also MEANFUNCTIONS.M and COVPREF.M.

if nargin<3, A = strrep(feval(mean{:}),'D','D/2'); return; end    % no of params

A =   feval(mean{:}, hyp, x(:,1      :end/2), varargin{:}) ...
    - feval(mean{:}, hyp, x(:,1+end/2:end  ), varargin{:});
