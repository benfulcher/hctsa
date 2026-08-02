function [m,dm] = meanPref(mean, hyp, x)

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
% Copyright (c) by Hannes Nickisch and Roman Garnett, 2016-04-16.
%
% See also meanFunctions.m and cov/covPref.m.

if nargin<3, m = strrep(feval(mean{:}),'D','D/2'); return; end    % no of params

[m1,dm1] = feval(mean{:}, hyp, x(:,1      :end/2));
[m2,dm2] = feval(mean{:}, hyp, x(:,1+end/2:end  ));
m = m1-m2; dm = @(q) dm1(q)-dm2(q);