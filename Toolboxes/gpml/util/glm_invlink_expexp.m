% Compute the log intensity for the inverse link function g(f) = exp(-exp(-f)).
%
% The function is used in GLM likelihoods such as likPoisson, likGamma, likBeta
% and likInvGauss.
%
% Copyright (c) by Hannes Nickisch, 2013-10-16.

function [lg,dlg,d2lg,d3lg] = glm_invlink_expexp(f)
  lg = -exp(-f);
  if nargout>1
    dlg = -lg;
    if nargout>2
      d2lg = lg;
      if nargout>2
        d3lg = -lg;
      end
    end
  end