% Compute the log intensity for the inverse link function g(f) = exp(f).
% Output range: 0 <= g(f).
%
% The function can be used in GLM likelihoods such as likPoisson, likGamma, and
% likInvGauss.
%
% Copyright (c) by Hannes Nickisch, 2016-10-04.

function [lg,dlg,d2lg,d3lg] = glm_invlink_exp(f)
  lg = f;
  if nargout>1
    dlg = ones(size(f));
    if nargout>2
      d2lg = zeros(size(f));
      if nargout>2
        d3lg = zeros(size(f));
      end
    end
  end