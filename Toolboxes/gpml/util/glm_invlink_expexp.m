% Compute the log intensity for the inverse link function g(f) = exp(-exp(-f)).
% Output range: 0 <= g(f) <= 1.
%
% The function can be used in GLM likelihoods such as likBeta.
%
% Copyright (c) by Hannes Nickisch, 2016-10-04.

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