% Compute the log intensity for the inverse link function g(f) = f^2+a.
% Output range: 0 <= g(f).
%
% The function can be used in GLM likelihoods such as likPoisson, likGamma, and
% likInvGauss.
%
% Copyright (c) by Hannes Nickisch, 2021-02-16.

function [lg,dlg,d2lg,d3lg] = glm_invlink_square(a,f)
  f2a = f.*f+a;
  lg = log(f2a);
  if nargout>1
    dlg = 2*f./f2a;
    if nargout>2
      d2lg = 2*(a-f.*f) ./ f2a.^2;
      if nargout>2
        d3lg = 4*f.*(f.*f-3*a) ./ f2a.^3;
      end
    end
  end