% Compute the log intensity for the inverse link function g(f) = log(1+exp(f))).
% Output range: 0 <= g(f).
%
% The function can be used in GLM likelihoods such as likPoisson, likGamma, and
% likInvGauss.
%
% Copyright (c) by Hannes Nickisch, 2016-10-04.

function [lg,dlg,d2lg,d3lg] = glm_invlink_logistic(f)
  l1pef = max(0,f) + log(1+exp(-abs(f)));         % safely compute log(1+exp(f))
  lg = log(l1pef); id = f<-15; lg(id) = f(id);   % fix log(log(1+exp(f))) limits
  if nargout>1
    sm = 1./(1+exp(-f));
    dlg = sm./l1pef; dlg(f<-15) = 1;
    if nargout>2
      sp = 1./(1+exp(f));
      d2lg = dlg.*(sp-dlg);
      if nargout>2
        d3lg = d2lg.*(sp-2*dlg) - dlg.*sp.*sm;
      end
    end
  end