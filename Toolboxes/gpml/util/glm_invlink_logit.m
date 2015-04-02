% Compute the log intensity for the inverse link function g(f) = 1/(1+exp(-f)).
%
% The function is used in GLM likelihoods such as likPoisson, likGamma, likBeta
% and likInvGauss.
%
% Copyright (c) by Hannes Nickisch, 2013-10-16.

function varargout = glm_invlink_logit(f)
  varargout = cell(nargout, 1);  % allocate the right number of output arguments
  [varargout{:}] = glm_invlink_logistic(f);
  if nargout>0
    elg = exp(varargout{1});
    varargout{1} = f - elg;
    if nargout>1
      dlg = varargout{2};
      varargout{2} = 1 - elg.*dlg;
      if nargout>2
        d2lg = varargout{3};
        varargout{3} = -elg.*(dlg.^2+d2lg);
        if nargout>3
          d3lg = varargout{4};
          varargout{4} = -elg.*(dlg.^3+3*d2lg.*dlg+d3lg);
        end
      end
    end
  end