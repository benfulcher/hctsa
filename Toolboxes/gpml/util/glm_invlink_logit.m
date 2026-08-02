% Compute the log intensity for the inverse link function g(f) = 1/(1+exp(-f)).
% Output range: 0 <= g(f) <= 1.
%
% The function can be used in GLM likelihoods such as likBeta.
%
% Copyright (c) by Hannes Nickisch, 2016-10-04.

function varargout = glm_invlink_logit(f)
  varargout = cell(nargout, 1);  % allocate the right number of output arguments
  [varargout{:}] = glm_invlink_logistic(f);
  if nargout>0
    elg = exp(varargout{1});
    varargout{1} = min(f-elg,0);                            % upper bound g by 1
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