% Compute the log intensity for the inverse link function (twice logistic)
% g(f) = h(f*(1+a*h(f))), where is the logistic h(f) = log(1+exp(f))).
% Output range: 0 <= g(f).
%
% The function can be used in GLM likelihoods such as likPoisson, likGamma, and
% likInvGauss.
%
% See Seeger et al., Bayesian Intermittent Demand Forecasting for Large
% Inventories, NIPS, 2016.
%
% Copyright (c) by Hannes Nickisch, 2016-10-04.

function [lg,dlg,d2lg,d3lg] = glm_invlink_logistic2(a,f)
  [lh,dlh,d2lh,d3lh] = glm_invlink_logistic(f); h = exp(lh);
  ft = f + a*f.*h;
  dft = 1 + a*h.*(1 + f.*dlh); w = a*h.*(dlh.^2 + d2lh);
  d2ft = 2*a*h.*dlh + f.*w;
  d3ft = 3*w + a*f.*h.*(dlh.^3+3*d2lh.*dlh+d3lh);
  [lgt,dlgt,d2lgt,d3lgt] = glm_invlink_logistic(ft);
  lg = lgt; dlg = dlgt.*dft;
  d2lg = d2lgt.*dft.^2 + dlgt.*d2ft;
  d3lg = d3lgt.*dft.^3 + 3*d2lgt.*dft.*d2ft + dlgt.*d3ft;