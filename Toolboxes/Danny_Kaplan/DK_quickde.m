function de = DK_quickde(ts,dim,lag,nmin)
% Does a quick-and-dirty characterization of determinism using de
% ts -- the time series
% dim -- the embedding dimension
% lag -- the embedding lag
% nmin -- optional: number of points to use for delta-eps fitting
%           default value: 500
% Copyright (c) 1996 by D. Kaplan, All Rights Reserved
% Tweaked ever so slightly by B. D. Fulcher

if nargin < 4
    nmin = 500;
end

% Embed the data
xx = DK_lagembed(ts,dim,lag);
% Get the 'forecasting' pre-image and image
[pre, post] = DK_getimage(xx,1);
% Find a 'typical' distance to use for delta-epsilon
% We want to have at least nmin points, so we pick

perc = (2*nmin)/(length(ts).^2);
deltamax = DK_disttyp(pre,perc);

[delta, epsilon] = DK_deltaeps(pre,post);
[de, b] = DK_defit(delta,epsilon,deltamax);

end