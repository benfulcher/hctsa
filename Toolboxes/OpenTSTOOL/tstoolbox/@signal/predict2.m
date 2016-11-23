function rs = predict2(s, dim, delay, len, nnr, mode, metric_lambda)
     
%tstoolbox/@signal/predict2
%   Syntax:
%     * rs = predict2(s, len, nnr, step, mode)
%
%   Input arguments:
%     * len - length of prediction (number of output values)
%     * nnr - number of nearest neighbors to use (default is one)
%     * step - stepsize (in samples) (default is one)
%     * mode:
%          + 0 = Output vectors are the mean of the images of the nearest
%            neighbors
%          + 1 = Output vectors are the distance weighted mean of the
%            images of the nearest neighbors
%          + 2 = Output vectors are calculated based on the local flow
%            using the mean of the images of the neighbors
%          + 3 = Output vectors are calculated based on the local flow
%            using the weighted mean of the images of the neighbors
%
%   Local constant iterative prediction for phase space data (e.g. data
%   stemming from a time delay reconstruction of a scalar time series),
%   using fast nearest neighbor search. Four methods of computing the
%   prediction output are possible.
%
% Copyright 1997-2001 DPI Goettingen, License http://www.physik3.gwdg.de/tstool/gpl.txt

narginchk(4,6);

if nargin < 5, nnr = 1; end
if nargin < 6, mode = 0; end
if nargin < 7, metric_lambda = 1; end

c = predict2(data(s), dim, delay, len, nnr, mode, metric_lambda); 	% do the prediction in phase space

c = c(:,1);

rs = signal(core(c), s);    
a = getaxis(rs, 1);
rs = setaxis(rs, 1, a);
rs = cut(rs, 1, dlens(s,1)+1, dlens(rs, 1));                    
rs = addhistory(rs,  'Local constant prediction for scalar time series');
rs = addcommandlines(rs, 's = predict2(s', dim, delay, len, nnr, mode);

					 
	
