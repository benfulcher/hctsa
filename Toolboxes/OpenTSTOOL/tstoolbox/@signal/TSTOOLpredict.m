function rs = TSTOOLpredict(s, len, nnr, step, mode)
     
%tstoolbox/@signal/predict
%   Syntax:
%     * rs = predict(s, dim, delay, len) => nnr=1
%     * rs = predict(s, dim, delay, len, nnr) => mode=0
%     * rs = predict(s, dim, delay, len, nnr, mode)
%
%   Input arguments:
%     * dim - dimension for time-delay reconstruction
%     * delay - delay time (in samples) for time-delay reconstruction
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
%   Local constant iterative prediction for scalar data, using fast
%   nearest neighbor search. Four methods of computing the prediction
%   output are possible.
%
% Copyright 1997-2001 DPI Goettingen, License http://www.physik3.gwdg.de/tstool/gpl.txt

narginchk(2,5);

if nargin < 3, nnr = 1; end
if nargin < 4, step = 1; end
if nargin < 5, mode = 0; end

N = dlens(s,1); 		% number of input samples

c = core(predict(data(s), len, nnr, step, mode));
rs = signal(c, s);    

a = getaxis(rs, 1);
a = setfirst(a, first(a) + N * delta(a));
rs = setaxis(rs, 1, a);                     
rs = addhistory(rs,  'Local constant prediction');
rs = addcommandlines(rs, 's = predict(s', len, nnr, step, mode);	 

					 
					 
