function varargout = covPPiso(v,varargin)

% Wrapper for Piecewise Polyonomial covariance function covPP.m.
%
% Piecewise polynomial covariance function with compact support, v = 0,1,2,3.
% The covariance functions are 2v times contin. diff'ble and the corresponding
% processes are hence v times  mean-square diffble. The covariance function is:
%
% k(x,z) = sf^2 * max(1-r,0)^(j+v) * f(r,j) with j = floor(D/2)+v+1
%
% where r is the distance sqrt((x-z)'*inv(P)*(x-z)), P is ell^2 times
% the unit matrix and sf2 is the signal variance. The hyperparameters are:
%
% hyp = [ log(ell)
%         log(sf) ]
%
% Copyright (c) by Carl Edward Rasmussen and Hannes Nickisch, 2016-04-27.
%
% See also covPP.m.

varargout = cell(max(1,nargout),1);
[varargout{:}] = covScale({'covPP','iso',[],v},varargin{:});