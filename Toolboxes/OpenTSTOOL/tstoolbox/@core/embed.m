function cout = embed(cin, dim, delay, shift, windowtype)

%tstoolbox/@core/embed
%   Syntax:
%     * cout = embed(cin, dim, delay, shift, windowtype)
%
%   Input Arguments:
%     * cin - core object
%     * dim - embed dimension
%     * delay - delay time in samples for time delay vectors
%     * shift - shift in samples for two sequent time delay vectors
%     * windowtype - type of window
%
%   Create time delay vectors with dimension dim, delay is measured in
%   samples
%   The input must be a scalar time series
%   The result is a n by dim array, each row contains the coordinates of
%   one point
%
% Copyright 1997-2001 DPI Goettingen, License http://www.physik3.gwdg.de/tstool/gpl.txt


N = dlens(cin,1);
M = floor((N-1-(dim-1)*delay)/shift)+1;

if M < 1
  error('time series to short for chosen embedding parameters')
end

d = zeros(M, dim);
len = (M-1)*shift+1;

for i=1:dim
	start = (i-1)*delay;
	d(:,i) = cin.data(start+1:shift:start+len);
end

if ~strcmp(windowtype, 'Rect')
	d = d .* repmat(window(dim, windowtype)', M, 1);
end

cout = core(d);

% l = dlens(cin,1);	% length in first dimension
% n = l - delay*(dim-1)-1;   % Berechnung der Anzahl Delayvektoren
% 
% if n < 1
% 	error('time series to short for chosen embedding parameters')
% end
% 
% d = cin.data(end-n:end);
% 
% for i=2:dim
% 	k = (i-1)*delay;
% 	d = [d cin.data(end-n-k:end-k)];
% end
% 
% cout = core(d);
