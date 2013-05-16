function [c_out, l_out, t_out] = corrgram(varargin)
%CORRGRAM Calculate windowed cross correlation between two signals.
%   C = CORRGRAM(A,B,MAXLAG,WINDOW,NOVERLAP) calculates the windowed cross   
%   correlation between the signals in vector A and vector B. CORRGRAM splits  
%   the signals into overlapping segments and forms the columns of C with 
%   their cross correlation values up to maximum lag specified by scalar 
%   MAXLAG. Each column of C contains the cross correlation function between 
%   the short-term, time-localized signals A and B. Time increases linearly 
%   across the columns of C, from left to right.  Lag increases linearly down 
%   the rows, starting at -MAXLAG. If lengths of A and B differ, the shorter
%   signal is filled with zeros. If N is the length of the signals, C is
%   a matrix with 2*MAXLAG+1 rows and 
%         k = fix((N-NOVERLAP)/(WINDOW-NOVERLAP)) 
%   columns.
%
%   [C,L,T] = CORRGRAM(...) returns a column of lag L and one of time T
%   at which the correlation coefficients are computed. L has length equal 
%   to the number of rows of C, T has length k.
%
%   C = CORRGRAM(A,B) calculates windowed cross correlation using defeault
%   settings; the defeaults are MAXLAG = floor(0.1*N), WINDOW = floor(0.1*N)
%   and NOVERLAP = 0. You can tell CORRGRAM to use the defeault for any 
%   parameter by leaving it off or using [] for that parameter, e.g. 
%   CORRGRAM(A,B,[],1000).
%
%   CORRGRAM(A,B) with no output arguments plots the windowed cross 
%   correlation using the current figure.
%
%   Example
%       x = cos(0:.01:10*pi)';
%       y = sin(0:.01:10*pi)' + .5 * randn(length(x),1);
%       corrgram(x,y)
%
%   See also CORRCOEF, CORR, XCORR, MIGRAM.

% Copyright (c) 2007 by AMRON
% Norbert Marwan, Potsdam University, Germany
% http://www.agnld.uni-potsdam.de
%
% $Date: 2007/12/20 16:26:05 $
% $Revision: 5.2 $


error(nargchk(2,5,nargin))
verbose = 1;

x = varargin{1}; y = varargin{2};
x = x(:); y = y(:);


% check input and inital setting of parameters

nx = length(x); ny = length(y);
if nx < ny    % zero-pad x if it has length less than y
    x(ny) = 0; nx = ny;
end

if ny < nx    % zero-pad y if it has length less than x
    y(nx) = 0;
end

maxlag = floor(nx/10); 
window = floor(nx/10);
noverlap = 0;

if length(varargin) > 2 & ~isempty(varargin{3})
    maxlag = varargin{3};
    if maxlag < 0, error('Requires positive integer value for maximum lag.'), end
    if length(maxlag) > 1, error('Requires MAXLAG to be a scalar.'), end
end

if length(varargin) > 3 & ~isempty(varargin{4})
    window = varargin{4};
    if window < 0, error('Requires positive integer value for window length.'), end
    if length(window) > 1, error('Requires WINDOW to be a scalar.'), end
end

if length(varargin) > 4 & ~isempty(varargin{5})
    noverlap = varargin{5};
    if noverlap < 0, error('Requires positive integer value for NOVERLAP.'), end
    if length(noverlap) > 1, error('Requires NOVERLAP to be a scalar.'), end
    if noverlap >= window, error('Requires NOVERLAP to be strictly less than the window length.'), end
end


% prepare time delayed signals
X = buffer(x,maxlag+1,maxlag)';
Y = fliplr(buffer(y,maxlag+1,maxlag)');

% divide the delayed signals into overlapping windows
% and compute the correlation coefficient
cnt = 1;

warning off
C = zeros(2*maxlag+1, fix((nx-noverlap)/(window-noverlap)));
if verbose, h = waitbar(0,'Compute cross correlation'); end

% -MAXLAG:0
[Yi dummy] = buffer(Y(:,1),window,noverlap,'nodelay'); 
Yi = normalize(Yi);
for i = 1:size(X,2), if verbose, waitbar(cnt/(2*size(X,2))), end
    [Xi dummy] = buffer(X(:,i),window,noverlap,'nodelay');
    Xi = normalize(Xi);
    C(cnt,:) = mean(Xi .* Yi);
    cnt = cnt + 1;
end

% 0:MAXLAG
[Xi dummy] = buffer(X(:,end),window,noverlap,'nodelay');
Xi = normalize(Xi);
for i = 2:size(Y,2), if verbose, waitbar(cnt/(2*size(X,2))), end
    [Yi dummy] = buffer(Y(:,i),window,noverlap,'nodelay');
    Yi = normalize(Yi);
    C(cnt,:) = mean(Xi .* Yi);
    cnt = cnt + 1;
end
if verbose, delete(h), end

warning on

% create time scale for the windows
t = (1:nx/size(Xi,2):nx)';
l = (-maxlag:maxlag)';

% display and output result
if nargout == 0
    newplot
    imagesc(t, l, C)
    xlabel('Time'), ylabel('Lag'), axis xy
    title('Windowed cross correlation', 'fontweight', 'bold')
    colorbar
elseif nargout == 1,
    c_out = C;
elseif nargout == 2,
    c_out = C;
    l_out = l;
elseif nargout == 3,
    c_out = C;
    l_out = l;
    t_out = t;
end


function Y = normalize(X)
    Y = X - repmat(mean(X), size(X,1), 1);
    s = sqrt( sum(conj(Y) .* Y) / (size(Y,1) - 1) );
    Y = Y ./ repmat(s, size(Y,1), 1);
