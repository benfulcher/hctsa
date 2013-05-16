function [c_out, l_out, t_out] = migram(varargin)
%MIGRAM Calculate windowed mutual information between two signals.
%   I = MIGRAM(A,B,MAXLAG,WINDOW,NOVERLAP) calculates the windowed mutual   
%   information between the signals in vector A and vector B. MIGRAM splits  
%   the signals into overlapping segments and forms the columns of I with 
%   their mutual information values up to maximum lag specified by scalar 
%   MAXLAG. Each column of I contains the mutual information function 
%   between the short-term, time-localized signals A and B.  Time increases
%   linearly across the columns of I, from left to right.  Lag increases 
%   linearly down the rows, starting at -MAXLAG. If lengths of A and B 
%   differ, the shorter signal is filled with zeros. If N is the length of 
%   the signals, I is a matrix with 2*MAXLAG+1 rows and 
%         k = fix((N-NOVERLAP)/(WINDOW-NOVERLAP)) 
%   columns.
%
%   I = MIGRAM(A,B,MAXLAG,WINDOW,NOVERLAP,NBINS) calculates the mutual
%   information based on histograms with the number of bins NBINS.
%
%   I = MIGRAM(...,'norm') calculates the renormalised mutual
%   information, which is I/log(NBINS) and ensures a value range [0 1].
%
%   [I,L,T] = MIGRAM(...) returns a column of lag L and one of time T
%   at which the mutual information is computed. L has length equal 
%   to the number of rows of I, T has length k.
%
%   I = MIGRAM(A,B) calculates windowed mutual information using defeault
%   settings; the defeaults are MAXLAG = floor(0.1*N), WINDOW = floor(0.1*N),
%   NOVERLAP = 0 and NBINS = 10. You can tell MIGRAM to use the defeault 
%   for any parameter by leaving it off or using [] for that parameter, e.g. 
%   MIGRAM(A,B,[],1000).
%
%   MIGRAM(A,B) with no output arguments plots the mutual information 
%   using the current figure.
%
%   Remark
%   Please note that the mutual information derived with MI slightly 
%   differs from the results derived with MIGRAM. The reason is that
%   MI also considers estimation errors. 
%
%   Example
%       x = cos(0:.01:10*pi)';
%       y = sin(0:.01:10*pi)' + .5 * randn(length(x),1);
%       migram(x,y)
%
%   See also MI, CORRGRAM.

% Copyright (c) 2008-2009
% Norbert Marwan, Potsdam Institute for Climate Impact Research, Germany
% http://www.pik-potsdam.de
%
% Copyright (c) 2007-2008
% Norbert Marwan, Potsdam University, Germany
% http://www.agnld.uni-potsdam.de
%
% $Date: 2009/03/24 08:30:55 $
% $Revision: 5.7 $


error(nargchk(2,7,nargin))
verbose = 0;

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
nbins = 10;
norm = 0;

i_num = find(cellfun('isclass',varargin,'double'));
i_char = find(cellfun('isclass',varargin,'char'));

if length(i_num) > 2 && ~isempty(varargin{i_num(3)})
    maxlag = varargin{i_num(3)};
    if maxlag < 0, error('Requires positive integer value for maximum lag.'), end
    if length(maxlag) > 1, error('Requires MAXLAG to be a scalar.'), end
end

if length(i_num) > 3 && ~isempty(varargin{i_num(4)})
    window = varargin{i_num(4)};
    if window <= 0, error('Requires positive integer value for window length.'), end
    if length(window) > 1, error('Requires WINDOW to be a scalar.'), end
end

if length(i_num) > 4 && ~isempty(varargin{i_num(5)})
    noverlap = varargin{i_num(5)};
    if noverlap < 0, error('Requires positive integer value for NOVERLAP.'), end
    if length(noverlap) > 1, error('Requires NOVERLAP to be a scalar.'), end
    if noverlap >= window, error('Requires NOVERLAP to be strictly less than the window length.'), end
end

if length(i_num) > 5 && ~isempty(varargin{i_num(6)})
    nbins = varargin{i_num(6)};
    if nbins <= 0, error('Requires positive integer value for NBINS.'), end
    if length(nbins) > 1, error('Requires NBINS to be a scalar.'), end
end

% normalise the result
for i = 1:length(i_char)
    if strcmpi(varargin(i_char(i)), 'norm'), norm = 1; end
end

% prepare time delayed signals
X = buffer(x,maxlag+1,maxlag)';
Y = fliplr(buffer(y,maxlag+1,maxlag)');

% divide the delayed signals into overlapping windows
% and compute the correlation coefficient
cnt = 1;

warning('off','MATLAB:divideByZero')

C = zeros(2*maxlag+1, fix((nx-noverlap)/(window-noverlap)));
if verbose, h = waitbar(0,'Compute mutual information'); end

% -MAXLAG:0
[Yi dummy] = buffer(Y(:,1),window,noverlap,'nodelay');
if exist('accumarray','builtin') == 5
    for i = 1:size(X,2), if verbose, waitbar(i/(2*size(X,2))), end
        [Xi dummy] = buffer(X(:,i),window,noverlap,'nodelay');
        C(i,:) = MI6(Xi, Yi, nbins);
        %cnt = cnt + 1;
    end
else
    for i = 1:size(X,2), if verbose, waitbar(cnt/(2*size(X,2))), end
        [Xi dummy] = buffer(X(:,i),window,noverlap,'nodelay');
        C(cnt,:) = MI5(Xi, Yi, nbins);
        cnt = cnt + 1;
    end
end
cnt = size(X,2)-1;
% 0:MAXLAG
[Xi dummy] = buffer(X(:,end),window,noverlap,'nodelay');
if exist('accumarray','builtin') == 5
    for i = 2:size(Y,2), if verbose, waitbar(cnt/(2*size(X,2))), end
        [Yi dummy] = buffer(Y(:,i),window,noverlap,'nodelay');
        C(i+cnt,:) = MI6(Xi, Yi, nbins);
       % cnt = cnt + 1;
    end
else
    for i = 2:size(Y,2), if verbose, waitbar(cnt/(2*size(X,2))), end
        [Yi dummy] = buffer(Y(:,i),window,noverlap,'nodelay');
        C(cnt,:) = MI5(Xi, Yi, nbins);
        cnt = cnt + 1;
    end
end

if verbose, delete(h), end

warning('on','MATLAB:divideByZero')

% create time scale for the windows
t = (1:nx/size(Xi,2):nx)';
l = (-maxlag:maxlag)';

% if result has to be normalised
if norm
    C = C / log(nbins);
end


% display and output result
if nargout == 0
    newplot
    imagesc(t, l, C)
    xlabel('Time'), ylabel('Lag'), axis xy
    title('Windowed mutual information', 'fontweight', 'bold')
    colorbar
elseif nargout == 1,
    c_out = C;
elseif nargout == 2,
    c_out = C;
    l_out = l;
elseif nargout == 3,
    c_out = C;
    t_out = t;
    l_out = l;
end


% mutual information for Matlab version >= 6
function Z = MI6(x, y, nbins)
    
    % normalise the data and replace the values with integers
    % in the range [1 nbins]
    x = x - repmat(min(x), size(x,1), 1); 
    y = y - repmat(min(y), size(y,1), 1);
    x = x ./ repmat(max(x) + eps, size(x,1), 1); 
    y = y ./ repmat(max(y) + eps, size(y,1), 1);
    
    x = floor(x * nbins) + 1;
    y = floor(y * nbins) + 1;
            
    % compute probabilities
    Z = zeros(1,size(x,2));
    for i = 1:size(x,2)
    
        Pxy = accumarray([x(:,i) y(:,i)] + 1, 1);
        Px = sum(Pxy,1);
        Py = sum(Pxy,2);

        Pxy = Pxy / sum(Pxy(:));
        Px = Px / sum(Px(:));
        Py = Py / sum(Py(:));

        % entropies
        Ix = -sum((Px(Px ~= 0)) .* log(Px(Px ~= 0)));
        Iy = -sum((Py(Py ~= 0)) .* log(Py(Py ~= 0)));
        Ixy = -sum(Pxy(Pxy ~= 0) .* log(Pxy(Pxy ~= 0)));                                                                    
        
        % mutual information
        Z(i) = Ix + Iy - Ixy;
    end    
    
    
% mutual information for Matlab version < 6
function Z = MI5(x, y, nbins)   
    
    % normalise the data and replace the values with integers
    % in the range [1 nbins]
    x = x - repmat(min(x), size(x,1), 1); 
    y = y - repmat(min(y), size(y,1), 1);
    x = x ./ repmat(max(x) + eps, size(x,1), 1); 
    y = y ./ repmat(max(y) + eps, size(y,1), 1);
    
    x = floor(x * nbins) + 1;
    y = floor(y * nbins) + 1;
    
    % compute probabilities
    Z = zeros(1,size(x,2));
    for i = 1:size(x,2)
    
        Pxy = full(sparse(x(:,i) + 1, y(:,i) + 1, 1));
        Px = sum(Pxy,1);
        Py = sum(Pxy,2);

        Pxy = Pxy / sum(Pxy(:));
        Px = Px / sum(Px(:));
        Py = Py / sum(Py(:));

        % entropies
        Ix = -sum((Px(Px ~= 0)) .* log(Px(Px ~= 0)));
        Iy = -sum((Py(Py ~= 0)) .* log(Py(Py ~= 0)));
        Ixy = -sum(Pxy(Pxy ~= 0) .* log(Pxy(Pxy ~= 0)));                                                                    
        
        % mutual information
        Z(i) = Ix + Iy - Ixy;
    end    
    
