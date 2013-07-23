% NL_LZcomplexity
% 
% Calculates the Lempel-Ziv complexity of a n-bit encoding of the time
% series using Michael Small's complexity code (renamed MS_complexity here),
% 
% cf. M. Small, Applied Nonlinear Time Series Analysis: Applications in Physics,
% Physiology, and Finance (book) World Scientific, Nonlinear Science Series A,
% Vol. 52 (2005)
% Code available at http://small.eie.polyu.edu.hk/matlab/
% 
% The code is a wrapper for Michael Small's original code and uses the
% associated mex file compiled from complexitybs.c (renamed MS_complexitybs.c
% here).
% 
% INPUTS:
% y, the input time series
% n, the (integer) number of bits to encode the data into
% preproc [opt], first apply a given preprocessing to the time series. For now,
%               just 'diff' is implemented, which zscores incremental
%               differences and then applies the complexity method.
% 
% The function has a single output: the normalized Lempel-Ziv complexity: i.e.,
% the number of distinct symbol sequences in the time series divided by the
% expected number of distinct symbols for a noise sequence.

function out = NL_LZcomplexity(y,n,preproc)
% Ben Fulcher, 19/2/2010

if nargin < 2 || isempty(n)
    n = 2; % n-bit encoding
end
if nargin < 3
    preproc = []; % no preprocessing
end

% Apply some pre-processing to the time series before performing the analysis
if ischar(preproc)
    switch preproc
    case 'diff'
        y = BF_zscore(diff(y));
    otherwise
        error('Unknown preprocessing setting ''%s''', preproc);
    end
end

% Run Michael Small's (mexed) code for calcaulting the Lempel-Ziv complexity:
out = MS_complexity(y,n);

% %_______________________________________________________________
% % cmp = complexity(x,n);
% %
% % calculate the Lempel-Ziv complexity of the n-bit encoding of x. 
% %
% % cmp is the normalised complexity, that is the number of distinct
% % symbol sequences in x, divided by the expected number of distinct 
% % symbols for a noise sequence.
% %
% % Algorithm is implemented in MS_complexitybs.c
% %
% % M. Small
% % ensmall@polyu.edu.hk
% % 7/10/04
% 
% if length(n) > 1
%    for ni = 1:length(n),
%         cmp(ni) = MS_complexitybs(y,n(ni));
%    end
% else
%     %do the binning, with equi-probably bins
%     y = y(:);
%     nx = length(y);
%     [xn, xi] = sort(y);
% %     y = zeros(nx,1);
%     y = 1:nx;
%     y = floor(y.*(n/(nx+1)));
%     y(xi) = y;
% 
%     % compute complexity with MS_complexitybs
%     cmp = MS_complexitybs(y+eps);
% end
% 
% %______________________________________________________________
% % Ok, so we have cmp, the complexity
% out = cmp;

end