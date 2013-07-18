function out = NL_LZcomplexity(y,n,preproc)
% Wrapper code around Michael Small's original 'complexity.m' code. Uses
% the mex file MS_complexitybs.
% http://small.eie.polyu.edu.hk/matlab/
% Ben Fulcher 19/2/2010

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