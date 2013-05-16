% Implements the Chung-Kennedy sliding window nonlinear step filter. This
% filter is similar to a centred moving average filter of length K, but the
% centre sample in the window is replaced by a weighted sum of forward and
% backward moving average filters. The weights are inversely proportional
% to the one-step-ahead error of the forward/backward filters at predicting
% the samples in the sliding window. The prediction error is estimated over
% 2M samples either side of the centre sample in the window. In this
% implementation, forwards/backwards filters of all lengths up to K are
% used. Larger p gives more sensitivity of the forward/backward filter
% weights to the filter prediction error.
%
% Usage:
%  x = ckfilter(y, K, M, p)
%
% Inputs
%  y - Input signal
%  K - Maximum forward/backward moving average filter length (samples)
%  M - Prediction error analysis window size (samples)
%  p - Positive scaling of prediction error
%
% Outputs
%  x - Step-filtered output signal
%
% (c) Max Little, 2010. The algorithm is described in:
% S.H. Chung, R.A. Kennedy (1991), "Forward-backward non-linear filtering
% technique for extracting small biological signals from noise",
% J. Neurosci. Methods. 40(1):71-86.
% If you use this code for your research, please cite:
% "Steps and bumps: precision extraction of discrete states of molecular
% machines using physically-based, high-throughput time series analysis"
% Max A. Little et al., 2010, arXiv:1004.1234v1 [q-bio.QM]

function x = ckfilter(y, K, M, p)

error(nargchk(4,4,nargin));
y = y(:);

N = length(y);
prior = 1/(2*K);

xf = zeros(N,K);
xb = zeros(N,K);

% Create forward and backwards filters of length L = 1,2..K
L = (1:K)';
for k = (K+1):(N-K-1)
    xf(k,:) = cumsum(y(k+1:k+K))./L;
    xb(k,:) = cumsum(y(k-1:-1:k-K))./L;
end

x = y;

for k = (max(K,M)+1):(N-max(K,M)-1)
    
    yb = repmat(y(k-M+1:k),1,K);
    yf = repmat(y(k:k+M-1),1,K);
    
    f = prior*sum((yb-xf(k-M+1:k,:)).^2).^(-p);
	b = prior*sum((yf-xb(k:k+M-1,:)).^2).^(-p);
   
    if (any(isinf(f)))
        f(isinf(f)) = 1;
        f(~isinf(f)) = 0;
    end
    if (any(isinf(b)))
        f(isinf(b)) = 1;
        f(~isinf(b)) = 0;
    end
    
    fb = sum(f+b);
    f = f/fb;
    b = b/fb;
    
    x(k) = f*xf(k,:)'+b*xb(k,:)';
end
