% function [lik,likv] = hmm_cl(X,T,K,Mu,Cov,P,Pi);
% 
% Calculate Likelihood for Hidden Markov Model 
%
% X - N x p data matrix
% T - length of each sequence (N must evenly divide by T, default T=N)
% K - number of states
% Mu - mean vectors
% Cov - output covariance matrix (full, tied across states)
% P - state transition matrix
% Pi - priors
%
% lik - log likelihood of X 
% likv - vector of log likelihoods of each sequence
%
% If 0 or 1 output arguments requested, lik is returned. If 2 output
% arguments requested, [lik likv] is returned.
% 
% Modified from the original: emission probabilities are rescaled in the log domain
% at each time step so that observations far from every state do not underflow.
%
% Machine Learning Toolbox
% Version 1.0  01-Apr-96
% Copyright (c) by Zoubin Ghahramani
% http://mlg.eng.cam.ac.uk/zoubin/software.html
%
% ------------------------------------------------------------------------------
% The MIT License (MIT)
% 
% Copyright (c) 1996, Zoubin Ghahramani
% 
% Permission is hereby granted, free of charge, to any person obtaining a copy
% of this software and associated documentation files (the "Software"), to deal
% in the Software without restriction, including without limitation the rights
% to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
% copies of the Software, and to permit persons to whom the Software is
% furnished to do so, subject to the following conditions:
% 
% The above copyright notice and this permission notice shall be included in
% all copies or substantial portions of the Software.
% 
% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
% IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
% FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
% AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
% LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
% OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
% THE SOFTWARE.
% ------------------------------------------------------------------------------

function [lik,likv] = ZG_hmm_cl(X,T,K,Mu,Cov,P,Pi)

p = length(X(1,:));
N = length(X(:,1));
tiny = exp(-700);

if (rem(N,T) ~= 0)
  error('Data matrix length must be multiple of sequence length T');
end;
N = N/T;

alpha = zeros(T,K);
B = zeros(T,K);      % P( output | s_i) 
	
k1 = (2*pi)^(-p/2);

Scale = zeros(T,1);

likv = zeros(1,N);

for n = 1:N
  
  B = zeros(T,K); 
  iCov = inv(Cov);      
  logk2 = log(k1)-0.5*log(det(Cov));
  logB = zeros(T,K);
  for l = 1:K
    d = X((n-1)*T+1:n*T,:)-ones(T,1)*Mu(l,:);
    logB(:,l) = logk2-0.5*sum((d*iCov).*d,2);
  end;
  % Rescale each row so that its largest emission probability is 1 (no underflow);
  % the scale factors are added back to the log likelihood
  shift = max(logB,[],2);
  B = exp(logB-shift*ones(1,K));
  
  scale = zeros(T,1);
  alpha(1,:) = Pi(:)'.*B(1,:);
  scale(1) = sum(alpha(1,:)); 
  alpha(1,:) = alpha(1,:)/(scale(1)+tiny);
  for i = 2:T
    alpha(i,:) = (alpha(i-1,:)*P).*B(i,:); 
    scale(i) = sum(alpha(i,:));
    alpha(i,:) = alpha(i,:)/(scale(i)+tiny);
  end;

  likv(n) = sum(log(scale+(scale == 0)*tiny))+sum(shift);
  Scale = Scale+log(scale+(scale == 0)*tiny)+shift;
end;

lik = sum(Scale);

end