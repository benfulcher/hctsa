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

function [lik,likv] = ZG_hmm_cl(X,T,K,Mu,Cov,P,Pi)

p = length(X(1,:));
N = length(X(:,1));
tiny = exp(-700);

if (rem(N,T) ~= 0)
  disp('Error: Data matrix length must be multiple of sequence length T');
  return;
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
  k2 = k1/sqrt(det(Cov));
  for i = 1:T
    for l = 1:K
      d = Mu(l,:)-X((n-1)*T+i,:);
      B(i,l) = k2*exp(-0.5*d*iCov*d');
    end; 
  end; 
  
  scale = zeros(T,1);
  alpha(1,:) = Pi(:)'.*B(1,:);
  scale(1) = sum(alpha(1,:)); 
  alpha(1,:) = alpha(1,:)/(scale(1)+tiny);
  for i = 2:T
    alpha(i,:) = (alpha(i-1,:)*P).*B(i,:); 
    scale(i) = sum(alpha(i,:));
    alpha(i,:) = alpha(i,:)/(scale(i)+tiny);
  end;

  likv(n) = sum(log(scale+(scale == 0)*tiny));
  Scale = Scale+log(scale+(scale == 0)*tiny);
end;

lik = sum(Scale);

end