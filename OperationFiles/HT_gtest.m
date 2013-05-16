function H=gtest(x, alpha)
%   GTEST: Single sample Geary goodness-of-fit hypothesis test.
%   H=GTEST(X,ALPHA) performs the Geary test to determine whether the null hypothesis
%   of composite normality PDF is a reasonable assumption regarding the population 
%   distribution of a random sample X with the desired significance level ALPHA.
%
%   H indicates the result of the hypothesis test according to the MATLAB rules 
%   of conditional statements:
%   H=1 => Do not reject the null hypothesis at significance level ALPHA.
%   H=0 => Reject the null hypothesis at significance level ALPHA.
% 
%   The Geary's hypotheses and test statistic are:
% 
%   Null Hypothesis:        X is normal with unknown mean and variance.
%   Alternative Hypothesis: X is not normal.
%
%   The test evaluates the "good estimator of STD(X) for normal distribution"/"reasonable
%   estimator of STD(X)" ratio. If X taken from a non-normal distribution this ratio is
%   considerably different from 1.0.
%   The decision to reject the null hypothesis is taken when the P value is less than 
%   significance level ALPHA. 
%
%   X must be a row vector representing a random sample. ALPHA must be a scalar.
%   The function doesn't check the formats of X and ALPHA, as well as a number of the
%   input and output parameters.
%
% Author: G. Levin, May, 2003.
%
% References:
%   R. E. Walpole, R. H. Mayers, S. L. Mayers, K. Ye, K. Yee, "Probability and Statistics 
%   for Engineers and Scientists", 7-th edition, Upper Saddle River, N.J., Prentice Hall, 2002.

N=length(x);
m=mean(x);
u=sqrt(pi/2/N)*sum(abs(x-m))./sqrt(sum((x-m).^2)); %Geary statistic
z=(u-1)*sqrt(N)/.2661; %normalization
p=1-0.5*(erfc(-abs(z)/sqrt(2))-erfc(abs(z)/sqrt(2))); %p-value

H=(p>=alpha);