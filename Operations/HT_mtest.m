function H=mtest(x, alpha)
%   MTEST: Single sample Smirnov-Cramer-Von Mises goodness-of-fit hypothesis test.
%   H = MTEST(X,ALPHA) performs the particular case of Smirnov-Cramer-Von Mises
%   test to determine whether the null hypothesis of composite normality CDF is 
%   a reasonable assumption regarding the population distribution of a random sample X
%   with the desired significance level ALPHA. The Smirnov-Cramer-Von Mises test 
%   is based on interpolation procedure, so the significance level is restricted to 
%   0.001 <= ALPHA <= 0.10.
%
%   H indicates the result of the hypothesis test according to the MATLAB rules 
%   of conditional statements:
%   H=1 => Do not reject the null hypothesis at significance level ALPHA.
%   H=0 => Reject the null hypothesis at significance level ALPHA.
% 
%   Let S(x) be the empirical c.d.f. estimated from the sample vector X, F(x) 
%   be the corresponding true normal population c.d.f., and CDF be a
%   normal c.d.f. with zero mean and unit standard deviation. The 
%   Smirnov-Cramer-Von Mises hypotheses and test statistic in this particular case are:
% 
%   Null Hypothesis:        F(x) is normal with zero mean and unit variance.
%   Alternative Hypothesis: F(x) is not normal with zero mean and unit variance.
%
%   Test Statistic:         W^2 = integral from 0 to 1 of (S(x)-F(x))^2 dF(x)
%
%   The decision to reject the null hypothesis is taken when the test statistic
%   exceeds the critical value. 
%
%   X must be a row vector representing a random sample. ALPHA must be a scalar.
%   The function doesn't check the formats of X and ALPHA, as well as a number of the
%   input and output parameters.
%
%   The asymptotic limit of the Smirnov-Cramer-Von Mises is reached when LENGTH(X)>=3.
%
% Author: G. Levin, May, 2003
%
% References:
%   W. T. Eadie, D. Drijard, F. E. James, M Roos and B. Sadoulet, "Statistical Methods
%   in Experimental Physics", North-Holland, Sec. Reprint, 1982.

%Find the critical value which corresponds to the given ALPHA
N=length(x);
a=[0.001, 0.010, 0.050, 0.100];  %alphas
q=[1.168, 0.743, 0.461, 0.347]/N; %asymptotic quantiles
CV=interp1(a, q, alpha, 'linear'); %the critical value

%Build the S(x), the empirical c.d.f.
I=[1:N];
S=I/N; %order statistics

%Build the Null Hypothesis c.d.f F(x), the normal distribution with zero mean and unit variance.
x=sort(x);
F=0.5*erfc(-x/sqrt(2)); %normal cumulative distribution

% cdfplot(x);
hold on;
% % plot(x,F);

%Calculate the Smirnov-Cramer-Von Mises metric
W2=(1/12/N+sum((F-(2*I-1)/2/N).^2))/N; %see reference

%Decision
H=(W2<=CV);
