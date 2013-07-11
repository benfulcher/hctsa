function out = ST_transition(y,discmeth,ng,tau)
% Calculates the transition probabilities given a discretization of the
% data given in y
% INPUTS: discmeth: the method of discretization ('quantile')
%         ng: number of groups in the course-graining
%         tau: the lag; transition matricies corresponding to this lag. We
%         can either downsample the time series at this lag and then do the
%         discretization as normal, or do the discretization and then just
%         look at this dicrete lag. Here we do the former.
% OUTPUTS: a structure containing the probabilities, size as appropriate
%           for dparam
% Ben Fulcher August 2009

if ng < 2; return; end

if strcmp(tau,'ac') % determine tau from first zero of autocorrelation
    tau = CO_fzcac(y);
end

if tau > 1; % calculate transition matrix at non-unity lag
    % downsample at rate 1:tau
    y = resample(y,1,tau);
end

N = length(y); % length of the time series

%% (((1))) Discretize the time series
switch discmeth
    case 'quantile'
        % 1) discretize the time series into a number of groups given by
        %    dparam
        th = quantile(y,linspace(0,1,ng+1)); % thresholds for dividing the time series values
        th(1) = th(1)-1; % this ensures the first point is included
        % turn the time series into a set of numbers from 1:ng
        yth = zeros(N,1);
        for i = 1:ng
            yth(y > th(i) & y <= th(i+1)) = i;
        end
        if any(yth == 0)
            % error -- they should all be assigned to a group
            error('Some time-series values not assigned to a group')
        end
    otherwise
        error('Unknown discritization method ''%s''',discmeth)
end

% Ok, at this stage we should have:
% (*) yth: a thresholded y containing integers from 1:ng (the number of groups)



%% (((2))) find 1-time transition matrix
% probably implemented already, but I'll do it myself
T = zeros(ng);
for i = 1:ng
    ri = find(yth == i);
    if isempty(ri)
        T(i,:) = 0;
    else
        if ri(end) == N; ri = ri(1:end-1); end
        for j = 1:ng
            T(i,j) = length(find(yth(ri+1) == j)); % the next element is off this class
        end
    end
end

% normalize to probability:
T = T/(N-1); % N-1 is appropriate because it's a 1-time transition matrix

%% (((3))) output measures from the transition matrix
% (i) the raw values of the transition matrix
% this has to be done bulkily (only for ng = 2,3):
if ng == 2;% return all elements of T
    for i = 1:4
        eval(sprintf('out.T%u = T(%u);',i,i));
    end
    % out.T1=T(1); out.T2=T(2); out.T3=T(3); out.T4=T(4);
elseif ng == 3;% return all elements of T
    for i = 1:9
        eval(sprintf('out.T%u = T(%u);',i,i));
    end
    % out.T1 = T(1); out.T2 = T(2); out.T3 = T(3); out.T4 = T(4); out.T5 = T(5); out.T6 = T(6);
    % out.T7 = T(7); out.T8 = T(8); out.T9 = T(9);
elseif ng == 4; % return diagonal elements of T
    for i = 1:4
        eval(sprintf('out.TD%u = T(%u,%u);',i,i,i));
    end
    % out.TD1 = T(1,1); out.TD2 = T(2,2); out.TD3 = T(3,3); out.TD4 = T(4,4);
elseif ng == 5; % return diagonal elements of T
    for i = 1:5
        eval(sprintf('out.TD%u = T(%u,%u);',i,i,i));
    end
    % out.TD1 = T(1,1); out.TD2 = T(2,2); out.TD3 = T(3,3); out.TD4 = T(4,4); out.TD5 = T(5,5);
elseif ng == 6; % return diagonal elements of T
    for i = 1:6
        eval(sprintf('out.TD%u = T(%u,%u);',i,i,i));
    end
    % out.TD1 = T(1,1); out.TD2 = T(2,2); out.TD3 = T(3,3); out.TD4 = T(4,4); out.TD5 = T(5,5); out.TD6 = T(6,6);
end
% ng>6, don't return elements of T.

% (ii) measures on the diagonal
out.ondiag = sum(diag(T)); % trace
out.stddiag = std(diag(T)); % std of diagonal elements

% (iii) measures of symmetry:
out.symdiff = sum(sum(abs((T-T')))); % sum of differences of individual elements
out.symsumdiff = sum(sum(tril(T,-1)))-sum(sum(triu(T,+1))); % difference in sums of upper and lower 
                                                          % triangular parts of T

% (iv) measures from covariance matrix:
out.sumdiagcov = sum(diag(cov(T))); % trace of covariance matrix

% (v) measures from eigenvalues of T
eigT = eig(T);
out.stdeig = std(eigT); % std of eigenvalues
out.maxeig = max(real(eigT)); % maximum eigenvalue
out.mineig = min(real(eigT)); % minimum eigenvalue
out.meaneig = mean(real(eigT)); % mean eigenvalue
out.maximeig = max(imag(eigT)); % maximum imaginary part of eigenavlues

% (vi) measures from eigenvalues of covariance matrix:
eigcovT = eig(cov(T));
out.stdeigcov = std(eigcovT); % std of eigenvalues of covariance matrix
out.maxeigcov = max(eigcovT); % max eigenvalue of covariance matrix
out.mineigcov = min(eigcovT); % min eigenvalue of covariance matrix
out.meaneigcov = mean(eigcovT); % mean eigenvalue of covariance matrix

end