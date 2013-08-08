% EN_ApEn
% 
% Estimates the Approximate Entropy of the time series, ApEn(m,r).
% 
% cf. "Approximate entropy as a measure of system complexity" Pincus, S. M., P.
% Natl. Acad. Sci. USA, 88(6) 2297 (1991)
%
% cf. for more information: http://physionet.org/physiotools/ApEn/
% 
% INPUTS:
% y, the input time series
% mnom, the embedding dimension
% rth, the threshold for judging closeness/similarity
%

function out = EN_ApEn(y,mnom,rth)
% Ben Fulcher cannot remember whether he wrote or obtained this code from
% somewhere. (Cannot find an original)

% Check inputs, set defaults:
if nargin < 2 || isempty(mnom)
    mnom = 1; % m = 1 (default)
end

if nargin < 3 || isempty(rth)
    rth = 0.2; % r = 0.2 (default)
end

r = rth*std(y); % threshold of similarity
N = length(y);% length of time series
phi = zeros(2,1);% phi(1)=phi_m, phi(2)=phi_{m+1}

for k = 1:2
    m = mnom+k-1; % pattern length
    C = zeros(N-m+1,1);
    % Define the matrix x, containing subsequences of u
    x = zeros(N-m+1,m);
    
    % Form vector sequences x from the time series y
    for i = 1:N-m+1
        x(i,:) = y(i:i+m-1);
    end
    
    ax = ones(N-m+1,m);
    for i = 1:N-m+1
        for j = 1:m
            ax(:,j) = x(i,j);
        end
        d = abs(x-ax);
        if m > 1 % Takes maximum distance
            d = max(d,[],2)';
        end
        dr = (d<=r);
        C(i) = sum(dr)/(N-m+1); % Number of x(j) within r of x(i)
    end
    phi(k) = mean(log(C));
end
out = phi(1)-phi(2);

end