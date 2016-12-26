function out = EN_ApEn(y,mnom,rth)
% EN_ApEn   Approximate Entropy of a time series
%
% ApEn(m,r).
%
% cf. S. M. Pincus, "Approximate entropy as a measure of system complexity",
% P. Natl. Acad. Sci. USA, 88(6) 2297 (1991)
%
% For more information, cf. http://physionet.org/physiotools/ApEn/
%
%---INPUTS:
% y, the input time series
% mnom, the embedding dimension
% rth, the threshold for judging closeness/similarity
%
%---NOTES:
% No record of where this was code was derived from :-/

% ------------------------------------------------------------------------------
% Copyright (C) 2016, Ben D. Fulcher <ben.d.fulcher@gmail.com>,
% <http://www.benfulcher.com>
%
% If you use this code for your research, please cite:
% B. D. Fulcher, M. A. Little, N. S. Jones, "Highly comparative time-series
% analysis: the empirical structure of time series and their methods",
% J. Roy. Soc. Interface 10(83) 20130048 (2013). DOI: 10.1098/rsif.2013.0048
%
% This function is free software: you can redistribute it and/or modify it under
% the terms of the GNU General Public License as published by the Free Software
% Foundation, either version 3 of the License, or (at your option) any later
% version.
%
% This program is distributed in the hope that it will be useful, but WITHOUT
% ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
% FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
% details.
%
% You should have received a copy of the GNU General Public License along with
% this program. If not, see <http://www.gnu.org/licenses/>.
% ------------------------------------------------------------------------------

% Check inputs, set defaults:
if nargin < 2 || isempty(mnom)
    mnom = 1; % m = 1 (default)
end

if nargin < 3 || isempty(rth)
    rth = 0.2; % r = 0.2 (default)
end

%-------------------------------------------------------------------------------

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
