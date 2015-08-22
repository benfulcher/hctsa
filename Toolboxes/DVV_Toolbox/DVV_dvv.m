function data = DVV_dvv(x, m, Nsub, nd, Ntv)
% Delay Vector Variance method for real and complex signals
%
%
% USAGE: C = dvv (X, m, Nsub, nd, Ntv)
%	X       original real-valued or complex time series
%	m       delay embedding dimension
%	Ntv     number of points on horizontal axes
%	Nsub	number of reference DVs to consider
%	nd      Span over which to perform DVV
%
%
%   A Delay Vector Variance (DVV) toolbox for MATLAB
%   (c) Copyright Danilo P. Mandic 2008
%   http://www.commsp.ee.ic.ac.uk/~mandic/dvv.htm
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   This program is free software; you can redistribute it and/or modify
%   it under the terms of the GNU General Public License as published by
%   the Free Software Foundation; either version 2 of the License, or
%   (at your option) any later version.
%
%   This program is distributed in the hope that it will be useful,
%   but WITHOUT ANY WARRANTY; without even the implied warranty of
%   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%   GNU General Public License for more details.
%
%   You can obtain a copy of the GNU General Public License from
%   http://www.gnu.org/copyleft/gpl.html or by writing to
%   Free Software Foundation, Inc.,675 Mass Ave, Cambridge, MA 02139, USA.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  


% ------------------------------------------------------------------------------
% Check parameters
% ------------------------------------------------------------------------------
if nargin < 1
	error('Not enough input arguments');
end
if nargin < 2 || isempty(m)
	m = 3;
end
if nargin < 3 || isempty(Nsub)
	Nsub = 200;
end
if nargin < 4 || isempty(nd)
	nd = 2.0;
end
if nargin < 5 || isempty(Ntv)
	Ntv = 25*nd;
end
if nargin < 6 || isempty(numSurr)
    numSurr = 10;
end

% ------------------------------------------------------------------------------
% Initial Conditions
% ------------------------------------------------------------------------------
N = length(x);              % Length of input vector
tau = 1;                    % Time delay parameter
d = zeros(N-m*tau, Nsub); 
y = zeros(Ntv,1);

% Make input vector x a column vector
if size(x,2) > size(x,1)
    x = x';
end

% Generate Nsub subset from existing DV's, randomly
temp = randperm (N - m*tau);
ref = temp(1:Nsub) + m*tau;

% Compute pairwise distances between reference DVs and all DVs
count = 0; 
acc = 0;
for i = 1:Nsub
    for j = m*tau+1:N
        d(j-m*tau,i) = norm (x(ref(i)-m*tau:tau:ref(i)-tau) - x(j-m*tau:tau:j-tau));
        if (ref(i) ~= j)
            acc = acc + d(j-m*tau,i); 
            count = count + 1;
        end
    end
end

% Mean and std variation calculation of input data
avg = acc/count;
count = 0;
acc = 0;
for i = 1:Nsub
    for j = m*tau + 1:N
        if (ref(i) ~= j)
            acc = acc + (d(j-m*tau,i)-avg).^2; 
            count = count + 1;
        end
    end
end
variance = sqrt(acc/(count-1));

% Calculates the range vector consisting of Ntv equally spaced regions
n = (1:Ntv)-1;
rd = avg-nd*variance + (2*nd*variance*n)/(Ntv-1);

% Creates sets of DV's, for each ref element of subset and value rd, which have norms closer than distance rd to ref 
for n = 1:length(rd)
    if rd(n) > 0
        tot = 0;
        count = 0;
        for k = 1:Nsub
            IND = find(d(:,k) <= rd(n)) + m*tau;
            IND = IND(IND~=k);
            % Only those variance values are considered for which the corresponding
            % sets have atleast 30 DVs
            if (length(IND) >= 30)
                tot = tot + var(x(IND));                 
                count = count+1;
            end
        end
        if (~count)
            y(n) = NaN;
        else
            y(n) = tot/(count*var(x));
        end 
    else
        y(n) = NaN;
    end    
end

% Horizontal axis 
T = (rd'-avg)/variance;

% DVV Output
data = [T,y];

end