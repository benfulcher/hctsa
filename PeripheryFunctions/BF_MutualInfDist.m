% ------------------------------------------------------------------------------
% BF_MutualInfDist
% ------------------------------------------------------------------------------
% 
% Converts a matrix of mutual informations (mis) into a distance matrix (D) by whatever
% method you select.
% See Dawy et al. (2005) for discussion of 'CL' and 'CR'
% 
%---INPUTS:
% v1, the first input vector
% v2, the second input vector
% r1, the bin-partitioning method for the first input vector, v1
% r2, the bin-partitioning method for the second input vector, v2
% nbins, the number of bins to partition each vector into.
% 
% NB: r1 and r2 can also be two-component vectors, that specify a custom range
%     for binning
% 
%---OUTPUT:
% mi, the mutual information computed between v1 and v2
% 
% ------------------------------------------------------------------------------
% Copyright (C) 2015, Ben D. Fulcher <ben.d.fulcher@gmail.com>,
% <http://www.benfulcher.com>
%
% This function is free software: you can redistribute it and/or modify it under
% the terms of the GNU General Public License as published by the Free Software
% Foundation, either version 3 of the License, or (at your option) any later
% version.
% 
% This program is distributed in the hope that it will be useful, but WITHOUT
% ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
% FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
% details.
% 
% You should have received a copy of the GNU General Public License along with
% this program. If not, see <http://www.gnu.org/licenses/>.
% ------------------------------------------------------------------------------

function D = BF_MutualInfDist(mis,mymethod,doasqf)

% ------------------------------------------------------------------------------
%% Check Inputs
% ------------------------------------------------------------------------------
if nargin < 2 || isempty(mymethod)
    % select method for converting MIs to distances
    mymethod = 'CL'; % 'CL', 'CR', 'mindiff', 'meandiff'
end
if nargin < 3 || isempty(doasqf)
    % return a vector of pairwise distances rather than the full matrix to
    % save memory
    doasqf = 0;
end

% ------------------------------------------------------------------------------
%% Do it
% ------------------------------------------------------------------------------
% difference between MI to itself and MI to this operation
D = zeros(size(mis));
switch mymethod
    case 'mindiff'
        % my naive way
        disp('converting mutual informations to distances by minimum difference')
        for j = 1:length(mis)
            D(j,:) = mis(j,j)-mis(j,:);
        end
        D = min(D,D'); % ** minimum MI difference
        D(D<0) = 0; % any small negative entries (numerical) should be set to zero
    case 'meandiff'
        % see Steuer et al. (2002)
        disp('converting mutual informations to distances by meandifference')
        for j = 1:length(mis)
            D(j,:) = mis(j,j)-mis(j,:); % H(j)-I(i,j)
        end
        D = D + D'; % ** take the sum => H(A) + H(B) - 2I(A,B)
        D(D<0) = 0; % any small negative entries (numerical) should be set to zero
    case 'norm' % normalized
        disp('Converting mutual informations to distances by simple normalization')
        for i = 1:length(mis)
            for j = i:length(mis)
                D(i,j) = mis(i,j)/(mis(i,i)+mis(j,j)-mis(i,j)); % H(j)-I(i,j)
            end
        end
        D = max(D,D'); % lower triangle and diagonal will be zeros since j = i:end
        D = 1-D;
        if any(D(:)<0)
            disp(['There were ' num2str(sum(D(:)<0)) ' negative entries in normalized MI matrix ---'...
                            ' all were set to 0'])
            D(D<0)=0;
        end
        % D should be in [0,1]
    case 'CL'
        % normalize by maximum entropy of both sources
        % this is a proper distance metric
        % Dawy et al. (2005)
        disp('converting mutual informations to distances by ''CL''')
        for j = 1:length(mis)
            D(j,:) = 1-mis(j,:)/mis(j,j);
        end
        D = max(D,D'); % ** normalize by maximum entropy of both sources
        D(D<0)=0; % any small negative entries (numerical) should be set to zero
    case 'CR'
        % normalize by minimum entropy of both sources
        % this is not a proper distance metric but is more 'robust'
        % Dawy et al. (2005)
        disp('converting mutual informations to distances by ''CR''')
        for j = 1:length(mis)
            D(j,:) = 1-mis(j,:)/mis(j,j);
        end
        D = min(D,D'); % ** normalize by the maximum possible MI (minimum individual entropy)
        D(D<0)=0; % any small negative entries (numerical) should be set to zero
end

% convert to a vector
if doasqf
    try
        D = squareform(D);
    catch
        disp('error squareforming D...')
        keyboard
    end
end


end