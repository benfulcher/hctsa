function A = EZ_VisibilityGraph(tsData,Nmax)
% EZ_VisibilityGraph    Enyu Zhuang implementation of visibility graph algorithm
%
%---INPUTS:
% tsData: time series
% Nmax: the maximum length of time series (number of nodes in vis graph)
%
%---OUTPUT:
% A: Adjacent matrix of the corresponding visibility graph
%
%-------------------------------------------------------------------------------
% Original code by Enyu ZHUANG (Zoey) 23/9/13
% Substantial modifications by Ben Fulcher, 26-8-2015
%-------------------------------------------------------------------------------

%-------------------------------------------------------------------------------
% Check inputs:
%-------------------------------------------------------------------------------
if nargin < 2
    Nmax = length(tsData);
end

% Make tsData a column vector
if size(tsData,1) < size(tsData,2)
    tsData = tsData';
end

% Initialize adjacency matrix:
A = zeros(Nmax,Nmax);
% Initialize gradient matrix:
S = zeros(Nmax,Nmax);

%-------------------------------------------------------------------------------
% Run the algorithm
%-------------------------------------------------------------------------------

for i = 1:(Nmax-1)
    % S(i,i+1:end) = (tsData(i+1:end)-tsData(i))./((i+1:Nmax)'-i);
    for j = i+1:Nmax
        % The gradient between the datapoints i,j
        S(i,j) = (tsData(j)-tsData(i))/(j-i);

        if j == i+1 % will always be visible to the next data point
            A(i,j) = 1;
        else
            % Check that the gradient to each intermediate point is lower than
            % the gradient to the point of interest

            % if all(S(i,i+1:j-1) < S(i,j))
                % A(i,j) = 1;
            % end

            % This loop implementation is faster than the vector one above:
            for k = i+1:j-1
                if S(i,k) >= S(i,j)
                    break;
                elseif k == (j-1)
                    A(i,j) = 1;
                end
            end
        end
    end
end

% Symmetrize:
At = A';
lowerT = logical(tril(ones(Nmax,Nmax)));
A(lowerT) = At(lowerT);


%-------------------------------------------------------------------------------
% My method is much slower:
%-------------------------------------------------------------------------------

% A = zeros(Nmax); % adjacency matrix
% for i = 1:Nmax-1
%     % compute all subsequent gradients
%     deltay = tsData(i+1:end) - tsData(i); % vector of deltay's
%     deltat = (1:Nmax-i)'; % time from current reference i
%     m = deltay./deltat; % gradients
%
%     % Compute maximum gradient so far:
%     cumMax = zeros(Nmax-i,1);
%     cumMax(1) = m(1);
%     if Nmax-i > 1
%         for j = 2:Nmax-i
%             cumMax(j) = max([m(j),cumMax(j-1)]);
%         end
%     end
%
%     % Links exist where gradient exceeds maximum so far
%     A(i,i+1:end) = (m >= cumMax);
% end

end
