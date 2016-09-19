function [e, p, A, B] = PN_sampenc(y,M,r,justM)
% PN_sampenc    Sample Entropy
%
%---INPUTS:
%
% y input time-series data
% M maximum template length (embedding dimension)
% r matching tolerance level
%
%---OUTPUTS
%
% e sample entropy estimates for m=0,1,...,M-1
% p probabilities
% A number of matches for m=1,...,M
% B number of matches for m=1,...,M excluding last point
%
% NOTE: The mexed implementation is much faster than this Matlab implementation
%  (sampen_mex.c)

%-------------------------------------------------------------------------------
% Modified by Ben Fulcher, from original code sampenc.m from
% http://physionet.org/physiotools/sampen/
% http://www.physionet.org/physiotools/sampen/matlab/1.1/sampenc.m
% Code by DK Lake (dlake@virginia.edu), JR Moorman and Cao Hanqing.
%
% BF added input checking, and altered a few minor things, including
% how the outputs are presented.
% BF also added input 'justM', to give e just for the given M, and not for
% all other m up to it
%
% This PhysioNet code, and the minor modifications to it, is shared under the
% GNU General Public License, cf.
% http://www.physionet.org/faq.shtml#license
%-------------------------------------------------------------------------------

%-------------------------------------------------------------------------------
% Check inputs
%-------------------------------------------------------------------------------
% Embedding dimension (template length)
if nargin < 2 || isempty(M)
    M = 1;
end
if nargin < 3 || isempty(r)
    r = 0.1*std(y);
end
if nargin < 4 || isempty(justM)
    justM = 0;
end

%-------------------------------------------------------------------------------
% Initialize variables:
%-------------------------------------------------------------------------------
N = length(y);
lastrun = zeros(1,N);
run = zeros(1,N);
A = zeros(M,1);
B = zeros(M,1);
p = zeros(M,1);
e = zeros(M,1);

%-------------------------------------------------------------------------------
% Get counting:
%-------------------------------------------------------------------------------
for i = 1:(N-1) % go through each point in the time series, counting matches
    y1 = y(i);
    for jj = 1:N-i % compare to points through the rest of the time series
        % Compare to future index, j:
        j = i + jj;
        % This future point, j, matches the time-series value at i:
        if abs(y(j)-y1) < r
            run(jj) = lastrun(jj) + 1; % increase run count for this lag
            M1 = min(M,run(jj));
            for m = 1:M1
                A(m) = A(m) + 1;
                if j < N
                    B(m) = B(m) + 1;
                end
            end
        else
            run(jj) = 0;
        end
    end
    for j = 1:N-i
        lastrun(j) = run(j);
    end
end

%-------------------------------------------------------------------------------
% Calculate for m = 1
NN = N*(N-1)/2;
p(1) = A(1)/NN;
e(1) = -log(p(1));

% Calculate for m > 1, up to M
for m = 2:M
    p(m) = A(m)/B(m-1);
    e(m) = -log(p(m));
end

% output vector e for m = 1, ..., M
% output vector p for m = 1, ..., M

%-------------------------------------------------------------------------------
% Flag to output the entropy and probability just at the maximum requested m
%-------------------------------------------------------------------------------
if justM
    e = e(end);
    p = p(end); % (just in case)
end

end
