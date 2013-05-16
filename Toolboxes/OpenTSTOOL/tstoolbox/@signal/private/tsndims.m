function n = tsndims(data)

% function s = tsndims(data)
% Replacement for ndims
% Returns one for a columns vector
% Two for a matrix/row vector

n = length(tssize(data));
