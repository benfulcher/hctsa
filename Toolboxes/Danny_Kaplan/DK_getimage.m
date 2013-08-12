function [data2, images] = DK_getimage(data, pred)
% [data2, images] = GETIMAGE(data,pred) finds the scalar images of
% the points in a time series <pred> time sets in the future
% data --- matrix of embedded data (from lagembed)
% pred --- look ahead time, default value 1
% Returns
% data2 --- a new embedded data matrix appropriately trimmed
% images --- the images (at time <pred>) of the points in data2
% This is a convenience program to trim an embedding appropriately.
% Copyright (c) 1996 by D. Kaplan, All Rights Reserved

if nargin < 2
  pred = 1;
end

images = data((1+pred):length(data),1);
data2 = data(1:(length(data)-pred),:);
