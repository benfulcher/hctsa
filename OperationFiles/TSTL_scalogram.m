
function out = TSTL_scalogram(y, scalemin, scalemax, scalestep, mlen)
% Uses TSTOOL code 'scalogram' to create a scalogram of signal using morlet
% wavelet. 
% INPUTS:
% none are explained in documentation
% Implemented in this form by Ben Fulcher 14/11/2009

% scalemin
if nargin<2 || isempty(scalemin)
    scalemin = 0.1;
end

% scalemax
if nargin<3 || isempty(scalemax)
    scalemax = 1;
end

% scalestep
if nargin<4 || isempty(scalestep)
    scalestep = 0.1;
end

% mlen
if nargin<5 || isempty(mlen)
   mlen = 10;
end


%% Signal the time series
s = signal(y);
% s = benembed(y,2,2,1);

%% Compute the scalogram using TSTOOL code
rs = scalogram(s);

% view(rs);
rss = data(rs);
keyboard




end