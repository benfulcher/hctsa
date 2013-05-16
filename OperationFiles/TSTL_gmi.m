function out = TSTL_gmi(y,D,eps,NNR,len,Nref)
% Uses TSTOOL code gmi
% Computes the generalized mutual information function for a scalar time
% series
% Documentation in TSTOOL is completely lacking for this function.
% Ben Fulcher November 2009

%% Preliminaries
% (1) D, embeds in this dimension at time lag 1
if nargin<2 || isempty(D)
    D=1;
end

% (2) eps
if nargin<3 || isempty(eps)
    eps=0.1;
end

% (3) Number of nearest neighbours to use, NNR
if nargin<4 || isempty(NNR)
    NNr = 5;
end

% (4) length, len
if nargin<5 || isempty(len)
    len=10;
end

% (5) Number of randomly-chosen reference points, Nref
if nargin<6 || isempty(Nref)
    Nref=100;
end
if Nref<1 && Nref>0
    Nref = round(N*Nref); % specify a proportion of time series length
end

N = length(y); % length of time series
s = signal(y);

%% Run
gmi(s,D,eps,NNR,len,Nref)

% view(rs);
keyboard

% Dq = data(rs);
% q = spacing(rs);
% plot(q,dq);
% keyboard

%% Get output stats




end