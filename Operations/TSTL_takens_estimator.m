function out = TSTL_takens_estimator(y, Nref, rad, past, embedparams)
% Takens estimator for correlation dimension
% Implements TSTOOL code takens_estimator
% Wrapped up in this form by Ben Fulcher 14/11/2009

%% Check inputs
N = length(y);


% 1) Nref
if nargin < 2 || isempty(Nref) 
    Nref = -1; % use all points
end

% 2) Maximum search radius (as proportion of attractor size)
if nargin < 3 || isempty(rad)
    rad = 0.05;
end

% 3) Theiler window
if nargin < 4 || isempty(past)
    past = 1; % just exclude current point
end
if past>0 && past<1
    past = floor(N*past); % specify a fraction of the time series length...
end

% 4) Embedding parameters
if nargin < 5 || isempty(embedparams)
    embedparams = {'ac','cao'};
    fprintf(1,'Using default time-delay embedding using autocorrelation and cao\n')
else
    if length(embedparams) ~= 2
        error('Embedding parameters are incorrectly formatted, we need {tau,m}')
    end
end

%% Embed the signal
% convert to embedded signal object for TSTOOL

s = BF_embed(y,embedparams{1},embedparams{2},1);

if ~strcmp(class(s),'signal') && isnan(s); % embedding failed
    error('Embedding failed')
    % out = NaN;
    % return
end


%% Run the estimation code
D2 = takens_estimator(s, Nref, rad, past);

out = D2;

end