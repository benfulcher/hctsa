function out = NL_embed_PCA(y,tau,m)
% Suggestion by Siddarth Arora to embed and look at PCA reduction
% Bioinformatics toolbox (pca code)
% Ben Fulcher 25/2/2010

if nargin < 2 || isempty(tau)
    tau = 'ac'; % embed by first zero-crossing of autocorrelation function
end

if nargin < 3 || isempty(m)
    m = 3; % three-dimensional embedding
end

% embed the signal via time-delay method
y_embed = benembed(y,tau,m,0);

if isnan(y_embed);
    % embedding parameters are unsuitable (likely that tau is too long...)
    out = NaN; return
end

% Do the pca using Bioinformatics toolbox routine 'princomp'
[pc, score, latent] = princomp(y_embed);

% perc=round(latent/sum(latent)*1000)/10; % percentage of variance explained (1 d.p.)
% sum(latent)
perc = latent/sum(latent); % percentage of variance explained

% plot(perc); ylim([0,1]);

%% Get statistics off the eigenvalue distribution
csperc = cumsum(perc);
out.std = std(perc);
out.range = max(perc) - min(perc);
out.max = max(perc);
out.min = min(perc);
out.top2 = sum(perc(1:2)); % variance explained in top two eigendirections
out.nto80 = find(csperc>0.8,1,'first'); % number of eigenvalues you need to reconstruct 80%
out.nto90 = find(csperc>0.9,1,'first'); % number of eigenvalues you need to reconstruct 90%

out.fb01 = find(perc < 0.1,1,'first'); % when perc goes below 0.1 for the first time
if isempty(out.fb01), out.fb01 = length(perc) + 1; end % could make it NaN...

out.fb001 = find(perc < 0.01,1,'first'); % when perc goes below 0.01 for the first time
if isempty(out.fb001), out.fb001 = length(perc)+1; end % could also make it NaN...


end