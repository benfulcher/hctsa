function out = WL_dwtcoeff(y,wname,level)
% Ben Fulcher January 2010

%% Check Inputs
N = length(y); % length of signal

if nargin < 2 || isempty(wname)
    wname = 'db3'; % Daubechies wavelet filter
end
if nargin < 3 || isempty(level)
    level = 3; % level of wavelet decomposition
end
if strcmp(level,'max')
    level = wmaxlev(N,wname);
end

maxlevelallowed = wmaxlev(N,wname);
if maxlevelallowed < level
    disp('Chosen level is too large for this wavelet on this signal. Will be returning some NaNs, methinks...');
end

%% Perform Wavelet Decomposition
% Computes the following:
%   (*) Wavelet decomposition vector c
%   (*) Bookkeeping vector l

if maxlevelallowed < level
    [c, l] = wavedec(y, maxlevelallowed, wname);
else
    [c, l] = wavedec(y, level, wname);
end

%% Expand DWT coefficients for visualization
% nbcol = 64; % color discretization steps
% 
% cfd = zeros(level,N); % detail coefficients
% for k = 1:level
%     d = detcoef(c,l,k);
%     d = d(:)';
%     d = d(ones(1,2^k),:);
%     cfd(k,:) = wkeep1(d(:)',N);
% end
% 
% cfd =  cfd(:);
% I = find(abs(cfd)<sqrt(eps));
% cfd(I) = zeros(size(I));
% cfd = reshape(cfd,level,N);
% cfd = wcodemat(cfd,nbcol,'row');

%% Do the plotting
% colormap(pink(nbcol));
% image(cfd);
% tics = 1:level;
% labs = int2str((1:level)');
% set(gca,'YTicklabelMode','manual','Ydir','normal', 'Box','On','Ytick',tics,'YTickLabel',labs);
% title('Discrete Wavelet Transform, Absolute Coefficients.');
% xlabel('Time (or Space)')
% ylabel('Level');

%% Get statistics on coefficients
for k = 1:level
    if k <= maxlevelallowed
        d = detcoef(c,l,k); % detail coefficients at level k
        % maximum coefficient at this level
        maxd = max(d);
        eval(sprintf('out.maxd_l%u = maxd;',k));
        % minimum coefficient at this level
        mind = min(d);
        eval(sprintf('out.mind_l%u = mind;',k));
        % std coefficients at this level
        stdd = std(d);
        eval(sprintf('out.stdd_l%u = stdd;',k));
        % 1-D noise coefficient estimate
        stddd = wnoisest(c,l,k);
        eval(sprintf('out.stddd_l%u = stddd;',k));
    else
        eval(sprintf('out.maxd_l%u = NaN;',k));
        eval(sprintf('out.mind_l%u = NaN;',k));
        eval(sprintf('out.stdd_l%u = NaN;',k));
        eval(sprintf('out.stddd_l%u = NaN;',k));
    end
end


% %% Compress Signal
% % Set approximation coefficients to zero
% % nc = wthcoef('a',c,l);
% NC = wthcoef('t',c,l,N,T,SORH);
% 
% 
% %% Single Level Reconstruction
% X = waverec(c,l,wname);
% plot(X);
% keyboard


%% Extract Approximation Coefficients from wavelet decomposition structure

% CA = appcoef(C,L,wname,level);

end