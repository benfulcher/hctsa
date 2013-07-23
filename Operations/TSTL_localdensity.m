% TSTL_localdensity
% 
% Uses TSTOOL code localdensity, which is very poorly documented in the TSTOOL
% documentation, but we can assume it returns local density estimates in the
% time-delay embedding space.
% 
% TSTOOL: http://www.physik3.gwdg.de/tstool/
% 
% INPUTS:
% 
% y, the time series as a column vector
% 
% NNR, number of nearest neighbours to compute
% 
% past, number of time-correlated points to discard (samples)
% 
% embedparams, the embedding parameters, inputs to BF_embed as {tau,m}, where
%               tau and m can be characters specifying a given automatic method
%               of determining tau and/or m (see BF_embed).
% 
% Outputs are various statistics on the local density estimates at each point in
% the time-delay embedding, including the minimum and maximum values, the range,
% the standard deviation, mean, median, and autocorrelation.
% 

function out = TSTL_localdensity(y,NNR,past,embedparams)
% Ben Fulcher, November 2009

%% Check inputs
if nargin < 2 || isempty(NNR)
    NNR = 3; % 3 nearest neighbours
end

if nargin < 3 || isempty(past)
    past = 40;
end

if nargin < 4 || isempty(embedparams)
    embedparams = {'ac','cao'};
    fprintf(1,'Using default embedding using autocorrelation and cao\n')
end

%% Embed the signal
s = BF_embed(y,embedparams{1},embedparams{2},1);

if ~strcmp(class(s),'signal') && isnan(s); % embedding failed
    error('Embedding failed.')
    % out = NaN;
    % return
end

%% Run the code
% try
rs = localdensity(s,NNR,past);
% catch emsg
%     if strcmp(emsg.message,'Fast nearest neighbour searcher : To many neighbors for each query point are requested')
%         out = NaN; return
%     end
% end
%% Convert output to data
locden = data(rs);
if all(locden == 0)
    out = NaN; return
end
% locden is a vector of length equal to the number of points in the
% embedding space (length of time series - m + 1), presumably the local
% at each point

out.minden = min(locden);
out.maxden = max(locden);
out.iqrden = iqr(locden);
out.rangeden = range(locden);
out.stdden = std(locden);
out.meanden = mean(locden);
out.medianden = median(locden);

F_acden = @(x) CO_autocorr(locden,x); % autocorrelation of locden
out.ac1den = F_acden(1);
out.ac2den = F_acden(2);
out.ac3den = F_acden(3);
out.ac4den = F_acden(4);
out.ac5den = F_acden(5);

% Estimates of correlation length:
out.tauacden = CO_fzcac(locden,'ac'); % first zero-crossing of autocorrelation function
out.taumiden = CO_FirstMin(locden,'mi'); % first minimum of automutual information function

end