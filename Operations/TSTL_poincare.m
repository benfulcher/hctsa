function out = TSTL_poincare(y,ref,embedparams)
% Uses TSTOOL code poincare
% looks at the poincare section and then tries to quantify structure in it
% ref: can be an absolute number (2 takes the second point in the (embedded)
% time series) or a string like 'max' or 'min' that takes the first
% maximum, ... in the (scalar) time series, ...
% embedparams: the usual thing to give BF_embed for the time-delay
% embedding. The second argument is 3 -- i.e., we embed in a 3 dimensional
% space so that the Poincare section is 2-dimensional
% Another thing that could be cool is to look for variation in the plots
% as ref changes... (not done here)
% Ben Fulcher November 2009


%% Check inputs
N = length(y); % length of the time series

if nargin < 2 || isempty(ref)
    ref = 'max'; % reference point is the first maximum of the time series
end

if nargin < 3 || isempty(embedparams)
    embedparams = {'mi', 3};
    fprintf(1,'Using default embedding settings: minimum of the automutual information for tau and m = 3')
end

if embedparams{2}~=3
    embedparams{2} = 3;
    fprintf(1,'Three-dimensional embedding')
end

% set ref
if ischar(ref)
    switch ref
        case 'max'
            % first local maximum
            dydt = diff(y);
            ref = find(dydt(1:end-1)>=0 & dydt(2:end)<0,1,'first')+1;
        case 'min'
            % first local minimum
            dydt = diff(y);
            ref = find(dydt(1:end-1)<=0 & dydt(2:end)>0,1,'first')+1;
        otherwise
            error('Unknown reference setting ''%s''',ref);
    end
end

if isempty(ref), out = NaN; return; end % ridiculous
if ref < 2, ref = 2; end % gives an error, uses previous value in algorithm

%% Do your magic, TSTOOL!:
% time-delay embed the signal:
s = BF_embed(y,embedparams{1},3,1);

if ~strcmp(class(s),'signal') && isnan(s); % embedding failed
    error('Embedding failed');
end

% run TSTOOL code:
try
    rs = poincare(s,ref);
catch me
    if strcmp(me.message,'No section points found') ...
            || strcmp(me.identifier,'MATLAB:badsubscript')
        fprintf(1,'No section points found to run TSTL_poincare\n');
        out = NaN; return
    else
        error(me.message)
    end
end

% convert to MATLAB forms
v = data(rs); % vectors on poincare surface
NN = length(v);
% labeling poincare surface plane x-y
x = (v(:,1));
y = (v(:,2)); 
% plot(x,y,'.k'); axis equal

%% Get statistics out

% basic stats
out.pcross = NN/N; % proportion of time series that crosses poincare surface

out.maxx = max(x);
out.minx = min(x);
out.stdx = std(x);
out.iqrx = iqr(x);
out.meanx = mean(x);
out.ac1x = CO_autocorr(x,1);
out.ac2x = CO_autocorr(x,2);
out.tauacx = CO_fzcac(x);

out.maxy = max(y);
out.miny = min(y);
out.stdy = std(y);
out.iqry = iqr(y);
out.meany = mean(y);
out.ac1y = CO_autocorr(y,1);
out.ac2y = CO_autocorr(y,2);
out.tauacy = CO_fzcac(y);

out.boxarea = range(x)*range(y);



% statistics on distance between adjacent points, ds
vdiff = v(2:end,:)-v(1:end-1,:);
ds = sqrt(vdiff(:,1).^2 + vdiff(:,2).^2);

% probability that next point in series is within radius r of current point
% in the poincare section
out.pwithinr01 = sum(ds<0.1)/(NN-1);
out.pwithin02 = sum(ds<0.2)/(NN-1);
out.pwithin03 = sum(ds<0.3)/(NN-1);
out.pwithin05 = sum(ds<0.5)/(NN-1);
out.pwithin1 = sum(ds<1)/(NN-1);
out.pwithin2 = sum(ds<2)/(NN-1);
out.meands = mean(ds);
out.maxds = max(ds);
out.minds = min(ds);
out.iqrds = iqr(ds);



% now normalize both axes and look for structure in the cloud of points
% don't normalize for standard deviation -- this probably reveals some
% structure...? But location is already noted.
x = x - mean(x);
y = y - mean(y);

% Statistics on distance on Poincare surface from (mean,mean)
D = sqrt(x.^2+y.^2); % distance from (mean,mean)
out.maxD = max(D);
out.minD = min(D);
out.stdD = std(D);
out.iqrD = iqr(D);
out.meanD = mean(D);
out.ac1D = CO_autocorr(D,1);
out.ac2D = CO_autocorr(D,2);
out.tauacD = CO_fzcac(D);


% entropy of boxed distribution
boxcounts = subcountboxes(x,y,10);% 10 partitions per axis
pbox = boxcounts/NN;


out.maxpbox10 = max(pbox(:));
out.minpbox10 = min(pbox(:));
out.zerospbox10 = sum(pbox(:) == 0);
out.meanpbox10 = mean(pbox(:));
out.rangepbox10 = range(pbox(:));
out.hboxcounts10 = -sum(pbox(pbox > 0).*log(pbox(pbox > 0)));
out.tracepbox10 = sum(diag(pbox)); % trace

boxcounts = subcountboxes(x,y,5);% 5 partitions per axis
pbox = boxcounts/NN;
out.maxpbox5 = max(pbox(:));
out.minpbox5 = min(pbox(:));
out.zerospbox5 = sum(pbox(:) == 0);
out.meanpbox5 = mean(pbox(:));
out.rangepbox5 = range(pbox(:));
out.hboxcounts5 = -sum(pbox(pbox>0).*log(pbox(pbox>0)));
out.tracepbox5 = sum(diag(pbox));


% can imagine doing many more things; like seeing different slabs of the
% space, or finding the line whos vicinity includes many points, etc. but I
% think this is enough for now.

    function boxcounts = subcountboxes(x,y,nbox)
        boxcounts = zeros(nbox);
        % boxes are quantiles along each axis
        xbox = quantile(x,linspace(0,1,nbox+1));
        ybox = quantile(y,linspace(0,1,nbox+1));
        xbox(end) = xbox(end)+1;
        ybox(end) = ybox(end)+1;
        
        for i = 1:nbox % x
            rx = (x >= xbox(i) & x < xbox(i+1)); % these x are in range
            for j = 1:nbox % y
                % only need to look at those ys for which the xs are in range
                boxcounts(i,j) = sum(y(rx) >= ybox(j) & y(rx) < ybox(j+1));
            end
        end     
    end








end