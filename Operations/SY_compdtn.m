function out = SY_compdtn(y,nseg,eachorpar,npoints)
% Compares the distribution of values in nseg consecutive partitions of the signal,
% returning the sum of differences of the kernel smoothed distributions
% eachorpar is either 'each': compares each subdistribution to each other
% subdistribtuion; or 'par': compares each subdistribtuion to the parent
% (full signal) distribution
% Ben Fulcher August 2009

% plot outputs
doplot = 0;

% Check inputs:
if nargin < 2 || isempty(nseg) % number of segments
    nseg = 5;
end
if nargin < 3 || isempty(eachorpar)
    eachorpar = 'par'; % compare each subsection to full (parent) distribution
end
if nargin < 4 || isempty(npoints)
    % number of points to compute the distribution across
    npoints = 200; % 200 by default
end

N = length(y); % number of samples in the time series
lseg = floor(N/nseg);
dns = zeros(npoints,nseg);
r = linspace(min(y),max(y),npoints); % make range of ksdensity uniform across all subsegments

% Compute the kernel-smoothed distribution in all nseg segments of the time series
for i = 1:nseg
    dns(:,i) = ksdensity(y((i-1)*lseg+1:i*lseg),r,'function','pdf');
end

if doplot
    figure('color','w')
    plot(dns,'k')
end

% Compare the local distributions
switch eachorpar
    case 'par'
        % Compares each subdistribtuion to the parent (full signal) distribution
        pardn = ksdensity(y,r,'function','pdf');
        divs = zeros(nseg,1);
        for i = 1:nseg
            divs(i) = sum(abs(dns(:,i)-pardn')); % each is just divergence to parent
        end
        if doplot
            hold on; plot(pardn,'r','LineWidth',2); hold off
        end
        
        % return same statistics as for the 'each' case
        out.meandiv = mean(divs);
        out.mediandiv = median(divs);
        out.mindiv = min(divs);
        out.maxdiv = max(divs);
        out.stddiv = std(divs);
        
    case 'each'
        % Compares each subdistribtuion to the parent (full signal) distribution
        if nseg == 2 % output is just an integer: only two distributions to compare
            out = sum(abs(dns(:,1)-dns(:,2)));
            return
        end
        
        % nseg > 2: need to compare a number of different distributions against each other
        diffmat = NaN * ones(nseg); % store pairwise differences
                                    % start as NaN to easily get upper triangle later
        for i = 1:nseg
            for j = 1:nseg
                if j > i
                    diffmat(i,j) = sum(abs(dns(:,i)-dns(:,j))); % store sum of absolute differences
                end
            end
        end
        
        divs = diffmat(~isnan(diffmat)); % (the upper triangle of diffmat)
                                         % set of divergences in all pairs of segments of the time series
        % divs = diffmat(diffmat > 0); % a set of non-zero divergences in all pairs of segments of the time series
        % if isempty(divs);
        %     fprintf(1,'That''s strange -- no changes in distribution??! This must be a really strange time series.\n');
        %     out = NaN; return
        % end
        
        % Return basic statistics on differences in distributions in different segments of the time series
        out.meandiv = mean(divs);
        out.mediandiv = median(divs);
        out.mindiv = min(divs);
        out.maxdiv = max(divs);
        out.stddiv = std(divs);
        
    otherwise
        error('Unknown method ''%s'', should be ''each'' or ''par''',eachorpar);
end

end