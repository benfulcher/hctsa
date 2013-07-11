function out = DN_kssimp(x,varargin)
% Outputs a range of simple statistics about the kernel-smoothed
% distribution of values in the time series
% Uses the ksdensity function in Matlab's Statistics Toolbox
% Ben Fulcher, 2009

% case 'numpeaks' % number of peaks
% case 'max'
% case 'entropy'
% case 'asym'
% case 'plsym'


% Define preliminaries
m = mean(x);

% First compute the smoothed empirical distribution of values in the time series
[f, xi] = ksdensity(x);

% 1. Number of peaks
df = diff(f);
ddf = diff(df);
sdsp = ddf(sgnchange(df));
out.npeaks = sum(sdsp < -0.0002); % 'large enough' maxima

% 2. Max
out.max = max(f); % maximum of the distribution

% 3. Entropy
out.entropy = - sum(f(f > 0).*log(f(f > 0))*(xi(2)-xi(1))); % entropy of the distribution

% 4. Assymetry
out1 = sum(f(xi > m).*(xi(2)-xi(1)));
out2 = sum(f(xi < m).*(xi(2)-xi(1)));
out.asym = out1/out2;

% 5. Plsym
out1 = sum(abs(diff(f(xi < m))).*(xi(2)-xi(1)));
out2 = sum(abs(diff(f(xi > m))).*(xi(2)-xi(1)));
out.plsym = out1/out2;

% 6. Numcross
% Specified in input
varhere = strcmp(varargin,'numcross');
if any(varhere)
    thresholds = varargin{find(varhere,1,'first') + 1};
    for i = 1:length(thresholds)
        ncrosses = length(sgnchange(f - thresholds(i)));
        outname = regexprep(sprintf('numcross_%.2f',thresholds(i)),'\.',''); % remove dots from 2-d.pl.
        eval(sprintf('out.%s = ncrosses;',outname));
    end
else
    % don't calculate numcross statistics
end

% 7. Area
% Specified in input
varhere = strcmp(varargin,'area');
if any(varhere)
    thresholds = varargin{find(varhere,1,'first') + 1};
    for i = 1:length(thresholds)
        areahere = sum(f(f < thresholds(i)).*(xi(2)-xi(1))); % integral under this portion
        outname = regexprep(sprintf('area_%.2f',thresholds(i)),'\.',''); % remove dots from 2-d.pl.
        eval(sprintf('out.%s = areahere;',outname));
    end
else
    % don't calculate area statistics
end

% 8. Arclength
% Specified in input
varhere = strcmp(varargin,'arclength');
if any(varhere)
    thresholds = varargin{find(varhere,1,'first') + 1};
    for i = 1:length(thresholds)
        % The integrand in the path length formula:
        fd = abs(diff(f(xi > m - thresholds(i) & xi < m + thresholds(i))));
        arclengthhere = sum(fd.*(xi(2)-xi(1)));
        outname = regexprep(sprintf('arclength_%.2f',thresholds(i)),'\.',''); % remove dots from 2-d.pl.
        eval(sprintf('out.%s = arclengthhere;',outname));
    end
else
    % don't calculate arc length statistics
end


% switch pinkbeanie
%     case 'numcross' % number of times crosses a given constant value
%         places = sgnchange(f-ange);
%         out = length(places);
%     case 'area'
%         if nargin < 3
%             error('DN_kssimp: Specify a number')
%         end
%         out = sum(f(f<ange).*(xi(2)-xi(1))); % integral under this portion
%     case 'arclength'  % path length
%         if nargin < 3
%             error('DN_kssimp: Specify a number')
%         end
%         m = mean(x);
%         fd = abs(diff(f(xi > m - ange & xi < m + ange))); % The integrand in the path length formula
%         out = sum(fd.*(xi(2)-xi(1)));
% end


end