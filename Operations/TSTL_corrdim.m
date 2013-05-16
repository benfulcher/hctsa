function out = TSTL_corrdim(y,nbins,embedparams)
% Uses TSTOOL code corrdim
% uses the box counting approach.
% y: column vector of time series data
% nbins: maximum number of partitions per axis
% embedparams [opt]: embedding parameters in 2-entry cell
% Ben Fulcher November 2009

%% Preliminaries
N = length(y); % length of time series

% (1) Maxmum number of partitions per axis, nbins
if nargin<2 || isempty(nbins)
    nbins = 100; % default number of bins per axis is 100
end

% (2) Set embedding parameters to defaults
if nargin<3 || isempty(embedparams)
    embedparams = {'ac','cao'};
else
    if length(embedparams)~=2
        disp('given embedding parameters incorrectly formatted -- need {tau,m}')
    end
end

%% Embed the signal
% convert to embedded signal object for TSTOOL
s = benembed(y,embedparams{1},embedparams{2},1);

if ~strcmp(class(s),'signal') && isnan(s); % embedding failed
	out = NaN; % set all outputs to NaN
    return
end

%% Run

rs = data(corrdim(s,nbins));
% Contains ldr as rows for embedding dimensions 1:m as columns;
% plot(rs);
% keyboard

%% Output Statistics
% Note: these aren't very well motivated.
m = size(rs,2); % number of embedding dimensions
ldr = size(rs,1); % I don't really know what this means; = 17
% keyboard
for i = 2:m
    eval(['out.meand' num2str(i) ' = mean(rs(:,' num2str(i) '));'])
    eval(['out.mediand' num2str(i) ' = median(rs(:,' num2str(i) '));'])
    eval(['out.mind' num2str(i) ' = min(rs(:,' num2str(i) '));'])
end

for i = 2:ldr
    eval(['out.meanr' num2str(i) ' = mean(rs(' num2str(i) ',:));'])
    eval(['out.medianr' num2str(i) ' = median(rs(' num2str(i) ',:));'])
    eval(['out.minr' num2str(i) ' = min(rs(' num2str(i) ',:));'])
    eval(['out.meanchr' num2str(i) ' = mean(diff(rs(' num2str(i) ',:)));'])
end

out.stdmean = std(mean(rs));
out.stdmedian = std(median(rs));

rsstretch = rs(:);
out.medianstretch = median(rsstretch);
out.minstretch = min(rsstretch);
out.iqrstretch = iqr(rsstretch);

	% used to be necessary when master function NaN was not yet implemented...
    % function out = SUB_allNaNs(m,ldr)
    %     
    %     for j=2:m
    %         eval(['out.meand' num2str(j) ' = NaN;'])
    %         eval(['out.mediand' num2str(j) ' = NaN;'])
    %         eval(['out.mind' num2str(j) ' = NaN;'])
    %     end
    %     
    %     for j=2:ldr
    %         eval(['out.meanr' num2str(i) ' = NaN;'])
    %         eval(['out.medianr' num2str(i) ' = NaN;'])
    %         eval(['out.minr' num2str(i) ' = NaN;'])
    %         eval(['out.meanchr' num2str(i) ' = NaN;'])
    %     end
    %     
    %     out.stdmean = NaN;
    %     out.stdmedian = NaN;
    %     out.medianstretch = NaN;
    %     out.minstretch = NaN;
    %     out.iqrstretch = NaN;
    % end


end