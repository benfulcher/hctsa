function out = TSTL_cao(y,maxdim,tau,NNR,Nref,justanum)
% Uses TSTOOL code cao
% INPUTS
% y: column vector of time series data
% maxdim: maximum embedding dimension to consider
% tau: time delay (can also be 'mi' or 'ac')
% NNR: number of nearest neighbours to use
% Nref: number of reference points (can also be a fraction; of data length)
% justanum [opt]: just return a number, the embedding dimension, based on
%                 the specified criterion
% OUTPUTS:
% a large structure with various statistics on the curves
% Ben Fulcher October 2009

%% Preliminaries
N = length(y); % length of time series
s = signal(y); % convert to signal object for TSTOOL

% (1) Maximum dimension, maxdim
if nargin < 2 || isempty(maxdim)
    maxdim = 10; % default maxdim is 10
end

% (2) Time delay, tau
if nargin < 3 || isempty(tau)
    tau = 'ac'; % choose from first zero crossing of ACF
end
if ischar(tau) % determine by some other method
    switch tau
        case 'mi'
            tau = CO_firstmin(y,'mi');
        case 'ac'
            tau = CO_fzcac(y);
        otherwise
            error('Unknown time-delay method ''%s''',tau);
    end
end

% (3) Number of nearest neighbours, NNR
if nargin < 4 || isempty(NNR)
    NNR = 3; % default to three nearest neighbours
end

% (4) Number of reference points, Nref
if nargin < 5 || isempty(Nref)
    Nref = -1; % default: use all points
end
if Nref < 1 && Nref > 0 % specify a fraction of data size
    Nref = round(N*Nref);
end

if nargin < 6
    justanum = [];
end

%% Run
try [caoo1 caoo2] = cao(s,maxdim,tau,NNR,Nref);
catch err
    % time series is too short for these embedding parameters; set all outputs to NaNs...
	if strcmp(err.message,'time series to short for chosen embedding parameters')
        if ~isempty(justanum) % just return a number for embedding dimension
            fprintf(1,'TSTL_cao: Time series too short for chosen embedding parameters. Using m = 10.\n');
            out = 10;
        else
            fprintf(1,'TSTL_cao: Time series too short for chosen embedding parameters. Returning NaNs.\n');
            out = NaN;
        end
        return
    elseif strcmp(err.identifier,'MATLAB:unassignedOutputs')
        error('TSTOOL''s ''cao'' function returned no output.')
    else
        error('Something weird happend with TSTOOL''s cao routine');
    end
end
caoo1 = data(caoo1);
caoo2 = data(caoo2);

if ~isempty(justanum) % JUST RETURN A NUMBER FOR EMBEDDING DIMENSION
    if ~iscell(justanum), justanum = {justanum}; end
    switch justanum{1}
        case 'thresh'
            % when caoo1 passes above a threshold
            % second element of array is the threshold; default is 10
            if length(justanum)<2
               th = 0.8;
            else
                th = justanum{2};
            end
            out = SUB_first(caoo1,'above',th,maxdim);
        case 'mthresh'
            % when gradient passes below a threshold (levels off)
            % second element of array is the threshold; default is 10
            if length(justanum)<2
                th = 0.1; % default
            else
                th = justanum{2};
            end
            m1 = diff(caoo1); % first differences
            out = SUB_first(m1,'below',th,maxdim-1)+1;
        case 'mmthresh'
            % when change in gradient passes above a threshold (levels off)
            % second element of justanum is the threshold; default is 10
            m1 = diff(caoo1); % first differences
            mm1 = abs(m1(1:end-1))./abs(m1(2:end));
            if length(justanum)<2
                th = 10;
            else
                th = justanum{2};
            end
            out = SUB_first(mm1,'above',th,maxdim-2)+1;
    end
else % RETURN STATISTICS ON CURVES
    
    %% Give output
    % plot(caoo1,'.-b'); hold on; plot(caoo2,'.-k'); hold off
    % (1) the raw values of each vector
    for i = 1:maxdim
        eval(sprintf('out.caoo1_%u = caoo1(%u);',i,i))
        eval(sprintf('out.caoo2_%u = caoo2(%u);',i,i))
    end
    
    % statistics on each vector
    out.min1 = min(caoo1);
    out.max1 = max(caoo1);
    out.min2 = min(caoo2);
    out.max2 = max(caoo2);
    out.median1 = median(caoo1);
    out.median2 = median(caoo2);
    out.std1 = std(caoo1);
    out.std2 = std(caoo2);
    
    %% When passes a threshold
    % if doesn't pass a threshold, assume it to be the next dimension
    % although this may not be the best way of dealing with this issue...
    % could just make NaN as 'unknown', but then would get alot of missing
    % values in the table...
    
    out.fp05_1 = SUB_first(caoo1,'above',0.5,maxdim);
    out.fp05_2 = SUB_first(caoo2,'above',0.5,maxdim);
    out.fp08_1 = SUB_first(caoo1,'above',0.8,maxdim);
    out.fp08_2 = SUB_first(caoo2,'above',0.8,maxdim);
    out.fp09_1 = SUB_first(caoo1,'above',0.9,maxdim);
    out.fp09_2 = SUB_first(caoo2,'above',0.9,maxdim);
    
    %% When gradient passes a threshold -- i.e., levels off
    m1 = diff(caoo1); % first differences
    m2 = diff(caoo2); % first differences (may not be appropriate for this component, but try anyway

    
    if all(isnan(m1))
        out.fm02_1 = NaN;
        out.fm01_1 = NaN;
		out.fm005_1 = NaN;
    else
        % first time drops below 0.1 (not the absolute value)
        % (+1s because of differencing)
        out.fm02_1 = SUB_first(m1,'below',0.2,maxdim-1)+1;
        out.fm01_1 = SUB_first(m1,'below',0.1,maxdim-1)+1;
        out.fm005_1 = SUB_first(m1,'below',0.05,maxdim-1)+1;
    end
    
    if all(isnan(m2))
        out.fm02_2 = NaN;
        out.fm01_2 = NaN;
        out.fm005_2 = NaN;
    else
        % first time drops below 0.1 (not the absolute value)
        % (+1s because of differencing)
        out.fm02_2 = SUB_first(m2,'below',0.2,maxdim-1)+1;
        out.fm01_2 = SUB_first(m2,'below',0.1,maxdim-1)+1;
        out.fm005_2 = SUB_first(m2,'below',0.05,maxdim-1)+1;    
    end



    % When magnitude of gradient decreases by some factor
    % (+1s because of differencing)
    mm2 = abs(m2(1:end-1))./abs(m2(2:end));
    mm1 = abs(m1(1:end-1))./abs(m1(2:end));


	if all(isnan(mm1))
		out.fmm10_1 = NaN;
		out.fmm20_1 = NaN;
		out.fmm40_1 = NaN;
		out.fmmmax_1 = NaN;		
	else
	    out.fmm10_1 = SUB_first(mm1,'above',10,maxdim-2)+1;
	    out.fmm20_1 = SUB_first(mm1,'above',20,maxdim-2)+1;
	    out.fmm40_1 = SUB_first(mm1,'above',40,maxdim-2)+1;
	    % where is maximum (+1s because of differencing
	    out.fmmmax_1 = find(mm1 == max(mm1),1,'first')+1;
	end


	if all(isnan(mm2))
		out.fmm10_2 = NaN;
		out.fmm20_2 = NaN;
		out.fmm40_2 = NaN;
		out.fmmmax_2 = NaN;
	else
		out.fmm10_2 = SUB_first(mm2,'above',10,maxdim-2)+1;
	    out.fmm20_2 = SUB_first(mm2,'above',20,maxdim-2)+1;
	    out.fmm40_2 = SUB_first(mm2,'above',40,maxdim-2)+1;
	    % where is maximum (+1s because of differencing)
	    out.fmmmax_2 = find(mm2 == max(mm2),1,'first')+1;
	end


end

% plot(boxdimo);

    function yep = SUB_first(x,ab,th,maxdim)
        % for input vector x, returns first index that it exceeds (ab = 'above') or
        % goes under ('below') the threshold th
        if strcmp(ab,'above')
            yep = find(x > th,1,'first');
        elseif strcmp(ab,'below')
            yep = find(x < th,1,'first');
        end
        if isempty(yep), yep = maxdim + 1; end 
        
    end

end