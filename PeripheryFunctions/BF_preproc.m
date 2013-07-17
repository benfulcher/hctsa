function [yp, best] = BF_preproc(y,choosebest,order,beatthis,dospectral)
% Returns a bunch of time series in the structure yp that have been
% preprocessed in a number of different ways.
% If second argument is specified, will choose amongst the preprocessings
% for the 'best' one according to the given criterion.
% Based on (really improvement/development of) MF_preproc.
% Ben Fulcher 25/2/2010

% I think a good way of 'normalizing' over autocorrelation should also be
% incorporated, but it's not obvious how, other than some downsampling...

% NOTE: yp is NOT z-scored -- needs to be z-scored after, if necessary.

%% Inputs
if nargin < 2
    choosebest = ''; % just return all the time series in structure yp.
end
if nargin < 3 || isempty(order);
   order = 2; % extra parameter for some choosebest settings.
end
if nargin < 4 || isempty(beatthis)
    beatthis = 0; % it has to beat doing nothing by this percentage. i.e., here it 
                % just has to beat doing nothing. Increasing will increase the
                % probablility of doing nothing.
end
if nargin < 5 || isempty(dospectral)
    dospectral = 1; % I did this because it's often worse to do a spectral
                    % method when another would do better. i.e., the remove
                    % around a peak can just overpower the structure in the
                    % time series, when indeed all is necessary is a linear
                    % detrending. Eventually we should have a complexity
                    % penalty.
end

%% Preliminaries
yp.nothing = y; % this *has* to be the first element of yp

%% 1) Differencing
yp.d1 = diff(y,1);
yp.d2 = diff(y,2);
yp.d3 = diff(y,3);

%% 2) Remove from power spectrum
% note that edge effects are dealt with simply by eliminating some small
% fraction from the start and end of the time series
if dospectral
    yp.lf_02 = SUB_remps(y,0.2,'lf');
    yp.lf_02_d1 = diff(yp.lf_02);
    yp.peaks_08 = SUB_remps(y,0.8,'biggest');
    yp.peaks_08_d1 = diff(yp.peaks_08);
end

%% 3) Remove piece-wise polynomials

yp.p1_5 = SUB_rempt(y,1,5);
yp.p1_10 = SUB_rempt(y,1,10);
yp.p1_20 = SUB_rempt(y,1,20);
yp.p1_40 = SUB_rempt(y,1,40);

yp.p2_5 = SUB_rempt(y,2,5);
yp.p2_10 = SUB_rempt(y,2,10);
yp.p2_20 = SUB_rempt(y,2,20);
yp.p2_40 = SUB_rempt(y,2,40);

%% Wavelet decomposition
% Maybe later - don't want to worry about toolboxes so much.

%% Rank map onto Gaussian distribution
x = sort(randn(N,1),'ascend'); % create the Normal distribution using random numbers
% could use mapping through linspace, but then a choice of the most
% extreme... I think it's statistically better to do this...? It makes it a
% stochastic algorithm, though...
[~, ix] = sort(y,'ascend');
XX = zeros(N,1);
XX(ix) = x; % the new mapped time series has entries in same rank ordering as y
            % but according to the new distribution defined by x.
yp.rmgd = XX;

%% Positive-only transformations
% log, log returns, sqrt, box-cox
if all(y > 0)
    yp.log = log(y); % log
    yp.logret = diff(log(y)); % log returns
    % Box-Cox
    yp.boxcox02 = SUB_boxcox(y,0.2);
    yp.boxcox05 = SUB_boxcox(y,0.5);
    yp.boxcox2 = SUB_boxcox(y,2);
    yp.boxcox20 = SUB_boxcox(y,20);
end
if all(y >= 0)
    yp.sqrt = sqrt(y);
end


%% -------------------------------------------------------------- %%
% Choose 'best' preprocessing according to some specified criterion, if
% 'choosebest' is specified
%  --------------------------------------------------------------

% Don't do any more (this saves us being indented in an if loop
%                       throughout the rest)
if isempty(choosebest)
    return
end

% We are now going to test based on the specified criterion the best
% preprocessing

% Preliminaries
fields = fieldnames(yp);
nfields = length(fields);

% Ideally would check in turn a bunch of things, like a strong linear
% trend, would check for and remove if necessary. I strong quadratic trend
% could also be checked (& higher polynomials in a hierarchy).
% Then could check seasonality, ...
% Anyway, this is something for now.
% In fact, can think of some regularization on 'complexity' akin to model
% selection -- more complex preprocessings may be penalized in favour of
% simpler ones (i.e., a linear detrend favoured over a quadratic detrend).
% Similarly for increasing levels of differencing. It's a compromise
% to 'lose' just enough of the original (trivial) structure and keep the
% 'interesting' structure...
% Could maybe also incorporate many measures that may be important for this
% 'prewhitening' procedure, rather than just specific choices implemented
% below. ++BF

switch choosebest
    case 'ar' % picks the *worst* fit to an AR(p) model
%         order = 2; % should be set from inputs/defaults
        rmserrs = zeros(nfields,1);
        
        for i = 1:nfields; % each preprocessing performed
            data = [];
            eval(sprintf('data = yp.%s;',fields{i}));
            data = BF_zscore(data);
            
            % (i) fit the model
            m = ar(data,order);
            
            % (ii) in-sample prediction error
            e = pe(m,data);
            rmserrs(i) = sqrt(mean(e.^2));
%             statstore.mabserr(i) = mean(abs(e));
%             statstore.ac1(i) = CO_autocorr(e,1);
        end

        % Now choose the one with the *most* error (!) -- i.e., the AR
        % model finds it hardest to fit
        if any(rmserrs > rmserrs(1)*(1+beatthis));
            best = fields{rmserrs == max(rmserrs)};
        else
            best = 'nothing';
        end

    case 'ac' % picks the *lowest* tau-step autocorrelated result
        acs = zeros(nfields,1);
        for i = 1:nfields
            data = [];
            eval(['data = yp.' fields{i} ';']);
            data = BF_zscore(data); % unnecessary for AC
            acs(i) = CO_autocorr(data,order);
        end
        
        % Now choose time series with lowest ac1
        if any(acs < acs(1)*(1-beatthis))
            best = fields{acs == min(acs)};
        else
            best = 'nothing';
        end
    otherwise
        error('Unknown method ''%s''',choosebest);
end

    function yboxcox = SUB_boxcox(x,lambda)
        yboxcox = (x.^lambda-1)/lambda;
    end

    function ydt =  SUB_remps(y,n,method)
        % Removes the first n (proportion) of power spectrum
        % Based on my deseasonalize1.m code
        

        %% Take the Fourier Transform

        Ny = length(y); % number of samples in y
%         t = linspace(0,1,Ny); % time vector
        NFFT = 2^nextpow2(Ny); % next power of 2
        Fy = fft(y,NFFT); % fast fourier transform of y
        Fy1 = Fy(1:NFFT/2+1);
%         f = 1/2*linspace(0,1,NFFT/2+1); % frequency vector

        %% Remove this range
        % set it to (mean of the rest) across this range
        switch method
            case 'lf'
                cullr = 1:floor(length(Fy1)*n);
            case 'biggest'
                cullr = find(abs(Fy1) > quantile(abs(Fy1),n));
            otherwise
                error('Unknown method ''%s''', method);
        end
            
        meanrest = mean(abs(Fy1(setxor(1:end,cullr))));
%         meanrest = 0;
        FyF = Fy;
        FyF(cullr) = meanrest;
        FyF(end-cullr+2) = meanrest;

        
        % PLOT
%         plot(abs(Fy)),hold on; plot(abs(FyF),'--r'); hold off
%         input('Here''s the filtered one...')
%         plot(abs(FyF),'k');
%         input('Again on it''s own...')

            
        %% Inverse Fourier Transform
        ydt = ifft(FyF,NFFT);
        ydt = BF_zscore(ydt(1:Ny)); % crop to desired length

        % CRUDE REMOVAL OF EDGE EFFECTS
        lose = min(5,floor(0.01*length(ydt))); % don't want to lose more than 1% of time series
        ydt = ydt(1+lose:end-lose);
        
        % PLOT
%         plot(zscore(ydt),'b'); hold on; plot(y,'r'); hold off;
%         input(['Mean difference is ' num2str(mean(y-ydt))])
    
    end

    function ydt = SUB_rempt(y,order,nbits)
        N = length(y);
        ydt = zeros(N,1);
        bits = round(linspace(0,N,nbits+1));
        for k = 1:nbits
            r = bits(k)+1 : bits(k+1); % range defined by adjacent 'bits'
            x = (1:length(r))'; % faux x-range
            ybit = y(r); % y-range
            p = polyfit(x,ybit,order);
            ydt(r) = ybit-polyval(p,x);
        end
        ydt = BF_zscore(ydt);
%         plot(y,'b'); hold on; plot(ydt,'r');
%         input('here we are')
    end



end