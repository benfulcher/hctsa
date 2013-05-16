function out = TSTL_boxdim(y,htembed,bins)
% Uses TSTOOL code boxdim
% y: column vector of time series data
% htembed: 'how to embed' -- (i) a scalar to embed in that dimension
% with unit delay, (ii) a vector to embed in that [dim,tau], (iii) a cell
% to embed with that method and that choice of time-delay, e.g.,
% {'fnn','mi'} will use false-nearest-neighbours code to determine the
% embedding dimension and first minimum of mutual information to determine
% the time-delay; {'fnn',5} will use a time-delay of 5. Etc.
% e.g., TSTL_boxdim(x,{'fnn','mi'},20)
% Ben Fulcher October 2009

%% Preliminaries
N = length(y); % length of time series
s = signal(y); % convert to signal object for TSTOOL

if nargin<3,
    disp('<<TSTL_boxdim>>: Using default of 100 bins')
    bins=100;
end

%% Embed the time series

% (1) What embedding dimension, m; time-delay, tau
if iscell(htembed)
    % set time-delay first; used later
    if length(htembed)>1
        if ischar(htembed{2}) % use a routine to inform tau
            switch htembed{2}
                case 'mi'
                    tau = CO_fmmi(y);
                case 'ac'
                    tau = CO_fzcac(y);
            end
        else
            tau=htembed{2};
        end
    else
        tau = 1; % default time-delay
        end
    
    if ischar(htembed{1}) % use a routine to inform m
        switch htembed{1}
            case 'fnn'
                % first time proportion of false nearest neighbours falls
                % below 10%:
                m = NL_fnnmar(y,10,2,tau,0.1);
        end
        
    else
        m = htembed{1}; % a fixed scalar is set
    end

else
    m = htembed(1); % a fixed scalar is set
    if length(htembed)>1
        tau = htembed{2};
    else
        tau = 1;
    end
end
m
tau
s_emb = embed(s,m,tau);

%% Run
boxdimo = data(boxdim(s_emb,bins));

%% Give output
% plot(boxdimo);
% keyboard

end