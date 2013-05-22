function out = EN_bennoisetit(y,testmeth,noisetype)
% Inspired by the noise titration paper: "Titration of chaos with added noise", Poon and Barahona, PNAS USA, Vol. 98, 7107 (2001)
% This is not that algorithm, I am just going to implement a somewhat similar idea in my own way
% (c.f., ST_benincdec)
% Ben Fulcher, September 2009

N = length(y);

% The same random noise signal will contaminate each column with increasing
% amplitude
switch noisetype
    case 'n'
        % Normally-distributed white noise
        n = randn(N,1);
    case 'u'
        % uniformly-distributed white noise
        error('Uniformly-distributed white noise not yet implemented')
    case 'corr'
        % Gaussian white noise with same linear correlation of data
        % like an MA model output
        [a e] = arcov(y,1);
        n = randn(N+1,1);
        n = a(1)*n(2:end) + a(2)*n(1:end-1);
end
noisecoeffs = (0:0.1:4);
N = n*noisecoeffs;

% Duplicate the time series along the columns in the matrix
Y = y*ones(1,length(noisecoeffs));

Y = Y + N; % signal plus noise
Y = zscore(Y);

% Now apply a test to each column, and look at variation
tests = zeros(length(noisecoeffs),1);
switch testmeth
    case 'BinEn_diff'
        % Discretize each column to 0,1 based on increase/decrease
        Y=((sign(diff(Y)))+1)/2;
        N=size(Y,1); % length of signal - 1 (difference operation)
        for i=1:length(tests)
            pup=length(find(Y(:,i)==1))/N;
            pdown=1-pup;
            p=[pup pdown];
            tests(i)=-sum(p.*log(p));
        end
    case 'BinEn_mean'
        % Discretize to 1 if above mean, zero otherwise
        Y = (sign(Y)+1)/2;
        N = size(Y,1); % length of signal - 1 (difference operation)
        for i = 1:length(tests)
            pup = length(find(Y(:,i)==1))/N;
            pdown = 1-pup;
            p = [pup, pdown];
            tests(i) = -sum(p.*log(p));
        end
    case 'BinEn_iqr'
        iqr = quantile(Y,[0.25 0.75]);
        for i = 1:size(Y,2);
            iniqr = (Y(:,i)>iqr(1,i) & Y(:,i)<=iqr(2,i));
            Y(:,i) = 0;
            Y(iniqr,i) = 1;
        end
    case 'permen2'
        for i = 1:length(tests)
            tests(i) = EN_permen(Y(:,i),2);
            % for some reason this doesn't work well...
        end
    case 'permen3'
        for i = 1:length(tests)
            tests(i) = EN_permen(Y(:,i),3);
        end
    case 'SampEn2_02'
        for i = 1:length(tests)
            tests(i) = EN_sampenc(Y(:,i),2,0.2);
        end
end


% plot(tests);



end