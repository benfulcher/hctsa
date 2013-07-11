function out = ST_momcorr(x,wl,olap,mom1,mom2,transf)
% Finds correlations between moments of the signal (x) [transformed according
% to transf], in windows (of length wl). x should be z-scored.
% Idea of Nick Jones.
% Ben Fulcher 5 July 2010

N = length(x); % number of samples in the input signal


% sliding window length (samples)
if nargin < 2
    wl = 0.02;
end
if wl < 1
    wl = ceil(N*wl);
end

% sliding window overlap length
if nargin < 3 || isempty(olap)
    olap = 1/5;
end
if olap < 1 % specify a fraction OF THE WINDOW LENGTH
    olap = floor(wl*olap);
end

if nargin < 4
    mom1 = 'mean';
end
if nargin < 5
    mom2 = 'std';
end

if nargin < 6
    transf = 'none';
end

% transform
switch transf
    case 'abs'
        x = abs(x);
    case 'sq'
        x = x.^2;
    case 'sqrt'
        x = sqrt(abs(x));
end

% ok, quit stuffing around

% create the sliding windows
x_buff = buffer(x,wl,olap);
Nw = (N/(wl-olap)); % number of windows

if size(x_buff,2)>Nw
    % fprintf(1,'Should have %u columns but we have %u: removing last one',Nw,size(x_buff,2))
    x_buff = x_buff(:,1:end-1);
end % lose last point


% ok, now we have the sliding window ('buffered') signal, x_buff
% first calculate the first moment in all the windows (each column is a
% 'window' of the signal

M1 = calcmemoments(x_buff,mom1);
M2 = calcmemoments(x_buff,mom2);

R = corrcoef(M1,M2);
out.R = R(2,1); % correlation coefficient
out.absR = abs(R(2,1)); % absolute value of correlation coefficient
out.density = range(M1)*range(M2)/N; % density of points in M1--M2 space
out.mi = benmi(M1,M2,[0,1],[0,1],floor(sqrt(N)));
% this is a poor choice of bin number -- M1 and M2 are not length N

% figure('color','w');
% plot(M1,M2,'.k');


function moms = calcmemoments(x_buff,momtype)
    switch momtype
        case 'mean'
            moms = mean(x_buff);
        case 'std'
            moms = std(x_buff);
        case 'median'
            moms = median(x_buff);
        case 'iqr'
            moms = iqr(x_buff);
    end        
end


end