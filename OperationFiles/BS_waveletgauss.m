% function out = BS_waveletgauss(y)

% % Load a fractal signal. 
% % load vonkoch 
% % vonkoch=vonkoch(1:510); 
% % lv = length(vonkoch);
% 
% % Load a fractal signal. 
% load vonkoch 
% % vonkoch=vonkoch(1:510); 
% lv = length(vonkoch);
% 
% subplot(311), plot(vonkoch);title('Analyzed signal.'); 
% % set(gca,'Xlim',[0 510])
% % Perform discrete wavelet transform at level 5 by sym2. 
% % Levels 1 to 5 correspond to scales 2, 4, 8, 16 and 32. 
% [c,l] = wavedec(vonkoch,5,'sym2');
% 
% % Expand discrete wavelet coefficients for plot. 
% % Levels 1 to 5 correspond to scales 2, 4, 8, 16 and 32. 
% cfd = zeros(5,lv); 
% for k = 1:5 
%     d = detcoef(c,l,k); 
%     d = d(ones(1,2^k),:); 
%     cfd(k,:) = wkeep(d(:)',lv); 
% end 
% 
% cfd = cfd(:); 
% I = find(abs(cfd)<sqrt(eps)); 
% cfd(I)=zeros(size(I)); 
% cfd = reshape(cfd,5,lv);
% 
% % Plot discrete coefficients. 
% subplot(312), colormap(pink(64)); 
% img = image(flipud(wcodemat(cfd,64,'row'))); 
% set(get(img,'parent'),'YtickLabel',[]); 
% title('Discrete Transform, absolute coefficients.') 
% ylabel('level')
% 
% % Perform continuous wavelet transform by sym2 at all integer 
% % scales from 1 to 32. 
% subplot(313)
% ccfs = cwt(vonkoch,1:32,'sym2','plot'); 
% title('Continuous Transform, absolute coefficients.') 
% colormap(pink(64)); 
% ylabel('Scale')




% Load signal. 
load noisdopp; x = noisdopp;

figure(1); subplot(211); 
plot(x); title('Original signal');

% Decompose x at depth 3 with db1 wavelet packets 
% using Shannon entropy. 
wpt = wpdec(x,3,'db1');

% Plot wavelet packet tree wpt. 
plot(wpt)
% Read packet (2,1) coefficients. 
cfs = wpcoef(wpt,[2 1]);

figure(1); subplot(212); 
plot(cfs); title('Packet (2,1) coefficients');

% % load vonkoch 
% % y=vonkoch(1:510);
% y=dlmread('AS_s5_f2_b8_l9910_45476.dat');
% y=y';
% N=length(y);
% 
% subplot(311), plot(y); title('Analyzed signal.'); 
% 
% % Perform discrete wavelet transform at level 5 by sym2. 
% % Levels 1 to 5 correspond to scales 2, 4, 8, 16 and 32. 
% [c,l] = wavedec(y,5,'sym2');
% 
% % Expand discrete wavelet coefficients for plot. 
% % Levels 1 to 5 correspond to scales 2, 4, 8, 16 and 32. 
% cfd = zeros(5,N); 
% for k = 1:5 
%     d = detcoef(c,l,k);
%     d = d(ones(1,2^k),:);
%     cfd(k,:) = wkeep(d(:)',N);
% end
% 
% cfd = cfd(:);
% I = find(abs(cfd)<sqrt(eps));
% cfd(I)=zeros(size(I));
% cfd = reshape(cfd,5,N);
% 
% % Plot discrete coefficients. 
% subplot(312), colormap(pink(64)); 
% img = image(flipud(wcodemat(cfd,64,'row')));
% set(get(img,'parent'),'YtickLabel',[]);
% title('Discrete Transform, absolute coefficients.')
% ylabel('level')
% 
% % Perform continuous wavelet transform by sym2 at all integer 
% % scales from 1 to 16.
% subplot(313)
% ccfs = cwt(y,1:16,'sym2','plot');
% title('Continuous Transform, absolute coefficients.')
% colormap(pink(64));
% ylabel('Scale')

% end