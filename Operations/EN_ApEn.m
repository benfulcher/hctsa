function out = EN_ApEn(u,mnom,rth)

% mnom = 1;
% rth=0.2;
r = rth*std(u); % threshold of similarity
N = length(u);% length of time series
phi = zeros(2,1);% phi(1)=phi_m, phi(2)=phi_{m+1}

% if matlabpool('size') == 0
%     matlabpool open
% else
%     disp(['matlab pool already open. Size: ' num2str(matlabpool('size'))])
% end
% tic
for k = 1:2
    m = mnom+k-1; % pattern length
    C = zeros(N-m+1,1);
    % Define the matrix x, containing subsequences of u
    x = zeros(N-m+1,m);
    
    % Form vector sequences x from the time series u
    for i=1:N-m+1
        x(i,:) = u(i:i+m-1);
    end
    
    ax=ones(N-m+1,m);
    for i = 1:N-m+1
        for j=1:m
            ax(:,j)=x(i,j);
        end
        d = abs(x-ax);
        if m > 1;% Takes maximum distance
            d = max(d,[],2)';
        end
        dr = (d<=r);
        C(i) = sum(dr)/(N-m+1); % Number of x(j) within r of x(i)
    end
    phi(k) = mean(log(C));
end
out = phi(1)-phi(2);
% toc

end