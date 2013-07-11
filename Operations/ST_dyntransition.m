function out=ST_dyntransition(y,ng,tau)
% Calculates the transition probabilities and measures of how they change
% as the number of groups; the lag time changes. Discretization is done by
% quantile separation.
% INPUTS: ng: number of groups in the course-graining (scalar for constant,
%         vector for range to look across)
%         tau: the lag; transition matricies corresponding to this lag. We
%         can either downsample the time series at this lag and then do the
%         discretization as normal, or do the discretization and then just
%         look at this dicrete lag. Here we do the former. (scalar for
%         constant tau, vector for range to vary across)
% OUTPUTS: a structure containing the probabilities, size as appropriate
%           for dparam
% Ben Fulcher August 2009


if strcmp(tau,'ac') % determine tau from first zero of autocorrelation
    tau = CO_fzcac(y);
    if tau>length(y)/50; % for highly-correlated signals
        tau = floor(length(y)/50); 
    end
end

nfeat = 9; % the number of features calculated at each point
if length(ng) == 1 && length(tau)>1 % vary tau
    if ng < 2; return; end % need more than 2 groups
    taur = tau; % the tau range
    store = zeros(length(taur),nfeat);
    
    for i = 1:length(taur)
        tau = taur(i);
        if tau>1; y = resample(y,1,tau); end % resample
        yth = loc_discretize(y,ng); % threshold
        store(i,:) = getmeasures(yth);
    end

    
elseif length(tau) == 1 && length(ng)>1 % vary ng
    if min(ng)<2; return; end % need more than 2 groups, always
    ngr = ng; % the ng range (ng is an input vector)
    store=zeros(length(ngr),nfeat);
    if tau>1; y=resample(y,1,tau); end % resample
    
    for i = 1:length(ngr)
        ng=ngr(i);
        yth=loc_discretize(y,ng); % thresholded data: yth
        store(i,:)=loc_getmeasures(yth,ng);
    end
    
    
    ngr = ngr'; % needs to be a column vector for the fitting routines
    
    % 1) mean of diagonal elements of the transition matrix: shows an exponential
    % decay to zero
    s = fitoptions('Method','NonlinearLeastSquares','StartPoint',[1 -0.2]);
    f = fittype('a*exp(b*x)','options',s);
    [c,gof] = fit(ngr,store(:,1),f);
    out.meandiagfexp_a = c.a;
    out.meandiagfexp_b = c.b;
    out.meandiagfexp_r2 = gof.rsquare;
    out.meandiagfexp_adjr2 = gof.adjrsquare;
    out.meandiagfexp_rmse = gof.rmse;
    
    % 2) maximum of diagonal elements of the transition matrix: shows an exponential
    % decay to zero
    s = fitoptions('Method','NonlinearLeastSquares','StartPoint',[1 -0.2]);
    f = fittype('a*exp(b*x)','options',s);
    [c,gof] = fit(ngr,store(:,2),f);
    out.maxdiagfexp_a = c.a;
    out.maxdiagfexp_b = c.b;
    out.maxdiagfexp_r2 = gof.rsquare;
    out.maxdiagfexp_adjr2 = gof.adjrsquare;
    out.maxdiagfexp_rmse = gof.rmse;
    
    % 3) trace of T
    % fit exponential
    s = fitoptions('Method','NonlinearLeastSquares','StartPoint',[1 -0.2]);
    f = fittype('a*exp(b*x)','options',s);
    [c,gof] = fit(ngr,store(:,3),f);
    out.trfexp_a = c.a;
    out.trfexp_b = c.b;
    out.trfexp_r2 = gof.rsquare;
    out.trfexp_adjr2 = gof.adjrsquare;
    out.trfexp_rmse = gof.rmse;
    
    % Also fit linear from the start to a fifth, a tenth of the starting
    % value
    s = fitoptions('Method','NonlinearLeastSquares','StartPoint',[-0.05 1]);
    f = fittype('a*x+b','options',s);
    
    r1 = find(store(:,3)>store(1,3)/5);
    if length(r1)>2
        [~,gof] = fit(ngr(r1),store(r1,3),f);
        out.trflin5_adjr2 = gof.adjrsquare;
    else
        out.trflin5_adjr2 = NaN;
    end
    
    r2 = find(store(:,3)>store(1,3)/10);
    if length(r2)>2
        [~,gof] = fit(ngr(r2),store(r2,3),f);
        out.trflin10adjr2 = gof.adjrsquare;
    else
        out.trflin10adjr2 = NaN;
    end

    % 4) Symmetry; differences in diagonal elements
    % return the slope
    s = fitoptions('Method','NonlinearLeastSquares','StartPoint',[0.1 0]);
    f = fittype('a*x+b','options',s);
    c = fit(ngr,store(:,4),f);
    out.symd_a = c.a;
    
    % return approximately when starts to rise; where means before and
    % after a moving dividing point are most different
    if all(store(:,4) == store(1,4)); % all the same
        out.symd_risept = NaN;
    else
        mba = zeros(length(ngr),2); % means before and after
        sba = zeros(length(ngr),2); % standard deviation before and after
        for i = 3:length(ngr-2)
            mba(i,1) = mean(store(1:i-1,4));
            sba(i,1) = std(store(1:i-1,4))/sqrt(i-1);
            mba(i,2) = mean(store(i+1:end,4));
            sba(i,2) = std(store(i+1:end,4))/sqrt(length(ngr)-i+1);
        end
        tstats = abs((mba(:,1)-mba(:,2))./sqrt(sba(:,1).^2 + sba(:,2).^2));
        out.symd_risept = find(tstats == max(tstats),1,'first');
    end
    
    % 5) trace of covariance matrix
    % check jump:
    out.trcov_jump = store(2,5)-store(1,5);
    if store(2,5)>store(1,5); r1 = 2:length(ngr); % jump
    else r1 = 1:length(ngr);
    end
    % fit exponential decay to range without possible first jump
    s = fitoptions('Method','NonlinearLeastSquares','StartPoint',[1 -0.5]);
    f = fittype('a*exp(b*x)','options',s);
    [c,gof] = fit(ngr(r1),store(r1,5),f);
    out.trcovfexp_a = c.a;
    out.trcovfexp_b = c.b;
    out.trcovfexp_r2 = gof.rsquare;
    out.trcovfexp_adjr2 = gof.adjrsquare;
    out.trcovfexp_rmse = gof.rmse;
    
    % 6) Standard deviation of eigenvalues of T
    % Fit an exponential decay
    s = fitoptions('Method','NonlinearLeastSquares','StartPoint',[1 -0.2]);
    f = fittype('a*exp(b*x)','options',s);
    [c,gof] = fit(ngr,store(:,6),f);
    out.stdeigfexp_a = c.a;
    out.stdeigfexp_b = c.b;
    out.stdeigfexp_r2 = gof.rsquare;
    out.stdeigfexp_adjr2 = gof.adjrsquare;
    out.stdeigfexp_rmse = gof.rmse;
    
    % 7) maximum (real) eigenvalue of T
    % Fit an exponential decay
    s = fitoptions('Method','NonlinearLeastSquares','StartPoint',[1 -0.2]);
    f = fittype('a*exp(b*x)','options',s);
    [c,gof] = fit(ngr,store(:,7),f);
    out.maxeig_fexpa = c.a;
    out.maxeig_fexpb = c.b;
    out.maxeig_fexpr2 = gof.rsquare;
    out.maxeig_fexpadjr2 = gof.adjrsquare;
    out.maxeig_fexprmse = gof.rmse;
    
    % 8) minimum (real) eigenvalue of T
    % Fit an exponential decay
    s = fitoptions('Method','NonlinearLeastSquares','StartPoint',[1 -0.2]);
    f = fittype('a*exp(b*x)','options',s);
    [c,gof] = fit(ngr,store(:,8),f);
    out.mineigfexp_a = c.a;
    out.mineigfexp_b = c.b;
    out.mineigfexp_r2 = gof.rsquare;
    out.mineigfexp_adjr2 = gof.adjrsquare;
    out.mineigfexp_rmse = gof.rmse;
    
    % 9) mean real eigenvalue of T
    % Fit an exponential decay
    s = fitoptions('Method','NonlinearLeastSquares','StartPoint',[1 -0.2]);
    f = fittype('a*exp(b*x)','options',s);
    [c,gof] = fit(ngr,store(:,9),f);
    out.meaneigfexp_a = c.a;
    out.meaneigfexp_b = c.b;
    out.meaneigfexp_r2 = gof.rsquare;
    out.meaneigfexp_adjr2 = gof.adjrsquare;
    out.meaneigfexp_rmse = gof.rmse;
    
    
end

%% %% %% NOW! Output features %% %% %%
% figure('color','w');

% out=store;


%% Subfunctions




    function yth=loc_discretize(y,ng)
        % 1) discretize the time series into a number of groups np
        th=quantile(y,linspace(0,1,ng+1)); % thresholds for dividing the time series values
        th(1)=th(1)-1; % this ensures the first point is included
        % turn the time series into a set of numbers from 1:ng
        yth=zeros(length(y),1);
        for li=1:ng
            yth(y>th(li) & y<=th(li+1))=li;
        end
        if any(yth == 0), yth=[]; return; end % error -- they should all be assigned to a group
        
    end

    function out=loc_getmeasures(yth,ng)
        % returns a bunch of metrics on the transition matrix
        N=length(yth);
%         ng=max(yth);
        
        % 1) Calculate the one-time transition matrix
        T=zeros(ng);
        for j=1:ng
            ri=find(yth == j);
            if isempty(ri) % yth is never j
                T(j,:)=0;
            else
                if ri(end) == N; ri=ri(1:end-1); end % looking at next element; remove last point
                for k=1:ng
                    T(j,k)=length(find(yth(ri+1) == k)); % the next element is of this class
                end
            end
        end
        T=T/(N-1); % N-1 is appropriate because it's a 1-time transition matrix
        
        % 2) return some quantities on the transition matrix, T
        %   (i) diagonal elements
        out(1)=mean(diag(T)); % mean of diagonal elements
        out(2)=max(diag(T)); % max of diagonal elements
        out(3)=sum(diag(T)); % sum of diagonal elements (trace)
%         out(4)=std(diag(T)); % std of diagonal elements

        %  (ii) measures of symmetry:
        out(4)=sum(sum(abs((T-T')))); % sum of differences of individual elements
%         out(5)=sum(sum(tril(T,-1)))-sum(sum(triu(T,+1))); % difference in sums of upper and lower 
                                                          % triangular parts of T

        % (iii) measures from covariance matrix:
        out(5)=sum(diag(cov(T))); % trace of covariance matrix
%         out(6)=max(max(cov(T))); % maximum value in the covariance matrix
%         out(7)=min(min(cov(T))); % minimum value in the covariance matrix

        % (iv) measures from eigenvalues of T
        eigT=eig(T);
        out(6)=std(eigT); % std of eigenvalues
        out(7)=max(real(eigT)); % maximum eigenvalue
        out(8)=min(real(eigT)); % minimum eigenvalue
        out(9)=mean(real(eigT)); % mean eigenvalue
%         out(14)=max(imag(eigT)); % maximum imaginary part of eigenavlues

        % % (vi) measures from eigenvalues of covariance matrix:
%         eigcovT=eig(cov(T));
%         out(10)=max(eigcovT); % max eigenvalue of covariance matrix
% %         out(16)=min(eigcovT); % min eigenvalue of covariance matrix
%         out(11)=std(eigcovT); % std of eigenvalues of covariance matrix
%         out(12)=mean(eigcovT); % mean eigenvalue of covariance matrix
        
    end


end