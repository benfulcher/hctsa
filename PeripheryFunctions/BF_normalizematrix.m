function F = BF_normalizemat(F,normopt,itrain)
% NaNs are ignored -- only real data is used for the normalization
% (assume NaNs are a minority of the data)
% Ben Fulcher 28/1/2011 -- Added this NaN capability 
% Ben Fulcher 12/9/2011 -- Added itrain input: obtain the trainsformation
% on this subset, apply it to all the data.

if nargin < 2 || isempty(normopt)
    fprintf(1,'We''re normalizing using sigmoid transform by default\n')
    normopt = 'sigmoid';
end

if nargin < 3
    itrain = [];
end

N2 = size(F,2);

if isempty(itrain)
    FT = F; % train the transformation on the full dataset
else
    FT = F(itrain,:); % train the transformation on the specified subset
    if ~strcmp(normopt,'scaledSQzscore')
        error('TRAINING SPECIFIER ONLY WORKS FOR ''scaledSQzscore''...');
    end
end

switch normopt
    case 'maxmin'
        % linear rescaling to the unit interval
        for i = 1:N2 % cycle through the metrics
            rr = ~isnan(F(:,i));
            kk = F(rr,i);
            if (max(kk)==min(kk)) % rescaling will blow up
                F(rr,i) = NaN;
            else
                F(rr,i) = (kk-min(kk))/(max(kk)-min(kk));
            end
        end
    case 'scaledSQzscore'
        % scaled sigmoided quantile zscore
        % problem is that if iqr=0, we're kind of fucked
        for i = 1:N2
            rr = ~isnan(F(:,i));
            rt = ~isnan(FT(:,i));
            FF = FT(rt,i); % good values in the training portion
            if iqr(FF)==0
                F(:,i) = NaN;
            else
                % sigmoid transformation (gets median and iqr only
                % from training data FT):
                F1 = (F(rr,i)-median(FF))/(iqr(FF)/1.35);
                kk = 1./(1+exp(-F1));
                % rescale to unit interval:
                F(rr,i) = (kk-min(kk))/(max(kk)-min(kk));
            end
        end
    case 'scaledsigmoid'
        % a sigmoid transform, then a rescaling to the unit interval
        for i = 1:N2 % cycle through the metrics
            rr = ~isnan(F(:,i));
            FF = F(rr,i);
            kk = 1./(1+exp(-zscore(FF)));
            if (max(kk)==min(kk)) % rescaling will blow up
                F(rr,i) = NaN;
            else
                F(rr,i) = (kk-min(kk))/(max(kk)-min(kk));
            end
        end
    case 'scaledsigmoid5q'
        % first caps at 5th and 95th quantile, then does scaled sigmoid
        for i = 1:N2 % cycle through the metrics
            rr = ~isnan(F(:,i));
            FF = F(rr,i);
            qs = quantile(FF,[0.05,0.95]);
            qr = (FF>=qs(1) & FF<=qs(2)); % quantile range
            % calculate mean and std based on quantile range only
            meanF = mean(FF(qr));
            stdF = std(FF(qr));
            if stdF==0
                F(rr,i) = NaN; % avoid +/- Infs
            else
%                 kk = 1./(1+exp(-zscore(FF)));
                kk = 1./(1+exp(-(FF-meanF)/stdF));
                F(rr,i) = (kk-min(kk))/(max(kk)-min(kk));
            end
        end
    case 'sigmoid'
        for i = 1:N2 % cycle through the metrics
            rr = ~isnan(F(:,i));
            FF = F(rr,i);
%             F(:,i) = 1./(1+exp(-pi*zscore(F(:,i))/sqrt(3)));
            F(rr,i) = 1./(1+exp(-zscore(FF)));
        end
% 	case 'maxmin'
% 		for i = 1:N2 % cycle through the metrics
%             rr = ~isnan(F(:,i));
%             FF = F(rr,i);
%             if range(FF)==0
%                 F(rr,i) = NaN;
%             else
%                 F(rr,i) = (FF-min(FF))/(max(FF)-min(FF));
%             end
%         end
    case 'zscore'
%         F = zscore(F);
        for i = 1:N2
            rr = ~isnan(F(:,i));
            F(rr,i) = zscore(F(rr,i));
        end
    case 'Qzscore'
        % quantile zscore
        % invented by me.
        for i = 1:N2
            rr = ~isnan(F(:,i));
            FF = F(rr,i);
            if iqr(FF)==0 % could get +/- Infs otherwise
                F(rr,i) = NaN;
            else
                F(rr,i) = (FF-median(FF))/(iqr(FF)/1.35);
            end
        end
    case 'SQzscore'
        % sigmoided quantile zscore
        for i = 1:N2
            rr = ~isnan(F(:,i));
            FF = F(rr,i);
            if iqr(FF)==0 % could get +/- Infs otherwise
                F(rr,i) = NaN;
            else
                F1 = (FF-median(FF))/(iqr(FF)/1.35);
                F(rr,i) = 1./(1+exp(-F1));
            end
        end
    case 'scaled2ways'
        for i = 1:N2
            rr = ~isnan(F(:,i));
            FF = F(rr,i);
            if iqr(FF)==0
                % then there's definitely no outlier problem: can safely do a
                % sigmoid
                kk = 1./(1+exp(-zscore(FF)));
                if (max(kk)==min(kk)) % rescaling will blow up
                    F(rr,i) = NaN;
                else
                    F(rr,i) = (kk-min(kk))/(max(kk)-min(kk));
                end
            else
                % wide distribution -- do a transformation that is not so
                % sensitive to outliers
                F1 = (FF-median(FF))/(iqr(FF)/1.35);
                kk = 1./(1+exp(-F1));
                F(rr,i) = (kk-min(kk))/(max(kk)-min(kk));
            end
        end
    case 'LDscaled'
        % (i) maxmin
        % linear rescaling to the unit interval
        for i = 1:N2 % cycle through the metrics
            rr = ~isnan(F(:,i));
            kk = F(rr,i);
            if (max(kk)==min(kk)) % rescaling will blow up
                F(rr,i) = 0; % constant column
            else
                F(rr,i) = (kk-min(kk))/(max(kk)-min(kk));
            end
        end
        
        
        % (ii) outliersigmoid or sigmoid or nothing
        for i = 1:N2
            rr = ~isnan(F(:,i));
            FF = F(rr,i);
            if iqr(FF)==0
                if std(FF)==0
                    F(rr,i) = 0;
                else
                    % then there's definitely no outlier problem: can safely do a
                    % normal sigmoid
                    kk = 1./(1+exp(-zscore(FF)));
                    if (max(kk)==min(kk)) % rescaling will blow up
                        F(rr,i) = NaN;
                    else
                        F(rr,i) = (kk-min(kk))/(max(kk)-min(kk));
                    end
                end
            else
                % wide distribution -- do a transformation that is not so
                % sensitive to outliers
                F1 = (FF-median(FF))/(iqr(FF)/1.35);
                kk = 1./(1+exp(-F1));
                F(rr,i) = (kk-min(kk))/(max(kk)-min(kk));
            end
        end
    otherwise
        error('Invalid normalization method')
end

end


% disp(['Normalization took ' num2str(toc) ' s']);

%% OLD METHOD: choose normalization method based on the statistic: mpostp
% for i=1:nm
%     nbr=find(isfinite(TS_loc(:,i))); % 'not bad' range -- not NaN or +/-Inf
% 	zeyngel=TS_loc(nbr,i);
%     if length(nbr)<nts; % some bad ones (i.e., NaNs or Infs) -- set them to NaNs in TS_loc_N
% 		nnbr=setxor(1:nts,nbr);
% 	    TS_loc_N(nnbr,i)=NaN;
% 	end
%     switch mpostp(i)
%         case 1 % min-max
%             mx=max(zeyngel); mn=min(zeyngel);
%             TS_loc_N(nbr,i)=(zeyngel-mn)/(mx-mn);
%         case 2 % 1.5 std
%             mu=mean(zeyngel); stu=std(zeyngel);
%             if stu==0; % all the same, make them all zeros
%                 disp(['problem with ' num2str(i)]); TS_loc_N(:,i)=0;
%             else
%                 mx=mu+1.5*stu; mn=mu-1.5*stu;
%                 TS_loc_N(nbr(zeyngel>mx),i)=1; TS_loc_N(nbr(zeyngel<mn),i)=0; % saturates
%                 r=nbr(zeyngel>=mn & zeyngel<=mx); % scales
%                 TS_loc_N(r,i)=(TS_loc(r,i)-mn)/(mx-mn);
%             end
%         case 3 % -1 => 0, 0 => 0.5, 1=>1
%             if any(abs(zeyngel))>1;
%                 disp(['error in ' num2str(i)])
%             end
%             TS_loc_N(nbr,i)=0.5*(zeyngel+1);
%         case 4 % 0 maps to 0.5 incorporates full range to either side (i.e. [0,infty) => [0.5,1))
%             TS_loc_N(nbr,i)=0.5*zeyngel/max(abs(zeyngel))+0.5;
%         case 5
%             p=quantile(zeyngel,[0.1 0.9]); % top and bottom 10% will saturate out of the plot
%             TS_loc_N(nbr(zeyngel>p(2)),i)=1; TS_loc_N(nbr(zeyngel<p(1)),i)=0; % saturates
%             r=nbr(zeyngel>=p(1) & zeyngel<=p(2)); % scales
%             TS_loc_N(r,i)=(TS_loc(r,i)-p(1))/(p(2)-p(1));
%         case 6 % for p-values: between 0 and 1
%             if any(abs(zeyngel-0.5))>0.5;
%                 disp(['error in plotting ' num2str(i) ' : range larger than [0,1]'])
%             end
%             TS_loc_N(nbr,i)=zeyngel;
%     end
% end