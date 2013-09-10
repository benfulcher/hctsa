% BF_NormalizeMatrix
% 
% Normalizes all columns of an input matrix.
% 
% INPUTS:
% F, the input matrix
% normopt, the normalization method to use (see body of the code for options)
% itrain, learn the normalization parameters just on these indices, then apply
%         it on the full dataset (required for training/testing procedures where
%         the testing data has to remain unseen).
% 
% OUTPUT:
% F, the normalized matrix
% 
% Note that NaNs are ignored -- only real data is used for the normalization
% (assume NaNs are a minority of the data)
% 
% ------------------------------------------------------------------------------
% Copyright (C) 2013,  Ben D. Fulcher <ben.d.fulcher@gmail.com>,
% <http://www.benfulcher.com>
%
% If you use this code for your research, please cite:
% B. D. Fulcher, M. A. Little, N. S. Jones., "Highly comparative time-series
% analysis: the empirical structure of time series and their methods",
% J. Roy. Soc. Interface 10(83) 20130048 (2010). DOI: 10.1098/rsif.2013.0048
%
% This function is free software: you can redistribute it and/or modify it under
% the terms of the GNU General Public License as published by the Free Software
% Foundation, either version 3 of the License, or (at your option) any later
% version.
% 
% This program is distributed in the hope that it will be useful, but WITHOUT
% ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
% FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
% details.
% 
% You should have received a copy of the GNU General Public License along with
% this program.  If not, see <http://www.gnu.org/licenses/>.
% ------------------------------------------------------------------------------
    
function F = BF_NormalizeMatrix(F,normopt,itrain)
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
        error('Invalid normalization method ''%s''', normopt)
end

end