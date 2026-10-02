function opIDs = TS_GiveMeFeatureSet(whatFeatureSet,Operations)
% TS_GiveMeFeatureSet Outputs a set of Operation IDs corresponding to a given set
%
% INPUTS:
% ---whatFeatureSet, the type of feature set to retrieve/filter.
% ---Operations, the Operations table to match to.

% ------------------------------------------------------------------------------
% Copyright (C) 2013-2026, Ben D. Fulcher <ben.d.fulcher@gmail.com>,
% <http://www.benfulcher.com>
%
% If you use this code for your research, please cite the following two papers:
%
% (1) B.D. Fulcher and N.S. Jones, "hctsa: A Computational Framework for Automated
% Time-Series Phenotyping Using Massive Feature Extraction", Cell Systems 5: 527 (2017).
% DOI: 10.1016/j.cels.2017.10.001
%
% (2) B.D. Fulcher, M.A. Little, N.S. Jones, "Highly comparative time-series
% analysis: the empirical structure of time series and their methods",
% J. Roy. Soc. Interface 10(83) 20130048 (2013).
% DOI: 10.1098/rsif.2013.0048
%
% This work is licensed under the Creative Commons
% Attribution-NonCommercial-ShareAlike 4.0 International License. To view a copy of
% this license, visit http://creativecommons.org/licenses/by-nc-sa/4.0/ or send
% a letter to Creative Commons, 444 Castro Street, Suite 900, Mountain View,
% California, 94041, USA.
% ------------------------------------------------------------------------------

%-------------------------------------------------------------------------------
if nargin < 1 || isempty(whatFeatureSet)
    whatFeatureSet = 'catch22';
end
if nargin < 2
    error('You must provide an Operations table to match features to by name');
end
%-------------------------------------------------------------------------------

switch whatFeatureSet
case 'noLengthLocationSpread'
    matchByName = false;
    % Remove length, location, spread-dependent features
    % Old datasets may still use the pre-rename keywords ('lengthdep', etc.);
    % detect which convention is present rather than assuming:
    doOld = isempty(TS_GetIDs('lengthDependent',Operations,'ops','Keywords'));
    if doOld
        lengthIDs = TS_GetIDs('lengthdep',Operations,'ops','Keywords');
        locIDs = TS_GetIDs('locdep',Operations,'ops','Keywords');
        spreadIDs = TS_GetIDs('spreaddep',Operations,'ops','Keywords');
    else
        lengthIDs = TS_GetIDs('lengthDependent',Operations,'ops','Keywords');
        locIDs = TS_GetIDs('locationDependent',Operations,'ops','Keywords');
        spreadIDs = TS_GetIDs('spreadDependent',Operations,'ops','Keywords');
    end
    depIDs = unique([lengthIDs; locIDs; spreadIDs]);
    % Exclude:
    opIDs = setxor(Operations.ID,depIDs);
case 'catch22'
    matchByName = false;
    % The catch22 feature set (EXCLUDES MEAN/SPREAD-DEPENDENT FEATURES),
    % cf. https://github.com/DynamicsAndNeuralSystems/catch22
    %
    % These are the features computed by catch22's own (independent, C)
    % implementation: the operations whose code strings start with 'catch22_',
    % as defined in FeatureSets/INP_ops_catch22.txt (computed by
    % TS_Init(...,'catch22') or 'catch24'). hctsa's own operations that catch22
    % was originally drawn from are NOT used: several have since been fixed or
    % redefined in hctsa (e.g., MD_hrv_classic, CO_Embed2_Dist, SC_FluctAnal,
    % FC_LocalSimple), so their values and names no longer match catch22.
    catch22Codes = {'catch22_DN_HistogramMode_5', 'catch22_DN_HistogramMode_10', ...
        'catch22_DN_OutlierInclude_p_001_mdrmd', 'catch22_DN_OutlierInclude_n_001_mdrmd', ...
        'catch22_CO_f1ecac', 'catch22_CO_FirstMin_ac', ...
        'catch22_SP_Summaries_welch_rect_area_5_1', 'catch22_SP_Summaries_welch_rect_centroid', ...
        'catch22_FC_LocalSimple_mean3_stderr', 'catch22_FC_LocalSimple_mean1_tauresrat', ...
        'catch22_CO_HistogramAMI_even_2_5', 'catch22_CO_trev_1_num', ...
        'catch22_MD_hrv_classic_pnn40', 'catch22_SB_BinaryStats_mean_longstretch1', ...
        'catch22_SB_BinaryStats_diff_longstretch0', 'catch22_SB_MotifThree_quantile_hh', ...
        'catch22_SB_TransitionMatrix_3ac_sumdiagcov', 'catch22_PD_PeriodicityWang_th0_01', ...
        'catch22_CO_Embed2_Dist_tau_d_expfit_meandiff', 'catch22_IN_AutoMutualInfoStats_40_gaussian_fmmi', ...
        'catch22_SC_FluctAnal_2_rsrangefit_50_1_logi_prop_r1', 'catch22_SC_FluctAnal_2_dfa_50_1_2_logi_prop_r1'};
    isCatch22 = ismember(Operations.CodeString,catch22Codes);
    if ~any(isCatch22)
        error(['This dataset contains no catch22 operations. The ''catch22'' feature set uses ', ...
                'catch22''s own implementation (operations ''catch22_*''); compute it with ', ...
                'TS_Init(INP_ts,''catch22'') (or ''catch24'') rather than selecting from the hctsa feature set.']);
    end
    opIDs = Operations.ID(isCatch22);
    fprintf(1,'Matched %u/%u catch22 features (catch22 implementation)\n',length(opIDs),length(catch22Codes));
    if length(opIDs) < length(catch22Codes)
        warning('%u catch22 feature(s) missing from this dataset: %s',length(catch22Codes)-length(opIDs), ...
                    strjoin(setdiff(catch22Codes,Operations.CodeString(isCatch22)),', '));
    end
case 'catchaMouse16'
    matchByName = true;
    % NOTE -- 'SC_FluctAnal_2_dfa_50_2_logi_r2_se2' keeps its upstream name here,
    % but its value diverges from upstream catchaMouse16 for the same reason
    % documented under the catch22 case above (shared split-point-search bug,
    % fixed in hctsa's SC_FluctAnal.m, left as-is in catch22's vendored C code).
    featureNames = {'SY_DriftingMean50_min',...
                    'MF_CompareAR_1_10_05_stddiff',...
                    'SC_FluctAnal_2_dfa_50_2_logi_r2_se2',...
                    'IN_AutoMutualInfoStats_diff_20_gaussian_ami8',...
                    'PH_Walker_momentum_5_w_propzcross',... % (was 'w_momentumzcross', a typo -- no such field ever existed)
                    {'MF_steps_ahead_arma_3_1_6_ac1_h6','MF_steps_ahead_arma_3_1_6_ac1_6'},...
                    {'CO_RemovePoints_absclose_05_remove_ac2rat','DN_RemovePoints_absclose_05_remove_ac2rat','DN_RemovePoints_absclose_05_ac2rat'},... % [new name first]
                    {'MF_steps_ahead_ar_2_6_stde_maxdiff','MF_steps_ahead_ar_2_6_maxdiffrms'},...
                    'SP_Summaries_fft_fpolysat_rmse',...
                    {'CO_HistogramAMI_even_2_3','CO_HistogramAMI_even_2bin_ami3'},...
                    'AC_nl_036',...
                    'AC_nl_112',...
                    'MF_StateSpace_n4sid_1_05_1_ac2',...
                    'ST_LocalExtrema_n100_diffmaxabsmin',...
                    'CO_TranslateShape_circle_35_pts_statav4_m',...
                    'CO_AddNoise_1_even_10_ami_at_10'};
case 'linearAutoCorrs'
    % Linear autocorrelations lags 1 through 20
    matchByName = true;
    featureNames = {'AC_1','AC_2','AC_3','AC_4','AC_5',...
                    'AC_6','AC_7','AC_8','AC_9','AC_10',...
                    'AC_11','AC_12','AC_13','AC_14','AC_15',...
                    'AC_16','AC_17','AC_18','AC_19','AC_20'};
case 'quantiles'
    % Quantiles of the distribution
    matchByName = true;
    featureNames = {'quantile_10',...
                    'quantile_20',...
                    'quantile_30',...
                    'quantile_40',...
                    'quantile_50',...
                    'quantile_60',...
                    'quantile_70',...
                    'quantile_80',...
                    'quantile_90'};
case 'linearAutoCorrsAndQuantiles'
    % Linear autocorrelations lags 1 through 20 plus quantiles
    matchByName = true;
    featureNames = {'AC_1','AC_2','AC_3','AC_4','AC_5',...
                    'AC_6','AC_7','AC_8','AC_9','AC_10',...
                    'AC_11','AC_12','AC_13','AC_14','AC_15',...
                    'AC_16','AC_17','AC_18','AC_19','AC_20',...
                    'quantile_10',...
                    'quantile_20',...
                    'quantile_30',...
                    'quantile_40',...
                    'quantile_50',...
                    'quantile_60',...
                    'quantile_70',...
                    'quantile_80',...
                    'quantile_90'};
otherwise
    error('Unknown feature set ''%s''',whatFeatureSet);
end

%-------------------------------------------------------------------------------
% Do the matching by feature name for feature sets that are lists of feature names
%-------------------------------------------------------------------------------
if matchByName
    isMatch = cellfun(@(x) find(ismember(Operations.Name,x)),featureNames,...
                                'UniformOutput',false);
    if any(cellfun(@isempty,isMatch))
        numMissing = sum(cellfun(@isempty,isMatch));
        warning('%u feature(s) were not found in the feature set ''%s''',...
                    numMissing,whatFeatureSet);
    end
    isMatch = [isMatch{:}];
    opIDs = Operations.ID(isMatch);
    fprintf(1,'Matched %u/%u features from %s!\n',length(opIDs),length(featureNames),whatFeatureSet);

    if length(opIDs) < length(featureNames)
        didNotMatch = find(~isMatch);
        for i = 1:length(didNotMatch)
            if iscell(featureNames{didNotMatch(i)})
                theFeatureName = featureNames{didNotMatch(i)}{1};
            else
                theFeatureName = featureNames{didNotMatch(i)};
            end
            fprintf(1,'''%s'' does not exist in this HCTSA dataset\n',theFeatureName);
        end
    end
end

end
