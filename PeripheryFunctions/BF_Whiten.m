function [y,best] = BF_Whiten(y,preProc,beVocal,randomSeed)
% BF_Whiten     Whiten a time series by comparing a range of preprocessings.

% ------------------------------------------------------------------------------
% Copyright (C) 2016, Ben D. Fulcher <ben.d.fulcher@gmail.com>,
% <http://www.benfulcher.com>
%
% If you use this code for your research, please cite:
% B. D. Fulcher, M. A. Little, N. S. Jones, "Highly comparative time-series
% analysis: the empirical structure of time series and their methods",
% J. Roy. Soc. Interface 10(83) 20130048 (2013). DOI: 10.1098/rsif.2013.0048
%
% This work is licensed under the Creative Commons
% Attribution-NonCommercial-ShareAlike 4.0 International License. To view a copy of
% this license, visit http://creativecommons.org/licenses/by-nc-sa/4.0/ or send
% a letter to Creative Commons, 444 Castro Street, Suite 900, Mountain View,
% California, 94041, USA.
% ------------------------------------------------------------------------------

% ------------------------------------------------------------------------------
% Check inputs:
% ------------------------------------------------------------------------------
if nargin < 3
    beVocal = 1;
end

% randomSeed: how to treat the randomization
if nargin < 4
    randomSeed = [];
end
% ------------------------------------------------------------------------------
best = ''; % in case not specified

if ismember(preProc,{'nothing','none'})
    % Don't do ANYTHING to the time series
    return
end

% Detrend it first:
y = detrend(y);

switch preProc
    case 'detrend'
        % Only detrend (already done)
        return
    case 'ar'
        % Apply the preprocessing that maximizes ar(2) whiteness
        % Apply a standard AR preprocessing to the data and return them in the
        % structure ypp. Also chooses the best preprocessing based on the worst fit
        % of an AR2 model to the processed time series.
        % has to beat doing nothing by 5% (error)
        % No spectral methods allowed...
        [ypp, best] = PP_PreProcess(y,'ar',2,0.05,0,randomSeed);
        y = ypp.(best);
        if beVocal
            fprintf(1,'Preprocessed the time series according to AR(2) criterion using %s.\n',best);
        end

    otherwise
        error('Unknwon preprocessing setting ''%s''',preProc);
end

end
