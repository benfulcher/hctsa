% BF_Whiten
% 
% Whitens a time series by comparing and testing a range of time-series
% preprocessings
% 
%---HISTORY:
% Ben Fulcher, 2015-03-20
% ------------------------------------------------------------------------------
% Copyright (C) 2013,  Ben D. Fulcher <ben.d.fulcher@gmail.com>,
% <http://www.benfulcher.com>
%
% If you use this code for your research, please cite:
% B. D. Fulcher, M. A. Little, N. S. Jones, "Highly comparative time-series
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

function [y,best] = BF_Whiten(y,preProc,beVocal,randomSeed)

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