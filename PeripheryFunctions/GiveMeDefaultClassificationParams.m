function params = GiveMeDefaultClassificationParams(TimeSeries,numClasses)
%-------------------------------------------------------------------------------
% Copyright (C) 2020, Ben D. Fulcher <ben.d.fulcher@gmail.com>,
% <http://www.benfulcher.com>
%
% If you use this code for your research, please cite these papers:
%
% (1) B.D. Fulcher and N.S. Jones, "hctsa: A Computational Framework for Automated
% Time-Series Phenotyping Using Massive Feature Extraction, Cell Systems 5: 527 (2017).
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
%-------------------------------------------------------------------------------

if nargin < 2
    numClasses = max(TimeSeries.Group);
    % Assumes group in form of integer class labels starting at 1
end
% Number of classes to classify
params.numClasses = numClasses;

% Set the classifier:
params.whatClassifier = 'svm_linear';

% Number of repeats of cross-validation:
% (reduce variance due to 'lucky splits')
params.numRepeats = 2;

% Number of folds:
params.numFolds = HowManyFolds(TimeSeries.Group,numClasses);

% Balance weighting
params.doReweight = true;

% Whether to output information about each fold, or average over folds
params.computePerFold = false;

end
