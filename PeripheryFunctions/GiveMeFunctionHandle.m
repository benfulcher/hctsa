function fn_handle = GiveMeFunctionHandle(cfnParams)
% GiveMeFunctionHandle  Returns a function handle for computing classification accuracy
%                       (fitting on training data and evaluating on test data)
%
%---INPUTS:
% cfnParams, parameters for the classification, e.g., set with
%           GiveMeDefaultClassificationParams.

% ------------------------------------------------------------------------------
% Copyright (C) 2013-2026, Ben D. Fulcher <ben.d.fulcher@gmail.com>,
% <http://www.benfulcher.com>
%
% If you use this code for your research, please cite the following two papers:
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
% ------------------------------------------------------------------------------

% The classifier is fit and evaluated by GiveMeCfn, with a single train/test
% split (no cross-validation), so the handle gives the same results as GiveMeCfn
% for every classifier. As a function of (XTrain,yTrain,XTest,yTest), it returns
% the performance measure (cfnParams.whatLoss) in the test set and, when a second
% output is requested, the trained model.
cfnParams.numFolds = 0;
fn_handle = @(XTrain,yTrain,XTest,yTest) GiveMeCfn(XTrain,yTrain,XTest,yTest,cfnParams,false);

end
