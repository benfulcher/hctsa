function out = EN_SampEn(y,M,r,preProcessHow)
% EN_SampEn     Sample Entropy of a time series
%
% SampEn(m,r), using code from PhysioNet.
% Uses a compiled C version of the code if available, otherwise uses a (slower)
% Matlab implementation (which can actually be faster for shorter time series
% due to overheads of reading/writing to disk)
%
% The publicly-available PhysioNet Matlab code, sampenc (renamed here to
% RN_sampenc) is available from:
% http://www.physionet.org/physiotools/sampen/matlab/1.1/sampenc.m
%
% cf. "Physiological time-series analysis using approximate entropy and sample
% entropy", J. S. Richman and J. R. Moorman, Am. J. Physiol. Heart Circ.
% Physiol., 278(6) H2039 (2000)
%
% This function can also calculate the SampEn of successive increments of time
% series, i.e., we using an incremental differencing pre-processing, as
% used in the so-called Control Entropy quantity:
%
% "Control Entropy: A complexity measure for nonstationary signals"
% E. M. Bollt and J. Skufca, Math. Biosci. Eng., 6(1) 1 (2009)
%
%---INPUTS:
% y, the input time series
% M, the embedding dimension
% r, the threshold
% preProcessHow [opt], (i) 'diff1', incremental differencing (as per 'Control Entropy').

% ------------------------------------------------------------------------------
% Copyright (C) 2015, Ben D. Fulcher <ben.d.fulcher@gmail.com>,
% <http://www.benfulcher.com>
%
% If you use this code for your research, please cite:
% B. D. Fulcher, M. A. Little, N. S. Jones, "Highly comparative time-series
% analysis: the empirical structure of time series and their methods",
% J. Roy. Soc. Interface 10(83) 20130048 (2013). DOI: 10.1098/rsif.2013.0048
%
% This function is free software: you can redistribute it and/or modify it under
% the terms of the GNU General Public License as published by the Free Software
% Foundation, either version 3 of the License, or (at your option) any later
% version.
%
% This program is distributed in the hope that it will be useful, but WITHOUT
% ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
% FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
% details.
%
% You should have received a copy of the GNU General Public License along with
% this program. If not, see <http://www.gnu.org/licenses/>.
% ------------------------------------------------------------------------------

% Check y is a column vector:
if size(y,1)==1
    y = y';
end

% Embedding dimension:
if nargin < 2
    M = 2;
end

% Tolerance:
if nargin < 3
    r = 0.1*std(y);
end

if nargin < 4
    preProcessHow = ''; % don't apply any preprocessing
end
%-------------------------------------------------------------------------------
% Can specify to first apply an incremental differencing of the time series
% thus yielding the 'Control Entropy':
% "Control Entropy: A complexity measure for nonstationary signals"
% E. M. Bollt and J. Skufca, Math. Biosci. Eng., 6(1) 1 (2009)
if ~isempty(preProcessHow)
    y = BF_preprocess(y,preProcessHow);
end

% ------------------------------------------------------------------------------
% Use the physionet code to calculate the Sample Entropy using these parameters:
% ------------------------------------------------------------------------------
% Check if a compiled C version exists:
% [a,b] = system('which sampen');

try
    sampEn = sampen_mex(y',M+1,r);
    sampEn = sampEn(1:M+1); % always that extra one for the M=0
catch
    warning('No mex file found: using a slower native Matlab implementation instead');
    % No mex version available; use (much slower) Matlab implementation
    % if isempty(b) || length(y) < 3000 % faster to run within Matlab
    % No compiled C version detected (for time series longer than 3000 samples):
    % run in Matlab (slower):
    % (length of 2000 because of read/write overhead)
    sampEn = PN_sampenc(y,M+1,r);
    % else
    %     fprintf('Using compiled C code~~~\n')
    %     % http://www.physionet.org/physiotools/sampen/c/
    %     % (use Makefile in Toolboxes/Physionet/ to run make, then make install)
    %
    %     % Run compiled C code:
    %     filePath = BF_WriteTempFile(y);
    %     command = sprintf('sampen -m %u -r %f < %s',M,r,filePath);
    %     [~,res] = system(command);
    %     fprintf(1,'%s\n',res)
    %     s = textscan(res,'%[^\n]'); s = s{1};
    %     sampEn = zeros(M,1);
    %     for i = 1:M
    %         [~,params] = regexp(s{i},'\((\S+)\)','tokens','match');
    %         params = regexp(params{1}(2:end-1),',','split');
    %         [~,result] = regexp(s{i},'= (\S+)','tokens','match');
    %         result = str2num(result{1}(3:end));
    %         sampEn(i) = result;
    %     end
    % end
end

% ------------------------------------------------------------------------------
% Compute outputs from the code
% ------------------------------------------------------------------------------
for i = 1:M+1
    % Sample entropy:
    out.(sprintf('sampen%u',i-1)) = sampEn(i);

    % Quadratic sample entropy (QSE), Lake (2006):
    % (allows better comparison across r values)
    out.(sprintf('quadSampEn%u',i-1)) = sampEn(i) + log(2*r);

    % COSEn (Lake and Moorman, 2011), doesn't really make sense in general;
    % especially for z-scored series!:
    % out.(sprintf('COSEn%u',i)) = sampEn(i) + log(2*r) - log(mean(y));
end

if M > 1
    out.meanchsampen = mean(diff(sampEn));
end

end
