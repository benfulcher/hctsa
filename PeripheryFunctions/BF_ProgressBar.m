function BF_ProgressBar(progressProp,PB_length,theChar,suffixText)
% Create/update a text progress bar.
%
%---INPUTS:
% progressProp:   EITHER Text string to initialize ('new') or terminate ('close')
%                     OR Proportion progress (numerical)
% PB_length:      Length of progress bar (in characters)
% theChar:        Character to use for progress bar
% suffixText:     [default: ''] Text to print after the bar (e.g., an ETA string
%                     like ' ~1.2min remaining'). Right-padded to a fixed field
%                     width (SUFFIX_WIDTH) so the backspace-erase below still
%                     covers the whole previous line even as the text shrinks
%                     (e.g. as an ETA counts down); text longer than that width
%                     is left unpadded and may leave stray characters behind on
%                     the next, shorter update.
%
% In a non-interactive (-batch) session, e.g. CI or a cluster job, the animated
% bar is replaced by a plain line at each 10% step: its backspace-redraws can't
% animate in a log file, and on GitHub Actions they swallow all later output.

%-------------------------------------------------------------------------------
%% Initialization
%-------------------------------------------------------------------------------
persistent progressBar; % Text progress bar
persistent lastSuffixText; % Suffix text last written (to detect a live-updating ETA even between bar ticks)
persistent lastLineStep; % (non-interactive mode) last 10% step written as a line

SUFFIX_WIDTH = 45; % fixed field width for suffixText -- see note above

% Vizualization parameters
if nargin < 2 || isempty(PB_length)
    PB_length = 40;   %   Length of progress bar
end
if nargin < 3 || isempty(theChar)
    theChar = ':';    %   Character to show progress
end
if nargin < 4 || isempty(suffixText)
    suffixText = '';
end
notTheChar = ' ';
lineMode = exist('batchStartupOptionUsed','builtin') && batchStartupOptionUsed; % non-interactive: no animation

%-------------------------------------------------------------------------------
%% Main
%-------------------------------------------------------------------------------
if isempty(progressBar) || (ischar(progressProp) && strcmp(progressProp,'new'))
    % Initialize progress bar
    progressBar = -1;
    lastSuffixText = '';
    lastLineStep = -1;
elseif ischar(progressProp) && strcmp(progressProp,'close')
    % Progress bar - termination (lines are already terminated in lineMode)
    if ~lineMode
        fprintf(1,'\n');
    end
    clear('progressBar','lastSuffixText','lastLineStep')
elseif isnumeric(progressProp) && lineMode
    assert(progressProp >= 0)
    assert(progressProp <= 1)

    % Write a full line only when progress reaches a new 10% step
    lineStep = floor(10*progressProp);
    if lineStep > lastLineStep
        progressBar = round(progressProp*PB_length);
        lastLineStep = lineStep;
        WriteProgressBar(true)
    end
elseif isnumeric(progressProp)
    assert(progressProp >= 0)
    assert(progressProp <= 1)

    % Check whether update to progress bar is required (bar has ticked over, or
    % the suffix text -- e.g. a live ETA -- has changed since the last redraw)
    newTick = round(progressProp*PB_length);
    if newTick > progressBar || progressBar == -1 || ~strcmp(suffixText,lastSuffixText)
        % clear current text:
        if progressBar~=-1
            fprintf(repmat('\b',1,PB_length+2+SUFFIX_WIDTH))
        end

        % Update and write:
        progressBar = newTick;
        lastSuffixText = suffixText;
        WriteProgressBar(false)
    end
else
    % Any other unexpected input
    warning('Unexpected input ''%s''.',progressProp);
    progressBar = [];
end

function WriteProgressBar(asLine)
    % Write progress bar text to commandline in a single call (rather than
    % looping fprintf per character); asLine: as a whole line, without padding
    barText = ['|',repmat(theChar,1,progressBar),repmat(notTheChar,1,PB_length-progressBar),'|'];
    if asLine
        fprintf(1,'%s%s\n',barText,suffixText);
    else
        fprintf(1,'%s%-*s',barText,SUFFIX_WIDTH,suffixText);
    end
end

end
