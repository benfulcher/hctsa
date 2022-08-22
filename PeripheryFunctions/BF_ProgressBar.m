function BF_ProgressBar(progressProp,PB_length,theChar)
% Create/update a text progress bar.
%
%---INPUTS:
% progressProp:   EITHER Text string to initialize ('new') or terminate ('close')
%                     OR Proportion progress (numerical)
% PB_length:      Length of progress bar (in characters)
% theChar:        Character to use for progress bar

%-------------------------------------------------------------------------------
%% Initialization
%-------------------------------------------------------------------------------
persistent progressBar; % Text progress bar

% Vizualization parameters
if nargin < 2
    PB_length = 40;   %   Length of progress bar
end
if nargin < 3
    theChar = ':';    %   Character to show progress
end
notTheChar = ' ';

%-------------------------------------------------------------------------------
%% Main
%-------------------------------------------------------------------------------
if isempty(progressBar) && ischar(progressProp) && strcmp(progressProp,'new')
    % Initialize progress bar
    progressBar = -1;
elseif ~isempty(progressBar) && ischar(progressProp) && strcmp(progressProp,'close')
    % Progress bar  - termination
    fprintf(1,'\n');
    clear('progressBar')
elseif isnumeric(progressProp)
    assert(progressProp >= 0)
    assert(progressProp <= 1)

    % Check whether update to progress bar is required
    if round(progressProp*PB_length) > progressBar
        % clear current text:
        if progressBar~=-1
            fprintf(repmat('\b',1,PB_length+2))
        end

        % Update and write:
        progressBar = round(progressProp*PB_length);
        WriteProgressBar
    end
else
    % Any other unexpected input
    warning('Unexpected input ''%s''',progressProp);
    progressBar = [];
end

function WriteProgressBar()
    % Write progress bar text to commandline
    fprintf('|');
    for i = 1:progressBar
        fprintf(theChar)
    end
    for i = 1:PB_length-progressBar
        fprintf(1,notTheChar);
    end
    fprintf('|');
end

end
