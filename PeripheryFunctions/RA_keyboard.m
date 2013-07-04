function RA_keyboard
% Romesh Abeysuriya's replacement of Matlab's 'keyboard' command
% Keyboard debug caller 
% Provides more information including the stack trace than simply using 
% 'keyboard'
% Use RA_keyboard() directly instead of keyboard()
% Note this function will not work if it is the last line in the program
% By Romesh Abeysuriya 15-11-12
    
[stack_trace] = dbstack;
fprintf(1,'\n')
for j = 2:length(stack_trace)
    fprintf('%s:%i << ',stack_trace(j).name,stack_trace(j).line);
end
fprintf('\b\b\b\b\n');
eval(sprintf('dbstop in %s at %u',stack_trace(2).name,stack_trace(2).line+1))

end