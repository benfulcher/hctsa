function p = HT_hyptests(x,thetest)
% Returns the p-value from one of a set of standard statistical hypothesis tests
% Ben Fulcher

switch thetest
    case 'signtest' % Statistics Toolbox
        [p, ~] = signtest(x);
        % for some reason this one has p-value as the first output
        
    case 'runstest' % Statistics Toolbox
        [~, p] = runstest(x);

    case 'vartest' % Statistics Toolbox
        [~, p] = vartest(x,1); % normal distribution of variance 1
        
    case 'ztest' % Statistics Toolbox
        [~, p] = ztest(x,0,1);
        
    case 'signrank' % Statistics Toolbox
        [p, ~] = signrank(x);
        
    case 'jb' % Statistics Toolbox
        [~, p] = jbtest(x);
        
    case 'lbq' % Econometrics Toolbox
        [~, p] = lbqtest(x);
end

end