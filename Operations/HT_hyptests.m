function p = HT_hyptests(x,ange)
% output is p-value from one of a set of standard statistical hypothesis tests
% Ben Fulcher

switch ange
    case 'signtest'
        [p h] = signtest(x);
        
    case 'runstest'
        [h p] = runstest(x);

    case 'vartest'
        [h p] = vartest(x,1); % normal distribution of variance 1
        
    case 'ztest'
        [h p] = ztest(x,0,1);
        
    case 'signrank'
        [p h] = signrank(x);
        
    case 'jb'
        [h p] = jbtest(x);
        
    case 'lbq'
        [h p] = lbqtest(x);
end

end