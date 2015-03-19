% ------------------------------------------------------------------------------
% crt_old_new_mets_man.m
% ------------------------------------------------------------------------------
% 
% Script to create a cell array containing manually adjusted old now map
% for dunctions. Necessary where function calls have changed, not only
% funtion names
% ------------------------------------------------------------------------------

function man_map = crt_old_new_mets_man()

man_map = { 'DV_dynsys_dblwell(y,[1,0.5,0.2]).mean','PH_ForcePotential(y,''dblwell'',[1,0.5,0.2])','.mean';
            'DV_dynsys_dblwell(y,[2,0.05,0.2]).ac50','PH_ForcePotential(y,''dblwell'',[2,0.05,0.2])','.ac50';
            'DV_dynsys_sine(y,[10,0.04,10]).mean','PH_ForcePotential(y,''sine'',[10,0.4,10])','.mean'; 
            'DV_dynsys_sine(y,[10,0.04,10]).std','PH_ForcePotential(y,''sine'',[10,0.4,10])','.std';
            'DF_mlefits(x,1,2)','DN_Fit_mle(y,''gaussian'')','.std';
            'SY_StatAv(y,50,''len'')','SY_StatAv(y,''len'',50)','';
            'SY_StatAv(y,100,''len'')','SY_StatAv(y,''len'',100)','';
            'ST_cumul(y,''mean'',''randcg'',500)','SY_LocalGlobal(y,''randcg'',500)','.mean' ;           
            'ST_cumul(y,''median'',''p'',0.01)','SY_LocalGlobal(y,''p'',0.01)','.median';      
            'ST_cumul(y,''iqr'',''l'',50)','SY_LocalGlobal(y,''l'',50)','.iqr';      
            'ST_cumul(y,''kurtosis'',''p'',0.1)','SY_LocalGlobal(y,''p'',0.1)','.kurtosis';             
            'ST_cumul(y,''AC1'',''l'',50)','SY_LocalGlobal(y,''l'',50)','.ac1';   
            'ST_cumul(y,''iqr'',''p'',0.005)','SY_LocalGlobal(y,''p'',0.005)','.iqr'; 
            'ST_cumul(y,''skewness'',''randcg'',100)','SY_LocalGlobal(y,''randcg'',100)','.skewness';    
            'SC_benfa(y,2,''range'',1).r2_se1','SC_FluctAnal(y,2,''range'',1,[],[],0)','.r2_se1';
            'SC_benfa(y,2,''rsrange'',1).r1_stats_coeffcorr','SC_FluctAnal(y,2,''rsrange'',1,[],[],0)','.stats_coeffcorr';
            'SC_benfa(y,2,''rsrangefit'',2,1).se2','SC_FluctAnal(y,2,''rsrangefit'',2,1,[],0)','.se2';
            'SC_benfa(y,2,''dfa'',2,3).resac1','SC_FluctAnal(y,2,''dfa'',2,3,[],0)','.resac1';
            };
end

