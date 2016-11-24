function [R, keepers] = BF_RemoveNaN_DistMat(R)
% BF_RemoveNaN_DistMat     Removes NaNs from an input distance matrix

% ------------------------------------------------------------------------------
% Copyright (C) 2016, Ben D. Fulcher <ben.d.fulcher@gmail.com>,
% <http://www.benfulcher.com>
%
% If you use this code for your research, please cite:
% B. D. Fulcher, M. A. Little, N. S. Jones, "Highly comparative time-series
% analysis: the empirical structure of time series and their methods",
% J. Roy. Soc. Interface 10(83) 20130048 (2013). DOI: 10.1098/rsif.2013.0048
%
% This work is licensed under the Creative Commons
% Attribution-NonCommercial-ShareAlike 4.0 International License. To view a copy of
% this license, visit http://creativecommons.org/licenses/by-nc-sa/4.0/ or send
% a letter to Creative Commons, 444 Castro Street, Suite 900, Mountain View,
% California, 94041, USA.
% ------------------------------------------------------------------------------

keepers = true(length(R),1);

if any(isnan(R(:)))
    areNaN = 1;

    while areNaN
        NumNaNs = sum(isnan(R(keepers,keepers)));
        [~,irem] = max(NumNaNs); % the index (of keepers==1) that has the most NaNs
        fkeep = find(keepers);
        keepers(fkeep(irem)) = 0; % remove this index from the full list, keepers
        R_keep_tmp = R(keepers,keepers);
        areNaN = any(isnan(R_keep_tmp(:))); % are there still more NaNs after removing this?
    end
    R = R(keepers,keepers);
end

end
