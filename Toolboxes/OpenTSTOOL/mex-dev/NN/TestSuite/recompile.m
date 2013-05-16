% recompile mex code to have actual version of the NN code

clear functions

eval('mex -I.. -I../.. ../fnearneigh.cpp -O -output fnearneigh');
eval('mex -I.. -I../.. ../fnearneigh.cpp -O -output fnearneigh_partial -DPARTIAL_SEARCH');
eval('mex -I.. -I../.. ../fnearneigh.cpp -O -DBRUTE -output nn_brute');
eval('mex -I.. -I../.. ../fnearneigh.cpp -O -DBRUTE -output nn_brute_partial -DPARTIAL_SEARCH');

eval('mex -I.. -I../.. ../nn_prepare.cpp -O -output nn_prepare');

eval('mex -I.. -I../.. ../nn_search.cpp -O -output nn_search');
eval('mex -I.. -I../.. ../nn_search.cpp -O -output nn_search_partial -DPARTIAL_SEARCH');

eval('mex -I.. -I../.. ../range_search.cpp -O -output range_search');
eval('mex -I.. -I../.. ../range_search.cpp -O -output range_search_partial -DPARTIAL_SEARCH');

eval('mex -I.. -I../.. ../corrsum.cpp -O -output corrsum');
eval('mex -I.. -I../.. ../corrsum.cpp -O -output corrsum_partial -DPARTIAL_SEARCH');
