! echo ""
! echo "**********************************************************************"
! echo "*  Correlation dimension                                             *"
! echo "**********************************************************************"
! echo ""
pause 1

set title "Correlation dimension"
set logscale x
set yran [0:3]
set title "local slopes"
plot '< henon -l5000 | d2 -d1 -t50; cat stdin.d2' with lines, 1.2
pause -1 "Press <return> when finished"
set title "local slopes, smoothed"
plot '< henon -l5000 | d2 -d1 -t50; av-d2 -a5 stdin.d2' with lines, 1.2
pause -1 "Press <return> when finished"
set title "Takens' estimator"
plot '< henon -l5000 | d2 -d1 -t50; c2t stdin.c2' with lines, 1.2
pause -1 "Press <return> when finished"
set title "Gaussian kernel"
plot '< henon -l5000 | d2 -d1 -t50; c2g stdin.c2' u 1:3 with lines, 1.2
pause -1 "Press <return> when finished"
! rm stdin.*
