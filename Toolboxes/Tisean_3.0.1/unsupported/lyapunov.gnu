! echo ""
! echo "**********************************************************************"
! echo "*  Lyapunov exponent                                                 *"
! echo "**********************************************************************"
! echo ""
pause 1

set title "Lyapunov exponent"
set data style lines
set logscale y
set yrange [0.001:1]
plot '< henon -l10000 | lyapunov -d1 -m2 -M3 -t50 -r0.005 -s20',  exp(0.42*x)/400
pause -1 "Press <return> when finished"
