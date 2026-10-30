! echo ""
! echo "**********************************************************************"
! echo "*  Information dimension                                             *"
! echo "**********************************************************************"
! echo ""
pause 1

set title "Information dimension"
set logscale x
set yran [0:3]
plot '< henon -l5000 | c1 -m1 -M6 -d1 -t50 -n500 ; c2d -a2 stdin_c1' with lines, 1.2
pause -1 "Press <return> when finished"
! rm stdin_*
