! echo ""
! echo "**********************************************************************"
! echo "*  Noise reduction. First: noisy data                                *"
! echo "**********************************************************************"
! echo ""
pause 1

set title "Noisy data"
plot '< henon -l5000 | makenoise -%2 | delay'with dots,'< henon -l5000 |delay' with dots 3
pause -1 "Press <return> to continue"

! echo ""
! echo "**********************************************************************"
! echo "*  Noise reduction                                                   *"
! echo "**********************************************************************"
! echo ""
pause 1

set title "Noise reduction"
plot '< henon -l5000 | makenoise -%2 | ghkss -m1,7 -q2 -r0.05 -k20 -i2 | delay' with dots,'< henon -l5000 |delay' with dots 3
pause -1 "Press <return> when finished"
! rm stdin.opt.[12]
