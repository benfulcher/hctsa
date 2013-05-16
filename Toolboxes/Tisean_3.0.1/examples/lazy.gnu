! echo ""
! echo "**********************************************************************"
! echo "*  Noise reduction. First: noisy data                                *"
! echo "**********************************************************************"
! echo ""
pause 1

set title "Noisy data"
plot '< henon -l5000 | makenoise -%2 | delay' with dots
pause -1 "Press <return> to continue"

! echo ""
! echo "**********************************************************************"
! echo "*  Simple noise reduction                                            *"
! echo "**********************************************************************"
! echo ""
pause 1
set title "Simple noise reduction"
plot '< henon -l5000 | makenoise -%2 | lazy -m7 -r0.06 -i3 | delay' with dots
pause -1 "Press <return> when finished"
