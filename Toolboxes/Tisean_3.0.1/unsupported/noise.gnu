! echo ""
! echo "**********************************************************************"
! echo "*  Noise reduction. First: noisy data                                *"
! echo "**********************************************************************"
! echo ""
pause 1

set title "Noisy data"
set data style dots
plot '< henon -l5000 | addnoise -v0.02 | embed -d1'
pause -1 "Press <return> to continue"

! echo ""
! echo "**********************************************************************"
! echo "*  Noise reduction in a data stream                                  *"
! echo "**********************************************************************"
! echo ""
pause 1

set title "Noise reduction in a data stream"
plot '< henon -l5000 | addnoise -v0.02 | noise -m7 -q2 -r0.3 -K500 | noise -m7 -q2 -r0.2 -K500 | embed -d1'
pause -1 "Press <return> when finished"
