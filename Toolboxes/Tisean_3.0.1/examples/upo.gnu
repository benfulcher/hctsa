! echo ""
! echo "**********************************************************************"
! echo "*  Orbits of period six (or subperiod)                               *"
! echo "**********************************************************************"
! echo ""
pause 1

set pointsize 3
set key default
set title "Orbits of period six (or subperiod)"
! henon -l1000 | makenoise -%10. | upo -p6 -m2 -v0.1 -n100 -V -o
plot '< henon -l1000 | makenoise -%10. | delay' with dots notitle,\
     '< cat stdin_upo_06 | upoembed -p1 -d1' with points title "fixed point",\
     '< cat stdin_upo_06 | upoembed -p2 -d1' with points title "period 2",\
     '< cat stdin_upo_06 | upoembed -p6 -d1' index 0 with points title "period 6",\
     '< cat stdin_upo_06 | upoembed -p6 -d1' index 1 with points title "period 6"
pause -1 "Press <return> when finished"
! rm stdin_*
