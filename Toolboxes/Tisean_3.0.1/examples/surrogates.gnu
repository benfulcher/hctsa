! echo ""
! echo "**********************************************************************"
! echo "*  Surrogate data                                                    *"
! echo "**********************************************************************"
! echo ""
pause 1

set title "Surrogate data"
set yran [-400:200]
set ytic (-50, 0, 50, "-50" -250, "0" -200, "50" -150)
plot 'spike.dat' with lines, '< surrogates spike.dat' u 0:($1-200) with lines
pause -1 "Press <return> to finish"
