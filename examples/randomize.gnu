! echo ""
! echo "**********************************************************************"
! echo "*  Constrained randomisation                                         *"
! echo "**********************************************************************"
! echo ""
pause 1

set title "Constrained randomisation"
set yran [-400:200]
set ytic (-50, 0, 50, "-50" -250, "0" -200, "50" -150)
! randomize_auto_exp_random -D50 -T0.01 -u0.95 -C0.001 spike.dat -o
plot 'spike.dat', 'spike.dat_rnd' u ($1-200) with lines
pause -1 "Press <return> to finish"
! rm -f spike.dat_rnd
