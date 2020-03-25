#set term pdfcairo
#set outp "vergleich.pdf"

plot "psi_lsg.dat" using 1 : 2 w l #, "psi_t_x.dat" using 1 : 2 w l

pause -1
#unset term
#unset outp