set term pdfcairo
set outp "Teilchent20.pdf"

set title "t=20"

plot "psi_lsg.dat" using 1 : 2 w l title "exakt" , "psi_t_x.dat" using 1 : 2 w l title "diskret"

unset term
unset outp