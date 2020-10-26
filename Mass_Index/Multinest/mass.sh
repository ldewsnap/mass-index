#!/bin/bash
rm resume.dat
./mass_index
./gr  > gr.out
tput clear

tput sgr0
tput rc

awk 'NR==1{printf ("tit = '"'"'s = %.3f @_{%.3f}^{+%.3f} '"'"'\n ", ($1-1), $2, $3)}' gr.out > plot_fgh

awk 'NR==1{a=$1;next}NR==2{b=$1;next}END{printf ("f(x) = %.3f * x + %.3f \n ", a,  b)}' gr.out >> plot_fgh

awk 'NR==1{printf ("%.3f,%.3f,%.3f", ($1-1), $2, $3)}' gr.out > results.out

tail plot_data -n1 | awk '{ print "a="$2}' >> plot_fgh
cat plot_temp >> plot_fgh

gnuplot plot_fgh

convert -quality 100 -density 150 -rotate 90 plot.eps plot.jpg
convert -quality 100 -density 150 -rotate 90 pgplot.ps pgplot.jpg

# evince plot.eps &
