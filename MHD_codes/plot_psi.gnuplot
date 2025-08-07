set title "Contours of Poloidal Flux Ψ(R,Z) "
set xlabel "R"
set ylabel "Z"

set size ratio -1
set view map
unset surface

set pm3d map
set pm3d implicit

set contour base
set cntrparam levels auto 20
unset clabel  # No contour labels

# Define white contour line style
set style line 1 lc rgb "white" lw 1.5

unset key

# Plot pm3d surface and contour lines separately
splot 'solovev1.dat' using 1:2:3 with pm3d notitle, \
      '' using 1:2:3 with lines ls 1 notitle