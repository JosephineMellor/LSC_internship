set terminal gif size 1000, 800 animate delay 50

set output 'MUSCL2.gif'

set pm3d
unset surface
set view map 
set cbrange[0:2.5]
set xrange[0:1]
set yrange[-1:1]
do for [idx =0:121]{\
splot 'MHD'.idx.'.dat' u 1:2:8 w l notitle}