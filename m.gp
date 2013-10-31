set terminal cga 
set output 
set noclip points
set clip one
set noclip two
set border
set dummy x,y
set format x "%g"
set format y "%g"
set format z "%g"
set nogrid
set nokey
set nolabel
set noarrow
set nologscale
set offsets 0, 0, 0, 0
set nopolar
set angles radians
set noparametric
set view 60, 30, 1, 1
set samples 100
set isosamples 10
set surface
set nocontour
set nohidden3d
set cntrparam order 4
set cntrparam linear
set cntrparam points 5
set size 1,1
set data style points
set function style lines
set tics in
set ticslevel 0.5
set xtics
set ytics
set ztics
set title "Ohybove momenty" 0,0
set notime
set rrange [-0 : 10]
set trange [-5 : 5]
set xlabel "" 0,0
set xrange [0 : 12]
set ylabel "" 0,0
set yrange [0 : 7.787]
set zlabel "" 0,0
set zrange [-10 : 10]
set autoscale r
set autoscale t
set autoscale xy
set autoscale z
set zero 1e-08
plot "c" using 2:3 w l,"c" using 8:9 w l
