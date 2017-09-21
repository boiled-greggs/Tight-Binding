#    G N U P L O T
#    Linux version 3.5 (pre 3.6)
#    patchlevel beta 340
#    last modified Tue Nov 25 22:57:44 GMT 1997
#
#    Copyright(C) 1986 - 1993, 1997
#    Thomas Williams, Colin Kelley and many others
#
#    Send comments and requests for help to info-gnuplot@dartmouth.edu
#    Send bugs, suggestions and mods to bug-gnuplot@dartmouth.edu
#
set terminal x11 0
set output
set noclip points
set clip one
set noclip two
set bar 1.000000
set border 31
set xdata
set ydata
set zdata
set x2data
set y2data
set boxwidth
set dummy x,y
set format x "%g"
set format y "%.2f"
set format x2 "%g"
set format y2 "%g"
set format z "%g"
set nogrid
set key title ""
set key right top Right noreverse nobox
set nolabel
set noarrow
set nolinestyle
set nologscale
set offsets 0, 0, 0, 0
set pointsize 1
set encoding default
set nopolar
set angles radians
set noparametric
set view 60, 30, 1, 1
set samples 100, 100
set isosamples 10, 10
set surface
set nocontour
set clabel '%8.3g'
set nohidden3d
set cntrparam order 4
set cntrparam linear
set cntrparam levels auto 5
set cntrparam points 5
set size ratio 0 1,1
set origin 0,0
set data style points
set function style lines
set tics in
set ticslevel 0.5
set tics scale 1
set mxtics default
set mytics default
set mx2tics default
set my2tics default
set xtics border mirror norotate 
set ytics border mirror norotate 
set ztics border nomirror norotate 
set nox2tics
set noy2tics
set title "TB Equation of State for Palladium" 0.000000,0.000000
set timestamp "" bottom norotate 0.000000,0.000000  ""
set rrange [ * : * ] noreverse nowriteback  # (currently [-0:10] )
set trange [ * : * ] noreverse nowriteback  # (currently [-5:5] )
set urange [ * : * ] noreverse nowriteback  # (currently [-5:5] )
set vrange [ * : * ] noreverse nowriteback  # (currently [-5:5] )
set xlabel "Volume/atom (a.u.^3)" 0.000000,0.000000
#set x2label "" 0.000000,0.000000  ""
set timefmt "%d/%m/%y\n%H:%M"
set xrange [ * : * ] noreverse nowriteback  # (currently [-10:10] )
set x2range [ * : * ] noreverse nowriteback  # (currently [-10:10] )
set ylabel "Energy/atom (Ry)" 0.000000,0.000000
#set y2label "" 0.000000,0.000000  ""
set yrange [ * : * ] noreverse nowriteback  # (currently [-10:10] )
set y2range [ * : * ] noreverse nowriteback  # (currently [-10:10] )
#set zlabel "" 0.000000,0.000000  ""
set zrange [ * : * ] noreverse nowriteback  # (currently [-10:10] )
set zero 1e-08
set lmargin -1
set bmargin -1
set rmargin -1
set tmargin -1
set locale ""
min(a,b) = a<b ? a:b
max(a,b) = a>b ? a:b
step(a) = a>0 ? 1:0
t(vo,v) = (vo/v)**.666666666666666667 - 1.0
e(eo,vo,ko,kop,v) = eo + 1.125*ko*vo*t(vo,v)*t(vo,v)* (1.0 + 0.5*(kop-4.0)*t(vo,v))
ef = .001
vf = 98
kf = .01
kfp = 4
fit e(ef,vf,kf,kfp,x) "< grep fcc SKENG | awk '{print $3, $5}'" \
via ef,vf,kf,kfp
eb = .015
vb = 100
kb = .01
kbf = 4
fit e(eb,vb,kb,kbp,x) "< grep bcc SKENG | awk '{print $3, $5}'" \
via eb,vb,kb,kbp
plot "< grep fcc SKENG | awk '{print $3, $5}'" t "fcc" w p 1 1,\
e(ef,vf,kf,kfp,x) t "Birch fit" w l 1,\
"< grep bcc SKENG | awk '{print $3, $5}'" t "bcc" w p 2 2,\
e(eb,vb,kb,kbp,x) t "Birch fit" w l 2
print "Results of 3rd order Birch fit:"
print "FCC Lattice:"
print "E_0 = ",ef," Ry"
print "V_0 = ",vf," Bohr**3"
af = (4.0*vf)**.333333333333333
print "a_0 = ",af," Bohr"
kfx = 14710.5164*kf
print "B_0 = ",kfx," GPa"
print "B_0'= ",kfp
print "BCC Lattice:"
print "E_0 = ",eb," Ry"
print "V_0 = ",vb," Bohr**3"
ab = (2.0*vb)**.333333333333333
print "a_0 = ",ab," Bohr"
kbx = 14710.5164*kb
print "B_0 = ",kbx," GPa"
print "B_0'= ",kbp
deltae = 13.6056981*(eb-ef)
print "FCC - BCC energy diference = ",deltae," eV"
save "pdplot.gnu"
pause -1
#    EOF
