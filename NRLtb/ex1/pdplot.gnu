#!/usr/local/Cellar/gnuplot/5.0.6/bin/gnuplot -persist
#
#    
#    	G N U P L O T
#    	Version 5.0 patchlevel 6    last modified 2017-03-18
#    
#    	Copyright (C) 1986-1993, 1998, 2004, 2007-2017
#    	Thomas Williams, Colin Kelley and many others
#    
#    	gnuplot home:     http://www.gnuplot.info
#    	faq, bugs, etc:   type "help FAQ"
#    	immediate help:   type "help"  (plot window: hit 'h')
# set terminal unknown 
# set output
unset clip points
set clip one
unset clip two
set bar 1.000000 front
set border 31 front lt black linewidth 1.000 dashtype solid
set zdata 
set ydata 
set xdata 
set y2data 
set x2data 
set boxwidth
set style fill  empty border
set style rectangle back fc  bgnd fillstyle   solid 1.00 border lt -1
set style circle radius graph 0.02, first 0.00000, 0.00000 
set style ellipse size graph 0.05, 0.03, first 0.00000 angle 0 units xy
set dummy x, y
set format x "%g" 
set format y "%.2f" 
set format x2 "%g" 
set format y2 "%g" 
set format z "%g" 
set format cb "% h" 
set format r "% h" 
set timefmt "%d/%m/%y
%H:%M"
set angles radians
set tics back
unset grid
set raxis
set style parallel front  lt black linewidth 2.000 dashtype solid
set key title "" center
set key inside right top vertical Right noreverse enhanced autotitle nobox
set key noinvert samplen 4 spacing 1 width 0 height 0 
set key maxcolumns 0 maxrows 0
set key noopaque
unset label
unset arrow
set style increment default
unset style line
unset style arrow
set style histogram clustered gap 2 title textcolor lt -1
unset object
set style textbox transparent margins  1.0,  1.0 border
unset logscale
set offsets 0, 0, 0, 0
set pointsize 1
set pointintervalbox 1
set encoding default
unset polar
unset parametric
unset decimalsign
unset micro
unset minussign
set view 60, 30, 1, 1
set samples 100, 100
set isosamples 10, 10
set surface 
unset contour
set cntrlabel  format '%8.3g' font '' start 5 interval 20
set mapping cartesian
set datafile separator whitespace
unset hidden3d
set cntrparam order 4
set cntrparam linear
set cntrparam levels auto 5
set cntrparam points 5
set size ratio 0 1,1
set origin 0,0
set style data points
set style function lines
unset xzeroaxis
unset yzeroaxis
unset zzeroaxis
unset x2zeroaxis
unset y2zeroaxis
set xyplane relative 0.5
set tics scale  1, 0.5, 1, 1, 1
set mxtics default
set mytics default
set mztics default
set mx2tics default
set my2tics default
set mcbtics default
set mrtics default
set xtics border in scale 1,0.5 mirror norotate  autojustify
set xtics  norangelimit autofreq 
set ytics border in scale 1,0.5 mirror norotate  autojustify
set ytics  norangelimit autofreq 
set ztics border in scale 1,0.5 nomirror norotate  autojustify
set ztics  norangelimit autofreq 
unset x2tics
unset y2tics
set cbtics border in scale 1,0.5 mirror norotate  autojustify
set cbtics  norangelimit autofreq 
set rtics axis in scale 1,0.5 nomirror norotate  autojustify
set rtics  norangelimit autofreq 
set paxis 1 tics border in scale 1,0.5 nomirror norotate  autojustify
set paxis 1 tics  rangelimit autofreq 
set paxis 2 tics border in scale 1,0.5 nomirror norotate  autojustify
set paxis 2 tics  rangelimit autofreq 
set paxis 3 tics border in scale 1,0.5 nomirror norotate  autojustify
set paxis 3 tics  rangelimit autofreq 
set paxis 4 tics border in scale 1,0.5 nomirror norotate  autojustify
set paxis 4 tics  rangelimit autofreq 
set paxis 5 tics border in scale 1,0.5 nomirror norotate  autojustify
set paxis 5 tics  rangelimit autofreq 
set paxis 6 tics border in scale 1,0.5 nomirror norotate  autojustify
set paxis 6 tics  rangelimit autofreq 
set paxis 7 tics border in scale 1,0.5 nomirror norotate  autojustify
set paxis 7 tics  rangelimit autofreq 
set title "TB Equation of State for Palladium" 
set title  font "" norotate
set timestamp bottom 
set timestamp "" 
set timestamp  font "" norotate
set rrange [ * : * ] noreverse nowriteback
set trange [ * : * ] noreverse nowriteback
set urange [ * : * ] noreverse nowriteback
set vrange [ * : * ] noreverse nowriteback
set xlabel "Volume/atom (a.u.^3)" 
set xlabel  font "" textcolor lt -1 norotate
set x2label "" 
set x2label  font "" textcolor lt -1 norotate
set xrange [ * : * ] noreverse nowriteback
set x2range [ * : * ] noreverse nowriteback
set ylabel "Energy/atom (Ry)" 
set ylabel  font "" textcolor lt -1 rotate by -270
set y2label "" 
set y2label  font "" textcolor lt -1 rotate by -270
set yrange [ * : * ] noreverse nowriteback
set y2range [ * : * ] noreverse nowriteback
set zlabel "" 
set zlabel  font "" textcolor lt -1 norotate
set zrange [ * : * ] noreverse nowriteback
set cblabel "" 
set cblabel  font "" textcolor lt -1 rotate by -270
set cbrange [ * : * ] noreverse nowriteback
set paxis 1 range [ * : * ] noreverse nowriteback
set paxis 2 range [ * : * ] noreverse nowriteback
set paxis 3 range [ * : * ] noreverse nowriteback
set paxis 4 range [ * : * ] noreverse nowriteback
set paxis 5 range [ * : * ] noreverse nowriteback
set paxis 6 range [ * : * ] noreverse nowriteback
set paxis 7 range [ * : * ] noreverse nowriteback
set zero 1e-08
set lmargin  -1
set bmargin  -1
set rmargin  -1
set tmargin  -1
set locale "en_US.UTF-8"
set pm3d explicit at s
set pm3d scansautomatic
set pm3d interpolate 1,1 flush begin noftriangles noborder corners2color mean
set palette positive nops_allcF maxcolors 0 gamma 1.5 color model RGB 
set palette rgbformulae 7, 5, 15
set colorbox default
set colorbox vertical origin screen 0.9, 0.2, 0 size screen 0.05, 0.6, 0 front  noinvert bdefault
set style boxplot candles range  1.50 outliers pt 7 separation 1 labels auto unsorted
set loadpath 
set fontpath 
set psdir
set fit brief errorvariables nocovariancevariables errorscaling prescale nowrap v5
min(a,b) = a<b ? a:b
max(a,b) = a>b ? a:b
step(a) = a>0 ? 1:0
t(vo,v) = (vo/v)**.666666666666666667 - 1.0
e(eo,vo,ko,kop,v) = eo + 1.125*ko*vo*t(vo,v)*t(vo,v)* (1.0 + 0.5*(kop-4.0)*t(vo,v))
GPFUN_min = "min(a,b) = a<b ? a:b"
GPFUN_max = "max(a,b) = a>b ? a:b"
GPFUN_step = "step(a) = a>0 ? 1:0"
GPFUN_t = "t(vo,v) = (vo/v)**.666666666666666667 - 1.0"
GPFUN_e = "e(eo,vo,ko,kop,v) = eo + 1.125*ko*vo*t(vo,v)*t(vo,v)* (1.0 + 0.5*(kop-4.0)*t(vo,v))"
ef = -0.00038649208543152
vf = 97.3474099529425
kf = 0.0157310242613339
kfp = 3.55235080916495
FIT_CONVERGED = 1
FIT_NDF = 5
FIT_STDFIT = 0.000466402372753846
FIT_WSSR = 1.08765586655209e-06
FIT_P = 1.0
FIT_NITER = 6
ef_err = 0.000130270865031668
vf_err = 0.0879944994028966
kf_err = 0.000292688624745649
kfp_err = 0.286769496149974
eb = 0.00848554892875125
vb = 99.0135973341695
kb = 0.0147522791128755
kbf = 4
kbp = 3.39317609034888
eb_err = 0.000272687997798211
vb_err = 0.139491694868377
kb_err = 0.00043362533500287
kbp_err = 0.279225281443218
af = 7.30233014881478
kfx = 231.411490385151
ab = 5.82874351192859
kbx = 217.013643827332
deltae = 0.120710311569789
## Last datafile plotted: "< grep bcc SKENG | awk '{print $3, $5}'"

## fit e(eb,vb,kb,kbp,x) "< grep bcc SKENG | awk '{print $3, $5}'" via eb,vb,kb,kbp
#    EOF
