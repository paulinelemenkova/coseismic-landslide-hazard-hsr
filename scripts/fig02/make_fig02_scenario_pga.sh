set -e

R=28.4/31.15/40.30/41.15

CONST=3.281234; SLOPE=-0.79016; B6SQ=61.88
B7=0.08753;     B8=0.01527

gmt grdmath -R$R -I0.005 -fg rupture.txt LDIST = rdist.nc
gmt grdmath rdist.nc SQR $B6SQ ADD SQRT LOG10 $SLOPE MUL $CONST ADD 10 EXCH POW 980.665 DIV = pga.nc

VS30=vs30_bestfit.tif
if [ -f "$VS30" ]; then
  echo "applying site amplification from $VS30"
  gmt grdsample "$VS30" -R$R -I0.005 -Gvs30r.nc
  gmt grdmath vs30r.nc 760 DENAN = vs30r.nc
  gmt grdmath vs30r.nc 360 LT = ss.nc
  gmt grdmath vs30r.nc 360 GE vs30r.nc 750 LT MUL = sa.nc
  gmt grdmath ss.nc $B7 MUL sa.nc $B8 MUL ADD 10 EXCH POW pga.nc MUL = pga.nc
fi

gmt math -T1/160/1 T SQR $B6SQ ADD SQRT LOG10 $SLOPE MUL $CONST ADD 10 EXCH POW 980.665 DIV = med.txt
awk '{print $1, $2*1.902}' med.txt > up.txt
awk '{print $1, $2/1.902}' med.txt > lo.txt
RA=1/160/0.02/0.6
JA=X5.6cl/3.4cl

gmt begin fig02_scenario_pga png E600
  gmt set FONT_TITLE=15p,Helvetica-Bold,black FONT_ANNOT_PRIMARY=9p,Helvetica,black \
          FONT_LABEL=10p,Helvetica,black MAP_FRAME_TYPE=plain MAP_FRAME_PEN=1p,black \
          FORMAT_GEO_MAP=ddd:mmF MAP_GRID_PEN_PRIMARY=0.25p,gray70 \
          MAP_TITLE_OFFSET=8p MAP_ANNOT_OFFSET_PRIMARY=3p PS_CHAR_ENCODING=ISOLatin1+

  gmt grdcut @earth_relief_03s -R$R -Grelief.nc
  gmt grdgradient relief.nc -A315 -Ne0.6 -Gint.nc
  gmt makecpt -Cgray -T-2500/2600 -H > topo.cpt
  gmt grdimage relief.nc -Iint.nc -Ctopo.cpt -Q -JM15c -R$R -Bpxa0.5f0.25g0.5 -Bpya0.25f0.125g0.25 -BWeSN
  gmt coast -Dh -Slightsteelblue2 -Clightsteelblue2 -W0.4p,gray45 -A8
  gmt makecpt -Cturbo -T0.02/0.40/0.005 -H > pga.cpt
  gmt grdimage pga.nc -Cpga.cpt -t45
  gmt grdcontour pga.nc -C0.05 -A0.1+f7p,Helvetica-Bold,black -Wa0.7p,black -Wc0.3p,black@40 -Gd5.5c

  gmt plot faults_nafz.gmt -W1.1p,gray25
  gmt plot rupture.txt -W3.4p,red2
  gmt plot hsr.gmt -W3.4p,gray10
  gmt plot hsr.gmt -W1.4p,white,-
  gmt plot seis.txt -i0,1,3+s0.026 -Sc -Gwhite@55 -W0.25p,gray25
  gmt plot epicentre.txt -Sa0.6c -Gyellow -W0.9p,black

  gmt plot -Sc0.13c -Gblack -W0.4p,white <<EOF
28.979 41.008
29.917 40.766
30.403 40.780
EOF
  gmt text -N -F+f9p,Helvetica-Bold,black+jLB -Gwhite@25 -W0.1p,gray70 <<EOF
29.005 41.020 Istanbul
29.955 40.775 Izmit
30.420 40.792 Adapazari
EOF
  gmt text -N -F+f8.5p,Helvetica-Bold,black+jRB -Gwhite@25 <<EOF
29.845 40.760 epicentre
EOF
  gmt text -N -F+f9p,Helvetica-BoldOblique,red3+jCB+a2 -Gwhite@25 <<EOF
30.42 40.672 1999 Izmit surface rupture (Mw 7.6)
EOF
  gmt text -N -F+f8.5p,Helvetica-BoldOblique,black+jCB+a-22 -Gwhite@20 <<EOF
29.31 40.905 HSR corridor
EOF
  gmt text -N -F+f11p,Helvetica-Oblique,navy+jCM <<EOF
28.70 40.767 Marmara Sea
EOF
  gmt text -N -F+f13p,Helvetica-Bold,black+jLT -Gwhite@15 -W0.2p,gray55 <<EOF
28.44 41.125 (a)
EOF

  gmt basemap -LjBL+w40k+o0.5c/0.55c+f+u --FONT_LABEL=8p --FONT_ANNOT_PRIMARY=7p
  gmt basemap -Tdg30.98/41.05+w0.7c+f2+l,,,N --FONT_TITLE=8p
  gmt basemap -B+t"Scenario ground motion \055 1999 Izmit (Kocaeli) Mw 7.6 rupture"

  gmt colorbar -Cpga.cpt -DJBC+jTC+o0c/0.85c+w11c/0.32c+h \
      -Bxa0.1f0.05+l"Median scenario PGA (g)" --FONT_LABEL=9p --FONT_ANNOT_PRIMARY=8p

  gmt text -N -F+f7.3p,Helvetica-Oblique,gray25+jTC <<EOF
29.775 40.020 PGA: median, Akkar & Bommer (2010) GMPE (rock, Vs30 = 760 m/s), from Joyner-Boore distance to the finite rupture.
29.775 39.980 Rupture after Barka et al. (2002); faults: GEM Global Active Faults; relief: SRTM 3 arc-sec; seismicity: IEB catalogue (M >= 4). Source: authors.
EOF

  gmt basemap -R$RA -J$JA -X16.4c -Y2.0c \
      -Bpxa1f3+l"R@-jb@- (km)" -Bpya1f3+l"PGA (g)" -BWSne \
      --FONT_LABEL=9p --FONT_ANNOT_PRIMARY=8p --MAP_FRAME_PEN=0.9p,gray30 \
      --MAP_GRID_PEN_PRIMARY=0.2p,gray88
  printf "1 0.02\n20 0.02\n20 0.6\n1 0.6\n" | gmt plot -R$RA -J$JA -Gwheat@55 -L
  gmt plot up.txt  -R$RA -J$JA -W0.9p,gray45,-
  gmt plot lo.txt  -R$RA -J$JA -W0.9p,gray45,-
  gmt plot med.txt -R$RA -J$JA -W2.6p,red2
  printf "4.6 0.5 HSR corridor\n"           | gmt text -R$RA -J$JA -N -F+f7p,Helvetica-Oblique,gray35+jCB
  printf "150 0.5 Akkar & Bommer (2010)\n"  | gmt text -R$RA -J$JA -N -F+f8p,Helvetica-Bold,black+jRT
  printf "150 0.36 Mw 7.6, rock\n"          | gmt text -R$RA -J$JA -N -F+f7.2p,Helvetica-Oblique,gray30+jRT
  printf "150 0.27 median +/- 1 s.d.\n"     | gmt text -R$RA -J$JA -N -F+f7.2p,Helvetica-Oblique,gray30+jRT
  printf "1 0.6 (b)  Ground-motion attenuation\n" | gmt text -R$RA -J$JA -N -F+f11p,Helvetica-Bold,black+jLB -Dj0c/0.35c
gmt end
echo "wrote fig02_scenario_pga.png"
