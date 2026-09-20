set -e

R=28/31/40.45/41.35
J=-JM21c

gmt begin fig01_study_area png E600
  gmt set FONT_TITLE=14p,Helvetica-Bold,black FONT_ANNOT_PRIMARY=9p,Helvetica,black \
          FONT_LABEL=10p,Helvetica,black MAP_FRAME_TYPE=plain MAP_FRAME_PEN=1p,black \
          FORMAT_GEO_MAP=ddd:mmF MAP_GRID_PEN_PRIMARY=0.4p,white \
          MAP_TITLE_OFFSET=9p MAP_ANNOT_OFFSET_PRIMARY=3p PS_CHAR_ENCODING=ISOLatin1+

  gmt grdcut @earth_relief_03s -R$R -Grelief.nc
  gmt grdgradient relief.nc -A315 -Ne0.6 -Gint.nc
  gmt makecpt -Cgeo -T-200/1600 -H > topo.cpt
  gmt grdimage relief.nc -Iint.nc -Ctopo.cpt $J -R$R -Bpxa0.5f0.25 -Bpya0.25f0.125 -BWESN

  gmt coast -Df -Slightsteelblue2 -Clightsteelblue2 -W0.5p,gray35 -A8

  gmt plot provinces.gmt   -W1p,gray25,-
  gmt plot hsr.gmt         -W4p,gray10
  gmt plot hsr.gmt         -W1.6p,white,-
  gmt plot faults_nafz.gmt -W2.2p,red2
  gmt plot cities.gmt      -Sc0.19c -Gblack -W0.6p,white

  gmt basemap -Bpxg0.5 -Bpyg0.25

  gmt text -N -F+f10p,Helvetica-Bold,black+jLB -Gwhite@30 -W0.1p,gray70 <<EOF
29.010 41.020 Istanbul
29.455 40.820 Gebze
29.945 40.784 Izmit
30.420 40.795 Adapazari
EOF
  gmt text -N -F+f13p,Helvetica-BoldOblique,navy+jCB <<EOF
29.35 41.28 BLACK SEA
EOF
  gmt text -N -F+f11p,Helvetica-BoldOblique,navy+jCB <<EOF
28.45 40.72 SEA OF MARMARA
EOF
  gmt text -N -F+f9p,Helvetica-Oblique,navy+jCB <<EOF
29.48 40.755 Gulf of Izmit
30.31 40.685 Sapanca L.
EOF
  gmt text -N -F+f15p,Helvetica-Bold,gray15@45+jCB <<EOF
28.556 41.226 ISTANBUL
30.056 40.916 KOCAELI
EOF
  gmt text -N -F+f15p,Helvetica-Bold,gold1+jCB <<EOF
28.55 41.23 ISTANBUL
30.05 40.92 KOCAELI
EOF
  gmt text -N -F+f9p,Helvetica-Bold,red3+jCB+a4  <<EOF
28.50 40.865 North Marmara Fault
EOF
  gmt text -N -F+f9p,Helvetica-Bold,red3+jCB+a2  <<EOF
29.70 40.702 Izmit Fault
EOF
  gmt text -N -F+f9p,Helvetica-Bold,red3+jCB+a-6 <<EOF
30.34 40.655 Sapanca Fault
EOF

  gmt basemap -Lg28.95/40.50+w40k+f+u --FONT_LABEL=9p --FONT_ANNOT_PRIMARY=8p
  gmt basemap -Tdg28.20/41.26+w0.9c+f2+l,,,N --FONT_TITLE=9p

  gmt inset begin -DjTR+w5.0c/2.31c+o0.12c/0.12c -F+gwhite+p0.8p,gray40
    gmt coast -R25.5/45/35.5/42.5 -JM5.0c -Gwheat -Slightsteelblue2 -Df -A200 \
        -N1/0.5p,gray45 -W0.2p,gray60
    gmt plot -R25.5/45/35.5/42.5 -JM5.0c -Sr+s -W1.1p,red2 <<EOF
28 40.45 31 41.35
EOF
  gmt inset end

  cat > leg.txt <<EOF
S 0.30c - 0.8c - 2.2p,red2 1.0c NAFZ active fault (GEM GAF)
S 0.30c - 0.8c - 4p,gray10 1.0c Ankara-Istanbul HSR corridor
S 0.36c c 0.19c black 0.6p,white 1.0c City / station
EOF
  gmt legend leg.txt -DJBL+jTL+o0c/0.9c+w7.6c -F+gwhite@6 --FONT_ANNOT_PRIMARY=9p
  gmt colorbar -Ctopo.cpt -DJBR+jTR+o0c/0.95c+w7.6c/0.32c+h -F+gwhite@8 \
      -Bxa400f200+l"Elevation (m)" -G0/1600 --FONT_LABEL=9p --FONT_ANNOT_PRIMARY=8p

  gmt basemap -B+t"Study area of the Istanbul-Kocaeli high-speed railway corridor"
  gmt text -N -F+f8p,Helvetica-Oblique,gray30+jLB <<EOF
28.0 40.14 Relief: SRTM 3 arc-sec (earth_relief_03s). Faults: GEM Global Active Faults. Provinces: Natural Earth. HSR alignment: OSM/TCDD (schematic). Projection: Mercator, WGS84.
EOF
gmt end show
