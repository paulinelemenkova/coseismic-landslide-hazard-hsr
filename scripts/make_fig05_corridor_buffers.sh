set -e
cd /tmp
python3 - <<'PY'
import urllib.request, urllib.parse, json, math, time
BBOX=(40.55,28.85,41.15,30.80)
EPS=["https://overpass-api.de/api/interpreter","https://overpass.kumi.systems/api/interpreter",
     "https://maps.mail.ru/osm/tools/overpass/api/interpreter","https://overpass.openstreetmap.ru/api/interpreter"]
def overpass(q):
    last=None
    for ep in EPS:
        for _ in range(2):
            try:
                req=urllib.request.Request(ep,data=urllib.parse.urlencode({"data":q}).encode(),headers={"User-Agent":"corridor-map/1.0"})
                return json.loads(urllib.request.urlopen(req,timeout=120).read().decode())
            except Exception as e: last=e; time.sleep(3)
    raise last
def rdp(p,eps):
    if len(p)<3: return p
    a,b=p[0],p[-1]; dm=0; k=0
    for i in range(1,len(p)-1):
        x0,y0=p[i]
        num=abs((b[1]-a[1])*x0-(b[0]-a[0])*y0+b[0]*a[1]-b[1]*a[0]); den=math.hypot(b[1]-a[1],b[0]-a[0]) or 1e-9
        d=num/den
        if d>dm: dm=d; k=i
    return rdp(p[:k+1],eps)[:-1]+rdp(p[k:],eps) if dm>eps else [a,b]
qr="""[out:json][timeout:90];( way["railway"="rail"]["highspeed"="yes"](%f,%f,%f,%f);
 way["railway"="rail"]["usage"="main"](%f,%f,%f,%f); ); out geom;"""%(BBOX+BBOX)
ways=[[(p["lon"],p["lat"]) for p in e["geometry"]] for e in overpass(qr)["elements"] if e.get("geometry") and len(e["geometry"])>1]
def key(p): return (round(p[0],5),round(p[1],5))
from collections import defaultdict
used=[False]*len(ways); ep=defaultdict(list)
for i,w in enumerate(ways): ep[key(w[0])].append(i); ep[key(w[-1])].append(i)
def build(s):
    path=list(ways[s]); used[s]=True; ch=True
    while ch:
        ch=False; t=key(path[-1])
        for j in ep[t]:
            if not used[j]:
                w=ways[j]; path+=(w[1:] if key(w[0])==t else w[::-1][1:]); used[j]=True; ch=True; break
    ch=True
    while ch:
        ch=False; h=key(path[0])
        for j in ep[h]:
            if not used[j]:
                w=ways[j]; path=(w[:-1] if key(w[-1])==h else w[::-1][:-1])+path; used[j]=True; ch=True; break
    return path
Ln=lambda p: sum(math.hypot((b[0]-a[0])*111.32*math.cos(math.radians((a[1]+b[1])/2)),(b[1]-a[1])*111.32) for a,b in zip(p,p[1:]))
hsr=rdp(sorted((build(i) for i in range(len(ways)) if not used[i]),key=Ln,reverse=True)[0],0.00035)
open("hsr.gmt","w").write("\n".join("%.6f %.6f"%(x,y) for x,y in hsr))
lat0=sum(p[1] for p in hsr)/len(hsr); mlon=111320*math.cos(math.radians(lat0)); mlat=111320; lon0=hsr[0][0]
xy=[((lo-lon0)*mlon,(la-lat0)*mlat) for lo,la in hsr]
def sn(a,b):
    dx,dy=b[0]-a[0],b[1]-a[1]; L=math.hypot(dx,dy) or 1e-9; return(-dy/L,dx/L)
n=len(xy); norm=[]
for i in range(n):
    if i==0: nv=sn(xy[0],xy[1])
    elif i==n-1: nv=sn(xy[n-2],xy[n-1])
    else:
        a=sn(xy[i-1],xy[i]); b=sn(xy[i],xy[i+1]); sx,sy=a[0]+b[0],a[1]+b[1]; L=math.hypot(sx,sy) or 1e-9; nv=(sx/L,sy/L)
    norm.append(nv)
def band(dd):
    lft=[(xy[i][0]+norm[i][0]*dd,xy[i][1]+norm[i][1]*dd) for i in range(n)]
    rgt=[(xy[i][0]-norm[i][0]*dd,xy[i][1]-norm[i][1]*dd) for i in range(n)]
    return [(lon0+x/mlon,lat0+y/mlat) for x,y in lft+rgt[::-1]+[lft[0]]]
for dd,nm in [(100,"buffer_100"),(250,"buffer_250"),(500,"buffer_500")]:
    open(nm+".gmt","w").write("\n".join("%.6f %.6f"%(x,y) for x,y in band(dd)))
try:
    roads=[rdp([(p["lon"],p["lat"]) for p in e["geometry"]],0.0009) for e in overpass('[out:json][timeout:90];way["highway"="motorway"](%f,%f,%f,%f);out geom;'%BBOX)["elements"] if e.get("geometry") and len(e["geometry"])>1]
    with open("roads.gmt","w") as f:
        for r in roads: f.write(">\n"+"\n".join("%.6f %.6f"%(x,y) for x,y in r)+"\n")
except Exception: open("roads.gmt","w").write("")
end=hsr[-1] if hsr[-1][0]>hsr[0][0] else hsr[0]
open("endlab.txt","w").write("%.4f %.4f  to Eskisehir / Ankara\n"%(end[0]+0.02,end[1]))
def merc(la): return math.log(math.tan(math.pi/4+math.radians(la)/2))
WA=15.0; WB=7.6; lonA=2.1; latA0,latA1=40.40,41.15
HA=WA*(merc(latA1)-merc(latA0))/(lonA*math.pi/180)
ci=min(range(n),key=lambda i:abs(hsr[i][0]-30.00)); cx,cy=hsr[ci]
lonspanB=0.19; targetdm=HA*(lonspanB*math.pi/180)/WB
lo,hi=0.0,0.5
for _ in range(60):
    mid=(lo+hi)/2
    if merc(cy+mid)-merc(cy-mid)<targetdm: lo=mid
    else: hi=mid
d=(lo+hi)/2
b=(round(cx-lonspanB/2,4),round(cx+lonspanB/2,4),round(cy-d,4),round(cy+d,4))
open("inset_box.gmt","w").write("\n".join("%.4f %.4f"%(a,c) for a,c in [(b[0],b[2]),(b[1],b[2]),(b[1],b[3]),(b[0],b[3]),(b[0],b[2])])+"\n")
open("insetR.txt","w").write("-R%f/%f/%f/%f\n"%b)
open("bscalepos.txt","w").write("%.4f/%.4f\n"%(b[0]+0.02,b[2]+0.012))
open("boxlab.txt","w").write("%.4f %.4f detail (b)\n"%(b[1],b[3]))
print("HA=%.3f cm  b-region"%HA,b)
PY

cat > rupture.txt <<'EOF'
29.40 40.715
29.65 40.715
29.83 40.720
30.05 40.722
30.26 40.720
30.50 40.720
30.62 40.725
30.80 40.735
30.95 40.755
31.10 40.775
EOF
cat > cities.txt <<'EOF'
28.979 41.008 Istanbul
29.431 40.803 Gebze
29.917 40.766 Izmit
30.403 40.780 Adapazari
EOF
cat > legB.txt <<'EOF'
N 4
S 0.3c - 0.8c - 1.7p,black 1.0c High-speed railway (YHT)
S 0.3c - 0.8c - 1.3p,goldenrod3 1.0c Motorway (OSM)
S 0.3c - 0.8c - 0.9p,gray15 1.0c NAFZ active fault
S 0.3c - 0.8c - 2.2p,red2 1.0c 1999 rupture (Mw 7.6)
G 0.16c
L 9p,Helvetica-Bold C Analysis buffer (to scale)
N 3
S 0.3c s 0.32c 230/85/13 0.2p,gray30 0.85c 100 m
S 0.3c s 0.32c 253/141/60 0.2p,gray30 0.85c 250 m
S 0.3c s 0.32c 253/208/162 0.2p,gray30 0.85c 500 m
EOF
cat > credit2.txt <<'EOF'
0.35 1.35 Relief: SRTM 3 arc-sec (jet, shaded); coastline: GSHHG; railway (Istanbul-Ankara YHT) and motorways: OpenStreetMap contributors (ODbL).
0.35 1.02 NAFZ faults: GEM Global Active Faults; 1999 rupture after Barka et al. (2002). Buffers (100/250/500 m) shown to scale in panel (b).
0.35 0.69 Buffers are geometric offsets of the OSM alignment (this study); not classified by hazard here. Source: authors.
EOF
printf '28.55 41.155 (a) Istanbul-Kocaeli HSR corridor, motorways and NAFZ setting\n' > titleA.txt
cat > faults_nafz.gmt <<'EOF'
>
29.58679 40.31422
29.55686 40.30000
>
29.81355 40.30000
29.88046 40.32964
29.96895 40.37431
>
30.41703 40.64534
30.53559 40.62257
30.65410 40.59981
30.77258 40.57705
30.89455 40.57154
31.01665 40.57484
31.13470 40.59531
>
29.32750 40.72590
29.44555 40.72737
29.56360 40.72884
29.68125 40.73646
29.79199 40.70546
29.90975 40.71155
30.02667 40.72394
30.14472 40.72527
30.26278 40.72537
30.37994 40.71440
30.49791 40.71115
30.61116 40.68600
30.72827 40.69699
30.83603 40.73351
30.92354 40.75529
>
29.05049 40.42447
29.28884 40.41470
29.52377 40.38294
29.76086 40.40155
29.99669 40.42953
30.22644 40.47573
30.45852 40.51905
30.61799 40.55624
30.68051 40.57070
>
28.05630 40.79680
28.35958 40.84621
28.66855 40.86934
28.85359 40.86610
28.96792 40.83304
29.12786 40.76085
29.22790 40.71570
EOF

gmt set FONT_ANNOT_PRIMARY 9p FONT_LABEL 10p MAP_FRAME_TYPE plain PS_CHAR_ENCODING ISOLatin1+ FORMAT_GEO_MAP ddd:mmF
gmt makecpt -Cjet -T-40/1700/25 -Z > relief.cpt
INSETR=$(cat insetR.txt); BSC=$(cat bscalepos.txt)

gmt begin fig05_corridor_buffers png E600
  R1=-R28.55/30.65/40.40/41.15
  gmt grdimage @earth_relief_03s $R1 -JM15c -Crelief.cpt -I+d
  gmt coast $R1 -JM15c -Slightsteelblue -Dh -W0.4p,gray40
  gmt plot faults_nafz.gmt $R1 -JM15c -W0.9p,gray15
  gmt plot roads.gmt  $R1 -JM15c -W1.2p,goldenrod3
  gmt plot rupture.txt $R1 -JM15c -W2.2p,red2
  gmt plot inset_box.gmt $R1 -JM15c -W2.4p,white
  gmt plot inset_box.gmt $R1 -JM15c -W1.2p,black,4_2
  gmt text boxlab.txt $R1 -JM15c -F+f8p,Helvetica-Bold+jBR -Gwhite@25 -Dj0.05c/0.06c
  gmt plot hsr.gmt $R1 -JM15c -W3.4p,white
  gmt plot hsr.gmt $R1 -JM15c -W1.7p,black
  gmt text endlab.txt $R1 -JM15c -F+f7.5p,black+jLM -Gwhite@25 -W0.2p,gray
  gmt plot cities.txt $R1 -JM15c -Sc0.15c -Gblack -W0.4p,white
  gmt text cities.txt $R1 -JM15c -F+f9p,Helvetica-Bold+jBL -Dj0.16c/0.10c -Gwhite@30 -W0.2p,gray -C10%
  gmt text titleA.txt $R1 -JM15c -F+f11p,Helvetica-Bold+jBL -N -Dj0/0.12c
  gmt basemap $R1 -JM15c -BWSne -Bxa0.5f0.25 -Bya0.25f0.125 -Lg29.05/40.455+w40k+f+u --FONT_ANNOT_PRIMARY=8p
  gmt basemap $R1 -JM15c -Tdg28.72/41.06+w0.8c+f2+l,,,N
  gmt grdimage @earth_relief_03s $INSETR -JM7.6c -Crelief.cpt -I+d -X17.0c -Y0
  gmt plot buffer_500.gmt $INSETR -JM7.6c -G253/208/162 -W0.25p,gray30
  gmt plot buffer_250.gmt $INSETR -JM7.6c -G253/141/60 -W0.25p,gray30
  gmt plot buffer_100.gmt $INSETR -JM7.6c -G230/85/13 -W0.25p,gray30
  gmt plot roads.gmt $INSETR -JM7.6c -W1.0p,goldenrod3
  gmt plot rupture.txt $INSETR -JM7.6c -W2.4p,red2
  gmt plot faults_nafz.gmt $INSETR -JM7.6c -W0.9p,gray15
  gmt plot hsr.gmt $INSETR -JM7.6c -W2.6p,white
  gmt plot hsr.gmt $INSETR -JM7.6c -W1.2p,black
  echo "29.917 40.766 Izmit" | gmt plot $INSETR -JM7.6c -Sc0.15c -Gblack -W0.4p,white
  echo "29.917 40.766 Izmit" | gmt text $INSETR -JM7.6c -F+f9p,Helvetica-Bold+jTL -Dj0.16c/0.10c -Gwhite@30 -W0.2p,gray -C10%
  echo "(b) Multi-buffer exposure zones (to scale)" | gmt text -F+cTC+f11p,Helvetica-Bold+jBC -N -D0/0.12c
  gmt basemap $INSETR -JM7.6c -BWSne -Bxa0.05f0.025 -Bya0.02f0.01 -Lg$BSC+w2k+f+u --FONT_ANNOT_PRIMARY=8p
  gmt colorbar -R0/24.6/0/5 -JX24.6c/5c -X-17.0c -Y-5.9c -Crelief.cpt \
       -Dx12.3/4.85+w13c/0.32c+h+e+jTC -Bxa500+l"Elevation (m)" --FONT_ANNOT_PRIMARY=10p --FONT_LABEL=12p
  gmt legend legB.txt -Dx12.3/3.62+w20c+jTC --FONT_ANNOT_PRIMARY=9p
  gmt text credit2.txt -F+f7.5p,gray30+jTL -N
gmt end

mkdir -p /sandbox/output && cp /tmp/fig05_corridor_buffers.png /sandbox/output/
echo DONE
