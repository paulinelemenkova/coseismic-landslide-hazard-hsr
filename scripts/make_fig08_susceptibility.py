import os, sys, math, json, zipfile, glob, shutil, numpy as np
import geopandas as gpd, pandas as pd, rasterio
from rasterio.features import rasterize
from rasterio.transform import from_origin
from scipy.spatial import cKDTree
from scipy.cluster.vq import kmeans2
from shapely.geometry import box, LineString, MultiLineString
from shapely.ops import unary_union, linemerge

UP="/mnt/user-data/uploads"; OUT="/mnt/user-data/outputs"
HERE=os.path.dirname(os.path.abspath(__file__))
W,S,E,N=28.0,40.35,30.9,41.35
WEIGHTS={'slope':0.312,'elev':0.094,'twi':0.112,'fault':0.241,'litho':0.178}
os.makedirs("/tmp/osm",exist_ok=True); os.makedirs("/tmp/geo",exist_ok=True)

if not glob.glob("/tmp/geo/**/Turkey 500k Lithology V1.TAB",recursive=True):
    zipfile.ZipFile(f"{UP}/Geological_map_data_of_Turkey.zip").extractall("/tmp/geo")
GEO=os.path.dirname(glob.glob("/tmp/geo/**/Turkey 500k Lithology V1.TAB",recursive=True)[0])
if not os.path.exists("/tmp/osm/corridor_osm.gpkg") and os.path.exists(f"{UP}/corridor_osm.zip"):
    zipfile.ZipFile(f"{UP}/corridor_osm.zip").extractall("/tmp/osm")
    _g=glob.glob("/tmp/osm/**/*.gpkg",recursive=True)
    if _g and _g[0]!="/tmp/osm/corridor_osm.gpkg": shutil.move(_g[0],"/tmp/osm/corridor_osm.gpkg")

NE="/tmp/ne_adm1.json"
if not os.path.exists(NE):
    import urllib.request
    urllib.request.urlretrieve("https://raw.githubusercontent.com/martynafford/natural-earth-geojson/"
        "master/10m/cultural/ne_10m_admin_1_states_provinces.json",NE)

FAC=next((q for q in [os.path.join(HERE,"fig06_factors_srtm30.tif"),
                      "/mnt/user-data/outputs/fig06_factors_srtm30.tif","/tmp/fac30.tif"]
          if os.path.exists(q)),None)
if FAC is None:
    sys.exit("Missing fig06_factors_srtm30.tif (SRTM 30 m factor grid). It ships with this "
             "script; regenerate with the GMT recipe in make_fig06_corridor_clh.py.")

d=rasterio.open(FAC); Fg=d.read().astype('float64'); trf=d.transform
r_slope,r_elev,r_twi=Fg[0],Fg[1],Fg[2]; ny,nx=r_slope.shape; rx=trf.a; ry=-trf.e
lons=trf.c+rx*(np.arange(nx)+0.5); lats=trf.f+trf.e*(np.arange(ny)+0.5)
LAT0=float(lats.mean()); kx=111320*math.cos(math.radians(LAT0)); ky=110570.0
land0=r_slope>0

def seg_read(fn):
    segs,s=[],[]
    for ln in open(fn):
        ln=ln.strip()
        if ln.startswith(">"):
            if len(s)>1: segs.append(s)
            s=[]
        elif ln and not ln.startswith("#"):
            p=ln.split()
            try: s.append((float(p[0]),float(p[1])))
            except Exception: pass
    if len(s)>1: segs.append(s)
    return segs
def densify(seg,st=0.004):
    o=[]
    for (x0,y0),(x1,y1) in zip(seg,seg[1:]):
        dd=math.hypot(x1-x0,y1-y0); n=max(1,int(dd/st))
        for t in range(n): o.append((x0+(x1-x0)*t/n,y0+(y1-y0)*t/n))
    o.append(seg[-1]); return o
faults=seg_read(f"{UP}/faults_nafz.gmt"); allpts=[p for s in faults for p in densify(s)]
LON,LAT=np.meshgrid(lons,lats); cell=np.column_stack([(LON*kx).ravel(),(LAT*ky).ravel()]); del LON,LAT
fdist=(cKDTree(np.array([[p[0]*kx,p[1]*ky] for p in allpts])).query(cell,workers=-1)[0]/1000).reshape(ny,nx); del cell
r_fault=(6-(np.digitize(fdist,[1,2,5,10])+1)).astype('float64')

def litho_rating(a,dd):
    a,dd=str(a),str(dd).lower()
    if a=="Water": return 0
    if "flysch" in dd: return 5
    if a in ("Unconsolidated and Semiconsolidated","Sedimentary","Volcaniclastic"): return 4
    if "Ultrabasic" in a or a=="Ophiolite": return 4
    if a=="Metamorphic": return 2 if ("marble" in dd or "crystalline limestone" in dd) else 3
    if "Volcanic" in a or "Hydrothermal" in a: return 3
    if "Intrusive" in a: return 2
    return 3
lit=gpd.clip(gpd.read_file(f"{GEO}/Turkey 500k Lithology V1.TAB"),box(W,S,E,N))
lit["rate"]=[litho_rating(a,dd) for a,dd in zip(lit["Lith_Association"],lit["Description"])]
r_litho=rasterize([(g,rt) for g,rt in zip(lit.geometry,lit["rate"]) if g is not None and rt>0],
                  out_shape=(ny,nx),transform=trf,fill=0,all_touched=True).astype('float64')
water=rasterize([(g,1) for g,rt in zip(lit.geometry,lit["rate"]) if rt==0],
                out_shape=(ny,nx),transform=trf,fill=0,all_touched=True)
r_litho[r_litho==0]=3.0; land=land0&(water==0)
tw=sum(WEIGHTS.values()); wn={k:v/tw for k,v in WEIGHTS.items()}
LSI=wn['slope']*r_slope+wn['elev']*r_elev+wn['twi']*r_twi+wn['fault']*r_fault+wn['litho']*r_litho
print("renorm weights:",{k:round(v,3) for k,v in wn.items()})

prov=gpd.read_file(NE)
prov=prov[(prov["iso_a2"]=="TR") & prov["name"].astype(str).str.contains("Istanbul|İstanbul|Kocaeli",case=False,na=False)]
provmask=rasterize([(g,1) for g in prov.geometry],out_shape=(ny,nx),transform=trf,fill=0,all_touched=True).astype(bool)
prov_area=float((prov.to_crs(6933).area/1e6).sum()); study=provmask&land
vals=LSI[study]; sub=np.random.default_rng(11).choice(vals,size=min(80000,vals.size),replace=False)
cen,_=kmeans2(sub.reshape(-1,1),5,minit="++",seed=11); cen=np.sort(cen.ravel()); bnd=(cen[:-1]+cen[1:])/2
cls=np.where(study,np.digitize(LSI,bnd)+1,np.nan)
dy=ry*110570.0; dxr=rx*111320.0*np.cos(np.radians(lats)); areacell=(dxr*dy/1e6)[:,None]*np.ones((1,nx))
names=["Very Low","Low","Moderate","High","Very High"]
rows=[float(np.nansum(areacell[cls==c])) for c in range(1,6)]; tk=sum(rows); pct=[100*r/tk for r in rows]
wnpanel={'slope':0.333,'fault':0.257,'litho':0.190,'twi':0.120,'elev':0.100}
print("=== SUSCEPTIBILITY CLASS AREAS (Istanbul+Kocaeli land, %.0f km2) ==="%tk)
for c in range(5): print("%-10s %7.0f km2  %4.1f%%"%(names[c],rows[c],pct[c]))
print("High+VeryHigh = %.1f%%"%(pct[3]+pct[4]))

HSR=None
try:
    if os.path.exists("/tmp/fig06_hsr_main.gpkg"):
        HSR=gpd.read_file("/tmp/fig06_hsr_main.gpkg")
    elif os.path.exists("/tmp/osm/corridor_osm.gpkg"):
        rail=gpd.read_file("/tmp/osm/corridor_osm.gpkg",layer="railways")
        hsf=rail[(rail['fclass']=='rail')&(rail['name'].astype(str).str.contains('hızlı',case=False,na=False))]
        hsf=gpd.clip(hsf,box(W,S,E,N)); mrg=linemerge(unary_union(hsf.to_crs(32635).geometry.values))
        allc=sorted([c for c in (list(mrg.geoms) if isinstance(mrg,MultiLineString) else [mrg]) if c.length>800],key=lambda c:c.bounds[0])
        cc=[]; cov=None
        for c in allc:
            x0,x1=c.bounds[0],c.bounds[2]
            if cov is not None and (max(0,min(x1,cov[1])-max(x0,cov[0]))/(x1-x0))>0.5: continue
            co=list(c.coords); co=co[::-1] if co[0][0]>co[-1][0] else co
            cc.append(LineString(co)); cov=(x0,x1) if cov is None else (min(cov[0],x0),max(cov[1],x1))
        HSR=gpd.GeoSeries(cc,crs=32635).to_crs(4326)
        HSR.to_file("/tmp/fig06_hsr_main.gpkg",driver="GPKG")
except Exception as ex: print("hsr context:",ex)

import warnings; warnings.filterwarnings("ignore")
import os, numpy as np, math, json, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.font_manager as fm
from matplotlib.colors import ListedColormap, BoundaryNorm
from matplotlib.patches import Patch, Rectangle
import matplotlib.patheffects as pe
import geopandas as gpd, rasterio

for _f in ["LiberationSans-Regular.ttf","LiberationSans-Bold.ttf","LiberationSans-Italic.ttf"]:
    try: fm.fontManager.addfont(f"/usr/share/fonts/truetype/liberation/{_f}")
    except Exception: pass
FAMILY="TeX Gyre Heros"
plt.rcParams.update({"font.family":"sans-serif",
    "font.sans-serif":[FAMILY,"Liberation Sans","Nimbus Sans","Arimo"],
    "mathtext.fontset":"stixsans","axes.unicode_minus":False,
    "pdf.fonttype":42,"ps.fonttype":42,"svg.fonttype":"none"})

W,E=float(lons[0]),float(lons[-1]); S,N=float(lats[-1]),float(lats[0]); ext=[W,E,S,N]
wn=wnpanel

COLS=["#1a9850","#91cf60","#fee08b","#fc8d59","#d73027"]
cmap=ListedColormap(COLS); norm=BoundaryNorm([0.5,1.5,2.5,3.5,4.5,5.5],cmap.N)

fig,ax=plt.subplots(figsize=(14.0,6.6),dpi=300); ax.set_facecolor("white")

_cand=["/tmp/hillshade_01s.tif","/mnt/user-data/outputs/fig08_hillshade_srtm90.tif",
       "/mnt/user-data/outputs/fig09_hillshade_srtm90.tif","/mnt/user-data/outputs/fig06_hillshade_srtm90.tif"]
HS=next((p for p in _cand if os.path.exists(p)),None)
if HS:
    _h=rasterio.open(HS); hs=_h.read(1).astype(float)
    hs_ext=[_h.bounds.left,_h.bounds.right,_h.bounds.bottom,_h.bounds.top]
    ax.imshow(np.where(hs>0,hs,np.nan),extent=hs_ext,origin="upper",cmap="gray",vmin=15,vmax=245,zorder=1)
ax.imshow(np.where(~land,1.0,np.nan),extent=ext,origin="upper",cmap=ListedColormap(["#cfe0ee"]),zorder=2)
ax.imshow(np.where(study,cls,np.nan),extent=ext,origin="upper",cmap=cmap,norm=norm,
          alpha=0.82,zorder=3,interpolation="nearest")

def seg_read(fn):
    segs=[];s=[]
    for ln in open(fn):
        ln=ln.strip()
        if ln.startswith(">"):
            if len(s)>1:segs.append(s)
            s=[]
        elif ln and not ln.startswith("#"):
            p=ln.split()
            try:s.append((float(p[0]),float(p[1])))
            except:pass
    if len(s)>1:segs.append(s)
    return segs

prov=gpd.read_file("/tmp/ne_adm1.json")
prov=prov[(prov["iso_a2"]=="TR") & prov["name"].astype(str).str.contains("Istanbul|İstanbul|Kocaeli",case=False,na=False)]
prov.boundary.plot(ax=ax,color="black",linewidth=1.3,zorder=6)

for s in seg_read(f"{UP}/faults_nafz.gmt"):
    xs,ys=zip(*s); ax.plot(xs,ys,color="#c1121f",lw=1.9,zorder=7,solid_capstyle="round")

try:
    for g in (HSR.geometry if HSR is not None else []):
        xs,ys=g.xy; ax.plot(xs,ys,color="white",lw=3.0,zorder=8,solid_capstyle="round")
        ax.plot(xs,ys,color="black",lw=1.4,zorder=9,solid_capstyle="round")
except Exception as ex: print("hsr",ex)

for nm,x,y in [("İstanbul",28.979,41.008),("Gebze",29.431,40.803),("İzmit",29.917,40.766),
               ("Adapazarı",30.403,40.780),("Sapanca",30.271,40.686)]:
    ax.plot(x,y,'o',ms=5,mfc="black",mec="white",mew=0.8,zorder=10)
    ax.text(x+0.02,y+0.015,nm,fontsize=9,fontweight="bold",zorder=11,
            path_effects=[pe.withStroke(linewidth=2.4,foreground="white")])
for nm,x,y in [("İSTANBUL",28.62,41.16),("KOCAELİ",30.05,40.98)]:
    ax.text(x,y,nm,fontsize=12,fontweight="bold",color="#333333",ha="center",zorder=6,
            path_effects=[pe.withStroke(linewidth=2.5,foreground="white")])
for nm,x,y in [("Black Sea",29.55,41.30),("Sea of Marmara",28.55,40.55),("İzmit Bay",29.72,40.745)]:
    ax.text(x,y,nm,fontsize=10,style="italic",color="#2c5d8f",ha="center",zorder=6)
ax.annotate("Fault-proximal flysch:\nHigh / Very High",xy=(29.95,40.85),xytext=(29.55,40.47),
            fontsize=9,fontweight="bold",ha="center",zorder=12,
            path_effects=[pe.withStroke(linewidth=2.6,foreground="white")],
            arrowprops=dict(arrowstyle="-|>",color="black",lw=1.1,shrinkA=0,shrinkB=2))
ax.text(30.45,40.66,"North Anatolian Fault Zone",fontsize=9.5,fontweight="bold",color="#c1121f",
        rotation=-6,zorder=11,path_effects=[pe.withStroke(linewidth=2.4,foreground="white")])

ax.set_xlim(W,E); ax.set_ylim(S,N); ax.set_aspect(1/math.cos(math.radians(LAT0)))
ax.set_xticks(np.arange(28.0,31.0,0.5)); ax.set_yticks(np.arange(40.4,41.4,0.2))
ax.xaxis.set_minor_locator(plt.MultipleLocator(0.1)); ax.yaxis.set_minor_locator(plt.MultipleLocator(0.05))
ax.set_xticklabels([f"{v:.1f}°E" for v in np.arange(28.0,31.0,0.5)],fontsize=8)
ax.set_yticklabels([f"{v:.1f}°N" for v in np.arange(40.4,41.4,0.2)],fontsize=8)
ax.tick_params(which="major",direction="out",length=4,top=True,right=True)
ax.tick_params(which="minor",direction="out",length=2.3,top=True,right=True)
for sp in ax.spines.values(): sp.set_linewidth(1.0)
ax.grid(which="major",color="white",linewidth=.4,alpha=.30,zorder=3.4)

km=40.0; dlon=km/(111.320*math.cos(math.radians(LAT0))); x0,y0=28.06,40.40
ax.add_patch(Rectangle((x0,y0),dlon,0.012,facecolor="black",edgecolor="black",zorder=13))
ax.add_patch(Rectangle((x0,y0),dlon/2,0.012,facecolor="white",edgecolor="black",zorder=13))
ax.text(x0,y0+0.03,"0",fontsize=7.5,ha="center",zorder=13)
ax.text(x0+dlon,y0+0.03,f"{km:.0f} km",fontsize=7.5,ha="center",zorder=13)
ax.annotate("",xy=(30.83,41.30),xytext=(30.83,41.19),arrowprops=dict(arrowstyle="-|>",color="black",lw=2.0),zorder=13)
ax.text(30.83,41.315,"N",fontsize=11,fontweight="bold",ha="center",zorder=13)

handles=[Patch(facecolor=COLS[4-i],edgecolor="gray",label=f"{names[4-i]}  ({pct[4-i]:.1f}%)") for i in range(5)]
leg=ax.legend(handles=handles,title="Landslide susceptibility",loc="lower left",
              bbox_to_anchor=(0.006,0.055),fontsize=8.5,title_fontsize=9.5,framealpha=0.92,
              borderpad=0.7,handlelength=1.3)
leg.get_frame().set_edgecolor("gray"); leg.set_zorder(14)

ax.set_title("Landslide Susceptibility of the İstanbul–Kocaeli High-Speed-Railway Corridor",
             fontsize=14,fontweight="bold",pad=16)
ax.text(0.5,1.006,"Expert-weighted AHP model — slope, elevation, wetness, fault proximity and MTA lithology — over the İstanbul + Kocaeli study area",
        transform=ax.transAxes,ha="center",fontsize=9,style="italic",color="#444444")

cred=("Susceptibility = AHP weighted linear combination of five reclassified factors: slope, elevation and topographic wetness index (SRTM 1 arc-sec, ~30 m; "
 "90 m analysis grid), Euclidean distance to the NAFZ (traces after Emre et al. 2013) and lithology (MTA 1:500 000 geology). AHP weights renormalised over "
 f"the five factors shown (slope {wn['slope']:.2f}, fault {wn['fault']:.2f}, lithology {wn['litho']:.2f}, TWI {wn['twi']:.2f}, elevation {wn['elev']:.2f}). "
 f"Five classes by natural breaks (k-means) over the study area ({tk:,.0f} km² of land within the İstanbul + Kocaeli provinces); sea and lakes masked. This is the same "
 "susceptibility surface (LSI) that the co-seismic hazard map (Fig. 6) scales by ground motion. Susceptibility from expert AHP weighting; classes are relative. "
 "Relief: SRTM 1 arc-sec hillshade. Railway: OpenStreetMap. Software: Python — Matplotlib, GeoPandas, rasterio, SciPy. Source: authors (P. Lemenkova & A. C. Zülfikar, İstanbul Technical University).")
fig.text(0.5,0.012,cred,ha="center",va="bottom",fontsize=6.1,color="#333333",wrap=True)

plt.subplots_adjust(left=0.045,right=0.99,top=0.90,bottom=0.11)
fig.savefig(f"{OUT}/fig08_susceptibility_map.png",dpi=300,facecolor="white")
fig.savefig(f"{OUT}/fig08_susceptibility_map.pdf",facecolor="white")
print("saved fig08 PNG + PDF (30 m)")
print("High+VH %.1f%%  study land %.0f km2"%(pct[3]+pct[4],tk))
