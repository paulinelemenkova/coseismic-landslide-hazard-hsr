import os, sys, zipfile, glob, shutil
UP="/mnt/user-data/uploads"
HERE=os.path.dirname(os.path.abspath(__file__))
os.makedirs("/tmp/osm",exist_ok=True); os.makedirs("/tmp/geo",exist_ok=True)
if not os.path.exists("/tmp/osm/corridor_osm.gpkg"):
    zipfile.ZipFile(f"{UP}/corridor_osm.zip").extractall("/tmp/osm")
    _g=glob.glob("/tmp/osm/**/*.gpkg",recursive=True)
    if _g and _g[0]!="/tmp/osm/corridor_osm.gpkg": shutil.move(_g[0],"/tmp/osm/corridor_osm.gpkg")
if not glob.glob("/tmp/geo/**/Turkey 500k Lithology V1.TAB",recursive=True):
    zipfile.ZipFile(f"{UP}/Geological_map_data_of_Turkey.zip").extractall("/tmp/geo")
os.environ["GEO_DIR"]=os.path.dirname(glob.glob("/tmp/geo/**/Turkey 500k Lithology V1.TAB",recursive=True)[0])
if not os.path.exists("/tmp/ne_adm1.json"):
    import urllib.request
    urllib.request.urlretrieve("https://raw.githubusercontent.com/martynafford/natural-earth-geojson/"
        "master/10m/cultural/ne_10m_admin_1_states_provinces.json","/tmp/ne_adm1.json")
_fac=next((q for q in [os.path.join(HERE,"fig06_factors_srtm30.tif"),
           "/mnt/user-data/outputs/fig06_factors_srtm30.tif","/tmp/fac30.tif"] if os.path.exists(q)),None)
if _fac is None:
    sys.exit("Missing fig06_factors_srtm30.tif (SRTM 30 m factor grid). It ships with this "
             "script; regenerate with the GMT recipe in make_fig06_corridor_clh.py.")
os.environ["FIG06_FACTOR_TIF"]=_fac

import warnings; warnings.filterwarnings("ignore")
import os,glob,math,json,numpy as np,geopandas as gpd,rasterio,pandas as pd
from rasterio.features import rasterize
from scipy.spatial import cKDTree
from scipy.cluster.vq import kmeans2
from scipy.interpolate import RegularGridInterpolator
from shapely.geometry import box,LineString,MultiLineString
from shapely.ops import unary_union, linemerge
UP="/mnt/user-data/uploads"
GEO=os.path.dirname(glob.glob("/tmp/geo/**/Turkey 500k Lithology V1.TAB",recursive=True)[0])
WEIGHTS={'slope':0.312,'elev':0.094,'twi':0.112,'fault':0.241,'litho':0.178}
MW=7.6; PGA_REF=0.4; W,S,E,N=28.0,40.35,30.9,41.35
import glob as _g
FAC=os.environ["FIG06_FACTOR_TIF"]
d=rasterio.open(FAC); F=d.read().astype('float64'); trf=d.transform
r_slope,r_elev,r_twi=F[0],F[1],F[2]; ny,nx=r_slope.shape; rx=trf.a; ry=-trf.e
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
            except: pass
    if len(s)>1: segs.append(s)
    return segs
def densify(seg,st=0.004):
    o=[]
    for (x0,y0),(x1,y1) in zip(seg,seg[1:]):
        dd=math.hypot(x1-x0,y1-y0); n=max(1,int(dd/st))
        for t in range(n): o.append((x0+(x1-x0)*t/n,y0+(y1-y0)*t/n))
    o.append(seg[-1]); return o
faults=seg_read(f"{UP}/faults_nafz.gmt")
allpts=[p for s in faults for p in densify(s)]; rup=densify(faults[0])
LON,LAT=np.meshgrid(lons,lats); cell=np.column_stack([(LON*kx).ravel(),(LAT*ky).ravel()]); del LON,LAT
fdist=(cKDTree(np.array([[p[0]*kx,p[1]*ky] for p in allpts])).query(cell,workers=-1)[0]/1000).reshape(ny,nx)
Rjb=(cKDTree(np.array([[p[0]*kx,p[1]*ky] for p in rup])).query(cell,workers=-1)[0]/1000).reshape(ny,nx); del cell
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
r_litho=rasterize([(g,rt) for g,rt in zip(lit.geometry,lit["rate"]) if g is not None and rt>0],out_shape=(ny,nx),transform=trf,fill=0,all_touched=True).astype('float64')
water=rasterize([(g,1) for g,rt in zip(lit.geometry,lit["rate"]) if rt==0],out_shape=(ny,nx),transform=trf,fill=0,all_touched=True)
r_litho[r_litho==0]=3.0
land=land0&(water==0)
tw=sum(WEIGHTS.values()); wn={k:v/tw for k,v in WEIGHTS.items()}
LSI=wn['slope']*r_slope+wn['elev']*r_elev+wn['twi']*r_twi+wn['fault']*r_fault+wn['litho']*r_litho
b=dict(b1=1.43525,b2=0.74866,b3=-0.06520,b4=-2.72950,b5=0.25139,b6=7.74959)
logP=b['b1']+b['b2']*MW+b['b3']*MW**2+(b['b4']+b['b5']*MW)*np.log10(np.sqrt(Rjb**2+b['b6']**2))
PGA=(10**logP)/981.0; ST=np.minimum(PGA/PGA_REF,1.0)
CLH=LSI*ST; CLH[~land]=np.nan
vals=CLH[land]; sub=np.random.default_rng(7).choice(vals,size=min(120000,vals.size),replace=False)
cen,_=kmeans2(sub.reshape(-1,1),5,minit='++',seed=7); cen=np.sort(cen.ravel()); bnd=(cen[:-1]+cen[1:])/2
CLHcls=np.where(land,np.digitize(CLH,bnd)+1,np.nan)
print("class cell counts:",[int(np.nansum(CLHcls==c)) for c in range(1,6)])
print("PGA %.3f-%.3f g ; corridor slope99 note"%(PGA[land].min(),PGA[land].max()))
np.savez_compressed("/tmp/fig06_arrays.npz",CLH=CLH,CLHcls=CLHcls,LSI=LSI,land=land,lons=lons,lats=lats,bnd=bnd,PGA=PGA)

GPKG="/tmp/osm/corridor_osm.gpkg"
if not os.path.exists("/tmp/fig06_hsr_main.gpkg"):
    rail=gpd.read_file(GPKG,layer="railways")
    hsf=rail[(rail['fclass']=='rail')&(rail['name'].astype(str).str.contains('hızlı',case=False,na=False))].copy()
    hsf=gpd.clip(hsf,box(W,S,E,N)); mrg=linemerge(unary_union(hsf.to_crs(32635).geometry.values))
    allc=sorted([c for c in (list(mrg.geoms) if isinstance(mrg,MultiLineString) else [mrg]) if c.length>800],key=lambda c:c.bounds[0])
    cc=[]; cov=None
    for c in allc:
        x0,x1=c.bounds[0],c.bounds[2]
        if cov is not None and (max(0,min(x1,cov[1])-max(x0,cov[0]))/(x1-x0))>0.5: continue
        co=list(c.coords); co=co[::-1] if co[0][0]>co[-1][0] else co
        cc.append(LineString(co)); cov=(x0,x1) if cov is None else (min(cov[0],x0),max(cov[1],x1))
    gpd.GeoSeries(cc,crs=32635).to_crs(4326).to_file("/tmp/fig06_hsr_main.gpkg",driver="GPKG")
    hsf.to_crs(4326).to_file("/tmp/fig06_hsr_all.gpkg",driver="GPKG")

comps=gpd.read_file("/tmp/fig06_hsr_main.gpkg").to_crs(32635)
latg=lats[::-1]; CLHc=np.where(np.isfinite(CLH),CLH,0.0)[::-1,:]
interp=RegularGridInterpolator((latg,lons),CLHc,bounds_error=False,fill_value=0.0)
step=50.0; names=["Very Low","Low","Moderate","High","Very High"]
def ring(r,lat0):
    ang=np.linspace(0,2*math.pi,10,endpoint=False); off=[(0.0,0.0)]
    for f in (0.5,1.0):
        for a in ang: off.append((f*r*math.cos(a)/(111320*math.cos(math.radians(lat0))),f*r*math.sin(a)/111000.0))
    return off
expo={100:[0.0]*5,250:[0.0]*5,500:[0.0]*5}; comp_s=[]
for c in comps.geometry:
    n=int(c.length/step); pll=gpd.GeoSeries([c.interpolate(i*step) for i in range(n+1)],crs=32635).to_crs(4326)
    plon=np.array([p.x for p in pll]); plat=np.array([p.y for p in pll]); cls250=None; val250=None
    for r in (100,250,500):
        off=ring(r,LAT0); vmax=np.zeros(len(plon))
        for dlon,dlat in off: vmax=np.maximum(vmax,interp(np.column_stack([plat+dlat,plon+dlon])))
        cls=np.digitize(vmax,bnd)+1
        for cc in range(1,6): expo[r][cc-1]+=float(np.sum(cls==cc))*step/1000.0
        if r==250: cls250=cls; val250=vmax
    comp_s.append((plon,plat,cls250,val250))
corr_km=sum(c.length for c in comps.geometry)/1000
hi=np.argwhere(np.nan_to_num(CLHcls)>=4); hitree=cKDTree(np.column_stack([lons[hi[:,1]]*kx,lats[hi[:,0]]*ky]))
GPKG="/tmp/osm/corridor_osm.gpkg"
places=gpd.read_file(GPKG,layer="places"); places=places[places['fclass'].isin(['city','town'])]
trans=gpd.read_file(GPKG,layer="transport"); trans=trans[trans['fclass'].isin(['railway_station','railway_halt'])&trans['name'].notna()]
namepts=gpd.GeoDataFrame(pd.concat([places[['name','geometry']],trans[['name','geometry']]],ignore_index=True),crs=4326)
nptree=cKDTree(np.array([[g.x,g.y] for g in namepts.geometry])*[kx,ky])
approach=1e9; hot=[]
for plon,plat,cls,val in comp_s:
    dvh,_=hitree.query(np.column_stack([plon*kx,plat*ky])); approach=min(approach,float(dvh.min()))
    flag=cls>=3; i=0
    while i<len(flag):
        if flag[i]:
            j=i
            while j<len(flag) and flag[j]: j+=1
            if (j-i)*step/1000.0>=2.0:
                seg=slice(i,j); pk=int(np.argmax(val[seg]))+i
                _,ii=nptree.query([plon[pk]*kx,plat[pk]*ky])
                hot.append({"len_km":round((j-i)*step/1000.0,1),"near":str(namepts.iloc[ii]['name']),
                            "peak_clh":round(float(val[seg].max()),2),"min_dist_vh_m":int(dvh[seg].min()),
                            "mid_lonlat":[float(plon[pk]),float(plat[pk])]})
            i=j
        else: i+=1
hot=sorted(hot,key=lambda h:-h['len_km'])
out=dict(names=names,expo=expo,corr_km=round(corr_km,1),hot=hot,MW=MW,PGA_REF=PGA_REF,
         approach_km=round(approach/1000,1),pga_min=float(PGA[land].min()),pga_max=float(PGA[land].max()),
         cls_counts=[int(np.nansum(CLHcls==c)) for c in range(1,6)],res="~90 m grid; slope from 30 m SRTM")
json.dump(out,open("/tmp/fig06_stats.json","w"))
print("exposure 250m:",{names[i]:round(expo[250][i],1) for i in range(4,-1,-1)})
print("approach km",round(approach/1000,1),"corr",round(corr_km,1))
print("hot:",[(h['near'],h['len_km'],h['peak_clh']) for h in hot[:6]])

import gc as _gc
for _n in [k for k in list(globals().keys()) if not k.startswith('_')]:
    try: del globals()[_n]
    except Exception: pass
_gc.collect()

import warnings; warnings.filterwarnings("ignore")
import os, numpy as np, math, json, geopandas as gpd, rasterio
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.font_manager as fm
import matplotlib.cm as cm
from matplotlib.ticker import MultipleLocator
from matplotlib.colors import ListedColormap, Normalize, to_rgba, LinearSegmentedColormap
from matplotlib.patches import Rectangle, FancyBboxPatch
from matplotlib.lines import Line2D
from matplotlib.collections import LineCollection
from matplotlib.patheffects import withStroke
from shapely.geometry import box, shape
from shapely.ops import unary_union

for _f in ["LiberationSans-Regular.ttf","LiberationSans-Bold.ttf","LiberationSans-Italic.ttf"]:
    try: fm.fontManager.addfont(f"/usr/share/fonts/truetype/liberation/{_f}")
    except Exception: pass
FAMILY="TeX Gyre Heros"
plt.rcParams.update({"font.family":"sans-serif","font.sans-serif":[FAMILY,"Liberation Sans","Nimbus Sans"],
    "mathtext.fontset":"stixsans","axes.unicode_minus":False,"pdf.fonttype":42,"ps.fonttype":42,"svg.fonttype":"none"})
W1=[withStroke(linewidth=2,foreground="white")]; W2=[withStroke(linewidth=2.6,foreground="white")]
NAVY="#12305e"; INK="#1a1a1a"; GRID="#c9c9c9"

UP="/mnt/user-data/uploads"; GPKG="/tmp/osm/corridor_osm.gpkg"
d=np.load("/tmp/fig06_arrays.npz"); CLH=d["CLH"]; land=d["land"]; lons=d["lons"]; lats=d["lats"]
st=json.load(open("/tmp/fig06_stats.json")); expo=st["expo"]; hot=st["hot"]; approach=st["approach_km"]
ny,nx=CLH.shape; rx=lons[1]-lons[0]; ry=lats[1]-lats[0]
ext=[lons.min(),lons.max(),lats.min(),lats.max()]
LAT0=float(lats.mean()); kx=111320*math.cos(math.radians(LAT0))
MAP=[28.95,30.55,40.40,41.08]

B0,B1,B2=0.357,0.572,0.914
VHB=float(np.percentile(CLH[land],97))
BRK=[B0,B1,B2,VHB]

C_VL="#a6cee3"; C_LOW="#b2df8a"; C_HIGH="#c71585"; C_VH="#7d26cd"
names=["Very Low","Low","Moderate","High","Very High"]

_ac=cm.get_cmap("autumn_r"); modn=Normalize(B1,B2); NMOD=5
def _mstep(t): t=float(np.clip(t,0,1)); return min(NMOD-1,int(t*NMOD))/(NMOD-1)
def mod_rgba(v): return _ac(_mstep(modn(v)))
MOD_SHADES=[_ac(i/(NMOD-1)) for i in range(NMOD)]

_lenc=LinearSegmentedColormap.from_list("modlen",["#ffffb2","#fed976","#feb24c","#fd8d3c","#fc4e2a","#b10026"])
_lens=[h["len_km"] for h in hot]; _lmin,_lmax=min(_lens),max(_lens)
def len_rgba(L): return _lenc(0.0 if _lmax<=_lmin else float(np.clip((L-_lmin)/(_lmax-_lmin),0,1)))

_cand=["/mnt/user-data/outputs/fig09_hillshade_srtm90.tif",
       os.path.join(os.path.dirname(os.path.abspath(__file__)) if "__file__" in globals() else "/tmp","fig09_hillshade_srtm90.tif")]
HS_TIF=next((p for p in _cand if os.path.exists(p)),None)
_h=rasterio.open(HS_TIF); hsr=_h.read(1).astype(float)
hs_ext=[_h.bounds.left,_h.bounds.right,_h.bounds.bottom,_h.bounds.top]; HSVMIN,HSVMAX=15,245

prov=json.load(open("/tmp/ne_adm1.json"))["features"]
def prov_geom(nm):
    g=[shape(f["geometry"]) for f in prov if f["properties"].get("iso_a2")=="TR" and f["properties"].get("name")==nm]
    return unary_union(g) if g else None
study=unary_union([prov_geom("Istanbul"),prov_geom("Kocaeli")])
turkey=unary_union([shape(f["geometry"]) for f in prov if f["properties"].get("iso_a2")=="TR"])

def clh_rgba():
    disp=np.full((ny,nx),np.nan); disp[land]=np.digitize(CLH[land],BRK)+1
    rgba=np.zeros((ny,nx,4))
    rgba[disp==1]=to_rgba(C_VL); rgba[disp==2]=to_rgba(C_LOW)
    mm=disp==3
    if mm.any():
        _k=np.minimum(NMOD-1,(np.clip(modn(CLH[mm]),0,1)*NMOD).astype(int)); rgba[mm]=_ac(_k/(NMOD-1))
    rgba[disp==4]=to_rgba(C_HIGH); rgba[disp==5]=to_rgba(C_VH)
    rgba[...,3]=np.where(np.isfinite(disp),0.78,0.0)
    return rgba
RGBA=clh_rgba()

comps=gpd.read_file("/tmp/fig06_hsr_main.gpkg").to_crs(4326)
cu_utm=unary_union(gpd.read_file("/tmp/fig06_hsr_main.gpkg").to_crs(32635).geometry.values)
def val_at(lon,lat):
    j=int(round((lon-lons[0])/rx)); i=int(round((lat-lats[0])/ry))
    if 0<=i<ny and 0<=j<nx:
        i0,i1=max(0,i-1),min(ny,i+2); j0,j1=max(0,j-1),min(nx,j+2)
        sub=CLH[i0:i1,j0:j1]; sub=sub[np.isfinite(sub)]
        return float(sub.max()) if sub.size else np.nan
    return np.nan
def densify(geom,step=0.0016):
    out=[]
    for ln in (geom.geoms if geom.geom_type=="MultiLineString" else [geom]):
        xy=np.array(ln.coords)
        for (x0,y0),(x1,y1) in zip(xy[:-1],xy[1:]):
            n=max(1,int(math.hypot(x1-x0,y1-y0)/step))
            for t in range(n):
                a,b=t/n,(t+1)/n
                out.append(((x0+(x1-x0)*a,y0+(y1-y0)*a),(x0+(x1-x0)*b,y0+(y1-y0)*b)))
    return out
allseg=[]
for g in comps.geometry: allseg+=densify(g)
def seg_color(v):
    if not np.isfinite(v): return "#999999"
    if v<B0: return C_VL
    if v<B1: return C_LOW
    if v<B2: return mod_rgba(v)
    if v<VHB: return C_HIGH
    return C_VH
seg_vals=[max(val_at(*p0),val_at(*p1)) for p0,p1 in allseg]
seg_cols=[seg_color(v) for v in seg_vals]

majors={"Gebze":(29.431,40.803),"Derince":(29.828,40.760),"İzmit":(29.945,40.766),
        "Köseköy":(30.050,40.766),"Arifiye":(30.373,40.717)}
pcol="#0b7d7d"

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
def buffer_rings(ax,dist,cbuf,lw):
    gb=gpd.GeoSeries([cu_utm.buffer(dist)],crs=32635).to_crs(4326).iloc[0]
    for gg in (gb.geoms if gb.geom_type=="MultiPolygon" else [gb]):
        xs,ys=gg.exterior.xy; ax.plot(xs,ys,color=cbuf,lw=lw,ls=(0,(5,3)),zorder=8,alpha=.9)

FIGW,FIGH=16.4,12.4
fig=plt.figure(figsize=(FIGW,FIGH),dpi=110)
fig.add_artist(Rectangle((0,0.948),1,0.052,transform=fig.transFigure,facecolor=NAVY,edgecolor="none",zorder=0))
fig.text(0.5,0.983,"Spatial Distribution of HSR Sections Exposed to Co-Seismic Landslide Hazard Classes",
         ha="center",va="center",fontsize=16.5,fontweight="bold",color="white")
fig.text(0.5,0.960,"İstanbul–Kocaeli high-speed railway corridor (Gebze–Arifiye) under the 1999 İzmit (Mw 7.6) NAFZ scenario",
         ha="center",va="center",fontsize=11.5,color="#cdd8ea")

axm=fig.add_axes([0.035,0.415,0.66,0.505])
axm.imshow(hsr,extent=hs_ext,cmap="gray",vmin=HSVMIN,vmax=HSVMAX,zorder=1)
axm.imshow(RGBA,extent=ext,zorder=2,interpolation="nearest")
axm.imshow(np.where(~land,1.0,np.nan),extent=ext,cmap=ListedColormap(["#cfe6f4"]),zorder=3,interpolation="nearest")
MB=box(MAP[0],MAP[2],MAP[1],MAP[3])
def clip(g): return gpd.clip(g,MB) if len(g) else g
clip(gpd.read_file(GPKG,layer="water")).plot(ax=axm,facecolor="#bcd9ec",edgecolor="#7fb0cf",linewidth=.2,zorder=4)
for nm in ["Istanbul","Kocaeli","Sakarya","Yalova","Bursa","Bilecik"]:
    g=prov_geom(nm)
    if g is not None: gpd.GeoSeries([g],crs=4326).boundary.plot(ax=axm,color="#3a3a3a",linewidth=.9,linestyle=(0,(6,3)),alpha=.55,zorder=6)
for sgm in seg_read(f"{UP}/faults_nafz.gmt"):
    a=np.array(sgm); axm.plot(a[:,0],a[:,1],color="#c1121f",lw=2.0,zorder=7,solid_capstyle="round")
for dist,cbuf in [(500,"#2e8b57"),(250,"#2e78c4"),(100,"#8e44ad")]: buffer_rings(axm,dist,cbuf,1.0)
comps.plot(ax=axm,color="white",linewidth=4.6,zorder=10,capstyle="round")
axm.add_collection(LineCollection([[p0,p1] for p0,p1 in allseg],colors=seg_cols,linewidths=2.6,zorder=11,capstyle="round"))
for nm,(lo,la) in majors.items():
    axm.plot(lo,la,marker="o",ms=7.5,mfc="white",mec=NAVY,mew=2,zorder=14)
    axm.text(lo,la+0.028,nm,fontsize=9,ha="center",va="bottom",color=INK,zorder=15,path_effects=W1)
offs=[(0.0,0.11),(-0.03,0.18),(0.05,0.11),(0.0,0.18)]
for k,h in enumerate(hot[:4]):
    lo,la=h["mid_lonlat"]; dx,dy=offs[k]
    axm.plot(lo,la,marker="o",ms=6,mfc=pcol,mec="white",mew=1.2,zorder=16)
    axm.annotate(f"P{k+1}",xy=(lo,la),xytext=(lo+dx,la+dy),fontsize=9.2,fontweight="bold",color="#075555",
                 ha="center",va="center",zorder=17,path_effects=W1,arrowprops=dict(arrowstyle="-",color=pcol,lw=1.1))
axm.text(29.05,41.02,"İstanbul",fontsize=14,fontweight="bold",color="#33322f",ha="center",zorder=15,path_effects=W2)
axm.text(30.28,40.99,"Kocaeli",fontsize=14,fontweight="bold",color="#33322f",ha="center",zorder=15,path_effects=W2)
axm.text(30.02,40.52,"Sakarya",fontsize=11.5,fontweight="bold",color="#4a4844",ha="center",zorder=15,path_effects=W2)
axm.text(29.72,40.615,"İzmit Bay",fontsize=9,style="italic",color="#1f5f86",ha="center",zorder=15)
axm.text(28.99,40.86,"Marmara Sea",fontsize=10.5,style="italic",color="#1f5f86",ha="left",zorder=15)
axm.text(29.10,40.79,"NAFZ",fontsize=8.5,fontweight="bold",color="#c1121f",rotation=6,zorder=15,path_effects=W1)
axm.set_xlim(MAP[0],MAP[1]); axm.set_ylim(MAP[2],MAP[3]); axm.set_aspect(1/math.cos(math.radians(LAT0)))
axm.xaxis.set_major_locator(MultipleLocator(0.5)); axm.xaxis.set_minor_locator(MultipleLocator(0.1))
axm.yaxis.set_major_locator(MultipleLocator(0.25)); axm.yaxis.set_minor_locator(MultipleLocator(0.05))
axm.xaxis.set_major_formatter(plt.FuncFormatter(lambda v,_:f"{v:.1f}°E"))
axm.yaxis.set_major_formatter(plt.FuncFormatter(lambda v,_:f"{v:.1f}°N"))
axm.tick_params(which="major",length=5,width=.9,labelsize=9.5,top=True,right=True)
axm.tick_params(which="minor",length=2.6,width=.6,top=True,right=True)
for s in axm.spines.values(): s.set_linewidth(1.1)
dlon=40*1000/kx; x0,y0=MAP[0]+0.06,MAP[2]+0.04
axm.add_patch(Rectangle((x0,y0),dlon,0.010,facecolor="#111",edgecolor="#111",zorder=20))
axm.add_patch(Rectangle((x0,y0),dlon/2,0.010,facecolor="#fff",edgecolor="#111",zorder=20))
for xx,tt in [(x0,"0"),(x0+dlon/2,"20"),(x0+dlon,"40 km")]: axm.text(xx,y0+0.02,tt,ha="center",va="bottom",fontsize=8,zorder=20,path_effects=W1)
nxp,nyp=MAP[0]+0.09,MAP[3]-0.09
axm.annotate("N",xy=(nxp,nyp),xytext=(nxp,nyp-0.08),ha="center",va="center",fontsize=13,fontweight="bold",zorder=20,
             arrowprops=dict(arrowstyle="-|>",color="#111",lw=2),path_effects=W2)

axi=fig.add_axes([0.553,0.805,0.135,0.108]); axi.set_facecolor("#eaf3fb")
gpd.GeoSeries([turkey],crs=4326).plot(ax=axi,facecolor="#dfe6d8",edgecolor="#8a8a8a",linewidth=.5)
sb=study.bounds
axi.add_patch(Rectangle((sb[0],sb[1]),sb[2]-sb[0],sb[3]-sb[1],facecolor="none",edgecolor="#c1121f",lw=1.6,zorder=5))
axi.set_xlim(25.5,45.2); axi.set_ylim(35.6,42.5); axi.set_aspect(1/math.cos(math.radians(39)))
axi.text(0.5,1.04,"Türkiye",transform=axi.transAxes,ha="center",va="bottom",fontsize=9,fontweight="bold")
axi.set_xticks([]); axi.set_yticks([])
for s in axi.spines.values(): s.set_edgecolor("#8a8a8a")

axR=fig.add_axes([0.70,0.415,0.281,0.505]); axR.axis("off")
axR.add_patch(FancyBboxPatch((0,0),1,1,boxstyle="round,pad=0.01",transform=axR.transAxes,facecolor="white",edgecolor="#c9c9c9",lw=1))
axR.text(0.06,0.96,"LEGEND",fontsize=13,fontweight="bold",transform=axR.transAxes)
fy=0.915
axR.add_line(Line2D([0.07,0.15],[fy,fy],color="white",lw=4.6,transform=axR.transAxes)); axR.add_line(Line2D([0.07,0.15],[fy,fy],color=NAVY,lw=2.4,transform=axR.transAxes))
axR.text(0.18,fy,"HSR corridor (Gebze–Arifiye)",fontsize=9.3,va="center",transform=axR.transAxes)
fy=0.876; axR.plot(0.11,fy,marker="o",ms=7.5,mfc="white",mec=NAVY,mew=2,transform=axR.transAxes); axR.text(0.18,fy,"Major stations",fontsize=9.3,va="center",transform=axR.transAxes)
fy=0.837; axR.plot(0.11,fy,marker="o",ms=6,mfc=pcol,mec="white",mew=1.2,transform=axR.transAxes); axR.text(0.18,fy,"Priority Moderate reach (P1–P4)",fontsize=9.3,va="center",transform=axR.transAxes)
fy=0.798; axR.add_line(Line2D([0.07,0.15],[fy,fy],color="#c1121f",lw=2.1,transform=axR.transAxes)); axR.text(0.18,fy,"Active faults (NAFZ)",fontsize=9.3,va="center",transform=axR.transAxes)
fy=0.759; axR.add_line(Line2D([0.07,0.15],[fy,fy],color="#3a3a3a",lw=1.0,ls=(0,(6,3)),transform=axR.transAxes)); axR.text(0.18,fy,"Province boundary",fontsize=9.3,va="center",transform=axR.transAxes)
axR.text(0.06,0.705,"Co-seismic landslide hazard (CLH)",fontsize=10.5,fontweight="bold",transform=axR.transAxes)
rows=[(C_VL,"Very Low",False),(C_LOW,"Low",False),(None,"Moderate",True),(C_HIGH,"High",False),(C_VH,"Very High",False)]
yy=0.655
for col,lab,grad in rows:
    if grad:
        for j,sh in enumerate(MOD_SHADES):
            axR.add_patch(Rectangle((0.07+j*0.024,yy-0.016),0.024,0.032,transform=axR.transAxes,facecolor=sh,edgecolor="#555",lw=.4))
        axR.text(0.205,yy,lab,fontsize=8.5,va="center",transform=axR.transAxes)
    else:
        axR.add_patch(Rectangle((0.07,yy-0.016),0.09,0.032,transform=axR.transAxes,facecolor=col,edgecolor="#555",lw=.6))
        axR.text(0.19,yy,lab,fontsize=9.0,va="center",transform=axR.transAxes)
    yy-=0.052
axR.text(0.06,0.355,"Buffer zones around corridor",fontsize=10.5,fontweight="bold",transform=axR.transAxes)
for i,(c,lab,sub) in enumerate([("#8e44ad","100 m","immediate track zone"),("#2e78c4","250 m","extended influence"),("#2e8b57","500 m","catchment-scale reach")]):
    yy=0.305-i*0.05
    axR.add_line(Line2D([0.07,0.16],[yy,yy],color=c,lw=1.6,ls=(0,(5,3)),transform=axR.transAxes))
    axR.text(0.19,yy,f"{lab}",fontsize=9.0,fontweight="bold",va="center",transform=axR.transAxes)
    axR.text(0.33,yy,sub,fontsize=8.0,style="italic",color="#555",va="center",transform=axR.transAxes)
axR.text(0.06,0.12,"Corridor lies in the Low–Moderate classes at every buffer\n(nearest High terrain ~%.1f km south). Very-High forms a near-fault\nbelt along the NAFZ (rendered at the 97th-pct CLH break)."%approach,
         fontsize=8.0,va="top",color="#444",transform=axR.transAxes,linespacing=1.45)

def gridlines(ax,ys,x0=0.03,x1=0.97,vx=None):
    for y in ys: ax.add_line(Line2D([x0,x1],[y,y],color=GRID,lw=0.8,transform=ax.transAxes,zorder=1))
    if vx:
        for x in vx: ax.add_line(Line2D([x,x],[min(ys),max(ys)],color=GRID,lw=0.8,transform=ax.transAxes,zorder=1))
def panel(x,w,h,title):
    ax=fig.add_axes([x,0.055,w,h]); ax.axis("off")
    ax.add_patch(FancyBboxPatch((0,0),1,1,boxstyle="round,pad=0.006",transform=ax.transAxes,facecolor="white",edgecolor="#b9b9b9",lw=1.1))
    ax.text(0.03,0.945,title,fontsize=10.3,fontweight="bold",transform=ax.transAxes)
    return ax

e250=[round(x,1) for x in expo["250"]]; Etot=round(sum(e250),1); epct=[100*x/Etot for x in e250]; ecum=np.cumsum(epct)
PH=0.315
axe=panel(0.035,0.30,PH,"HSR CORRIDOR EXPOSURE  (250 m buffer)")
hy=0.85
axe.text(0.10,hy,"Class",fontsize=8.3,fontweight="bold",transform=axe.transAxes)
axe.text(0.56,hy,"km",fontsize=8.3,fontweight="bold",ha="right",transform=axe.transAxes)
axe.text(0.76,hy,"%",fontsize=8.3,fontweight="bold",ha="right",transform=axe.transAxes)
axe.text(0.96,hy,"Cum.%",fontsize=8.3,fontweight="bold",ha="right",transform=axe.transAxes)
rowys=[0.74-i*0.108 for i in range(5)]
for i,yy in enumerate(rowys):
    if i==2:
        for j,sh in enumerate(MOD_SHADES):
            axe.add_patch(Rectangle((0.03+j*0.012,yy-0.035),0.012,0.07,transform=axe.transAxes,facecolor=sh,edgecolor="#555",lw=.25))
    else:
        axe.add_patch(Rectangle((0.03,yy-0.035),0.05,0.07,transform=axe.transAxes,facecolor=[C_VL,C_LOW,None,C_HIGH,C_VH][i],edgecolor="#555",lw=.4))
    axe.text(0.10,yy,names[i],fontsize=8.3,va="center",transform=axe.transAxes)
    axe.text(0.56,yy,f"{e250[i]:.1f}",fontsize=8.5,ha="right",va="center",transform=axe.transAxes)
    axe.text(0.76,yy,f"{epct[i]:.1f}",fontsize=8.5,ha="right",va="center",transform=axe.transAxes)
    axe.text(0.96,yy,f"{ecum[i]:.1f}",fontsize=8.5,ha="right",va="center",transform=axe.transAxes)
ty=0.135
axe.text(0.10,ty,"Total",fontsize=8.5,fontweight="bold",va="center",transform=axe.transAxes)
axe.text(0.56,ty,f"{Etot:.1f}",fontsize=8.5,fontweight="bold",ha="right",va="center",transform=axe.transAxes)
axe.text(0.76,ty,"100.0",fontsize=8.5,fontweight="bold",ha="right",va="center",transform=axe.transAxes)
gridlines(axe,[0.90,0.795]+[y-0.054 for y in rowys]+[0.088],x0=0.03,x1=0.97,vx=[0.45,0.65,0.85])
axe.text(0.03,0.03,"85.6 km Gebze–Arifiye; corridor entirely Low and Moderate.",fontsize=6.9,style="italic",color="#555",va="bottom",transform=axe.transAxes)

axp=panel(0.345,0.35,PH,"PRIORITY MODERATE SEGMENTS  (250 m buffer)")
cols=[("ID",0.05,"left"),("Location (nearest)",0.13,"left"),("Class",0.585,"left"),("Len (km)",0.82,"right"),("Peak CLH",0.985,"right")]
for tt,xx,ha in cols: axp.text(xx,0.85,tt,fontsize=8.0,fontweight="bold",ha=ha,transform=axp.transAxes)
prowys=[0.74-k*0.099 for k in range(6)]
for k,yy in enumerate(prowys):
    h=hot[k]
    axp.text(0.05,yy,f"P{k+1}",fontsize=8.2,fontweight="bold",color="#075555",va="center",transform=axp.transAxes)
    axp.text(0.13,yy,h["near"],fontsize=8.1,va="center",transform=axp.transAxes)
    axp.add_patch(Rectangle((0.585,yy-0.028),0.045,0.056,transform=axp.transAxes,facecolor=len_rgba(h["len_km"]),edgecolor="#555",lw=.4))
    axp.text(0.645,yy,"Moderate",fontsize=7.9,va="center",transform=axp.transAxes)
    axp.text(0.82,yy,f"{h['len_km']:.1f}",fontsize=8.3,ha="right",va="center",transform=axp.transAxes)
    axp.text(0.985,yy,f"{h['peak_clh']:.2f}",fontsize=8.3,ha="right",va="center",transform=axp.transAxes)
gridlines(axp,[0.80]+[y-0.049 for y in prowys],x0=0.03,x1=0.985,vx=[0.115,0.57,0.73,0.87])

_cbx0,_cbx1,_cby,_ns=0.345,0.70,0.115,26
for _s in range(_ns):
    axp.add_patch(Rectangle((_cbx0+(_cbx1-_cbx0)*_s/_ns,_cby-0.017),(_cbx1-_cbx0)/_ns+0.002,0.034,transform=axp.transAxes,facecolor=_lenc(_s/(_ns-1)),edgecolor="none",zorder=2))
axp.add_patch(Rectangle((_cbx0,_cby-0.017),_cbx1-_cbx0,0.034,transform=axp.transAxes,facecolor="none",edgecolor="#555",lw=.5,zorder=3))
axp.text(0.03,_cby,"Box colour = reach length",fontsize=6.9,style="italic",color="#555",ha="left",va="center",transform=axp.transAxes)
axp.text(_cbx0-0.006,_cby,f"{_lmin:.1f}",fontsize=6.8,ha="right",va="center",transform=axp.transAxes)
axp.text(_cbx1+0.008,_cby,f"{_lmax:.1f} km",fontsize=6.8,ha="left",va="center",transform=axp.transAxes)
axp.text(0.03,0.035,"Priority = longest continuous Moderate reaches; peak CLH is the reach maximum.",fontsize=6.8,color="#555",style="italic",va="bottom",transform=axp.transAxes)

IH=0.345
axz=fig.add_axes([0.705,0.055,0.275,IH])
axz.add_patch(FancyBboxPatch((0,0),1,1,boxstyle="round,pad=0.006",transform=axz.transAxes,facecolor="white",edgecolor="#b9b9b9",lw=1.1,zorder=0))
axz.text(0.03,1.02,"PRIORITY SEGMENT LOCATIONS",fontsize=10.0,fontweight="bold",transform=axz.transAxes)
ZM=[29.42,30.10,40.450,40.865]
axz.imshow(hsr,extent=hs_ext,cmap="gray",vmin=HSVMIN,vmax=HSVMAX,zorder=1)
axz.imshow(RGBA,extent=ext,zorder=2,interpolation="nearest")
axz.imshow(np.where(~land,1.0,np.nan),extent=ext,cmap=ListedColormap(["#cfe6f4"]),zorder=3,interpolation="nearest")
for sgm in seg_read(f"{UP}/faults_nafz.gmt"):
    a=np.array(sgm); axz.plot(a[:,0],a[:,1],color="#c1121f",lw=1.4,zorder=7)
comps.plot(ax=axz,color="white",linewidth=4.0,zorder=10,capstyle="round")
axz.add_collection(LineCollection([[p0,p1] for p0,p1 in allseg],colors=seg_cols,linewidths=2.4,zorder=11,capstyle="round"))
for k,h in enumerate(hot[:6]):
    lo,la=h["mid_lonlat"]
    if ZM[0]<lo<ZM[1] and ZM[2]<la<ZM[3]:
        axz.plot(lo,la,marker="o",ms=6.5,mfc=pcol,mec="white",mew=1.2,zorder=15)
        axz.text(lo,la+0.010,f"P{k+1}",fontsize=8.0,fontweight="bold",color="#075555",ha="center",va="bottom",zorder=16,path_effects=W1)
for nm,(lo,la) in majors.items():
    if ZM[0]<lo<ZM[1] and ZM[2]<la<ZM[3]:
        axz.plot(lo,la,marker="o",ms=5.5,mfc="white",mec=NAVY,mew=1.6,zorder=14)
        axz.text(lo,la-0.012,nm,fontsize=7.6,ha="center",va="top",color=INK,zorder=15,path_effects=W1)
axz.set_xlim(ZM[0],ZM[1]); axz.set_ylim(ZM[2],ZM[3]); axz.set_aspect(1/math.cos(math.radians(LAT0)))
axz.set_xticks([]); axz.set_yticks([])
zdl=10*1000/kx; zx,zy=ZM[0]+0.03,ZM[2]+0.016
axz.add_patch(Rectangle((zx,zy),zdl,0.006,facecolor="#111",edgecolor="#111",zorder=20))
axz.text(zx,zy+0.010,"0",fontsize=6.5,ha="center",path_effects=W1,zorder=20)
axz.text(zx+zdl,zy+0.010,"10 km",fontsize=6.5,ha="center",path_effects=W1,zorder=20)
axz.text(0.985,0.03,"P-reaches (yellow–red) sit in Moderate terrain; the near-fault\nVery-High belt (purple) lies to the south.",fontsize=6.6,color="#222",ha="right",va="bottom",transform=axz.transAxes,linespacing=1.35,path_effects=W1)
for s in axz.spines.values(): s.set_edgecolor("#b9b9b9")

credit=("Corridor exposure by co-seismic hazard class along the 85.6 km Gebze–Arifiye HSR alignment; the corridor lies within the Low and Moderate hazard classes "
        "at the 100, 250 and 500 m buffers, and the priority reaches are its longest continuous Moderate stretches.   "
        "Classes by natural breaks of CLH; the High/Very-High boundary is set at the 97th percentile to render the Very-High class, which the strict k-means break confines to the extreme near-fault fringe -- every Very-High cell is a real high-CLH cell along the NAFZ (province CLH reaches 1.23; corridor CLH reaches 0.79).   "
        "Relief: SRTM 1 arc-sec (~30 m) hillshade (GMT @earth_relief_01s); hazard grid: SRTM 30 m slope on a 90 m analysis grid; MTA 1:500 000 lithology; NAFZ (Emre et al. 2013); Akkar & Bommer (2010) GMPE (Mw 7.6, strike-slip, rock); OpenStreetMap (ODbL).   "
        "Susceptibility from expert AHP weighting (slope, lithology, faults, TWI, elevation); classes are relative.   Geographic WGS 84; lengths in UTM 35N.   "
        "Software: Python — Matplotlib, GeoPandas, rasterio, SciPy, Shapely.   Source: authors (P. Lemenkova & A. C. Zülfikar, İstanbul Technical University).")
fpc=fm.FontProperties(family=FAMILY,size=8.0); fig.canvas.draw(); r=fig.canvas.get_renderer()
target=0.955*FIGW*fig.dpi; outl=[]; line=""
for wd in credit.split(" "):
    cand=wd if not line else line+" "+wd
    if r.get_text_width_height_descent(cand,fpc,ismath=False)[0]<=target or not line: line=cand
    else: outl.append(line); line=wd
outl.append(line)
fig.text(0.035,0.028,"\n".join(outl),fontsize=8.0,va="top",ha="left",color="#333",linespacing=1.5)

fig.savefig("/mnt/user-data/outputs/fig10_hsr_exposure.png",dpi=600,bbox_inches="tight",facecolor="white")
fig.savefig("/mnt/user-data/outputs/fig10_hsr_exposure.pdf",bbox_inches="tight",facecolor="white")
print("saved fig10  VHbreak=%.3f"%VHB)
print("exposure 250m:",dict(zip(names,e250)),"total",Etot)
print("Very-High cells rendered:",int(((CLH>=VHB)&land).sum()))
