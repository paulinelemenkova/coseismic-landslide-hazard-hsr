import os, sys, zipfile, glob, shutil
UP = "/mnt/user-data/uploads"
HERE = os.path.dirname(os.path.abspath(__file__))
os.makedirs("/tmp/osm", exist_ok=True); os.makedirs("/tmp/geo", exist_ok=True)

if not os.path.exists("/tmp/osm/corridor_osm.gpkg"):
    zipfile.ZipFile(f"{UP}/corridor_osm.zip").extractall("/tmp/osm")
    _g = glob.glob("/tmp/osm/**/*.gpkg", recursive=True)
    if _g and _g[0] != "/tmp/osm/corridor_osm.gpkg":
        shutil.move(_g[0], "/tmp/osm/corridor_osm.gpkg")

if not glob.glob("/tmp/geo/**/Turkey 500k Lithology V1.TAB", recursive=True):
    zipfile.ZipFile(f"{UP}/Geological_map_data_of_Turkey.zip").extractall("/tmp/geo")
os.environ["GEO_DIR"] = os.path.dirname(
    glob.glob("/tmp/geo/**/Turkey 500k Lithology V1.TAB", recursive=True)[0])

_fac_cands = [os.path.join(HERE, "fig06_factors_srtm30.tif"),
              "/mnt/user-data/outputs/fig06_factors_srtm30.tif",
              "/tmp/fac30.tif"]
FACTOR_TIF = next((p for p in _fac_cands if os.path.exists(p)), None)
if FACTOR_TIF is None:
    sys.exit("Missing fig06_factors_srtm30.tif (SRTM 30 m factor grid). It ships "
             "with this script; regenerate it with the GMT recipe in the header.")
os.environ["FIG06_FACTOR_TIF"] = FACTOR_TIF

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

import warnings; warnings.filterwarnings("ignore")
import numpy as np, math, json, geopandas as gpd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.font_manager as fm
from matplotlib.ticker import MultipleLocator
from matplotlib.colors import ListedColormap, BoundaryNorm, LightSource
from matplotlib.patches import Rectangle
from matplotlib.lines import Line2D
from matplotlib.collections import LineCollection
from matplotlib.patheffects import withStroke
from shapely.geometry import box
from shapely.ops import unary_union
from scipy.interpolate import RegularGridInterpolator

for _f in ["LiberationSans-Regular.ttf","LiberationSans-Bold.ttf","LiberationSans-Italic.ttf"]:
    try: fm.fontManager.addfont(f"/usr/share/fonts/truetype/liberation/{_f}")
    except Exception: pass
FAMILY="TeX Gyre Heros"
plt.rcParams.update({
    "font.family":"sans-serif",
    "font.sans-serif":[FAMILY,"Liberation Sans","Nimbus Sans","Arimo"],
    "mathtext.fontset":"stixsans",
    "axes.unicode_minus":False,
    "pdf.fonttype":42,"ps.fonttype":42,"svg.fonttype":"none",
})
WHITE=[withStroke(linewidth=2,foreground="white")]

GPKG="/tmp/osm/corridor_osm.gpkg"; UP="/mnt/user-data/uploads"
import os, rasterio
d=np.load("/tmp/fig06_arrays.npz")
CLH=d["CLH"]; CLHcls=d["CLHcls"]; lons=d["lons"]; lats=d["lats"]; bnd=d["bnd"]; land=d["land"]
st=json.load(open("/tmp/fig06_stats.json")); hot=st["hot"]
W,E=lons.min(),lons.max(); S,N=lats.min(),lats.max(); ext=[W,E,S,N]
MAP=[29.20,30.42,40.47,40.99]
LAT0=float(lats.mean()); kx=111320*math.cos(math.radians(LAT0))

CLHcol=["#3288bd","#66c2a5","#f6d543","#f46d43","#c4001d"]
names=["Very Low","Low","Moderate","High","Very High"]
cmap=ListedColormap(CLHcol); norm=BoundaryNorm([.5,1.5,2.5,3.5,4.5,5.5],cmap.N)

_cand=["/tmp/hillshade_01s.tif","/mnt/user-data/outputs/fig09_hillshade_srtm90.tif","/mnt/user-data/outputs/fig06_hillshade_srtm90.tif"]
HS_TIF=next((p for p in _cand if os.path.exists(p)),None)
if HS_TIF:
    _h=rasterio.open(HS_TIF); hs=_h.read(1).astype(float)
    hs_ext=[_h.bounds.left,_h.bounds.right,_h.bounds.bottom,_h.bounds.top]; HSVMIN,HSVMAX=15,245
else:
    hs=np.zeros((10,10)); hs_ext=ext; HSVMIN,HSVMAX=0,1

FIGW,FIGH=16.6,13.9
fig=plt.figure(figsize=(FIGW,FIGH),dpi=110)
MAPL,MAPW=0.045,0.915
axm=fig.add_axes([MAPL,0.375,MAPW,0.560])
axL=fig.add_axes([MAPL,0.180,MAPW,0.175]); axL.axis("off")

axm.imshow(hs,extent=hs_ext,cmap="gray",vmin=HSVMIN,vmax=HSVMAX,zorder=1)
axm.imshow(np.where(land,CLHcls,np.nan),extent=ext,cmap=cmap,norm=norm,alpha=0.70,zorder=2,interpolation="nearest")
axm.imshow(np.where(~land,1.0,np.nan),extent=ext,cmap=ListedColormap(["#bfe0ef"]),zorder=3,interpolation="nearest")

MAPBOX=box(MAP[0],MAP[2],MAP[1],MAP[3])
def clip(g): return gpd.clip(g,MAPBOX) if len(g) else g
clip(gpd.read_file(GPKG,layer="water")).plot(ax=axm,facecolor="#aecfe0",edgecolor="#6fa8c8",linewidth=.2,zorder=4)
wl=gpd.read_file(GPKG,layer="waterways"); clip(wl[wl['fclass'].isin(['river','canal'])]).plot(ax=axm,color="#4a90c0",linewidth=.6,alpha=.75,zorder=5)
adm=gpd.read_file(GPKG,layer="adminareas"); clip(adm[adm['fclass']=='admin_level6']).boundary.plot(ax=axm,color="#6b6b6b",linewidth=.4,linestyle=(0,(5,4)),alpha=.35,zorder=5)
rd=gpd.read_file(GPKG,layer="roads"); clip(rd[rd['fclass'].isin(['motorway','motorway_link','trunk'])]).plot(ax=axm,color="#8a8a8a",linewidth=.7,alpha=.6,zorder=6)

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
for sgm in seg_read(f"{UP}/faults_nafz.gmt"):
    a=np.array(sgm); axm.plot(a[:,0],a[:,1],color="#6a0dad",lw=2.0,zorder=8,solid_capstyle="round")
    axm.plot(a[:,0],a[:,1],color="#6a0dad",lw=5,alpha=.18,zorder=7)

clip(gpd.read_file("/tmp/fig06_hsr_all.gpkg")).plot(ax=axm,color="#2b2b2b",linewidth=.5,alpha=.25,zorder=9)
comps=gpd.read_file("/tmp/fig06_hsr_main.gpkg").to_crs(32635)
u=unary_union(comps.geometry.values)
for r,al,ec in [(500,.10,"#333333"),(250,.13,None),(100,.18,None)]:
    gpd.GeoSeries([u.buffer(r)],crs=32635).to_crs(4326).plot(ax=axm,facecolor="#444444",alpha=al,edgecolor=ec,linewidth=.5,zorder=10)
CLHc=np.where(np.isfinite(CLH),CLH,0.0)[::-1,:]; interp=RegularGridInterpolator((lats[::-1],lons),CLHc,bounds_error=False,fill_value=0.0)
segsL=[]; segC=[]
for g in comps.geometry:
    n=max(2,int(g.length/120)); pts=[g.interpolate(i*g.length/n) for i in range(n+1)]
    pll=gpd.GeoSeries(pts,crs=32635).to_crs(4326); xy=np.array([[p.x,p.y] for p in pll])
    v=interp(np.column_stack([xy[:,1],xy[:,0]])); cl=np.clip(np.digitize(v,bnd),0,4)
    for i in range(len(xy)-1): segsL.append([xy[i],xy[i+1]]); segC.append(CLHcol[int(cl[i])])
axm.add_collection(LineCollection(segsL,colors="black",linewidths=4.4,zorder=11,capstyle="round"))
axm.add_collection(LineCollection(segsL,colors=segC,linewidths=2.6,zorder=12,capstyle="round"))

tr=gpd.read_file(GPKG,layer="transport"); tr=tr[tr['fclass'].isin(['railway_station','railway_halt']) & tr['name'].notna()]; tr=clip(tr)
cu=unary_union(comps.to_crs(4326).geometry.values); cu_u=gpd.GeoSeries([cu],crs=4326).to_crs(32635).iloc[0]
tr=tr[tr.geometry.apply(lambda g: gpd.GeoSeries([g],crs=4326).to_crs(32635).distance(cu_u).iloc[0]<1500)]
for _,rr in tr.iterrows(): axm.plot(rr.geometry.x,rr.geometry.y,marker="s",ms=4.5,mfc="#111",mec="white",mew=.6,zorder=14)

KEEP={"Gebze","Darıca","Çayırova","Dilovası","Körfez","Derince","İzmit","Kocaeli","Köseköy","Kartepe",
      "Başiskele","Gölcük","Karamürsel","Yalova","Sapanca","Arifiye","Adapazarı","Serdivan","Hendek","Kandıra"}
pl=gpd.read_file(GPKG,layer="places"); pl=clip(pl[pl['fclass'].isin(['city','town'])])
for _,rr in pl.iterrows():
    big=rr['fclass']=='city'
    axm.plot(rr.geometry.x,rr.geometry.y,marker="o",ms=6 if big else 3.4,mfc="#ffffff",mec="#111",mew=.8,zorder=15,clip_on=True)
    if big or rr['name'] in KEEP:
        axm.text(rr.geometry.x,rr.geometry.y+0.010,rr['name'],fontsize=9.5 if big else 8.0,ha="center",va="bottom",
                 fontweight="bold" if big else "normal",color="#111",zorder=16,clip_on=True,path_effects=WHITE)

for k,h in enumerate(hot[:4],1):
    lo,la=h["mid_lonlat"]
    axm.plot(lo,la,marker="o",ms=15,mfc="none",mec="#c4001d",mew=2.2,zorder=17)
    axm.text(lo,la-0.016,f"H{k}",fontsize=10.5,fontweight="bold",ha="center",va="top",color="#c4001d",zorder=18,clip_on=True,path_effects=WHITE)

axm.set_xlim(MAP[0],MAP[1]); axm.set_ylim(MAP[2],MAP[3]); axm.set_aspect(1/math.cos(math.radians(LAT0)))
axm.xaxis.set_major_locator(MultipleLocator(0.25)); axm.xaxis.set_minor_locator(MultipleLocator(0.05))
axm.yaxis.set_major_locator(MultipleLocator(0.10)); axm.yaxis.set_minor_locator(MultipleLocator(0.025))
axm.xaxis.set_major_formatter(plt.FuncFormatter(lambda v,_:f"{v:.2f}"))
axm.yaxis.set_major_formatter(plt.FuncFormatter(lambda v,_:f"{v:.2f}"))
axm.tick_params(which="major",length=5,width=.9,labelsize=10.5,direction="out",top=True,right=True)
axm.tick_params(which="minor",length=2.6,width=.6,direction="out",top=True,right=True)
axm.set_xlabel("Longitude (°E)",fontsize=12); axm.set_ylabel("Latitude (°N)",fontsize=12)
for s in axm.spines.values(): s.set_linewidth(1.0)
axm.grid(which="major",color="white",linewidth=.5,alpha=.30,zorder=3.4)

seg_km=20; dlon=seg_km*1000/kx; x0,y0=MAP[0]+0.06,MAP[2]+0.035
axm.add_patch(Rectangle((x0,y0),dlon,0.006,facecolor="#111",edgecolor="#111",zorder=20))
axm.add_patch(Rectangle((x0,y0),dlon/2,0.006,facecolor="#fff",edgecolor="#111",zorder=20))
axm.text(x0+dlon/2,y0+0.012,"20 km",ha="center",va="bottom",fontsize=8.5,fontweight="bold",zorder=20,path_effects=WHITE)
nxp,nyp=MAP[1]-0.045,MAP[3]-0.055
axm.annotate("N",xy=(nxp,nyp),xytext=(nxp,nyp-0.05),ha="center",va="center",fontsize=12.5,fontweight="bold",zorder=20,
             arrowprops=dict(arrowstyle="-|>",color="#111",lw=2.0),path_effects=WHITE)

axm.set_title("Co-Seismic Landslide Hazard (CLH) and Exposure of the İstanbul–Kocaeli High-Speed Railway Corridor\n"
              "1999 İzmit Scenario (Mw 7.6, North Anatolian Fault Zone)",fontsize=14.5,fontweight="bold",pad=10)

HEAD=13.0; ENTRY=12.0; TABLE=11.0
def sw(x,y,c): axL.add_patch(Rectangle((x,y),0.020,0.13,transform=axL.transAxes,facecolor=c,edgecolor="#333",lw=.6,zorder=3))

axL.text(0.0,0.95,"Co-seismic landslide hazard (CLH)",fontsize=HEAD,fontweight="bold",transform=axL.transAxes)
for i,(c,nm) in enumerate(zip(CLHcol,names)):
    x=i*0.098; sw(x,0.60,c); axL.text(x+0.026,0.665,nm,fontsize=ENTRY,va="center",transform=axL.transAxes)

hx=0.55
axL.text(hx,0.95,"Railway corridor & context",fontsize=HEAD,fontweight="bold",transform=axL.transAxes)
items=[("YHT corridor (colour = CLH class)","line","#f6d543"),
       ("100 / 250 / 500 m exposure buffers","band","#444444"),
       ("NAFZ active fault","line","#6a0dad"),
       ("Priority (Moderate) segment  (H)","hot","#c4001d"),
       ("Railway station","sta","#111"),
       ("Motorway","road","#8a8a8a"),
       ("River / canal","river","#4a90c0")]
for i,(lab,kind,c) in enumerate(items):
    yy=0.80-i*0.118
    if kind=="line" and c=="#f6d543":
        axL.add_line(Line2D([hx+0.004,hx+0.044],[yy,yy],color="black",lw=4.6,transform=axL.transAxes,zorder=1,solid_capstyle="round"))
        axL.add_line(Line2D([hx+0.004,hx+0.044],[yy,yy],color=c,lw=2.8,transform=axL.transAxes,zorder=2,solid_capstyle="round"))
    elif kind=="line":
        axL.add_line(Line2D([hx+0.004,hx+0.044],[yy,yy],color=c,lw=3.2,transform=axL.transAxes,solid_capstyle="round"))
    if kind=="band": axL.add_patch(Rectangle((hx+0.004,yy-0.03),0.040,0.06,transform=axL.transAxes,facecolor="#444",alpha=.3,edgecolor="#333",lw=.5))
    if kind=="hot": axL.plot(hx+0.024,yy,marker="o",ms=13,mfc="none",mec="#c4001d",mew=2.2,transform=axL.transAxes)
    if kind=="sta": axL.plot(hx+0.024,yy,marker="s",ms=7,mfc="#111",mec="white",mew=.7,transform=axL.transAxes)
    if kind=="road": axL.add_line(Line2D([hx+0.004,hx+0.044],[yy,yy],color="#8a8a8a",lw=2.2,transform=axL.transAxes,solid_capstyle="round"))
    if kind=="river": axL.add_line(Line2D([hx+0.004,hx+0.044],[yy,yy],color="#4a90c0",lw=2.2,transform=axL.transAxes,solid_capstyle="round"))
    axL.text(hx+0.052,yy,lab,fontsize=ENTRY,va="center",transform=axL.transAxes)

tx=0.82; axL.text(tx,0.95,"Corridor exposure  (track-km)",fontsize=HEAD,fontweight="bold",transform=axL.transAxes)
rows=[("Class","100 m","250 m","500 m")]+[(names[i],f"{st['expo']['100'][i]:.1f}",f"{st['expo']['250'][i]:.1f}",f"{st['expo']['500'][i]:.1f}") for i in range(4,-1,-1)]
xoff=[0.0,0.066,0.110,0.154]
for r,row in enumerate(rows):
    yy=0.80-r*0.132
    for cci,val in enumerate(row):
        axL.text(tx+xoff[cci],yy,val,fontsize=TABLE,va="center",transform=axL.transAxes,
                 fontweight="bold" if r==0 else "normal",
                 color=(CLHcol[4-(r-1)] if (r>=1 and cci==0) else "#111"))

CREDIT_PT=11.0
credit=("Relief: SRTM 1 arc-sec (~30 m) hillshade (GMT @earth_relief_01s).   Terrain factors (slope, TWI, elevation): SRTM 1 arc-sec, 90 m analysis grid.   Lithology: MTA 1:500 000 Geological Map of Türkiye.   "
        "Active faults: NAFZ, Emre et al. (2013).   Ground motion: Akkar & Bommer (2010) GMPE — İzmit 1999 scenario "
        "Mw 7.6, strike-slip, rock site; modelled PGA "
        f"{st['pga_min']:.2f}–{st['pga_max']:.2f} g.   CLH = LSI (5-factor expert-weighted AHP: slope, "
        "distance-to-fault, lithology, TWI, elevation) × ST, with ST = min(PGA / 0.4g, 1).   "
        "Railways (high-speed = OSM ‘yüksek hızlı demiryolu’), roads, water, places and districts: "
        "© OpenStreetMap contributors (ODbL).   "
        f"Corridor {st['corr_km']:.0f} km; peak corridor CLH upper-Moderate, lying in the Low–Moderate classes within 500 m; "
        f"closest approach to High/Very-High hazard terrain {st['approach_km']:.1f} km.   "
        "Susceptibility from expert AHP weighting; classes are relative.   "
        "Projection UTM 35N (EPSG:32635).   Software: Python — Matplotlib, GeoPandas, rasterio, SciPy, Shapely.   "
        "Source: authors (P. Lemenkova & A. C. Zülfikar, İstanbul Technical University).")

def wrap_to_width(fig,text,target_px,fontsize,family):
    fp=fm.FontProperties(family=family,size=fontsize)
    r=fig.canvas.get_renderer(); out=[]; line=""
    for w in text.split(" "):
        cand=w if not line else line+" "+w
        if r.get_text_width_height_descent(cand,fp,ismath=False)[0]<=target_px or not line: line=cand
        else: out.append(line); line=w
    out.append(line); return "\n".join(out)

fig.canvas.draw()
target_px=MAPW*FIGW*fig.dpi*0.995
credit_w=wrap_to_width(fig,credit,target_px,CREDIT_PT,FAMILY)
fig.text(MAPL,0.152,credit_w,fontsize=CREDIT_PT,va="top",ha="left",color="#222",linespacing=1.55,
         bbox=dict(boxstyle="round,pad=0.6",fc="#f4f4f2",ec="#cfcfcf",lw=.8))

fig.savefig("/mnt/user-data/outputs/fig06_corridor_clh.png",dpi=600,bbox_inches="tight",facecolor="white")
fig.savefig("/mnt/user-data/outputs/fig06_corridor_clh.pdf",bbox_inches="tight",facecolor="white")
print("saved fig06_corridor_clh.png / .pdf")
