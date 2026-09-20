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

import warnings; warnings.filterwarnings("ignore")
import os, numpy as np, math, json, geopandas as gpd, rasterio
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.font_manager as fm
from matplotlib.ticker import MultipleLocator
from matplotlib.colors import ListedColormap, BoundaryNorm, LightSource
from matplotlib.patches import Rectangle, FancyBboxPatch
from matplotlib.lines import Line2D
from matplotlib.patheffects import withStroke
from rasterio.features import rasterize
from rasterio.transform import Affine
from shapely.geometry import shape, box
from shapely.ops import unary_union

for _f in ["LiberationSans-Regular.ttf","LiberationSans-Bold.ttf","LiberationSans-Italic.ttf"]:
    try: fm.fontManager.addfont(f"/usr/share/fonts/truetype/liberation/{_f}")
    except Exception: pass
FAMILY="TeX Gyre Heros"
plt.rcParams.update({"font.family":"sans-serif","font.sans-serif":[FAMILY,"Liberation Sans","Nimbus Sans"],
    "mathtext.fontset":"stixsans","axes.unicode_minus":False,"pdf.fonttype":42,"ps.fonttype":42,"svg.fonttype":"none"})
W1=[withStroke(linewidth=2,foreground="white")]; W2=[withStroke(linewidth=2.6,foreground="white")]
NAVY="#12305e"; INK="#1a1a1a"

UP="/mnt/user-data/uploads"; GPKG="/tmp/osm/corridor_osm.gpkg"
d=np.load("/tmp/fig06_arrays.npz"); CLHcls=d["CLHcls"]; land=d["land"]; lons=d["lons"]; lats=d["lats"]
st=json.load(open("/tmp/fig06_stats.json")); expo=st["expo"]
ny,nx=CLHcls.shape; rx=lons[1]-lons[0]; ry=lats[1]-lats[0]
tr=Affine.translation(lons[0]-rx/2,lats[0]-ry/2)*Affine.scale(rx,ry)
ext=[lons.min(),lons.max(),lats.min(),lats.max()]
LAT0=float(lats.mean()); kx=111320*math.cos(math.radians(LAT0))
MAP=[28.75,30.90,40.40,41.30]

CLHcol=["#3288bd","#66c2a5","#f6d543","#f46d43","#c4001d"]; names=["Very Low","Low","Moderate","High","Very High"]
cmap=ListedColormap(CLHcol); norm=BoundaryNorm([.5,1.5,2.5,3.5,4.5,5.5],cmap.N)

prov=json.load(open("/tmp/ne_adm1.json"))["features"]
def prov_geom(nm):
    g=[shape(f["geometry"]) for f in prov if f["properties"].get("iso_a2")=="TR" and f["properties"].get("name")==nm]
    return unary_union(g) if g else None
study=unary_union([prov_geom("Istanbul"),prov_geom("Kocaeli")])
turkey=unary_union([shape(f["geometry"]) for f in prov if f["properties"].get("iso_a2")=="TR"])

mask=rasterize([(study,1)],out_shape=(ny,nx),transform=tr,fill=0,all_touched=True).astype(bool)
LATg=np.repeat(lats[:,None],nx,axis=1)
cell=(abs(rx)*111320*np.cos(np.radians(LATg)))*(abs(ry)*110570)/1e6
area=[float(cell[(CLHcls==c)&mask&land].sum()) for c in range(1,6)]; Atot=sum(area)
apct=[100*a/Atot for a in area]
e250=expo["250"]; Etot=sum(e250); epct=[100*x/Etot for x in e250]

_cand=["/tmp/hillshade_01s.tif","/mnt/user-data/outputs/fig09_hillshade_srtm90.tif",
       os.path.join(os.path.dirname(os.path.abspath(__file__)) if "__file__" in globals() else "/tmp","fig09_hillshade_srtm90.tif")]
HS_TIF=next((p for p in _cand if os.path.exists(p)),None)
if HS_TIF:
    _h=rasterio.open(HS_TIF); hsr=_h.read(1).astype(float)
    hs_ext=[_h.bounds.left,_h.bounds.right,_h.bounds.bottom,_h.bounds.top]
    HSVMIN,HSVMAX=15,245; RELIEF="SRTM 1 arc-sec (~30 m), GMT @earth_relief_01s"
else:
    hsr=np.full((10,10),128.0)
    hs_ext=ext; HSVMIN,HSVMAX=0,255; RELIEF="(relief backdrop unavailable)"

FIGW,FIGH=16.4,11.6
fig=plt.figure(figsize=(FIGW,FIGH),dpi=110)

fig.add_artist(Rectangle((0,0.945),1,0.055,transform=fig.transFigure,facecolor=NAVY,edgecolor="none",zorder=0))
fig.text(0.5,0.983,"Co-Seismic Landslide Hazard under the 1999 İzmit (Kocaeli) Scenario Earthquake",
         ha="center",va="center",fontsize=17,fontweight="bold",color="white")
fig.text(0.5,0.958,"Expert-weighted AHP susceptibility × Akkar–Bommer (2010) PGA triggering factor — İstanbul–Kocaeli",
         ha="center",va="center",fontsize=11.5,color="#cdd8ea")

axm=fig.add_axes([0.035,0.335,0.66,0.585])
axm.imshow(hsr,extent=hs_ext,cmap="gray",vmin=HSVMIN,vmax=HSVMAX,zorder=1)
axm.imshow(np.where(land,CLHcls,np.nan),extent=ext,cmap=cmap,norm=norm,alpha=0.72,zorder=2,interpolation="nearest")
axm.imshow(np.where(~land,1.0,np.nan),extent=ext,cmap=ListedColormap(["#cfe6f4"]),zorder=3,interpolation="nearest")

def clipL(g):
    try: return gpd.clip(g,box(*MAP[:1],*MAP[2:3],*MAP[1:2],*MAP[3:4])) if len(g) else g
    except Exception: return g
MB=box(MAP[0],MAP[2],MAP[1],MAP[3])
def clip(g): return gpd.clip(g,MB) if len(g) else g
clip(gpd.read_file(GPKG,layer="water")).plot(ax=axm,facecolor="#bcd9ec",edgecolor="#7fb0cf",linewidth=.2,zorder=4)
wl=gpd.read_file(GPKG,layer="waterways"); clip(wl[wl['fclass'].isin(['river','canal'])]).plot(ax=axm,color="#5a97c4",linewidth=.5,alpha=.7,zorder=5)

for nm in ["Istanbul","Kocaeli","Sakarya","Yalova","Bursa","Bilecik"]:
    g=prov_geom(nm)
    if g is not None:
        gpd.GeoSeries([g],crs=4326).boundary.plot(ax=axm,color="#3a3a3a",linewidth=.9,linestyle=(0,(6,3)),alpha=.55,zorder=6)

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
    a=np.array(sgm); axm.plot(a[:,0],a[:,1],color="#c1121f",lw=2.1,zorder=8,solid_capstyle="round")

comps=gpd.read_file("/tmp/fig06_hsr_main.gpkg")
comps.plot(ax=axm,color="white",linewidth=4.2,zorder=10,capstyle="round")
comps.plot(ax=axm,color=NAVY,linewidth=2.2,zorder=11,capstyle="round")

majors={"Gebze":(29.431,40.803),"Derince":(29.828,40.760),"İzmit":(29.945,40.766),
        "Köseköy":(30.050,40.766),"Arifiye":(30.373,40.717)}
cu=unary_union(comps.geometry.values)
for nm,(lo,la) in majors.items():
    axm.plot(lo,la,marker="o",ms=8,mfc="white",mec=NAVY,mew=2,zorder=14)
    axm.text(lo,la+0.028,nm,fontsize=9,ha="center",va="bottom",color=INK,zorder=15,path_effects=W1)

axm.text(29.05,41.13,"İstanbul",fontsize=15,fontweight="bold",color="#33322f",ha="center",zorder=15,path_effects=W2)
axm.text(30.28,40.98,"Kocaeli",fontsize=15,fontweight="bold",color="#33322f",ha="center",zorder=15,path_effects=W2)
axm.text(30.62,40.86,"Sakarya",fontsize=12.5,fontweight="bold",color="#4a4844",ha="center",zorder=15,path_effects=W2)
axm.text(29.45,41.22,"Black Sea",fontsize=12,style="italic",color="#1f5f86",ha="center",zorder=15)
axm.text(28.93,40.80,"Marmara Sea",fontsize=12,style="italic",color="#1f5f86",ha="center",zorder=15)
axm.text(29.72,40.66,"İzmit Bay",fontsize=8.5,style="italic",color="#1f5f86",ha="center",zorder=15)
axm.text(29.10,40.83,"NAFZ",fontsize=8.5,fontweight="bold",color="#c1121f",rotation=8,zorder=15,path_effects=W1)

axm.set_xlim(MAP[0],MAP[1]); axm.set_ylim(MAP[2],MAP[3]); axm.set_aspect(1/math.cos(math.radians(LAT0)))
axm.xaxis.set_major_locator(MultipleLocator(0.5)); axm.xaxis.set_minor_locator(MultipleLocator(0.1))
axm.yaxis.set_major_locator(MultipleLocator(0.25)); axm.yaxis.set_minor_locator(MultipleLocator(0.05))
axm.xaxis.set_major_formatter(plt.FuncFormatter(lambda v,_:f"{v:.1f}°E"))
axm.yaxis.set_major_formatter(plt.FuncFormatter(lambda v,_:f"{v:.1f}°N"))
axm.tick_params(which="major",length=5,width=.9,labelsize=9.5,top=True,right=True)
axm.tick_params(which="minor",length=2.6,width=.6,top=True,right=True)
for s in axm.spines.values(): s.set_linewidth(1.1)

seg=40; dlon=seg*1000/kx; x0,y0=MAP[0]+0.06,MAP[2]+0.05
axm.add_patch(Rectangle((x0,y0),dlon,0.010,facecolor="#111",edgecolor="#111",zorder=20))
axm.add_patch(Rectangle((x0,y0),dlon/2,0.010,facecolor="#fff",edgecolor="#111",zorder=20))
axm.text(x0,y0+0.02,"0",ha="center",va="bottom",fontsize=8,zorder=20,path_effects=W1)
axm.text(x0+dlon/2,y0+0.02,"20",ha="center",va="bottom",fontsize=8,zorder=20,path_effects=W1)
axm.text(x0+dlon,y0+0.02,"40 km",ha="center",va="bottom",fontsize=8,zorder=20,path_effects=W1)
nxp,nyp=MAP[0]+0.10,MAP[3]-0.10
axm.annotate("N",xy=(nxp,nyp),xytext=(nxp,nyp-0.10),ha="center",va="center",fontsize=13,fontweight="bold",zorder=20,
             arrowprops=dict(arrowstyle="-|>",color="#111",lw=2),path_effects=W2)

axi=fig.add_axes([0.545,0.775,0.15,0.135]); axi.set_facecolor("#eaf3fb")
gpd.GeoSeries([turkey],crs=4326).plot(ax=axi,facecolor="#dfe6d8",edgecolor="#8a8a8a",linewidth=.5)
sb=study.bounds
axi.add_patch(Rectangle((sb[0],sb[1]),sb[2]-sb[0],sb[3]-sb[1],facecolor="none",edgecolor="#c1121f",lw=1.6,zorder=5))
axi.set_xlim(25.5,45.2); axi.set_ylim(35.6,42.5); axi.set_aspect(1/math.cos(math.radians(39)))
axi.text(0.5,1.04,"Türkiye",transform=axi.transAxes,ha="center",va="bottom",fontsize=9,fontweight="bold")
axi.set_xticks([]); axi.set_yticks([])
for s in axi.spines.values(): s.set_edgecolor("#8a8a8a")

axR=fig.add_axes([0.70,0.335,0.281,0.585]); axR.axis("off")
axR.add_patch(FancyBboxPatch((0.0,0.0),1,1,boxstyle="round,pad=0.01",transform=axR.transAxes,
              facecolor="white",edgecolor="#c9c9c9",lw=1))
axR.text(0.06,0.955,"LEGEND",fontsize=13,fontweight="bold",transform=axR.transAxes)

fy=0.905
axR.add_line(Line2D([0.07,0.15],[fy,fy],color="white",lw=4.2,transform=axR.transAxes)); axR.add_line(Line2D([0.07,0.15],[fy,fy],color=NAVY,lw=2.2,transform=axR.transAxes))
axR.text(0.18,fy,"HSR corridor (Ankara–İstanbul YHT)",fontsize=9.5,va="center",transform=axR.transAxes)
fy=0.86; axR.plot(0.11,fy,marker="o",ms=8,mfc="white",mec=NAVY,mew=2,transform=axR.transAxes); axR.text(0.18,fy,"Major stations",fontsize=9.5,va="center",transform=axR.transAxes)
fy=0.815; axR.add_line(Line2D([0.07,0.15],[fy,fy],color="#c1121f",lw=2.2,transform=axR.transAxes)); axR.text(0.18,fy,"Active faults (NAFZ)",fontsize=9.5,va="center",transform=axR.transAxes)
fy=0.77; axR.add_line(Line2D([0.07,0.15],[fy,fy],color="#3a3a3a",lw=1.1,ls=(0,(6,3)),transform=axR.transAxes)); axR.text(0.18,fy,"Province boundary",fontsize=9.5,va="center",transform=axR.transAxes)

axR.text(0.06,0.70,"Co-seismic landslide hazard",fontsize=11,fontweight="bold",transform=axR.transAxes)
axR.text(0.06,0.665,"(CLH class)",fontsize=9.5,transform=axR.transAxes)
for i,(c,nm) in enumerate(zip(CLHcol,names)):
    yy=0.60-i*0.052
    axR.add_patch(Rectangle((0.07,yy-0.017),0.09,0.034,transform=axR.transAxes,facecolor=c,edgecolor="#555",lw=.6))
    axR.text(0.19,yy,nm,fontsize=9.5,va="center",transform=axR.transAxes)

axR.text(0.06,0.30,"Hazard calculation",fontsize=11,fontweight="bold",transform=axR.transAxes)
note=("CLH = AHP susceptibility (slope, distance-to-fault, lithology, TWI, "
      "elevation) × ground-motion scaling ST, with ST = min(PGA / 0.4g, 1). "
      "PGA from the Akkar & Bommer (2010) GMPE for a 1999 İzmit-type rupture "
      "(Mw 7.6, strike-slip, rock). Five classes by natural breaks (Jenks/k-means).")
def wrapax(ax,text,x,y,w_frac,fs):
    fp=fm.FontProperties(family=FAMILY,size=fs); r=fig.canvas.get_renderer()
    tw=w_frac*ax.get_window_extent().width; out=[]; line=""
    for wd in text.split(" "):
        cand=wd if not line else line+" "+wd
        if r.get_text_width_height_descent(cand,fp,ismath=False)[0]<=tw or not line: line=cand
        else: out.append(line); line=wd
    out.append(line); ax.text(x,y,"\n".join(out),fontsize=fs,va="top",transform=ax.transAxes,linespacing=1.5)
fig.canvas.draw()
wrapax(axR,note,0.06,0.265,0.88,9.0)

def panel(x,w,title):
    ax=fig.add_axes([x,0.055,w,0.245]); ax.axis("off")
    ax.add_patch(FancyBboxPatch((0,0),1,1,boxstyle="round,pad=0.008",transform=ax.transAxes,facecolor="white",edgecolor="#c9c9c9",lw=1))
    ax.text(0.03,0.93,title,fontsize=10.5,fontweight="bold",transform=ax.transAxes)
    return ax

axd=panel(0.035,0.30,"HAZARD CLASS — DESCRIPTION")
desc=["Very low probability of co-seismic landsliding","Low probability of co-seismic landsliding",
      "Moderate probability of co-seismic landsliding","High probability of co-seismic landsliding",
      "Very high probability of co-seismic landsliding"]
for i,(c,nm,ds) in enumerate(zip(CLHcol,names,desc)):
    yy=0.79-i*0.15
    axd.add_patch(Rectangle((0.03,yy-0.055),0.20,0.11,transform=axd.transAxes,facecolor=c,edgecolor="#555",lw=.5))
    axd.text(0.13,yy,nm,fontsize=8.3,ha="center",va="center",fontweight="bold",transform=axd.transAxes,color="#111",path_effects=W1)
    axd.text(0.27,yy,ds,fontsize=8.6,va="center",transform=axd.transAxes)

axs=panel(0.345,0.185,"HAZARD-CLASS AREA")
axs.text(0.30,0.80,"km²",fontsize=8.6,fontweight="bold",ha="right",transform=axs.transAxes)
axs.text(0.62,0.80,"%",fontsize=8.6,fontweight="bold",ha="right",transform=axs.transAxes)
for i in range(5):
    yy=0.66-i*0.115
    axs.add_patch(Rectangle((0.03,yy-0.045),0.055,0.09,transform=axs.transAxes,facecolor=CLHcol[i],edgecolor="#555",lw=.4))
    axs.text(0.11,yy,names[i],fontsize=8.3,va="center",transform=axs.transAxes)
    axs.text(0.62,yy,f"{area[i]:,.0f}",fontsize=8.6,ha="right",va="center",transform=axs.transAxes)
    axs.text(0.95,yy,f"{apct[i]:.1f}",fontsize=8.6,ha="right",va="center",transform=axs.transAxes)
axs.text(0.11,0.075,"Total",fontsize=8.6,fontweight="bold",va="center",transform=axs.transAxes)
axs.text(0.62,0.075,f"{Atot:,.0f}",fontsize=8.6,fontweight="bold",ha="right",va="center",transform=axs.transAxes)
axs.text(0.95,0.075,"100.0",fontsize=8.6,fontweight="bold",ha="right",va="center",transform=axs.transAxes)
axs.text(0.03,0.02,"İstanbul + Kocaeli provinces only",fontsize=7.3,style="italic",color="#555",va="bottom",transform=axs.transAxes)

axe=panel(0.545,0.205,"CORRIDOR EXPOSURE  (250 m buffer)")
axe.text(0.34,0.80,"km",fontsize=8.6,fontweight="bold",ha="right",transform=axe.transAxes)
axe.text(0.62,0.80,"%",fontsize=8.6,fontweight="bold",ha="right",transform=axe.transAxes)
for i in range(5):
    yy=0.66-i*0.115
    axe.add_patch(Rectangle((0.03,yy-0.045),0.05,0.09,transform=axe.transAxes,facecolor=CLHcol[i],edgecolor="#555",lw=.4))
    axe.text(0.10,yy,names[i],fontsize=8.3,va="center",transform=axe.transAxes)
    axe.text(0.60,yy,f"{e250[i]:.1f}",fontsize=8.6,ha="right",va="center",transform=axe.transAxes)
    axe.text(0.92,yy,f"{epct[i]:.1f}",fontsize=8.6,ha="right",va="center",transform=axe.transAxes)
axe.text(0.10,0.075,"Total",fontsize=8.6,fontweight="bold",va="center",transform=axe.transAxes)
axe.text(0.60,0.075,f"{Etot:.1f}",fontsize=8.6,fontweight="bold",ha="right",va="center",transform=axe.transAxes)
axe.text(0.92,0.075,"100.0",fontsize=8.6,fontweight="bold",ha="right",va="center",transform=axe.transAxes)
axe.text(0.03,0.02,"Gebze–Arifiye alignment; corridor lies in Low–Moderate within 500 m",fontsize=7.0,style="italic",color="#555",va="bottom",transform=axe.transAxes)

axb=panel(0.755,0.21,"CORRIDOR BUFFERS & ELEVATION")
for i,(c,lab,sub) in enumerate([("#8e44ad","100 m","immediate track zone"),("#2e78c4","250 m","extended influence"),("#2e8b57","500 m","catchment-scale reach")]):
    yy=0.80-i*0.13
    axb.add_line(Line2D([0.04,0.15],[yy,yy],color=c,lw=1.8,ls=(0,(5,3)),transform=axb.transAxes))
    axb.text(0.19,yy,lab,fontsize=8.8,fontweight="bold",va="center",transform=axb.transAxes)
    axb.text(0.34,yy,sub,fontsize=8.0,style="italic",color="#555",va="center",transform=axb.transAxes)
axb.text(0.03,0.36,"Base relief & projection",fontsize=9,fontweight="bold",transform=axb.transAxes)
axb.text(0.03,0.28,"Relief: SRTM 1 arc-sec (~30 m) hillshade\n(GMT @earth_relief_01s). Hazard grid: SRTM\n30 m slope, 90 m analysis grid.",fontsize=7.6,va="top",color="#444",transform=axb.transAxes,linespacing=1.4)
axb.text(0.03,0.13,"WGS 84 geographic; UTM 35N for lengths.\nScale ≈ 1:400 000.",fontsize=7.8,va="top",color="#444",transform=axb.transAxes,linespacing=1.4)

credit=("Area statistics cover the İstanbul + Kocaeli provinces; higher-hazard terrain continues south of the study area (dashed boundaries).   "
        "Relief: SRTM 1 arc-sec hillshade (GMT @earth_relief_01s).   Hazard grid & terrain factors: SRTM 1 arc-sec (~30 m) slope on a 90 m analysis grid; MTA 1:500 000 lithology; NAFZ (Emre et al. 2013); Akkar & Bommer (2010) GMPE; OpenStreetMap (ODbL).   "
        "Susceptibility from expert AHP weighting; classes are relative.   "
        "Software: Python — Matplotlib, GeoPandas, rasterio, SciPy, Shapely.   Source: authors (P. Lemenkova & A. C. Zülfikar, İstanbul Technical University).")
fpc=fm.FontProperties(family=FAMILY,size=8.2); r=fig.canvas.get_renderer()
target=0.955*FIGW*fig.dpi; out=[]; line=""
for wd in credit.split(" "):
    cand=wd if not line else line+" "+wd
    if r.get_text_width_height_descent(cand,fpc,ismath=False)[0]<=target or not line: line=cand
    else: out.append(line); line=wd
out.append(line)
fig.text(0.035,0.028,"\n".join(out),fontsize=8.2,va="top",ha="left",color="#333",linespacing=1.5)

fig.savefig("/mnt/user-data/outputs/fig09_cosei_hazard.png",dpi=600,bbox_inches="tight",facecolor="white")
fig.savefig("/mnt/user-data/outputs/fig09_cosei_hazard.pdf",bbox_inches="tight",facecolor="white")
print("saved fig09")
