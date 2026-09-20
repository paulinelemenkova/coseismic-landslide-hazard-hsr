# Co-seismic landslide hazard — İstanbul–Kocaeli high-speed-railway corridor

Figure-generation code for a study of co-seismic landslide hazard along the
İstanbul–Kocaeli segment of the Ankara–İstanbul high-speed railway (HSR), under a
1999 İzmit (Mw 7.6) North Anatolian Fault Zone (NAFZ) scenario.

The workflow couples a knowledge-driven Analytic Hierarchy Process (AHP)
susceptibility index (LSI) with a deterministic ground-motion triggering factor
(Akkar & Bommer, 2010 GMPE), giving a relative co-seismic landslide hazard (CLH)
surface that is then intersected with the railway corridor and its buffers.

## Figures

| Script | Output | Tool |
|--------|--------|------|
| `scripts/fig01/make_fig01_study_area.sh` | Study-area / locator map | GMT |
| `scripts/fig02/make_fig02_scenario_pga.sh` | Scenario peak ground acceleration (PGA) | GMT |
| `scripts/make_fig03_workflow.py` | Methodology workflow diagram | Matplotlib |
| `scripts/make_fig04_ahp_sensitivity.py` | AHP weight sensitivity | Matplotlib |
| `scripts/make_fig05_corridor_buffers.sh` | Corridor + multi-buffer exposure map | GMT + OSM |
| `scripts/make_fig06_corridor_clh.py` | Corridor co-seismic hazard (CLH) detail | Matplotlib |
| `scripts/make_fig07_ew_framework.py` | Cascading early-warning framework | Matplotlib |
| `scripts/make_fig08_susceptibility.py` | Landslide susceptibility (LSI) map | Matplotlib |
| `scripts/make_fig09_cosei_hazard.py` | Regional co-seismic hazard map | Matplotlib |
| `scripts/make_fig10_hsr_exposure.py` | Corridor exposure by hazard class | Matplotlib |

## Requirements

- **Python 3.11+** with the packages in `requirements.txt`
  (numpy, scipy, pandas, matplotlib, geopandas, shapely, rasterio, pyogrio).
- **GMT 6.x** for the shell scripts (`fig01`, `fig02`, `fig05`); `fig01`/`fig02`
  need internet access to the GMT data server on first run, and `fig05` queries
  the OpenStreetMap Overpass API.
- A sans-serif font (Helvetica / Arial / Nimbus Sans / TeX Gyre Heros) for the
  Matplotlib figures.

```bash
python -m venv .venv && source .venv/bin/activate
pip install -r requirements.txt
```

## Input data

The scripts read the following datasets (not redistributed here; obtain from the
providers). The SRTM 1 arc-sec hillshade used as a relief backdrop
(`scripts/fig09_hillshade_srtm90.tif`) is included.

- **SRTM 1 arc-sec (~30 m)** elevation — terrain factors (slope, TWI, elevation)
  and hillshade (via GMT `@earth_relief_01s`).
- **MTA 1:500 000 Geological Map of Türkiye** — lithology.
- **NAFZ active-fault traces** (Emre et al., 2013).
- **OpenStreetMap** (ODbL) — railway, roads, water, places.
- **Natural Earth admin-1** — province polygons (downloaded automatically if absent).

The Python hazard scripts (`fig06`, `fig08`, `fig09`, `fig10`) additionally read a
pre-computed factor grid, `fig06_factors_srtm30.tif` — a 3-band `uint8` GeoTIFF of
slope / elevation / TWI reclassified to 1–5 on a 90 m analysis grid. It is derived
from the SRTM tile with GMT, e.g.:

```
gmt grdcut @earth_relief_01s -R28/30.9/40.35/41.35 -Gdem.nc
# slope at 1", aggregate 3x3 -> 90 m mean; elevation and TWI at 90 m;
# reclass slope [5,12,20,30]° and elevation/TWI by 20/40/60/80th land percentiles
# to classes 1-5; stack as a 3-band uint8 GeoTIFF -> fig06_factors_srtm30.tif
```

## Method notes

- Susceptibility (LSI) is an expert-weighted AHP linear combination of five
  reclassified factors (slope, distance-to-fault, lithology, TWI, elevation);
  the hazard classes are **relative**.
- The triggering factor is `ST = min(PGA / 0.4g, 1)`; `CLH = LSI × ST`.
- Coordinates are WGS 84 geographic; lengths are measured in UTM 35N (EPSG:32635).

## Paths

Input and output locations are set as absolute paths near the top of each script
(inputs directory, outputs directory). Edit those paths, or place the input files
accordingly, before running.

## Authors

P. Lemenkova and A. C. Zülfikar, İstanbul Technical University.

## License

Released under the MIT License (see `LICENSE`). Input datasets remain under their
own licenses (OpenStreetMap ODbL; MTA; SRTM; Natural Earth; etc.).
