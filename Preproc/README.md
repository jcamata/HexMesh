# Preproc — geometry preparation for HexMesh (Python)

Builds the two GTS surfaces HexMesh needs — topography and bathymetry (the "interface") —
from a global elevation grid plus a vector coastline. It is the Python replacement for the
Matlab workflow in `../matlab/mainSRTM.m`, and reproduces its `bathy1 + bathy2 + vertical
wall` construction.

## Install

```sh
pip install -r requirements.txt
```

`scipy` (Delaunay) and `netCDF4`/`xarray` are the only hard requirements; `geopandas`/`shapely`
are needed for vector coastlines, `trimesh` only for STL export (currently commented out in
`mesh_builder.py`).

## Run

```sh
python3 preproc.py mauna_loa.input      # one config
./run_all.sh                            # every *.input that has an nc_file key
./run_all.sh hawaii.input kefalonia.input
```

`run_all.sh` writes `log/<name>.log` per config and returns non-zero if any run failed.
It selects configs by `grep -l '^\s*nc_file'`, which excludes the `*_material.input` files
(those are *outputs*, not inputs).

Move the resulting `.gts` files to `../input/` and point `topo` / `inter` in
`../HexMesh.input` at them.

## Pipeline

`preproc.py` runs, in order:

1. `downloader.py` — fetch the NetCDF grid and the coastline GeoJSON if missing
   (`{lat_min}`, `{lat_max}`, `{lon_min}`, `{lon_max}` are substituted into `gebco_url`).
2. `elevation.py` — crop the grid to the bounding box, flip so **row 0 is `lat_min`**,
   then Gaussian-smooth with `sigma = (h_min_feature / pixel_res) / 2` to remove features
   smaller than `h_min_feature`.
3. `coastline.py` — read the coastline vectors, clip to the box and simplify (Douglas–Peucker,
   `coastline_tolerance` in metres). Falls back to the DEM zero-crossing if the file is absent.
   The bounding box is passed to `gpd.read_file(bbox=...)` so OGR's spatial index does the
   filtering — reading a global full-resolution shapefile costs seconds, not minutes.
4. `mesh_builder.py` — build the two surfaces (below) and write GTS/VTK.
5. `geology.py` — write a PREM velocity profile to `material_output_file`.

## Configuration keys

Parsed by `config.py`; unknown keys are ignored, `#` and `;` start a comment.

| Key | Meaning |
|---|---|
| `lat_min` `lat_max` `lon_min` `lon_max` | bounding box, decimal degrees |
| `h_min_feature` | smallest terrain feature kept, in metres (Gaussian smoothing of the DEM raster). Below one pixel (~464 m for ETOPO 15") it is a no-op; `0` disables it explicitly |
| `coastline_tolerance` | Douglas–Peucker tolerance on the coastline *polyline*, in metres. Independent of `h_min_feature` — different data, different filter. `0` keeps every vertex |
| `z_min` | `minwater` of `mainSRTM.m`: DEM values ≥ 0 that fall in water are pushed to this depth. Only has an effect if negative. |
| `nc_file` `gebco_url` | elevation grid and where to download it |
| `coastline_file` `coastline_url` | coastline vectors — see *Coastline source* below. A `.zip` URL is downloaded and unpacked automatically; `coastline_file` names a file inside it |
| `include_coastline` | insert the coastline into the triangulation and build the wall |
| `coastline_vtk` | debug output: the clipped input coastlines **before** simplification, as VTK polylines at z = 0 (empty ⇒ not written) |
| `geology_model` `material_output_file` | `PREM` velocity profile output |
| `bathymetry_gts/stl/vtk`, `topography_gts/stl/vtk` | output file names (empty ⇒ not written) |

## Coastline source

The bundled configs use **GSHHG full resolution** (Global Self-consistent, Hierarchical,
High-resolution Geography — Wessel & Smith), level 1, the land/ocean boundary:

```
coastline_file = gshhg/GSHHS_shp/f/GSHHS_f_L1.shp
coastline_url  = https://www.soest.hawaii.edu/pwessel/gshhg/gshhg-shp-2.3.7.zip
```

The archive is 142 MB and unpacks to ~400 MB in `Preproc/gshhg/` (gitignored); the download
happens once, on the first run that needs it. Levels 2–4 in the same directory are lakes,
islands in lakes and ponds in islands — not used yet, but they are what `include_lakes` /
`include_rivers` would need.

This replaced Natural Earth `ne_10m_coastline`, which is 1:10 000 000 cartography — generalised
for world maps, with errors of hundreds of metres to kilometres on indented coasts. Measured on
Kefalonia small1, same bounding box:

| | polylines | raw vertices | after 100 m simplification |
|---|---|---|---|
| Natural Earth 10m | 2 | 153 | 142 |
| GSHHG full | 35 | 2668 | 768 |

Neither is what the Matlab pipeline used: `mainSRTM.m` reads **SWBD** (SRTM Water Body Data,
3 arc-sec, 60°S–60°N) through `swbd_shore.m`, which is also where its Ocean / Land / River /
Lake / Isle classes come from. GSHHG is global, has no latitude limit, and is comparable in
resolution (~100 m).

Anything GDAL can read works — shapefile, GeoJSON, GPKG — so switching source is a config
change, not a code change. Other options: OSM coastline extracts (most detailed, quality varies
by region), or EMODnet/national hydrography for a specific area.

### Measured across the bundled cases

All 12 configs, GSHHG full resolution, run through `./run_all.sh`:

| case | polylines | raw pts | after simplification | wall triangles |
|---|---|---|---|---|
| belle_ile | 41 | 3792 | 1155 | 2802 |
| cadarache | 1 | 5 | 5 | 0 |
| hawaii | 69 | 11438 | 1448 | 9018 |
| hyeres | 18 | 1441 | 1441 | 2250 |
| kashiwazaki | 1 | 495 | 58 | 568 |
| kefalonia | 102 | 6081 | 1811 | 4876 |
| kefalonia_small1 | 35 | 2668 | 768 | 2458 |
| kefalonia_small2 | 35 | 2668 | 768 | 2416 |
| mauna_loa | 24 | 6190 | 716 | 4826 |
| mauna_loa_small | 7 | 761 | 206 | 950 |
| preproc_catalog_example | 69 | 11438 | 1448 | 9018 |
| test_preproc | 727 | 46125 | 12350 | 294 |

Reading the odd rows:

- **cadarache** is inland: one 5-point fragment, no land/water interface, 0 wall triangles.
- **hyeres** has `coastline_tolerance = 1.0`, below one vertex spacing, so nothing is removed —
  see the note on disabling smoothing under *Configuration keys*.
- **kefalonia_small1 and small2** report identical coastline counts although their boxes differ:
  the same 35 polygons lie entirely inside both, so the clip is a no-op for both. Their wall
  counts differ (2458 vs 2416) because the DEM extent does.
- **test_preproc** covers 10°×10° of the Atlantic but has no NetCDF file, so `elevation.py`
  falls back to a synthetic plane. The 46125 coastline points are real, the terrain is not,
  hence only 294 wall triangles. It is a plumbing test, not a physical case.

## Coordinate projection

Spherical equirectangular, R = 6371 km, `nm = π·R/180 = 111194.9 m/degree`, centred on the
bounding box (`mesh_builder.project_to_meters`):

```
x = (lon - lon_c) · π/180 · R · cos(lat_c)
y = (lat - lat_c) · π/180 · R
```

**This is not identical to the Matlab `lonlat2m.m`.** That function writes
`cos(ym/(180*pi))` where it means `cos(ym*pi/180)` — it divides by 180π instead of
multiplying, so the cosine is ≈ 1 and the meridian convergence is essentially not applied:

| | latitude | Matlab m/degree in x | Python m/degree in x | ratio |
|---|---|---|---|---|
| Mauna Loa | 19.75° | 111127 | 104654 | 1.062 |
| Kefalonia | 38.27° | 110940 | 87303 | 1.271 |
| Kashiwazaki | 37.4° | 110952 | 88335 | 1.256 |

Matlab meshes therefore come out stretched in longitude — 6 % at Hawaii, 27 % at Kefalonia.
The Python projection is the physically correct one. Matlab also shifts the origin to the SW
corner (`x - min(x)`) while Python centres on the box; that is irrelevant to HexMesh, which
takes the bounding box from the GTS surface itself (`intercept_surface.cpp:321`). What does
matter is that topography and bathymetry share one origin — they do, both derive
`center_lat/lon` from the same config bounds.

## Bathymetry surface — why it looks like that

`Apply_material` (`../src/aplly_material.cpp:137`) classifies each hexahedron by shooting a
ray **straight up** from the centroid of its top face to `2·z2` of the bathymetry bounding
box: if the ray hits the surface the element is solid (`n_mat = 0`), otherwise it is fluid.
The surface is built to make that test give the right answer everywhere:

- **bathy1 — water sheet.** The sea floor at its DEM depth. An element below the floor is hit
  from above ⇒ solid; one in the water column is not ⇒ fluid.
- **bathy2 — land plateau.** A *flat lid* at `z_top = 2·max(z)`, on a duplicated set of nodes,
  covering the island. Every element of a land column, including the topmost, hits it ⇒ the
  whole column is solid. This is what makes the raytrace see "solid ground" over land.
- **wall.** The nodes on the water/land seam sit at `z = 0` and are duplicated up to `z_top`;
  two triangles per seam edge close the two sheets with a vertical curtain. The wall is what
  `intercept_surface` cuts the hexahedra against, so it is what produces a sharp, vertical
  coastline in the hex mesh instead of a slanted DEM ramp.

Implementation: a single Delaunay over DEM grid points **plus** the inserted coastline points;
seam edges are the edges shared by a land and a water triangle (`_seam_edges`), which makes the
wall watertight by construction. Nodes on the seam get `z = 0` on the water copy and `z = z_top`
on the land copy.

### Land or water: point-in-coastline

Each triangle is classified at its centroid by **`shapely.contains_xy` against the coastline
polygons** (`_is_land`) — the criterion `mainSRTM.m` used through `distanceCoastline`, where
`dist > 0` is land. The polygons come from the same simplified geometry that is inserted into
the triangulation, so the seam and the test cannot disagree.

The sign of the DEM is only the fallback, used when the coastline source has no polygons (open
polylines such as Natural Earth) or when there is no coastline file at all. It is a proxy and it
is wrong in two ways that matter:

- a crater, polder or lagoon **below sea level but inside the coastline** is called water, and
  gets a sea floor plus a wall around it;
- a **shallow bank above sea level but offshore** is called land, gets a plateau at `z_top`, and
  `Apply_material` then marks that whole water column solid.

Both are what `z_min` is for: a DEM node reading ≥ 0 that lies outside the coastline is pushed
to `z_min` (`mainSRTM.m`'s `ind = dist<0 & z>=0; z(ind) = minwater`). Under the old DEM-sign
test that rule could only ever touch seam nodes, which are forced to 0 anyway — it did nothing.
The count of clamped nodes is printed per run; it is 189–1755 across the bundled cases.

The seam still follows the coastline closely even though the Delaunay is unconstrained, because
the inserted coastline vertices are denser than the DEM. Measured maximum deviation of a seam
node from the polygon, against a 464–469 m pixel:

| case | mean | median | p95 | max | above 1 pixel |
|---|---|---|---|---|---|
| kefalonia_small1 | 48 m | 0 m | 208 m | 343 m | 0 % |
| mauna_loa | 93 m | 72 m | 273 m | 342 m | 0 % |
| belle_ile | 40 m | 0 m | 194 m | 344 m | 0 % |
| hyeres | 6 m | 0 m | 34 m | 278 m | 0 % |
| kashiwazaki | 9 m | 0 m | 73 m | 307 m | 0 % |

If a coarser `coastline_tolerance` ever pushes this past a pixel, the cheap fix is
`shapely.segmentize` on the simplified coastline before insertion; the real one is a constrained
Delaunay (`triangle`, or CGAL). Neither is needed today.

Differences remaining from `mainSRTM.m`: it used `1.0·max(z)` for the wall top and `1.5·max(z)`
for the plateau, leaving a gap between them; here both are `2·max(z)`.

## Topography surface

Following `mainSRTM.m`, with the same point-in-coastline test:

- everything **outside** the coastline is flattened to `z = 0`;
- land within `coastline_tolerance` of the shore with `z < 0` is flattened to 0 — that band is
  where the simplified coastline and the DEM disagree, and a pit right at the shore is an
  artefact, not terrain (the band is applied by eroding the land polygon, which is steadier
  than a distance transform);
- land beyond that band **keeps its elevation, negative included** — a closed depression below
  sea level is real terrain, not sea.

Coastline points are merged into the surface at `z = 0` through the same Delaunay. No wall — the
topography is a single-valued height field. Without polygons it falls back to `z < 0 ⇒ 0`
everywhere, the previous behaviour.

## Elevation data sources

The bundled configs use ETOPO 2022 (15 arc-sec) served by the NOAA ERDDAP, which is a combined
topo+bathy grid. Better products exist, per category.

### Topography

| Source | Resolution | Coverage | Note |
|---|---|---|---|
| National LiDAR (IGN RGE ALTI, USGS 3DEP, GSI Japan) | 1–5 m | national | the real top of the range where it exists |
| ArcticDEM / REMA (PGC) | 2 m | polar | |
| ALOS AW3D | 5 m (paid) / 30 m (free) | global | |
| **FABDEM** | 30 m | global | Copernicus with forests and buildings removed — bare earth, best global choice for wave modelling |
| Copernicus DEM GLO-30 | 30 m | global | DSM (canopy/roof top), the basis of FABDEM |
| NASADEM | 30 m | ±60° | reprocessed SRTM; replaces the SRTM used by `mainSRTM.m` |
| ASTER GDEM v3 | 30 m | global | noisy, avoid |

### Bathymetry

| Source | Resolution | Coverage | Note |
|---|---|---|---|
| NOAA CUDEM | 1/9 arc-sec (~3 m) | US coast | already-merged topobathy |
| SHOM Litto3D | 1 m | French coast | bathymetric LiDAR |
| **EMODnet Bathymetry** | 1/16 arc-min (~115 m) | European seas | best option for Kefalonia, Hyères, Belle-Île |
| National multibeam (NCEI, JAMSTEC/JODC) | 10–100 m | surveyed areas | |
| GEBCO 2024 | 15 arc-sec (~450 m) | global | the global standard |
| SRTM15+ V2.6 | 15 arc-sec | global | satellite altimetry (Sandwell/Smith) — what `bathymetrySRTM.m` fetched from topex.ucsd.edu |
| ETOPO 2022 | 15 arc-sec | global | combined topo+bathy; what the bundled configs use |

Suggested pairings: Kefalonia / Hyères / Belle-Île → EMODnet + FABDEM; Japan → JODC/JAMSTEC +
GSI 5 m; Hawaii → CUDEM + Copernicus GLO-30.

### Using two different sources

Not supported yet: `preproc.py` extracts **one** grid and passes it to both mesh builders. It
would take a second config key, a second `extract_elevation_region` call, a resample of the
coarser grid onto the finer one, and `elev = np.where(topo_grid > 0, topo_grid, bathy_grid)`.

Three things to get right when doing it:

1. **Vertical datum.** GEBCO/ETOPO are referenced to mean sea level, EMODnet to LAT (lowest
   astronomical tide), CUDEM to NAVD88. Without conversion a 1–3 m step appears exactly at the
   shoreline. Negligible at 500 m mesh resolution, but real.
2. **The two shorelines disagree.** The `z = 0` contour of the topography grid and of the
   bathymetry grid are offset by tens to hundreds of metres. Since `build_bathymetry_mesh`
   classifies land/water by the sign of the DEM, the classification must come from a single
   criterion. The `np.where(topo_grid > 0, ...)` merge above gives that for free: the sign is
   always the topography's.
3. **Format.** `extract_elevation_region` reads a single NetCDF or GeoTIFF file. Copernicus and
   FABDEM ship 1°×1° GeoTIFF tiles and would need mosaicking first (`gdalbuildvrt` or
   `rioxarray.merge`). EMODnet serves NetCDF over WCS, which loads directly.

## Tests

```sh
python3 test_mesh_builder.py
```

Builds a synthetic conical island and asserts: the plateau exists and is flat at `2·max(z)`,
every `z = 0` seam node has a duplicate at `z_top`, vertical wall triangles spanning the full
range exist, the sea floor is not flattened, an upward ray at the island centre hits the
plateau, and one offshore does not hit anything above sea level.

Then it adds a crater reaching −300 m *inside* the coastline and a bank reaching +80 m *outside*
it — the two cases where the sign of the DEM and the polygon disagree — and asserts the crater
is plateau, the bank is water clamped to `z_min`, the seam stays within one DEM pixel of the
polygon, the topography keeps the crater below sea level, and the no-polygon fallback still
behaves the old way.
