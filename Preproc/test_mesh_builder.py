#!/usr/bin/env python3
"""Self-check for the bathymetry wall (mainSRTM.m bathy1+bathy2 spec). Run: python3 test_mesh_builder.py"""
import numpy as np
from mesh_builder import build_bathymetry_mesh

# synthetic island: cone peaking at +200 m in the middle of a -1000 m ocean
n = 60
lat = np.linspace(0.0, 1.0, n)[:, None]
lon = np.linspace(0.0, 1.0, n)[None, :]
r = np.sqrt((lat - 0.5) ** 2 + (lon - 0.5) ** 2)
elev = np.where(r < 0.25, 200.0 * (1 - r / 0.25), -1000.0 * (r - 0.25) / 0.25 - 1.0)

verts, tris = build_bathymetry_mesh(elev, 0.0, 1.0, 0.0, 1.0, [], z_min_bathymetry=-100.0)
z_top = 2.0 * elev.max()

# 1. the land plateau exists and is flat at 2*max(z)
plateau = np.isclose(verts[:, 2], z_top)
assert plateau.sum() > 100, f"no land plateau at z_top={z_top}: {plateau.sum()} nodes"

# 2. the seam nodes at z=0 are duplicated up to z_top (same x,y, both z present)
seam = np.isclose(verts[:, 2], 0.0)
assert seam.sum() > 20, f"no coastline nodes at z=0: {seam.sum()}"
xy_seam = {tuple(np.round(p, 6)) for p in verts[seam][:, :2]}
xy_top = {tuple(np.round(p, 6)) for p in verts[plateau][:, :2]}
assert xy_seam <= xy_top, "coastline nodes at z=0 were not duplicated to z_top"

# 3. vertical wall triangles exist (span the full 0..z_top range)
tz = verts[tris][:, :, 2]
wall = (tz.min(axis=1) < 1e-9) & (tz.max(axis=1) > z_top - 1e-9)
assert wall.sum() > 20, f"no vertical wall triangles: {wall.sum()}"

# 4. nothing above the plateau, sea floor kept negative
assert verts[:, 2].max() <= z_top + 1e-9
assert verts[:, 2].min() < -500.0, "sea floor was flattened"

# 5. an upward ray at the island centre hits the plateau (Apply_material -> solid),
#    a ray in the open ocean above the sea floor does not.
def hit_z(x, y):
    """z of the triangles whose 2D projection contains (x, y)."""
    p = verts[tris][:, :, :2]
    d = lambda a, b, c: (b[:, 0] - a[:, 0]) * (c[1] - a[:, 1]) - (b[:, 1] - a[:, 1]) * (c[0] - a[:, 0])
    s1, s2, s3 = d(p[:, 0], p[:, 1], (x, y)), d(p[:, 1], p[:, 2], (x, y)), d(p[:, 2], p[:, 0], (x, y))
    inside = ((s1 >= 0) & (s2 >= 0) & (s3 >= 0)) | ((s1 <= 0) & (s2 <= 0) & (s3 <= 0))
    return verts[tris][inside][:, :, 2]

centre = hit_z(0.0, 0.0)
assert centre.size and np.isclose(centre.max(), z_top), "island centre is not covered by the plateau"
offshore = hit_z(verts[:, 0].max() * 0.9, 0.0)
assert offshore.size and offshore.max() < 0.0, "open ocean is covered by a lid above sea level"

print("OK: plateau at 2*max(z), coastline nodes duplicated z=0 -> z_top, vertical wall present.")

# ---------------------------------------------------------------------------
# point-in-coastline classification: built so the DEM sign and the polygon disagree
# ---------------------------------------------------------------------------
import shapely
from mesh_builder import build_topography_mesh, meters_to_lonlat

# same island, plus two features the sign of the DEM gets wrong:
#   - a crater inside the island reaching -300 m  -> is land, must be plateau, no wall around it
#   - a shallow bank offshore reaching +80 m      -> is water, must be clamped, no plateau
crater = np.sqrt((lat - 0.42) ** 2 + (lon - 0.42) ** 2) < 0.06
bank = np.sqrt((lat - 0.80) ** 2 + (lon - 0.80) ** 2) < 0.05
elev2 = np.where(crater, -300.0, elev)
elev2 = np.where(bank, 80.0, elev2)
island = shapely.Point(0.5, 0.5).buffer(0.25)  # (lon, lat) == the r < 0.25 cone footprint

v2, t2 = build_bathymetry_mesh(elev2, 0.0, 1.0, 0.0, 1.0, [], z_min_bathymetry=-100.0, land=island)
z_top2 = 2.0 * elev2.max()

def z_over(lon_, lat_, verts, tris):
    """z of the triangles covering (lon, lat), given in degrees."""
    x = (lon_ - 0.5) * (np.pi / 180) * 6371000.0 * np.cos(np.radians(0.5))
    y = (lat_ - 0.5) * (np.pi / 180) * 6371000.0
    p = verts[tris][:, :, :2]
    d = lambda a, b, c: (b[:, 0] - a[:, 0]) * (c[1] - a[:, 1]) - (b[:, 1] - a[:, 1]) * (c[0] - a[:, 0])
    s1, s2, s3 = d(p[:, 0], p[:, 1], (x, y)), d(p[:, 1], p[:, 2], (x, y)), d(p[:, 2], p[:, 0], (x, y))
    inside = ((s1 >= 0) & (s2 >= 0) & (s3 >= 0)) | ((s1 <= 0) & (s2 <= 0) & (s3 <= 0))
    return verts[tris][inside][:, :, 2]

# the crater is inside the coastline -> plateau, despite the DEM reading -300 m
zc = z_over(0.42, 0.42, v2, t2)
assert zc.size and np.isclose(zc.max(), z_top2), "crater below sea level was classified as water"
# the bank is outside the coastline -> water, despite the DEM reading +80 m
zb = z_over(0.80, 0.80, v2, t2)
assert zb.size and zb.max() < 0.0, "shallow bank above sea level was classified as land"
assert np.isclose(zb.max(), -100.0, atol=1e-6), "bank was not clamped to minwater"

# the wall follows the polygon: every seam node within one DEM pixel of the coastline
seam2 = v2[np.isclose(v2[:, 2], 0.0)]
lon_s, lat_s = meters_to_lonlat(seam2[:, 0], seam2[:, 1], 0.5, 0.5)
d_m = shapely.distance(shapely.points(lon_s, lat_s), island.boundary) * 111194.9
pixel_m = (1.0 / (n - 1)) * 111194.9
assert seam2.size and d_m.max() < pixel_m, f"seam strays {d_m.max():.0f} m from the coastline (pixel {pixel_m:.0f} m)"

# topography: land below sea level survives away from the shore, water is flat at 0
tv, tt = build_topography_mesh(elev2, 0.0, 1.0, 0.0, 1.0, [], land=island, coast_band_m=1000.0)
lon_t, lat_t = meters_to_lonlat(tv[:, 0], tv[:, 1], 0.5, 0.5)
on_land = shapely.contains_xy(island, lon_t, lat_t)
assert tv[on_land, 2].min() < -200.0, "the crater was flattened although it is inland"
assert np.allclose(tv[~on_land, 2], 0.0), "water is not flat at z = 0 in the topography"

# and the fallback (no polygon) must still behave as before
v3, t3 = build_bathymetry_mesh(elev2, 0.0, 1.0, 0.0, 1.0, [], z_min_bathymetry=-100.0, land=None)
zc3 = z_over(0.42, 0.42, v3, t3)
assert zc3.max() < 0.0, "without a polygon the crater should fall back to being water"

print("OK: point-in-coastline classifies the inland crater as land and the offshore bank as water;")
print("    seam within one DEM pixel of the polygon; below-sea-level land kept in the topography;")
print("    DEM-sign fallback unchanged.")
