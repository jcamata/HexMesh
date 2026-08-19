import numpy as np
import os

EARTH_RADIUS_M = 6371000.0

def project_to_meters(lat, lon, z, center_lat, center_lon):
    """
    Projects (lat, lon, z) to local metric grid in meters
    """
    x = (lon - center_lon) * (np.pi / 180.0) * EARTH_RADIUS_M * np.cos(np.radians(center_lat))
    y = (lat - center_lat) * (np.pi / 180.0) * EARTH_RADIUS_M
    return x, y, z

def meters_to_lonlat(x, y, center_lat, center_lon):
    """Inverse of project_to_meters (the projection is affine, so this is exact)."""
    lat = y / ((np.pi / 180.0) * EARTH_RADIUS_M) + center_lat
    lon = x / ((np.pi / 180.0) * EARTH_RADIUS_M * np.cos(np.radians(center_lat))) + center_lon
    return lon, lat

def _is_land(x, y, land, elevations, center_lat, center_lon, x0, y0, dx, dy, width, height):
    """
    Land/water for a set of points given in metres.

    With a polygon coastline this is the inside/outside test -- the criterion mainSRTM.m used
    through distanceCoastline (dist > 0 is land). Without one (open polylines, or the DEM
    zero-crossing fallback) it uses the sign of the DEM, which is only a proxy: it calls land
    below sea level water, and shallow banks above sea level land.
    """
    if land is not None:
        import shapely
        lon, lat = meters_to_lonlat(x, y, center_lat, center_lon)
        return shapely.contains_xy(land, lon, lat)

    ci = np.clip(np.rint((x - x0) / dx).astype(int), 0, width - 1)
    ri = np.clip(np.rint((y - y0) / dy).astype(int), 0, height - 1)
    return elevations[ri, ci] > 0

def _build_grid_triangles(height: int, width: int):
    triangles = []
    for r in range(height - 1):
        for c in range(width - 1):
            tl = r * width + c
            tr = tl + 1
            bl = (r + 1) * width + c
            br = bl + 1
            triangles.append([tl, bl, tr])
            triangles.append([tr, bl, br])
    return np.array(triangles, dtype=np.int32)

def _seam_edges(water_tris: np.ndarray, land_tris: np.ndarray):
    """Edges shared by a water and a land triangle: the coastline of the triangulation."""
    def edge_set(T):
        e = np.concatenate([T[:, [0, 1]], T[:, [1, 2]], T[:, [2, 0]]])
        return set(map(tuple, np.sort(e, axis=1)))
    return edge_set(water_tris) & edge_set(land_tris)

def write_coastline_vtk(path: str, polylines: list, lat_min: float, lat_max: float, lon_min: float, lon_max: float):
    """
    Debug output: the coastlines alone, as VTK polylines at z = 0, in the same metric
    frame as the surface meshes. 'polylines' is a list of [(lat, lon), ...] in degrees.
    """
    center_lat = (lat_min + lat_max) / 2.0
    center_lon = (lon_min + lon_max) / 2.0

    cells = []
    pts = []
    for line in polylines:
        ids = []
        for lat, lon in line:
            px, py, _ = project_to_meters(lat, lon, 0.0, center_lat, center_lon)
            ids.append(len(pts))
            pts.append((px, py))
        if len(ids) >= 2:
            cells.append(ids)
        else:
            del pts[len(pts) - len(ids):]

    size = sum(len(c) + 1 for c in cells)
    print(f"Writing coastline debug VTK: '{path}' ({len(cells)} polylines, {len(pts)} points)")
    with open(path, 'w', encoding='utf-8') as f:
        f.write("# vtk DataFile Version 3.0\nRaw input coastlines\nASCII\nDATASET POLYDATA\n")
        f.write(f"POINTS {len(pts)} double\n")
        for px, py in pts:
            f.write(f"{px:.4f} {py:.4f} 0.0000\n")
        f.write(f"\nLINES {len(cells)} {size}\n")
        for c in cells:
            f.write(f"{len(c)} " + " ".join(str(i) for i in c) + "\n")
        f.write(f"\nCELL_DATA {len(cells)}\nSCALARS polyline_id int 1\nLOOKUP_TABLE default\n")
        for i in range(len(cells)):
            f.write(f"{i}\n")

def _triangulate_2d(pts: np.ndarray, segments: list = None, use_cdt: bool = True):
    """
    Deduplicates points and cleans segments, then triangulates 2D points.
    If use_cdt=True and segments are provided, uses 'triangle' library for CDT.
    Returns (cleaned_pts, triangles, inv_map).
    """
    if pts is None or len(pts) == 0:
        return pts, np.zeros((0, 3), dtype=np.int32), np.array([], dtype=np.int32)

    # 1. Round to 1 mm (1e-3 m) to merge coincident/near-coincident vertices and prevent CDT segfaults
    rounded_pts = np.round(pts, decimals=3)
    cleaned_pts, inv_map = np.unique(rounded_pts, axis=0, return_inverse=True)

    cleaned_segs = None
    if segments is not None and len(segments) > 0:
        segs_arr = np.array(segments, dtype=np.int32)
        mapped_segs = inv_map[segs_arr]
        # remove self loops (start == end)
        valid = mapped_segs[:, 0] != mapped_segs[:, 1]
        mapped_segs = mapped_segs[valid]
        if len(mapped_segs) > 0:
            sorted_segs = np.sort(mapped_segs, axis=1)
            cleaned_segs = np.unique(sorted_segs, axis=0)

    if use_cdt and cleaned_segs is not None and len(cleaned_segs) > 0:
        try:
            import triangle as tr
            data = {
                'vertices': np.ascontiguousarray(cleaned_pts, dtype=np.float64),
                'segments': np.ascontiguousarray(cleaned_segs, dtype=np.int32)
            }
            res = tr.triangulate(data, 'pc')
            print(f"  CDT: Constrained Delaunay Triangulation performed with {len(cleaned_segs)} segments "
                  f"via 'triangle' package ({len(res['triangles'])} triangles).")
            return cleaned_pts, res['triangles'].astype(np.int32), inv_map
        except Exception as e:
            print(f"  Warning: CDT via 'triangle' failed ({e}), falling back to scipy.spatial.Delaunay.")

    from scipy.spatial import Delaunay
    tris = Delaunay(cleaned_pts).simplices.astype(np.int32)
    return cleaned_pts, tris, inv_map

def build_bathymetry_mesh(elevations: np.ndarray, lat_min: float, lat_max: float, lon_min: float, lon_max: float, coastlines: list, z_min_bathymetry: float = None, land=None, use_cdt: bool = True):
    """
    Replicates the mainSRTM.m bathymetry (bathy1 + bathy2 + vertical coastline wall).

    One Delaunay (or CDT) over DEM grid points + inserted coastline points; each triangle is
    classified land/water by the DEM sign at its centroid, then:
      - bathy1 (water): sea floor at the DEM depth; nodes on the coastline seam forced to z = 0
      - bathy2 (land) : flat plateau at z_top = 2*max(z), on a duplicated set of nodes, so the
                        upward ray of Apply_material always hits it -> the whole column is solid
      - wall          : the seam nodes at z = 0 duplicated up to z_top, closing the two sheets
                        with a vertical curtain -> that is what cuts the hexahedra at the coast.
    """
    height, width = elevations.shape
    center_lat = (lat_min + lat_max) / 2.0
    center_lon = (lon_min + lon_max) / 2.0

    max_z = float(np.max(elevations))
    if max_z <= 0:
        max_z = 100.0
    z_top = 2.0 * max_z

    # mainSRTM.m 'minwater': positive DEM values that fall in water are pushed down
    minwater = float(z_min_bathymetry) if (z_min_bathymetry is not None and z_min_bathymetry < 0) else 0.0

    print(f"Bathymetry Mesh (mainSRTM.m spec): Max Z = {max_z:.1f} m, land plateau / wall top Z_top = {z_top:.1f} m")

    lats = np.linspace(lat_min, lat_max, height)
    lons = np.linspace(lon_min, lon_max, width)
    lon_grid, lat_grid = np.meshgrid(lons, lats)

    x_grid = (lon_grid - center_lon) * (np.pi / 180.0) * EARTH_RADIUS_M * np.cos(np.radians(center_lat))
    y_grid = (lat_grid - center_lat) * (np.pi / 180.0) * EARTH_RADIUS_M

    dem_pts = np.column_stack([x_grid.ravel(), y_grid.ravel()])
    n_dem = len(dem_pts)

    coast_pts = []
    segments = []
    for coastline in coastlines or []:
        line_indices = []
        for lat, lon in coastline:
            if lat_min <= lat <= lat_max and lon_min <= lon <= lon_max:
                cx, cy, _ = project_to_meters(lat, lon, 0.0, center_lat, center_lon)
                idx = n_dem + len(coast_pts)
                coast_pts.append([cx, cy])
                line_indices.append(idx)
        for i in range(len(line_indices) - 1):
            segments.append([line_indices[i], line_indices[i + 1]])

    try:
        raw_pts = np.vstack([dem_pts, np.array(coast_pts)]) if coast_pts else dem_pts
    except Exception:
        raw_pts = dem_pts

    print(f"Inserting {len(coast_pts)} coastline points into the bathymetry triangulation...")
    pts, tris, inv_map = _triangulate_2d(raw_pts, segments, use_cdt=use_cdt)

    # land / water per triangle, at its centroid
    cent = pts[tris].mean(axis=1)
    dx = x_grid[0, 1] - x_grid[0, 0] if width > 1 else 1.0
    dy = y_grid[1, 0] - y_grid[0, 0] if height > 1 else 1.0
    is_land = _is_land(cent[:, 0], cent[:, 1], land, elevations, center_lat, center_lon,
                       x_grid[0, 0], y_grid[0, 0], dx, dy, width, height)
    print(f"Classification: {'point-in-coastline' if land is not None else 'sign of the DEM'} "
          f"-- {is_land.sum()} land / {(~is_land).sum()} water triangles.")

    land_tris = tris[is_land]
    water_tris = tris[~is_land]
    if len(land_tris) == 0 or len(water_tris) == 0:
        print("Warning: no land/water interface found - no coastline wall will be generated.")

    # z of the water sheet: sea floor, coastline seam at z = 0.
    # mainSRTM.m: 'ind = dist<0 & z>=0; z(ind) = minwater' -- a DEM cell that reads at or above
    # sea level but lies outside the coastline is a shallow bank, pushed down to minwater.
    z_water = np.zeros(len(pts))
    e = elevations.ravel()
    z_water[inv_map[:n_dem]] = np.where(e >= 0, minwater, e)
    if land is not None:
        n_bank = int(np.sum((e >= 0) & ~_is_land(dem_pts[:, 0], dem_pts[:, 1], land, elevations,
                                                 center_lat, center_lon, x_grid[0, 0], y_grid[0, 0],
                                                 dx, dy, width, height)))
        if n_bank:
            print(f"  {n_bank} DEM nodes at/above sea level but outside the coastline "
                  f"(shallow banks) clamped to {minwater:.1f} m.")

    seam = _seam_edges(water_tris, land_tris)
    seam_nodes = np.array(sorted({n for edge in seam for n in edge}), dtype=np.int32)
    z_water[seam_nodes] = 0.0

    # duplicate nodes: one copy for the water sheet, one for the land plateau at z_top
    w_used = np.unique(water_tris)
    l_used = np.unique(land_tris)
    w_map = np.full(len(pts), -1, dtype=np.int32)
    l_map = np.full(len(pts), -1, dtype=np.int32)
    w_map[w_used] = np.arange(len(w_used))
    l_map[l_used] = np.arange(len(l_used)) + len(w_used)

    verts = np.vstack([
        np.column_stack([pts[w_used], z_water[w_used]]),
        np.column_stack([pts[l_used], np.full(len(l_used), z_top)]),
    ])

    # vertical wall: seam node at z = 0 -> its duplicate at z_top
    wall = [[w_map[a], w_map[b], l_map[b]] for a, b in seam] + \
           [[w_map[a], l_map[b], l_map[a]] for a, b in seam]

    tri_all = np.vstack([w_map[water_tris], l_map[land_tris]] + ([np.array(wall, dtype=np.int32)] if wall else []))

    print(f"Generated Bathymetry Mesh: {len(verts)} vertices, {len(tri_all)} triangles "
          f"({len(water_tris)} water, {len(land_tris)} land plateau, {len(wall)} wall).")
    return verts, tri_all.astype(np.int32)

def build_topography_mesh(elevations: np.ndarray, lat_min: float, lat_max: float, lon_min: float, lon_max: float, coastlines: list = None, land=None, coast_band_m: float = 0.0, use_cdt: bool = True):
    """
    Generates Topography surface mesh, following mainSRTM.m:
    - everything outside the coastline (water) is flattened to z = 0
    - land within `coast_band_m` of the coast with z < 0 is flattened to z = 0 (that band is
      where the simplified coastline and the DEM disagree, and a pit right at the shore is an
      artefact, not terrain)
    - land beyond that band keeps its elevation, negative included -- polders and closed
      depressions are real terrain below sea level
    - Coastline points (z = 0) are merged into the surface via 2D Delaunay triangulation

    Without a polygon coastline (land is None) it falls back to 'z < 0 becomes 0' everywhere.
    """
    height, width = elevations.shape
    center_lat = (lat_min + lat_max) / 2.0
    center_lon = (lon_min + lon_max) / 2.0

    lats = np.linspace(lat_min, lat_max, height)
    lons = np.linspace(lon_min, lon_max, width)
    lon_grid, lat_grid = np.meshgrid(lons, lats)

    x_grid = (lon_grid - center_lon) * (np.pi / 180.0) * EARTH_RADIUS_M * np.cos(np.radians(center_lat))
    y_grid = (lat_grid - center_lat) * (np.pi / 180.0) * EARTH_RADIUS_M

    if land is None:
        z_final = np.where(elevations < 0, 0.0, elevations)
    else:
        import shapely
        on_land = shapely.contains_xy(land, lon_grid, lat_grid)
        z_final = np.where(on_land, elevations, 0.0)

        # land within coast_band_m of the shore: eroding the polygon is cheaper and steadier
        # than a distance transform, and empty erosion just means "the whole island is coastal"
        if coast_band_m > 0:
            inner = land.buffer(-coast_band_m / 111194.9)
            near_coast = on_land & ~shapely.contains_xy(inner, lon_grid, lat_grid) if not inner.is_empty else on_land
            z_final = np.where(near_coast & (z_final < 0), 0.0, z_final)

        n_deep = int(np.sum(on_land & (z_final < 0)))
        print(f"Topography: {on_land.sum()} land / {(~on_land).sum()} water DEM nodes"
              + (f"; {n_deep} land nodes kept below sea level." if n_deep else "."))

    dem_pts_2d = np.column_stack([x_grid.ravel(), y_grid.ravel()])
    n_dem = len(dem_pts_2d)

    coast_pts_2d = []
    segments = []
    if coastlines:
        for coastline in coastlines:
            line_indices = []
            for lat, lon in coastline:
                if lat_min <= lat <= lat_max and lon_min <= lon <= lon_max:
                    cx, cy, _ = project_to_meters(lat, lon, 0.0, center_lat, center_lon)
                    idx = n_dem + len(coast_pts_2d)
                    coast_pts_2d.append([cx, cy])
                    line_indices.append(idx)
            for i in range(len(line_indices) - 1):
                segments.append([line_indices[i], line_indices[i + 1]])

    if len(coast_pts_2d) > 0:
        print(f"Merging {len(coast_pts_2d)} coastline points into DEM topography surface mesh via Delaunay triangulation...")
        all_pts_2d = np.vstack([dem_pts_2d, np.array(coast_pts_2d)])
        cleaned_pts_2d, triangles, inv_map = _triangulate_2d(all_pts_2d, segments, use_cdt=use_cdt)
        
        pts_3d_z = np.zeros(len(cleaned_pts_2d))
        pts_3d_z[inv_map[:n_dem]] = z_final.ravel()
        pts_3d_z[inv_map[n_dem:]] = 0.0
        all_pts_3d = np.column_stack([cleaned_pts_2d, pts_3d_z])
    else:
        all_pts_3d = np.column_stack([dem_pts_2d, z_final.ravel()])
        triangles = _build_grid_triangles(height, width)

    print(f"Generated Topography Mesh: {len(all_pts_3d)} vertices, {len(triangles)} triangles.")
    return all_pts_3d, np.array(triangles, dtype=np.int32)

def export_mesh_files(vertices: np.ndarray, triangles: np.ndarray, stl_file: str, gts_file: str, vtk_file: str):
    """
    Exports mesh to STL, GTS, and VTK formats
    """
    # 1. Export GTS (GNU Triangulated Surface)
    if gts_file:
        print(f"Writing GTS file: '{gts_file}'")
        num_edges = len(triangles) * 3
        with open(gts_file, 'w', encoding='utf-8') as f:
            f.write(f"{len(vertices)} {num_edges} {len(triangles)}\n")
            for v in vertices:
                f.write(f"{v[0]:.4f} {v[1]:.4f} {v[2]:.4f}\n")

            for i, t in enumerate(triangles):
                f.write(f"{t[0] + 1} {t[1] + 1}\n")
                f.write(f"{t[1] + 1} {t[2] + 1}\n")
                f.write(f"{t[2] + 1} {t[0] + 1}\n")

            for i in range(len(triangles)):
                base = i * 3 + 1
                f.write(f"{base} {base + 1} {base + 2}\n")

    # 2. Export STL (Binary STL Format)
    if stl_file:
        print(f"Writing STL file: '{stl_file}'")
        try:
            import trimesh
            mesh = trimesh.Trimesh(vertices=vertices, faces=triangles)
            mesh.export(stl_file)
        except Exception:
            import struct
            with open(stl_file, 'wb') as f:
                header = b"GEBCO Python Preproc Metric STL".ljust(80, b'\x00')
                f.write(header)
                f.write(struct.pack('<I', len(triangles)))

                for t in triangles:
                    v0, v1, v2 = vertices[t[0]], vertices[t[1]], vertices[t[2]]
                    u = v1 - v0
                    v = v2 - v0
                    n = np.cross(u, v)
                    norm = np.linalg.norm(n)
                    n = n / norm if norm > 0 else np.array([0, 0, 1])

                    f.write(struct.pack('<fff', *n))
                    f.write(struct.pack('<fff', *v0))
                    f.write(struct.pack('<fff', *v1))
                    f.write(struct.pack('<fff', *v2))
                    f.write(struct.pack('<H', 0))

    # 3. Export VTK (ASCII PolyData Surface Format)
    if vtk_file:
        print(f"Writing VTK file: '{vtk_file}'")
        with open(vtk_file, 'w', encoding='utf-8') as f:
            f.write("# vtk DataFile Version 3.0\n")
            f.write("Preproc Surface Mesh\n")
            f.write("ASCII\n")
            f.write("DATASET POLYDATA\n")
            f.write(f"POINTS {len(vertices)} double\n")
            for v in vertices:
                f.write(f"{v[0]:.4f} {v[1]:.4f} {v[2]:.4f}\n")

            f.write(f"\nPOLYGONS {len(triangles)} {len(triangles) * 4}\n")
            for t in triangles:
                f.write(f"3 {t[0]} {t[1]} {t[2]}\n")
