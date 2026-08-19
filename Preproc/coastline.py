import os
import numpy as np

def _as_polylines(geom):
    """Every geometry type the coastline file may hold, as lists of (lat, lon)."""
    if geom.geom_type == 'LineString':
        return [[(lat, lon) for lon, lat in geom.coords]]
    if geom.geom_type == 'MultiLineString':
        return [[(lat, lon) for lon, lat in line.coords] for line in geom.geoms]
    if geom.geom_type in ('Polygon', 'MultiPolygon'):
        return _as_polylines(geom.boundary)
    return []

def load_and_smooth_coastlines(filepath: str, lat_min: float, lat_max: float, lon_min: float, lon_max: float, tolerance_m: float, raw_vtk: str = None, max_segment_m: float = None):
    """
    Loads coastline / hydrography vectors using geopandas & shapely, clips to bounding box,
    and applies Douglas-Peucker polyline simplification (shapely.simplify).

    Returns (polylines, land):
      polylines: list of list of (lat, lon), the simplified coastline, for insertion into
                 the surface triangulations;
      land:      shapely geometry (lon, lat) covering the land, for the inside/outside test
                 that classifies triangles -- or None when the source holds no polygons
                 (Natural Earth ships open lines), in which case callers fall back to the
                 sign of the DEM.

    raw_vtk: optional path for a debug VTK holding the clipped input coastlines *before*
    simplification -- the raw file data, to compare against what the mesh ends up using.
    """
    polylines = []
    raw = []
    land = None

    if os.path.exists(filepath):
        try:
            import geopandas as gpd
            from shapely.geometry import box

            print(f"Loading coastline vector file '{filepath}' via geopandas...")
            # bbox= pushes the filter down to OGR's spatial index: reading a global
            # full-resolution shapefile (GSHHG L1 is 161 MB) costs seconds, not minutes.
            gdf = gpd.read_file(filepath, bbox=(lon_min, lat_min, lon_max, lat_max))

            # Clip to bounding box
            bbox = box(lon_min, lat_min, lon_max, lat_max)
            clipped = gdf.clip(bbox)

            # Convert tolerance in meters to degrees approx
            tol_deg = tolerance_m / 111320.0

            land_parts = []
            for geom in clipped.geometry:
                if geom.is_empty:
                    continue

                # The land test must use the *same* geometry that gets inserted into the
                # triangulation, otherwise the seam and the inside/outside test disagree
                # by up to `tolerance`.
                simplified = geom.simplify(tolerance=tol_deg, preserve_topology=True)
                # Segmentize to prevent long straight segments from causing Delaunay zigzags
                if max_segment_m and max_segment_m > 0 and hasattr(simplified, 'segmentize'):
                    import shapely
                    seg_deg = max_segment_m / 111320.0
                    simplified = shapely.segmentize(simplified, max_segment_length=seg_deg)

                raw.extend(_as_polylines(geom))
                polylines.extend(_as_polylines(simplified))
                if simplified.geom_type in ('Polygon', 'MultiPolygon'):
                    land_parts.append(simplified)

            n_raw = sum(len(l) for l in raw)
            n_smooth = sum(len(l) for l in polylines)
            print(f"Loaded {len(raw)} coastline polylines, {n_raw} points; "
                  f"smoothed to {n_smooth} points (tolerance {tolerance_m} m).")

            if land_parts:
                import shapely
                land = shapely.union_all(land_parts)
                if not land.is_valid:
                    land = land.buffer(0)
                print(f"Land polygons available: classification will use point-in-coastline.")
            else:
                print("Coastline source has no polygons (open lines): "
                      "classification falls back to the sign of the DEM.")

            if raw_vtk:
                from mesh_builder import write_coastline_vtk
                write_coastline_vtk(raw_vtk, raw, lat_min, lat_max, lon_min, lon_max)

            return polylines, land
        except Exception as e:
            print(f"Could not load vector file '{filepath}': {e}. Extracting zero-crossing from grid.")

    return [], None

def extract_coastline_from_grid(elevations: np.ndarray, lat_min: float, lat_max: float, lon_min: float, lon_max: float):
    """
    Extracts zero-crossing boundary contour directly from elevation matrix (Marching Squares subset)
    """
    height, width = elevations.shape
    d_lon = (lon_max - lon_min) / max(width - 1, 1)
    d_lat = (lat_max - lat_min) / max(height - 1, 1)

    coastline_pts = []

    for r in range(height - 1):
        for c in range(width - 1):
            cell = elevations[r:r+2, c:c+2]
            if np.any(cell < 0) and np.any(cell >= 0):
                lat = lat_min + (r + 0.5) * d_lat  # row 0 is lat_min (see extract_elevation_region)
                lon = lon_min + (c + 0.5) * d_lon
                coastline_pts.append((lat, lon))

    if coastline_pts:
        print(f"Extracted {len(coastline_pts)} coastline boundary points from grid zero-crossings.")
        return [coastline_pts]
    return []
