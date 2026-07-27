#!/usr/bin/env python3
"""
Preproc Pipeline in Python
Processes GEBCO NetCDF bathymetry/topography, coastlines, and PREM geological velocity models.
"""

import sys
import os

from config import parse_config_file
from downloader import download_file_if_missing
from elevation import extract_elevation_region, apply_smoothing_filter
from coastline import load_and_smooth_coastlines, extract_coastline_from_grid
from geology import get_prem_profile, export_material_file
from mesh_builder import build_bathymetry_mesh, build_topography_mesh, export_mesh_files

def main():
    config_file = sys.argv[1] if len(sys.argv) > 1 else "preproc.input"
    print(f"Loading configuration from: {config_file}")
    config = parse_config_file(config_file)

    print("\n--- Preproc Python Configuration ---")
    print(f"Bounding Box: Lat [{config.lat_min}, {config.lat_max}], Lon [{config.lon_min}, {config.lon_max}]")
    print(f"Feature Smoothing h_min: {config.h_min_feature} m")
    print(f"Coastline Tolerance: {config.coastline_tolerance} m")
    print(f"NetCDF File: {config.nc_file}")
    print(f"Coastline File: {config.coastline_file}")
    print(f"Geology Model: {config.geology_model}")
    print("------------------------------------\n")

    # 1. Download files if missing
    if config.gebco_url:
        download_file_if_missing(config.gebco_url, config.nc_file, config=config)

    if config.include_coastline and config.coastline_url:
        download_file_if_missing(config.coastline_url, config.coastline_file, config=config)

    # 2. Extract elevation region
    elevations, width, height, pixel_res_m = extract_elevation_region(
        config.nc_file, config.lat_min, config.lat_max, config.lon_min, config.lon_max
    )

    # 3. Apply smoothing filter to terrain
    elevations = apply_smoothing_filter(elevations, config.h_min_feature, pixel_res_m)

    # 4. Load or extract coastlines & waterbodies
    coastlines = []
    land = None
    if config.include_coastline:
        coastlines, land = load_and_smooth_coastlines(
            config.coastline_file, config.lat_min, config.lat_max, config.lon_min, config.lon_max, config.coastline_tolerance,
            raw_vtk=config.coastline_vtk
        )
        if not coastlines:
            coastlines = extract_coastline_from_grid(
                elevations, config.lat_min, config.lat_max, config.lon_min, config.lon_max
            )

    # 5. Build Bathymetry Mesh (z < 0 preserved, z > 0 plateau z_top = 2*max(z), duplicated boundary nodes)
    print("\nBuilding Bathymetry Mesh...")
    bathy_verts, bathy_tris = build_bathymetry_mesh(
        elevations, config.lat_min, config.lat_max, config.lon_min, config.lon_max, coastlines,
        z_min_bathymetry=config.z_min_bathymetry, land=land
    )

    # 6. Build Topography Mesh (z < 0 flattened to 0, z >= 0 preserved)
    print("\nBuilding Topography Mesh...")
    topo_verts, topo_tris = build_topography_mesh(
        elevations, config.lat_min, config.lat_max, config.lon_min, config.lon_max, coastlines,
        land=land, coast_band_m=config.coastline_tolerance
    )

    # 7. Export surface meshes
    print("\nExporting Mesh Files...")
    export_mesh_files(bathy_verts, bathy_tris, config.bathymetry_stl, config.bathymetry_gts, config.bathymetry_vtk)
    export_mesh_files(topo_verts, topo_tris, config.topography_stl, config.topography_gts, config.topography_vtk)

    # 8. Export Geological Velocity Profile (PREM)
    if config.geology_model.upper() == "PREM":
        print("\nGenerating PREM Geological Velocity Profile...")
        profile = get_prem_profile(max_depth_km=100.0, step_km=1.0)
        export_material_file(config.material_output_file, profile)

    print("\n--- Preproc Pipeline Completed Successfully ---")

if __name__ == "__main__":
    main()
