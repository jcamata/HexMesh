from dataclasses import dataclass, field
import os

@dataclass
class PreprocConfig:
    lat_min: float = 35.0
    lat_max: float = 45.0
    lon_min: float = -10.0
    lon_max: float = 0.0

    h_min_feature: float = 5000.0
    coastline_tolerance: float = 500.0
    z_min_bathymetry: float = 0.0  # Se setado (ex: z > z_min -> z = z_min)

    nc_file: str = "GEBCO_2023.nc"
    gebco_url: str = "https://www.bodc.ac.uk/data/open_download/gebco/gebco_2023/zip/"
    coastline_file: str = "coastline.geojson"
    coastline_url: str = "https://www.naturalearthdata.com/http//www.naturalearthdata.com/download/10m/physical/ne_10m_coastline.zip"

    include_coastline: bool = True
    include_lakes: bool = False
    include_rivers: bool = False
    include_dams: bool = False

    geology_model: str = "PREM"
    material_output_file: str = "material.input"

    bathymetry_gts: str = "output_bathymetry.gts"
    bathymetry_stl: str = "output_bathymetry.stl"
    bathymetry_vtk: str = "output_bathymetry.vtk"

    coastline_vtk: str = ""  # debug: the coastlines alone, as VTK lines (empty = not written)

    topography_gts: str = "output_topography.gts"
    topography_stl: str = "output_topography.stl"
    topography_vtk: str = "output_topography.vtk"

def parse_config_file(filepath: str) -> PreprocConfig:
    config = PreprocConfig()
    if not os.path.exists(filepath):
        print(f"Warning: Configuration file '{filepath}' not found. Using default parameters.")
        return config

    with open(filepath, 'r', encoding='utf-8') as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#') or line.startswith(';'):
                continue
            if '=' not in line:
                continue

            key, val = line.split('=', 1)
            key = key.strip()
            val = val.strip()

            if key == "z_min":
                key = "z_min_bathymetry"

            if hasattr(config, key):
                attr_type = type(getattr(config, key))
                if attr_type == bool:
                    setattr(config, key, val.lower() in ('true', '1', 'yes'))
                elif attr_type == float:
                    setattr(config, key, float(val))
                elif attr_type == int:
                    setattr(config, key, int(val))
                else:
                    setattr(config, key, val)

    return config
