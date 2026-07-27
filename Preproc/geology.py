import numpy as np

EARTH_RADIUS_KM = 6371.0

def get_prem_profile(max_depth_km: float = 100.0, step_km: float = 1.0):
    """
    Computes 1D PREM (Preliminary Reference Earth Model) velocity profile (Vp, Vs, Density) by depth
    """
    profile = []

    for depth in np.arange(0.0, max_depth_km + step_km, step_km):
        r = EARTH_RADIUS_KM - depth
        x = r / EARTH_RADIUS_KM

        if depth < 3.0:
            # Ocean layer
            vp = 1.45
            vs = 0.0
            density = 1.02
        elif depth < 15.0:
            # Upper Crust
            vp = 5.80
            vs = 3.20
            density = 2.60
        elif depth < 24.4:
            # Lower Crust
            vp = 6.80
            vs = 3.90
            density = 2.90
        else:
            # Upper Mantle / LID (PREM polynomials)
            density = 2.6910 + 0.6924 * x
            vp = 0.8393 + 11.7156 * x - 5.5763 * x * x
            vs = 0.3808 + 7.4251 * x - 3.4270 * x * x

        profile.append({
            'depth_km': depth,
            'vp_kms': vp,
            'vs_kms': vs,
            'density': density
        })

    return profile

def export_material_file(filepath: str, profile: list) -> bool:
    """
    Exports PREM velocity profile to HexMesh material.input file
    """
    try:
        with open(filepath, 'w', encoding='utf-8') as f:
            f.write("# PREM Geological Velocity Profile (Depth_km, Vp_kms, Vs_kms, Density_g_cm3)\n")
            f.write("# Depth(km)\tVp(km/s)\tVs(km/s)\tDensity(g/cm3)\n")
            for p in profile:
                f.write(f"{p['depth_km']:.2f}\t{p['vp_kms']:.4f}\t{p['vs_kms']:.4f}\t{p['density']:.4f}\n")

        print(f"Exported geological velocity profile to '{filepath}'")
        return True
    except Exception as e:
        print(f"Error writing material file '{filepath}': {e}")
        return False
