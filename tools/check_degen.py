#!/usr/bin/env python3
"""Check mesh_0001_0000.h5 for invalid hexes (the ones ParaView culls).

Primary criterion: minimum corner Jacobian <= 0 (matches ParaView/SEM3D).
Also reports the 6-tet signed volume for reference. The h5 connectivity is
already in the output order, so no reordering is needed here.

Usage: python3 tools/check_degen.py [mesh_0001_0000.h5]
"""
import sys
import h5py
import numpy as np

path = sys.argv[1] if len(sys.argv) > 1 else 'mesh_0001_0000.h5'
f = h5py.File(path)
X = f['Nodes'][:]
H = f['Sem3D/Hexa8'][:]
M = f['Sem3D/Mat'][:]
P = f['Sem3D/PillowType'][:] if 'Sem3D/PillowType' in f else None

c = X[H]  # (ne, 8, 3)

# corner Jacobians (VTK hex corner -> its 3 edge neighbors)
nbr = {0: (1, 3, 4), 1: (2, 0, 5), 2: (3, 1, 6), 3: (0, 2, 7),
       4: (7, 5, 0), 5: (4, 6, 1), 6: (5, 7, 2), 7: (6, 4, 3)}
minJ = np.full(len(H), np.inf)
for k, (a, b, d) in nbr.items():
    A = c[:, a] - c[:, k]
    B = c[:, b] - c[:, k]
    C = c[:, d] - c[:, k]
    J = np.einsum('ij,ij->i', np.cross(A, B), C)
    minJ = np.minimum(minJ, J)

# 6-tet signed volume, for reference
T = [(0, 1, 3, 7), (0, 1, 7, 4), (1, 2, 3, 7), (1, 7, 2, 6), (1, 5, 7, 4), (1, 6, 7, 5)]
v = np.zeros(len(H))
for a, b, cc, d in T:
    A = c[:, b] - c[:, a]
    B = c[:, cc] - c[:, a]
    C = c[:, d] - c[:, a]
    v += np.einsum('ij,ij->i', np.cross(A, B), C) / 6.0

bad = np.where(minJ <= 0)[0]
print(f"elements: {len(H)}, nodes: {len(X)}")
print(f"median vol: {np.median(v):.3e}")
print(f"INVALID (min corner Jacobian <= 0): {len(bad)}")
print(f"negative 6-tet volume (reference): {(v <= 0).sum()}")
if len(bad):
    print("\ninvalid by mat:", dict(zip(*map(list, np.unique(M[bad], return_counts=True)))))
    if P is not None:
        print("invalid by pillowtype:", dict(zip(*map(list, np.unique(P[bad], return_counts=True)))))
    zc = c[:, :, 2].mean(axis=1)
    print("\nworst 10 (idx, minJ, mat, pil, z-centroid):")
    for i in bad[np.argsort(minJ[bad])][:10]:
        pt = P[i] if P is not None else -1
        print(f"  {i}  {minJ[i]:.3e}  mat={M[i]}  pil={pt}  zc={zc[i]:.1f}")
sys.exit(1 if len(bad) else 0)
