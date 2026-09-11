#!/usr/bin/env python3
"""Characterise the inverted elements a pillow run leaves behind.

Reads the ASCII VTU of a run plus its run.log (for the pillow element/node
counts) and answers two questions:

  1. WHAT is inverted -- pillow-layer hexes vs host hexes whose corners were
     remapped onto buffer nodes, how many corners go negative, on which side of
     the layer, how thin/anisotropic the elements are.
  2. WHY -- inversion rate against the two geometric signatures of the layer:
     the worst angle between the four corner offset directions of a hex, and the
     thickness ratio between its thickest and thinnest corner. Rates are taken
     over ALL pillow hexes, so the valid ones are the control group.

Usage:  python3 tools/pillow_report.py out/Argostoli_ref4 [out/hawaii_ref4 ...]

Needs the VTU written in ascii (the default) and a run.log containing the
"N new pillow elements added" / "N new nodes created" summary lines.
"""
import glob
import os
import re
import sys

import numpy as np

# VTK hex corner -> the three corners sharing an edge with it
NBR = {0: (1, 3, 4), 1: (2, 0, 5), 2: (3, 1, 6), 3: (0, 2, 7),
       4: (7, 5, 0), 5: (4, 6, 1), 6: (5, 7, 2), 7: (6, 4, 3)}
EDGES = [(0, 1), (1, 2), (2, 3), (3, 0), (4, 5), (5, 6), (6, 7), (7, 4),
         (0, 4), (1, 5), (2, 6), (3, 7)]


def data_array(txt, name, dtype=float):
    m = re.search(r'<DataArray[^>]*Name="%s"[^>]*>(.*?)</DataArray>' % name, txt, re.S)
    if m is None:
        raise SystemExit("no DataArray named %s (is the VTU ascii?)" % name)
    return np.fromstring(m.group(1), sep=' ', dtype=dtype)


def corner_jacobians(c):
    """(ne,8,3) corner coords -> (ne,8) corner Jacobians. <= 0 means inverted."""
    J = np.zeros(c.shape[:2])
    for k, (a, b, d) in NBR.items():
        J[:, k] = np.einsum('ij,ij->i', np.cross(c[:, a] - c[:, k], c[:, b] - c[:, k]),
                            c[:, d] - c[:, k])
    return J


def layer_metrics(P, conn, isbuf):
    """Per pillow hex: worst angle between the 4 offset directions (deg), and the
    thickness ratio thickest/thinnest corner. Each buffer node is paired with its
    nearest interface node -- they are the two ends of the same lattice corner."""
    ang = np.zeros(len(conn))
    ratio = np.zeros(len(conn))
    thick = np.zeros((len(conn), 4))
    for i, (row, m) in enumerate(zip(conn, isbuf)):
        i_n, b_n = row[~m], row[m]
        D = np.linalg.norm(P[b_n][:, None, :] - P[i_n][None, :, :], axis=2)
        v = P[b_n] - P[i_n][D.argmin(1)]
        n = np.linalg.norm(v, axis=1)
        u = v / np.maximum(n, 1e-12)[:, None]
        ang[i] = np.degrees(np.arccos(np.clip((u @ u.T).min(), -1.0, 1.0)))
        ratio[i] = n.max() / max(n.min(), 1e-9)
        thick[i] = n
    return ang, ratio, thick


def rate_table(label, x, inv, bins):
    print("\n  %s -> inversion rate" % label)
    for lo, hi in bins:
        s = (x >= lo) & (x < hi)
        if s.sum():
            print("    %8.1f - %-8.1f : %5d / %7d  = %6.2f %%"
                  % (lo, hi, inv[s].sum(), s.sum(), 100 * inv[s].mean()))


def report(run_dir):
    vtus = [f for f in glob.glob(os.path.join(run_dir, "*_1_0.vtu")) if "lattice" not in f]
    if not vtus:
        raise SystemExit("no mesh VTU in %s" % run_dir)
    txt = open(vtus[0]).read()
    log = open(os.path.join(run_dir, "run.log")).read()
    n_new_e = int(re.search(r'(\d+) new pillow elements added', log).group(1))
    n_new_n = int(re.search(r'(\d+) new nodes created', log).group(1))

    P = data_array(txt, 'Position').reshape(-1, 3)
    H = data_array(txt, 'connectivity', np.int64).reshape(-1, 8)
    mat = data_array(txt, 'ElemType', np.int64)
    fold = data_array(txt, 'NeighborFold', np.int64)
    n_host_e, n_orig_n = len(H) - n_new_e, len(P) - n_new_n

    J = corner_jacobians(P[H])
    nneg = (J <= 0).sum(1)
    bad = np.where(J.min(1) <= 0)[0]

    print("\n===== %s: %d inverted of %d (%d host + %d pillow) ====="
          % (run_dir, len(bad), len(H), n_host_e, n_new_e))
    if len(bad) == 0:
        return
    is_pil = bad >= n_host_e
    print("  pillow-layer hexes %d   host hexes (corner remapped to a buffer node) %d"
          % (is_pil.sum(), (~is_pil).sum()))
    print("  negative corners of 8: " + "  ".join(
        "%d/8:%d" % (n, (nneg[bad] == n).sum()) for n in range(1, 9) if (nneg[bad] == n).sum()))
    print("  material: " + "  ".join("mat%d:%d" % (m, (mat[bad] == m).sum())
                                     for m in np.unique(mat[bad])))
    print("  also neighbour-folded (real 3D overlap): %d" % (fold[bad] != 0).sum())

    c = P[H[bad]]
    ed = np.stack([np.linalg.norm(c[:, a] - c[:, b], axis=1) for a, b in EDGES], 1)
    print("  edge length: shortest %.2f m   median shortest %.1f m   median longest %.1f m"
          % (ed.min(), np.median(ed.min(1)), np.median(ed.max(1))))
    print("  aspect (longest/shortest edge): median %.1f   max %.1f"
          % (np.median(ed.max(1) / ed.min(1)), (ed.max(1) / ed.min(1)).max()))

    host = bad[~is_pil]
    if len(host):
        nbuf = (H[host] >= n_orig_n).sum(1)
        print("  host hexes by number of buffer corners: " + "  ".join(
            "%d:%d" % (k, (nbuf == k).sum()) for k in range(9) if (nbuf == k).sum()))

    # --- control group: every pillow hex, valid or not -----------------------
    pil = np.arange(n_host_e, len(H))
    isbuf = H[pil] >= n_orig_n
    clean = isbuf.sum(1) == 4          # a well-formed 4 interface + 4 buffer hex
    pil, isbuf = pil[clean], isbuf[clean]
    if len(pil) == 0:
        return
    ang, ratio, thick = layer_metrics(P, H[pil], isbuf)
    inv = J[pil].min(1) <= 0
    print("\n  --- against all %d pillow hexes (%d inverted) ---" % (len(pil), inv.sum()))
    rate_table("worst angle between the 4 corner offset directions (deg)", ang, inv,
               [(0, 30), (30, 60), (60, 90), (90, 120), (120, 180.01)])
    rate_table("thickness ratio within a hex (thickest corner / thinnest)", ratio, inv,
               [(1, 1.5), (1.5, 2), (2, 3), (3, 5), (5, 1e9)])
    print("\n  layer thickness on the inverted hexes: min %.2f m  median %.1f m  max %.1f m"
          % (thick[inv].min(), np.median(thick[inv]), thick[inv].max()))
    negb = np.array([(J[e][m] <= 0).sum() for e, m in zip(pil[inv], isbuf[inv])])
    negi = np.array([(J[e][~m] <= 0).sum() for e, m in zip(pil[inv], isbuf[inv])])
    print("  negative corner sits on: buffer side only %d   interface side only %d   both %d"
          % (((negb > 0) & (negi == 0)).sum(), ((negi > 0) & (negb == 0)).sum(),
             ((negb > 0) & (negi > 0)).sum()))


if __name__ == '__main__':
    if len(sys.argv) < 2:
        raise SystemExit(__doc__)
    for d in sys.argv[1:]:
        report(d.rstrip('/'))
