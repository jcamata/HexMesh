#ifndef PILLOW_UNTANGLE_H
#define PILLOW_UNTANGLE_H

#include <vector>
#include <unordered_map>
#include "hexa.h"

// Where a buffer node came from and how far it may stray: its interface node, the extrusion
// direction (unit), and the thickness it was placed at. Feeds the feasible-set clamp that keeps
// the optimiser from smoothing the layer away.
struct BufferAnchor {
	double orig_xyz[3];
	double u[3];
	double t_nom;
};

// Patch-local elliptic untangling: one small solve per connected cluster of inverted elements,
// free nodes = the buffer nodes of that cluster, everything else frozen. Returns how many
// inverted elements were removed; *n_patches_out (optional) gets the number of clusters tried.
int EllipticPatchUntangle(hexa_tree_t *mesh, std::vector<double> &coords, int ref,
                          const std::unordered_map<int, BufferAnchor> &anchors,
                          double alpha_min, double angle_cos, int *n_patches_out);

#endif
