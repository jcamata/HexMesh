#ifndef PILLOW_UNTANGLE_H
#define PILLOW_UNTANGLE_H

#include <vector>
#include <unordered_map>
#include <cstdint>
#include "hexa.h"

// Where a buffer node came from and how far it may stray: its interface node, the extrusion
// direction (unit), and the thickness it was placed at. Feeds the feasible-set clamp that keeps
// the optimiser from smoothing the layer away.
struct BufferAnchor {
	int orig_id;
	double orig_xyz[3];
	double u[3];
	double t_nom;
};

// Interface node metadata for tangential sliding on the GTS surface inside an elliptic patch.
struct InterfaceAnchor {
	int gts_id;
	double orig_xyz[3];
	double max_disp;
	uint8_t lock_mask;
	int buf_count;
	int buf_ids[4];
};

// Patch-local elliptic untangling: one small solve per connected cluster of inverted elements,
// free nodes = the buffer nodes of that cluster (+ interface nodes of inverted elements if sliding
// is enabled), everything else frozen. Returns how many inverted elements were removed;
// *n_patches_out (optional) gets the number of clusters tried.
int EllipticPatchUntangle(hexa_tree_t *mesh, std::vector<double> &coords, int ref,
                          const std::unordered_map<int, BufferAnchor> &anchors,
                          const std::unordered_map<int, InterfaceAnchor> &iface_anchors,
                          double alpha_min, double angle_cos, int *n_patches_out);

#endif
