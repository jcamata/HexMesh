#ifndef OPTIMIZE_MESH_H
#define OPTIMIZE_MESH_H

#include <cstdint>
#include <vector>
#include "hexa.h"

// Per-node freedom for the untangler. A node with all three LOCK bits set never
// moves; gts_surface_id >= 0 makes the untangler slide it ALONG that surface
// instead of freely in 3D, so an interface node stays on the interface.
struct NodeConstraint {
	uint8_t lock_mask;  // LOCK_X | LOCK_Y | LOCK_Z for boundary planes/lines/corners
	int gts_surface_id; // -1 if not on GTS surface, 0..N-1 = gdata_vec index, 1000 = topography (tdata)
};

std::vector<NodeConstraint> classify_node_constraints(hexa_tree_t *mesh,
                                                      const std::vector<double> &coords,
                                                      const std::vector<int> &nodes_b_mat,
                                                      std::vector<uint8_t> *wall_out);

// Relax nodes until no element is inverted or folded. Only nodes the `cons`
// entry leaves free are touched. Returns how many inverted elements remain.
int untangle_inversions(hexa_tree_t *mesh, std::vector<double> &coords,
                        const std::vector<NodeConstraint> &cons,
                        const std::vector<uint8_t> &wall_lock, int ref);

#endif // OPTIMIZE_MESH_H
