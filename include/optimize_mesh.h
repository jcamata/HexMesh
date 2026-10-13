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

// Level 2, interface step: tangential Laplacian of the material-interface nodes,
// re-projected vertically onto their GTS surface. Before the pillow pass -1; after the
// pillow pass the node count before it, so the pillow inner nodes move with the interface.
void SmoothInterfaceSurface(hexa_tree_t *mesh, std::vector<double> &coords, int first_pillow_node = -1);

// Level 2, top step (before the pillow): smart Laplacian on the free surface -- water lid
// in (x, y) at z = 0, land nodes re-projected vertically onto the topography.
void SmoothTopSurface(hexa_tree_t *mesh, std::vector<double> &coords);

// Level 2, wall step (before the pillow): smart Laplacian in the plane of the 4 side walls
// and the bottom, with box edges, surfaces and zcuts transitions fixed.
void SmoothWallSurfaces(hexa_tree_t *mesh, std::vector<double> &coords);

// Level 1, before the pillow: 1D Laplacian on the 12 edges of the domain box.
void SmoothBoxEdges(hexa_tree_t *mesh, std::vector<double> &coords, const std::vector<int> &nodes_b_mat);

// Level 3, before the pillow: smart Laplacian on the interior nodes, guarding the thickness
// of the hosts the pillow will split.
void SmoothVolume(hexa_tree_t *mesh, std::vector<double> &coords);

// After the pillow: patch-wise pattern search on the worst-dt elements' nodes, maximizing
// the local minimum critical time step (serial only for now).
void OptimizeDtPatches(hexa_tree_t *mesh, std::vector<double> &coords);

// After the pillow: the same goal by nonlinear CG on a smooth barrier of the sigma_min(J)/vp
// dt proxy, no eigen-dt in the loop (serial only for now).
void OptimizeSizeGradient(hexa_tree_t *mesh, std::vector<double> &coords);

#endif // OPTIMIZE_MESH_H
