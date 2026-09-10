#ifndef VERIFY_MESH_H
#define VERIFY_MESH_H

#include <vector>
#include "hexa.h"

// Orientation-aware verification for the (z-inverted) hex mesh. "Inverted" is
// judged against the mesh-dominant reference sign, not a fixed positive sign.

struct MeshAnalysis {
	int reference_sign;              // +1 or -1
	int n_inverted;
	std::vector<int> inverted_ids;
	double h_min;                    // global minimum edge length
	double global_min_sj;
	double global_max_sj;
};

// Load one element's 8 corner coords in reordered (h5) slot order.
void load_elem_xyz(hexa_tree_t *mesh, const std::vector<double> &coords,
                   int iel, double X[8], double Y[8], double Z[8]);

// Full analysis (no printing).
MeshAnalysis analyze_mesh(hexa_tree_t *mesh, const std::vector<double> &coords);

// Backwards-compatible entry: prints the report, returns inverted count.
int VerifyMeshInversion(hexa_tree_t *mesh, const std::vector<double> *coords);

// Inversion map: classifies every inverted element (parent-octree cut pattern,
// lattice extent, node displacement, warp column state) and prints cross-tabs
// plus per-class inversion rates; writes invmap_<tag>.csv in the cwd. Pass an
// empty `prev` when there is no previous-stage snapshot to compare against.
void DumpInversionMap(hexa_tree_t *mesh, const std::vector<double> &coords,
                      const std::vector<double> &prev, const char *tag);

// Face planarity: prints the warp report (max gap / mean edge over the 6 faces
// of every element), returns how many elements exceed tol.
int VerifyFacePlanarity(hexa_tree_t *mesh, const std::vector<double> &coords, double tol = 1e-3);

#endif // VERIFY_MESH_H
