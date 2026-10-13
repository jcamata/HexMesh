#include <numeric>
#include <cstdio>
#include <cmath>
#include <vector>
#include <array>
#include <algorithm>
#include <unordered_map>
#include <unordered_set>
#include <map>
#include <set>
#include <mpi.h>
#include <sc.h>
#include <sc_containers.h>
#include "hexa.h"
#include "verify_mesh.h"
#include <functional>
#include "mesh_geom.h"
#include "pillow_geometry.h"
#include "optimize_mesh.h"
#include "stability.h"
#ifdef _OPENMP
#include <omp.h>
#endif

/**
 * =========================================================================================
 *   HEXAHEDRAL MESH OPTIMIZATION & UNTANGLING ENGINE FOR SPECTRAL ELEMENT WAVE PROPAGATION
 * =========================================================================================
 *
 * This module implements a comprehensive two-pass optimization framework designed to
 * guarantee valid, conforming, high-order Spectral Element (SEM / GLL) meshes for 3D seismic
 * wave propagation (e.g. SPECFEM3D, SES3D) over complex realistic topography and bathymetry.
 *
 * -----------------------------------------------------------------------------------------
 *   METHODOLOGY & SCIENTIFIC LITERATURE REFERENCES:
 * -----------------------------------------------------------------------------------------
 *
 * [1] REGULARIZED BARRIER UNTANGLING & ALGEBRAIC OBJECTIVES:
 *     - Escobar, J. M., Rodríguez, E., Montenegro, R., Montero, G., & González-Yuste, J. M.
 *       "Simultaneous untangling and smoothing of tetrahedral and hexahedral meshes."
 *       Computer Methods in Applied Mechanics and Engineering, 192(25-26), 2775-2787 (2003).
 *     - Escobar, J. M., Montero, G., Montenegro, R., & Rodríguez, E.
 *       "An algebraic method for smoothing surface triangulations on Riemannian manifolds."
 *       International Journal for Numerical Methods in Engineering, 66(13), 2040-2060 (2006).
 *     - Knupp, P. M.
 *       "Algebraic mesh quality metrics for unstructured initial meshes."
 *       SIAM Journal on Scientific Computing, 23(1), 193-218 (2001).
 *
 * [2] PILLOWING, SHEET OPERATIONS & CONFORMAL COLLAPSE:
 *     - Mitchell, S. A., & Tautges, T. J.
 *       "Pillowing doublets: refining hexahedral meshes via sheet insertion."
 *       In Proceedings of the 4th International Meshing Roundtable, pp. 147-158 (1995).
 *     - Shepherd, J. F., & Johnson, C. R.
 *       "Hexahedral mesh generation constraints."
 *       Engineering with Computers, 24(3), 195-213 (2008).
 *
 * [3] SPECTRAL ELEMENT GLL QUADRATURE & STABILITY ANALYSIS (CFL):
 *     - Komatitsch, D., & Tromp, J.
 *       "Introduction to the spectral element method for three-dimensional seismic wave propagation."
 *       Geophysical Journal International, 139(3), 806-822 (1999).
 *     - De Basabe, J. D., & Sen, M. K.
 *       "Grid dispersion and stability criteria of some common finite-difference and spectral element methods."
 *       Geophysics, 72(6), T81-T95 (2007).
 *
 * [4] TARGET-MATRIX OPTIMIZATION & SMART LAPLACIAN QUALITY ENHANCEMENT:
 *     - Freitag, L. A.
 *       "An adaptive strategy for mesh improvement using Laplacian and optimization-based techniques."
 *       Computational Geometry, Theoretical and Practical (1997).
 *     - Garimella, R. V., Knupp, P. M., & Shashkov, M. J.
 *       "Triangular and quadrilateral surface mesh quality optimization using local subproblem formulations."
 *       Computer Methods in Applied Mechanics and Engineering, 193(9-11), 913-928 (2004).
 * =========================================================================================
 */

// ---- Tunables -------------------------------------------------------------
static const int    MAX_UNTANGLE_ITERS = 2000;
static const int    STALL_PATIENCE     = 100;   // stalled iters to tolerate before giving up
static const double UNTANGLE_STEP0     = 1.5;   // initial line-search step = this * shortest incident edge
static const int    PHASE_A_SWEEPS     = 20;
static const int    PHASE_B_ROUNDS     = 10;
static const double PHASE_B_FRACTION   = 0.05;  // fraction of smallest elems attacked
static const double ESCALATION_CAP     = 0.15;  // max interface move = 15% shortest edge
// Off: measured on 4 cases it moved GTS surface nodes up to 3450 m off their surface and
// removed exactly one inversion in total. The residual is tangled slivers, not nodes that
// need more freedom, so this only costs surface fidelity. The fix is upstream, in the
// projection limiter in moving_nodes.cpp.
static const bool   ESCALATION_ENABLED = false;
static const double UNTANGLE_CAP       = 1.0;   // max total displacement per relax_node call, in local edges
static const double H_FLOOR_FRAC       = 0.05;  // a local edge below this fraction of the mesh scale is not trusted

// ---- Mesh length scale ----------------------------------------------------
// Every step size in this file is a multiple of either the node's own shortest
// incident edge or of this. Absolute lengths in metres (the old 50 / 150 / 250 /
// 350 / 500 literals) are a bug: they assume a ~500 m element. On belle_ile the
// shortest edges are ~1 m and on HexMesh_plane.input the whole domain is under
// 1 unit, so a "small" 50 m step teleported the node several elements away --
// which is how the untangler produced the folds it was supposed to remove.
// ponytail: one file-scope value, one writer (set_mesh_scale), single-threaded.
static double g_h_ref = 1.0;

// Local edge length, floored against a pathological (near-zero) local edge.
static inline double h_floor(double he) {
	double lo = H_FLOOR_FRAC * g_h_ref;
	return (he > lo) ? he : lo;
}

// Median element characteristic size |vol|^(1/3). Independent of
// shortest_incident_edge, which needs g_h_ref for its own fallback.
static void set_mesh_scale(hexa_tree_t *mesh, const std::vector<double> &coords) {
	std::vector<double> h;
	h.reserve(mesh->elements.elem_count);
	for (int iel = 0; iel < mesh->elements.elem_count; iel++) {
		double X[8], Y[8], Z[8];
		load_elem_xyz(mesh, coords, iel, X, Y, Z);
		double v = std::fabs(mgeom::hex_signed_volume(X, Y, Z));
		if (v > 0.0) h.push_back(std::cbrt(v));
	}
	if (h.empty()) { g_h_ref = 1.0; return; }
	std::nth_element(h.begin(), h.begin() + h.size()/2, h.end());
	g_h_ref = h[h.size()/2];
}

// ---- Topology helpers -----------------------------------------------------
struct face_key {
	int n[4];
	bool operator==(const face_key &o) const {
		return n[0]==o.n[0] && n[1]==o.n[1] && n[2]==o.n[2] && n[3]==o.n[3];
	}
};
struct face_key_hash {
	size_t operator()(const face_key &f) const {
		size_t h = 1469598103934665603ULL;
		for (int i = 0; i < 4; i++) { h ^= (size_t)f.n[i]; h *= 1099511628211ULL; }
		return h;
	}
};
static face_key make_face_key(int a, int b, int c, int d) {
	int t[4] = {a,b,c,d};
	std::sort(t, t+4);
	return face_key{ {t[0],t[1],t[2],t[3]} };
}

// eval_gts_height lives in intercept_surface.cpp (declared in hexa.h) so
// GetMeshFromSurface can share it too.

std::vector<NodeConstraint> classify_node_constraints(hexa_tree_t *mesh,
                                                      const std::vector<double> &coords,
                                                      const std::vector<int> &nodes_b_mat,
                                                      std::vector<uint8_t> *wall_out) {
	int nn = coords.size() / 3;
	std::vector<NodeConstraint> cons(nn);
	std::vector<uint8_t> lock(nn, 0);

	for (int i = 0; i < nn; i++) {
		cons[i].lock_mask = 0;
		cons[i].gts_surface_id = -1;
	}

	// External faces: a face appearing exactly once across all elements.
	// Store one representative (elem,iface) per key; single-count keys are walls.
	std::unordered_map<face_key, std::pair<int,int>, face_key_hash> seen;
	std::unordered_set<face_key, face_key_hash> shared;
	seen.reserve(mesh->elements.elem_count * 6);

	for (int iel = 0; iel < mesh->elements.elem_count; iel++) {
		octant_t *e = (octant_t *) sc_array_index(&mesh->elements, iel);
		for (int f = 0; f < 6; f++) {
			int id0 = e->nodes[FaceNodesMap[f][0]].id;
			int id1 = e->nodes[FaceNodesMap[f][1]].id;
			int id2 = e->nodes[FaceNodesMap[f][2]].id;
			int id3 = e->nodes[FaceNodesMap[f][3]].id;
			face_key key = make_face_key(id0, id1, id2, id3);
			auto it = seen.find(key);
			if (it == seen.end()) seen[key] = {iel, f};
			else shared.insert(key);
		}
	}

	double dom_min_x = 1e300, dom_max_x = -1e300;
	double dom_min_y = 1e300, dom_max_y = -1e300;
	double dom_min_z = 1e300, dom_max_z = -1e300;
	for (int i = 0; i < nn; i++) {
		double x = coords[3*i+0], y = coords[3*i+1], z = coords[3*i+2];
		if (x < dom_min_x) dom_min_x = x;
		if (x > dom_max_x) dom_max_x = x;
		if (y < dom_min_y) dom_min_y = y;
		if (y > dom_max_y) dom_max_y = y;
		if (z < dom_min_z) dom_min_z = z;
		if (z > dom_max_z) dom_max_z = z;
	}

	// Wall-plane tolerance. Must scale with the domain: a fixed 1.0 m is wider than the
	// WHOLE domain on a unit-scale model (HexMesh_plane.input spans 0.96), where it locked
	// every node in X and Y and left the untangler unable to move anything at all.
	// min(1.0, 1e-3*L) keeps metre-scale models byte-identical to the old behaviour.
	const double wall_x = std::min(1.0, 1e-3 * (dom_max_x - dom_min_x));
	const double wall_y = std::min(1.0, 1e-3 * (dom_max_y - dom_min_y));
	const double wall_z = std::min(1.0, 1e-3 * (dom_max_z - dom_min_z));

	for (auto &kv : seen) {
		if (shared.count(kv.first)) continue; // internal face, skip
		int iel = kv.second.first, f = kv.second.second;
		octant_t *e = (octant_t *) sc_array_index(&mesh->elements, iel);
		int nid[4];
		for (int k = 0; k < 4; k++) {
			nid[k] = e->nodes[FaceNodesMap[f][k]].id;
		}

		uint8_t m = 0;
		if (nid[0] >= 0 && nid[0] < nn && nid[1] >= 0 && nid[1] < nn &&
		    nid[2] >= 0 && nid[2] < nn && nid[3] >= 0 && nid[3] < nn) {
			bool cx = (std::fabs(coords[3*nid[0]+0] - coords[3*nid[1]+0]) < 1e-4 &&
			           std::fabs(coords[3*nid[0]+0] - coords[3*nid[2]+0]) < 1e-4 &&
			           std::fabs(coords[3*nid[0]+0] - coords[3*nid[3]+0]) < 1e-4);
			bool cy = (std::fabs(coords[3*nid[0]+1] - coords[3*nid[1]+1]) < 1e-4 &&
			           std::fabs(coords[3*nid[0]+1] - coords[3*nid[2]+1]) < 1e-4 &&
			           std::fabs(coords[3*nid[0]+1] - coords[3*nid[3]+1]) < 1e-4);
			bool cz = (std::fabs(coords[3*nid[0]+2] - coords[3*nid[1]+2]) < 1e-4 &&
			           std::fabs(coords[3*nid[0]+2] - coords[3*nid[2]+2]) < 1e-4 &&
			           std::fabs(coords[3*nid[0]+2] - coords[3*nid[3]+2]) < 1e-4);

			if (cx) {
				double xval = coords[3*nid[0]+0];
				if (std::fabs(xval - dom_min_x) < wall_x || std::fabs(xval - dom_max_x) < wall_x)
					m |= mgeom::LOCK_X;
			}
			if (cy) {
				double yval = coords[3*nid[0]+1];
				if (std::fabs(yval - dom_min_y) < wall_y || std::fabs(yval - dom_max_y) < wall_y)
					m |= mgeom::LOCK_Y;
			}
			if (cz) {
				double zval = coords[3*nid[0]+2];
				if (std::fabs(zval - dom_min_z) < wall_z)
					m |= mgeom::LOCK_Z;
			}
		}

		for (int k = 0; k < 4; k++) {
			if (nid[k] >= 0 && nid[k] < nn) {
				lock[nid[k]] |= m;
				cons[nid[k]].lock_mask |= m;
			}
		}
	}

	// Also directly lock any node sitting on the global bounding planes (X-, X+, Y-, Y+, Z-)
	for (int i = 0; i < nn; i++) {
		if (std::fabs(coords[3*i+0] - dom_min_x) < wall_x || std::fabs(coords[3*i+0] - dom_max_x) < wall_x) {
			lock[i] |= mgeom::LOCK_X;
			cons[i].lock_mask |= mgeom::LOCK_X;
		}
		if (std::fabs(coords[3*i+1] - dom_min_y) < wall_y || std::fabs(coords[3*i+1] - dom_max_y) < wall_y) {
			lock[i] |= mgeom::LOCK_Y;
			cons[i].lock_mask |= mgeom::LOCK_Y;
		}
		if (std::fabs(coords[3*i+2] - dom_min_z) < wall_z) {
			lock[i] |= mgeom::LOCK_Z;
			cons[i].lock_mask |= mgeom::LOCK_Z;
		}
	}

	// Wall-only mask
	if (wall_out) *wall_out = lock;

	// Interface nodes: associate each with its closest matching GTS surface id
	for (int nid : nodes_b_mat) {
		if (nid >= 0 && nid < nn) {
			int best_k = -1;
			double best_dz = 1e300;
			double nx = coords[3*nid+0], ny = coords[3*nid+1], nz = coords[3*nid+2];
			for (size_t k = 0; k < mesh->gdata_vec.size(); k++) {
				double z_eval = 0.0;
				if (eval_gts_height(mesh, (int)k, nx, ny, z_eval)) {
					double dz = std::fabs(nz - z_eval);
					if (dz < best_dz) {
						best_dz = dz;
						best_k = (int)k;
					}
				}
			}
			cons[nid].gts_surface_id = best_k;
		}
	}

	return cons;
}

// ---- Face adjacency & folding --------------------------------------------
// (elem,face) -> neighbour packed as other*8 + other_face, or -1 for a boundary face.
// ponytail: one global, one writer (build_face_neighbours), single-threaded -- same
// reason as g_h_ref: node_state has ~20 call sites and does not need a new parameter.
static std::vector<int> g_face_nbr;

static void build_face_neighbours(hexa_tree_t *mesh) {
	int ne = mesh->elements.elem_count;
	g_face_nbr.assign((size_t)6*ne, -1);
	std::unordered_map<face_key, int, face_key_hash> owner;
	owner.reserve((size_t)ne*6);
	for (int iel = 0; iel < ne; iel++) {
		octant_t *e = (octant_t *) sc_array_index(&mesh->elements, iel);
		for (int f = 0; f < 6; f++) {
			face_key k = make_face_key(e->nodes[FaceNodesMap[f][0]].id, e->nodes[FaceNodesMap[f][1]].id,
			                           e->nodes[FaceNodesMap[f][2]].id, e->nodes[FaceNodesMap[f][3]].id);
			auto it = owner.find(k);
			if (it == owner.end()) { owner[k] = iel*8 + f; continue; }
			int oe = it->second / 8, of = it->second % 8;
			g_face_nbr[6*iel + f] = oe*8 + of;
			g_face_nbr[6*oe  + of] = iel*8 + f;
			owner.erase(it);
		}
	}
}

static inline void face_complement_centroid(hexa_tree_t *mesh, const std::vector<double> &coords,
                                            int iel, int f, double c[3]) {
	octant_t *e = (octant_t *) sc_array_index(&mesh->elements, iel);
	bool on_face[8] = {false};
	for (int k = 0; k < 4; k++) on_face[FaceNodesMap[f][k]] = true;
	c[0] = c[1] = c[2] = 0.0;
	for (int k = 0; k < 8; k++) {
		if (on_face[k]) continue;
		int id = e->nodes[k].id;
		for (int d = 0; d < 3; d++) c[d] += 0.25 * coords[3*id+d];
	}
}

// How many of this element's faces it has folded through into the neighbour.
static int elem_fold_count(hexa_tree_t *mesh, const std::vector<double> &coords, int iel) {
	if (g_face_nbr.empty()) return 0;
	octant_t *e = (octant_t *) sc_array_index(&mesh->elements, iel);
	int n = 0;
	for (int f = 0; f < 6; f++) {
		int packed = g_face_nbr[6*iel + f];
		if (packed < 0) continue;
		double q[4][3], ca[3], cb[3];
		for (int k = 0; k < 4; k++) {
			int id = e->nodes[FaceNodesMap[f][k]].id;
			for (int d = 0; d < 3; d++) q[k][d] = coords[3*id+d];
		}
		face_complement_centroid(mesh, coords, iel, f, ca);
		face_complement_centroid(mesh, coords, packed / 8, packed % 8, cb);
		if (mgeom::faces_folded(q, ca, cb)) n++;
	}
	return n;
}

// Mesh-level badness, same lexicographic order as better_state: inversions then folds.
struct MeshBad { int inv, fold; };
static inline bool bad_worse(const MeshBad &b, const MeshBad &a) {
	return (b.inv != a.inv) ? (b.inv > a.inv) : (b.fold > a.fold);
}

static int count_folded_elements(hexa_tree_t *mesh, const std::vector<double> &coords) {
	int n = 0;
	for (int iel = 0; iel < mesh->elements.elem_count; iel++)
		if (elem_fold_count(mesh, coords, iel) > 0) n++;
	return n;
}

// node -> list of incident element ids
static std::vector<std::vector<int>> build_incidence(hexa_tree_t *mesh, int n_nodes) {
	std::vector<std::vector<int>> inc(n_nodes);
	for (int iel = 0; iel < mesh->elements.elem_count; iel++) {
		octant_t *e = (octant_t *) sc_array_index(&mesh->elements, iel);
		for (int ino = 0; ino < 8; ino++) {
			int id = e->nodes[ino].id;
			if (id >= 0 && id < n_nodes) inc[id].push_back(iel);
		}
	}
	return inc;
}

// shortest edge incident to a node (physical)
static double shortest_incident_edge(hexa_tree_t *mesh, const std::vector<double> &coords,
                                     const std::vector<int> &inc, int node) {
	static const int E[12][2] = {
		{0,1},{1,2},{2,3},{3,0},{4,5},{5,6},{6,7},{7,4},{0,4},{1,5},{2,6},{3,7}
	};
	double h = 1e300;
	for (int iel : inc) {
		octant_t *e = (octant_t *) sc_array_index(&mesh->elements, iel);
		for (int k = 0; k < 12; k++) {
			int a = e->nodes[E[k][0]].id, b = e->nodes[E[k][1]].id;
			if (a != node && b != node) continue;
			if (a == b) continue;
			double dx=coords[3*a]-coords[3*b], dy=coords[3*a+1]-coords[3*b+1], dz=coords[3*a+2]-coords[3*b+2];
			double len = std::sqrt(dx*dx+dy*dy+dz*dz);
			if (len > 1e-9 * g_h_ref && len < h) h = len;
		}
	}
	return (h < 1e299) ? h : g_h_ref;
}

static const double SEARCH_DIRS[26][3] = {
	{ 1, 0, 0}, {-1, 0, 0}, { 0, 1, 0}, { 0,-1, 0}, { 0, 0, 1}, { 0, 0,-1},
	{ 1, 1, 0}, {-1, 1, 0}, { 1,-1, 0}, {-1,-1, 0},
	{ 1, 0, 1}, {-1, 0, 1}, { 1, 0,-1}, {-1, 0,-1},
	{ 0, 1, 1}, { 0,-1, 1}, { 0, 1,-1}, { 0,-1,-1},
	{ 1, 1, 1}, {-1, 1, 1}, { 1,-1, 1}, {-1,-1, 1},
	{ 1, 1,-1}, {-1, 1,-1}, { 1,-1,-1}, {-1,-1,-1}
};

// Local state of a node's incident elements: how many are inverted, how many have
// folded through a face into a neighbour, worst scaled Jacobian, sum of negative
// Jacobians, and face warp.
//
// The old z_disparity / xy_disparity terms are gone: penalising a node for sitting at a
// different height from its horizontal neighbours is a penalty on TOPOGRAPHY, and it
// carried weight 0.50 against min_sj -- the optimizer was being paid to flatten relief.
struct NodeState {
	int n_inv;
	int n_fold;              // face-folds among incident elements (3D interpenetration)
	double min_sj;
	double sum_neg_sj;
	double max_face_warp;    // Max face non-planarity ratio among incident elements
	double sum_sj;
};

static NodeState node_state(hexa_tree_t *mesh, const std::vector<double> &coords,
                            const std::vector<std::vector<int>> &adj,
                            const std::vector<int> &inc, int node, int ref) {
	NodeState s;
	s.n_inv = 0;
	s.n_fold = 0;
	s.min_sj = 1e300;
	s.sum_neg_sj = 0.0;
	s.sum_sj = 0.0;
	s.max_face_warp = 0.0;

	for (int iel : inc) {
		double X[8], Y[8], Z[8];
		load_elem_xyz(mesh, coords, iel, X, Y, Z);
		double vol = mgeom::hex_signed_volume(X, Y, Z);
		double sj  = mgeom::hex_min_corner_sj(X, Y, Z);
		double ssj = sj * ref;
		if (mgeom::is_inverted(vol, sj, ref)) {
			s.n_inv++;
			s.sum_neg_sj += (ssj < 0.0) ? ssj : -1e-4;
		}
		s.sum_sj += ssj;
		if (ssj < s.min_sj) s.min_sj = ssj;

		double warp = mgeom::hex_max_face_warp(X, Y, Z);
		if (warp > s.max_face_warp) s.max_face_warp = warp;

		s.n_fold += elem_fold_count(mesh, coords, iel);
	}

	(void)adj;   // kept in the signature only to avoid touching ~20 call sites
	return s;
}

// Ordering is LEXICOGRAPHIC and the first key is the inversion count, alone.
// An inverted element stops the wave-propagation solver from starting at all, so it is a
// hard constraint: no amount of fold reduction may ever pay for re-introducing one.
// Folds rank immediately below -- they are unusable geometry too, but a mesh carrying
// them still runs, so they are only decided once the inversion counts tie.
static bool better_state(const NodeState &b, const NodeState &a, double h_ref = 0.0) {
	(void)h_ref;
	if (b.n_inv != a.n_inv) return b.n_inv < a.n_inv;

	if (b.n_inv > 0) {
		// Still untangling. Drive the worst corner Jacobian up; folds are a tiebreak here
		// and must not block a move that is making progress toward removing an inversion.
		if (b.min_sj > a.min_sj + 1e-6) return true;
		if (b.min_sj < a.min_sj - 1e-6) return false;
		if (b.sum_neg_sj > a.sum_neg_sj + 1e-6) return true;
		if (b.sum_neg_sj < a.sum_neg_sj - 1e-6) return false;
		if (b.n_fold != a.n_fold) return b.n_fold < a.n_fold;
		return b.sum_sj > a.sum_sj + 1e-6;
	}

	// No inversions on either side: folds become the primary target.
	if (b.n_fold != a.n_fold) return b.n_fold < a.n_fold;

	// Valid mesh phase: worst corner Jacobian, with face warp as a mild shape penalty.
	double score_a = a.min_sj - 0.20 * a.max_face_warp;
	double score_b = b.min_sj - 0.20 * b.max_face_warp;
	if (score_b > score_a + 1e-5) return true;
	if (score_b < score_a - 1e-5) return false;
	return b.sum_sj > a.sum_sj + 1e-12;
}

// Perform surface-constrained tangential relaxation for a node on a GTS surface
static bool relax_gts_surface_node(hexa_tree_t *mesh, std::vector<double> &coords,
                                   const std::vector<std::vector<int>> &adj,
                                   const std::vector<int> &inc, int node, int ref,
                                   const NodeConstraint &cons, double he = 0.0) {
	if (inc.empty() || adj[node].empty() || cons.gts_surface_id < 0) return false;
	double px = coords[3*node+0], py = coords[3*node+1], pz = coords[3*node+2];
	NodeState s0 = node_state(mesh, coords, adj, inc, node, ref);
	bool moved = false;

	// 1. Try moving toward the horizontal centroid of incident neighbors on the surface
	double cx = 0.0, cy = 0.0;
	int cnt = 0;
	for (int nb : adj[node]) {
		if (nb >= 0 && nb < (int)coords.size()/3) {
			cx += coords[3*nb+0]; cy += coords[3*nb+1]; cnt++;
		}
	}
	if (cnt > 0) {
		cx /= cnt; cy /= cnt;
		double dx = cx - px, dy = cy - py, dz = 0.0;
		mgeom::apply_lock(cons.lock_mask, dx, dy, dz);
		for (double alpha : {1.0, 0.75, 0.5, 0.25, 0.1, 0.05, 0.01}) {
			double cand_x = px + alpha * dx, cand_y = py + alpha * dy, cand_z = pz;
			if (eval_gts_height(mesh, cons.gts_surface_id, cand_x, cand_y, cand_z)) {
				coords[3*node+0] = cand_x; coords[3*node+1] = cand_y; coords[3*node+2] = cand_z;
				NodeState s = node_state(mesh, coords, adj, inc, node, ref);
				if (better_state(s, s0, he)) {
					s0 = s; moved = true; break;
				}
				coords[3*node+0] = px; coords[3*node+1] = py; coords[3*node+2] = pz;
			}
		}
	}

	// 2. Try 16 tangential directions along the GTS surface with scale-relative steps
	double h_surf = h_floor(he);
	double SURF_STEPS[] = {0.8*h_surf, 0.5*h_surf, 0.25*h_surf, 0.1*h_surf, 0.05*h_surf,
	                       0.01*h_surf, 0.005*h_surf, 0.002*h_surf, 0.001*h_surf, 0.0005*h_surf};
	for (int dir = 0; dir < 16; dir++) {
		double theta = dir * (2.0 * M_PI / 16.0);
		double vx = std::cos(theta), vy = std::sin(theta), vz = 0.0;
		mgeom::apply_lock(cons.lock_mask, vx, vy, vz);
		double vlen = std::sqrt(vx*vx + vy*vy);
		if (vlen < 1e-6) continue;
		vx /= vlen; vy /= vlen;

		for (double step : SURF_STEPS) {
			double cand_x = coords[3*node+0] + vx * step;
			double cand_y = coords[3*node+1] + vy * step;
			double cand_z = coords[3*node+2];

			if (eval_gts_height(mesh, cons.gts_surface_id, cand_x, cand_y, cand_z)) {
				double cur_x = coords[3*node+0], cur_y = coords[3*node+1], cur_z = coords[3*node+2];
				coords[3*node+0] = cand_x; coords[3*node+1] = cand_y; coords[3*node+2] = cand_z;
				NodeState s = node_state(mesh, coords, adj, inc, node, ref);
				if (better_state(s, s0, he)) {
					s0 = s; moved = true; break;
				}
				coords[3*node+0] = cur_x; coords[3*node+1] = cur_y; coords[3*node+2] = cur_z;
			}
		}
	}

	return moved;
}

// Direct volume-gradient ascent: moves node in the direction that maximizes incident elements' volume / jacobian
static bool relax_volume_gradient(hexa_tree_t *mesh, std::vector<double> &coords,
                                   const std::vector<std::vector<int>> &adj,
                                   const std::vector<int> &inc, int node, int ref,
                                   const NodeConstraint &cons, double he = 0.0) {
	if (cons.gts_surface_id >= 0 || inc.empty()) return false;
	double px = coords[3*node+0], py = coords[3*node+1], pz = coords[3*node+2];
	NodeState s0 = node_state(mesh, coords, adj, inc, node, ref);

	double eps = 1.0;
	coords[3*node+0] = px + eps;
	NodeState sx_p = node_state(mesh, coords, adj, inc, node, ref);
	coords[3*node+0] = px - eps;
	NodeState sx_m = node_state(mesh, coords, adj, inc, node, ref);
	coords[3*node+0] = px;
	double gx = ((sx_p.sum_neg_sj - sx_m.sum_neg_sj) + 2.0 * (sx_p.min_sj - sx_m.min_sj)) / (2.0 * eps);

	coords[3*node+1] = py + eps;
	NodeState sy_p = node_state(mesh, coords, adj, inc, node, ref);
	coords[3*node+1] = py - eps;
	NodeState sy_m = node_state(mesh, coords, adj, inc, node, ref);
	coords[3*node+1] = py;
	double gy = ((sy_p.sum_neg_sj - sy_m.sum_neg_sj) + 2.0 * (sy_p.min_sj - sy_m.min_sj)) / (2.0 * eps);

	coords[3*node+2] = pz + eps;
	NodeState sz_p = node_state(mesh, coords, adj, inc, node, ref);
	coords[3*node+2] = pz - eps;
	NodeState sz_m = node_state(mesh, coords, adj, inc, node, ref);
	coords[3*node+2] = pz;
	double gz = ((sz_p.sum_neg_sj - sz_m.sum_neg_sj) + 2.0 * (sz_p.min_sj - sz_m.min_sj)) / (2.0 * eps);

	mgeom::apply_lock(cons.lock_mask, gx, gy, gz);
	double glen = std::sqrt(gx*gx + gy*gy + gz*gz);
	if (glen < 1e-12) return false;
	gx /= glen; gy /= glen; gz /= glen;

	bool moved = false;
	double h_scale = h_floor(he);
	double GRAD_STEPS[] = {1.0*h_scale, 0.6*h_scale, 0.3*h_scale, 0.15*h_scale, 0.05*h_scale,
	                       0.01*h_scale, 0.005*h_scale, 0.002*h_scale, 0.001*h_scale, 0.0005*h_scale};
	for (double sign : {1.0, -1.0}) {
		for (double step : GRAD_STEPS) {
			coords[3*node+0] = px + sign * gx * step;
			coords[3*node+1] = py + sign * gy * step;
			coords[3*node+2] = pz + sign * gz * step;

			NodeState s = node_state(mesh, coords, adj, inc, node, ref);
			if (better_state(s, s0, he)) {
				s0 = s; moved = true; break;
			}
			coords[3*node+0] = px; coords[3*node+1] = py; coords[3*node+2] = pz;
		}
		if (moved) break;
	}
	return moved;
}

// Relax a pillow layer node along the segment connecting an interior neighbor and a surface neighbor
static bool relax_pillow_node(hexa_tree_t *mesh, std::vector<double> &coords,
                              const std::vector<std::vector<int>> &adj,
                              const std::vector<int> &inc, int node, int ref,
                              const std::vector<NodeConstraint> &cons, double he = 0.0) {
	if (cons[node].gts_surface_id >= 0) return false; // purely internal/pillow node
	int n_surf = -1, n_oct = -1;
	for (int nb : adj[node]) {
		if (nb < 0 || nb >= (int)coords.size()/3) continue;
		if (cons[nb].gts_surface_id >= 0 && n_surf < 0) n_surf = nb;
		else if (cons[nb].gts_surface_id < 0 && n_oct < 0) n_oct = nb;
	}
	if (n_surf < 0 || n_oct < 0) return false;

	double px = coords[3*node+0], py = coords[3*node+1], pz = coords[3*node+2];
	double sx = coords[3*n_surf+0], sy = coords[3*n_surf+1], sz = coords[3*n_surf+2];
	double ox = coords[3*n_oct+0], oy = coords[3*n_oct+1], oz = coords[3*n_oct+2];

	NodeState s0 = node_state(mesh, coords, adj, inc, node, ref);
	bool moved = false;

	for (double alpha : {0.5, 0.35, 0.65, 0.2, 0.8, 0.45, 0.55, 0.1, 0.9}) {
		double cand_x = (1.0 - alpha) * ox + alpha * sx;
		double cand_y = (1.0 - alpha) * oy + alpha * sy;
		double cand_z = (1.0 - alpha) * oz + alpha * sz;

		double dx = cand_x - px, dy = cand_y - py, dz = cand_z - pz;
		mgeom::apply_lock(cons[node].lock_mask, dx, dy, dz);
		coords[3*node+0] = px + dx;
		coords[3*node+1] = py + dy;
		coords[3*node+2] = pz + dz;

		NodeState s = node_state(mesh, coords, adj, inc, node, ref);
		if (better_state(s, s0, he)) {
			s0 = s;
			moved = true;
			break;
		}
		coords[3*node+0] = px; coords[3*node+1] = py; coords[3*node+2] = pz;
	}
	return moved;
}

// 3D Directional line-search that improves the node's incident state monotonically.
static bool relax_node(hexa_tree_t *mesh, std::vector<double> &coords,
                       const std::vector<std::vector<int>> &adj,
                       const std::vector<int> &inc, int node, int ref,
                       uint8_t eff_mask, double step0, double cap, double he = 0.0) {
	double sx = coords[3*node+0], sy = coords[3*node+1], sz = coords[3*node+2];
	NodeState s0 = node_state(mesh, coords, adj, inc, node, ref);
	bool moved = false;

	// step0 is the caller's budget (UNTANGLE_STEP0 local edges). It used to be
	// accepted and then ignored in favour of a fixed metre ladder.
	double h = h_floor(he);
	if (step0 <= 0.0) step0 = UNTANGLE_STEP0 * h;
	double max_z_step = 0.7 * h;
	const double FRAC[] = {1.0, 0.6, 0.4, 0.2, 0.1, 0.04, 0.02, 0.008, 0.004, 0.002, 0.0004};
	double STEPS[sizeof(FRAC)/sizeof(FRAC[0])];
	for (size_t i = 0; i < sizeof(FRAC)/sizeof(FRAC[0]); i++) STEPS[i] = FRAC[i] * step0;

	for (int d = 0; d < 26; d++) {
		double vx = SEARCH_DIRS[d][0], vy = SEARCH_DIRS[d][1], vz = SEARCH_DIRS[d][2];
		mgeom::apply_lock(eff_mask, vx, vy, vz);
		double vlen = std::sqrt(vx*vx + vy*vy + vz*vz);
		if (vlen < 1e-6) continue;
		vx /= vlen; vy /= vlen; vz /= vlen;

		for (double step : STEPS) {
			if (std::fabs(vz * step) > max_z_step) continue;

			double cand_x = coords[3*node+0] + vx * step;
			double cand_y = coords[3*node+1] + vy * step;
			double cand_z = coords[3*node+2] + vz * step;

			if (cap > 0.0) {
				double dx = cand_x - sx, dy = cand_y - sy, dz = cand_z - sz;
				double dist = std::sqrt(dx*dx + dy*dy + dz*dz);
				if (dist > cap) {
					cand_x = sx + dx * (cap / dist);
					cand_y = sy + dy * (cap / dist);
					cand_z = sz + dz * (cap / dist);
				}
			}

			double px = coords[3*node+0], py = coords[3*node+1], pz = coords[3*node+2];
			coords[3*node+0] = cand_x; coords[3*node+1] = cand_y; coords[3*node+2] = cand_z;
			NodeState s = node_state(mesh, coords, adj, inc, node, ref);
			if (better_state(s, s0, he)) {
				s0 = s;
				moved = true;
				break;
			}
			coords[3*node+0] = px; coords[3*node+1] = py; coords[3*node+2] = pz; // revert
		}
	}
	return moved;
}

// Exhaustive 3D box probe: tests a grid of 3D displacements scaled to local element size he
static bool probe_node_3d_box(hexa_tree_t *mesh, std::vector<double> &coords,
                              const std::vector<std::vector<int>> &adj,
                              const std::vector<int> &inc, int node, int ref,
                              const NodeConstraint &cons, double he = 0.0) {
	if (cons.gts_surface_id >= 0 || inc.empty()) return false;
	double px = coords[3*node+0], py = coords[3*node+1], pz = coords[3*node+2];
	NodeState s0 = node_state(mesh, coords, adj, inc, node, ref);
	NodeState best_s = s0;
	double best_x = px, best_y = py, best_z = pz;
	bool found = false;

	double scale = h_floor(he);
	static const double FX[] = {-0.6, -0.3, -0.1, 0.0, 0.1, 0.3, 0.6};
	static const double FY[] = {-0.6, -0.3, -0.1, 0.0, 0.1, 0.3, 0.6};
	static const double FZ[] = {-0.8, -0.5, -0.3, -0.1, 0.0, 0.1, 0.3, 0.5, 0.8};

	for (double fx : FX) {
		for (double fy : FY) {
			for (double fz : FZ) {
				if (fx == 0.0 && fy == 0.0 && fz == 0.0) continue;
				double mx = fx * scale, my = fy * scale, mz = fz * scale;
				mgeom::apply_lock(cons.lock_mask, mx, my, mz);
				if (std::fabs(mx) < 1e-6 && std::fabs(my) < 1e-6 && std::fabs(mz) < 1e-6) continue;

				coords[3*node+0] = px + mx;
				coords[3*node+1] = py + my;
				coords[3*node+2] = pz + mz;

				NodeState s = node_state(mesh, coords, adj, inc, node, ref);
				if (better_state(s, best_s, scale)) {
					best_s = s;
					best_x = coords[3*node+0];
					best_y = coords[3*node+1];
					best_z = coords[3*node+2];
					found = true;
				}
			}
		}
	}

	if (found) {
		coords[3*node+0] = best_x;
		coords[3*node+1] = best_y;
		coords[3*node+2] = best_z;
		return true;
	}
	coords[3*node+0] = px; coords[3*node+1] = py; coords[3*node+2] = pz;
	return false;
}

// Escobar-Montenegro Regularized Barrier Inversion-Free Energy:
// For corner Jacobian J, regularized h(J) = 0.5 * (J + sqrt(J^2 + 4*eps^2))
// Energy eta_c = ||J_c||_F^2 / (h(J_c))^(2/3)
static double eval_escobar_node_energy(hexa_tree_t *mesh, const std::vector<double> &coords,
                                      const std::vector<int> &inc, int node, int ref, double eps) {
	double total_energy = 0.0;
	for (int ie : inc) {
		double X[8], Y[8], Z[8];
		load_elem_xyz(mesh, coords, ie, X, Y, Z);
		for (int c = 0; c < 8; c++) {
			int i0 = c, i1 = mgeom::CORNER_NB[c][0], i2 = mgeom::CORNER_NB[c][1], i3 = mgeom::CORNER_NB[c][2];
			double e1x = X[i1] - X[i0], e1y = Y[i1] - Y[i0], e1z = Z[i1] - Z[i0];
			double e2x = X[i2] - X[i0], e2y = Y[i2] - Y[i0], e2z = Z[i2] - Z[i0];
			double e3x = X[i3] - X[i0], e3y = Y[i3] - Y[i0], e3z = Z[i3] - Z[i0];

			double J_c = (double)ref * (e1x * (e2y * e3z - e2z * e3y)
			                          - e1y * (e2x * e3z - e2z * e3x)
			                          + e1z * (e2x * e3y - e2y * e3x));

			double frob2 = (e1x*e1x + e1y*e1y + e1z*e1z)
			             + (e2x*e2x + e2y*e2y + e2z*e2z)
			             + (e3x*e3x + e3y*e3y + e3z*e3z);

			double h_J = 0.5 * (J_c + std::sqrt(J_c * J_c + 4.0 * eps * eps));
			double denom = std::pow(std::max(h_J, 1e-18), 2.0 / 3.0);
			total_energy += frob2 / denom;
		}
	}
	return total_energy;
}

// Relax a node using analytical gradient descent on the Escobar Regularized Inversion-Free Barrier
static bool relax_node_escobar_barrier(hexa_tree_t *mesh, std::vector<double> &coords,
                                       const std::vector<std::vector<int>> &adj,
                                       const std::vector<int> &inc, int node, int ref,
                                       const NodeConstraint &cons, double he = 0.0) {
	if (inc.empty()) return false;
	if (cons.gts_surface_id >= 0 || cons.lock_mask == (mgeom::LOCK_X|mgeom::LOCK_Y|mgeom::LOCK_Z)) return false;

	double px = coords[3*node+0], py = coords[3*node+1], pz = coords[3*node+2];
	NodeState s0 = node_state(mesh, coords, adj, inc, node, ref);

	double h_scale = h_floor(he);
	double eps = 1e-3 * h_scale * h_scale * h_scale; // scale-invariant regularization parameter

	double E0 = eval_escobar_node_energy(mesh, coords, inc, node, ref, eps);

	// Compute gradient of Escobar energy via central differences
	double delta = 1e-3 * h_scale;
	coords[3*node+0] = px + delta;
	double Ex_p = eval_escobar_node_energy(mesh, coords, inc, node, ref, eps);
	coords[3*node+0] = px - delta;
	double Ex_m = eval_escobar_node_energy(mesh, coords, inc, node, ref, eps);
	coords[3*node+0] = px;

	coords[3*node+1] = py + delta;
	double Ey_p = eval_escobar_node_energy(mesh, coords, inc, node, ref, eps);
	coords[3*node+1] = py - delta;
	double Ey_m = eval_escobar_node_energy(mesh, coords, inc, node, ref, eps);
	coords[3*node+1] = py;

	coords[3*node+2] = pz + delta;
	double Ez_p = eval_escobar_node_energy(mesh, coords, inc, node, ref, eps);
	coords[3*node+2] = pz - delta;
	double Ez_m = eval_escobar_node_energy(mesh, coords, inc, node, ref, eps);
	coords[3*node+2] = pz;

	double gx = (Ex_p - Ex_m) / (2.0 * delta);
	double gy = (Ey_p - Ey_m) / (2.0 * delta);
	double gz = (Ez_p - Ez_m) / (2.0 * delta);

	mgeom::apply_lock(cons.lock_mask, gx, gy, gz);
	double glen = std::sqrt(gx*gx + gy*gy + gz*gz);
	if (glen < 1e-12) return false;
	gx /= glen; gy /= glen; gz /= glen;

	// Backtracking line search along negative gradient direction
	bool moved = false;
	double STEPS[] = {1.5*h_scale, 1.0*h_scale, 0.6*h_scale, 0.3*h_scale, 0.15*h_scale, 0.05*h_scale, 0.01*h_scale};
	for (double step : STEPS) {
		coords[3*node+0] = px - gx * step;
		coords[3*node+1] = py - gy * step;
		coords[3*node+2] = pz - gz * step;

		NodeState s = node_state(mesh, coords, adj, inc, node, ref);
		double E = eval_escobar_node_energy(mesh, coords, inc, node, ref, eps);
		if ((s.n_inv < s0.n_inv) || (s.n_inv == s0.n_inv && E < E0 - 1e-6) || better_state(s, s0, he)) {
			s0 = s; E0 = E; moved = true; break;
		}
		coords[3*node+0] = px; coords[3*node+1] = py; coords[3*node+2] = pz;
	}

	return moved;
}

// Relax an interior node to the midpoint of its topological vertical column neighbors
static bool relax_vertical_column_midpoint(hexa_tree_t *mesh, std::vector<double> &coords,
                                           const std::vector<std::vector<int>> &adj,
                                           const std::vector<int> &inc, int node, int ref,
                                           const NodeConstraint &cons) {
	if (inc.empty() || cons.gts_surface_id >= 0 || (cons.lock_mask & mgeom::LOCK_Z)) return false;

	static const int assign[8] = {4, 5, 6, 7, 0, 1, 2, 3}; // 0..3 bot, 4..7 top
	int n_above = -1, n_below = -1;

	for (int iel : inc) {
		octant_t *e = (octant_t *) sc_array_index(&mesh->elements, iel);
		for (int ino = 0; ino < 8; ino++) {
			if (e->nodes[assign[ino]].id == node) {
				if (ino < 4) {
					int top_id = e->nodes[assign[ino + 4]].id;
					if (top_id >= 0 && top_id < (int)coords.size()/3 && top_id != node) {
						n_above = top_id;
					}
				} else {
					int bot_id = e->nodes[assign[ino - 4]].id;
					if (bot_id >= 0 && bot_id < (int)coords.size()/3 && bot_id != node) {
						n_below = bot_id;
					}
				}
			}
		}
	}

	double px = coords[3*node+0], py = coords[3*node+1], pz = coords[3*node+2];
	double target_z = pz;
	if (n_above >= 0 && n_below >= 0) {
		target_z = 0.5 * (coords[3*n_above+2] + coords[3*n_below+2]);
	} else if (n_above >= 0) {
		double he = shortest_incident_edge(mesh, coords, inc, node);
		target_z = coords[3*n_above+2] - h_floor(he);
	} else if (n_below >= 0) {
		double he = shortest_incident_edge(mesh, coords, inc, node);
		target_z = coords[3*n_below+2] + h_floor(he);
	} else {
		return false;
	}

	double dz = target_z - pz;
	if (std::fabs(dz) < 1e-6 * g_h_ref) return false;

	NodeState s0 = node_state(mesh, coords, adj, inc, node, ref);
	bool moved = false;

	for (double alpha : {1.0, 0.75, 0.5, 0.25, 0.1, 0.05}) {
		coords[3*node+2] = pz + alpha * dz;
		NodeState s = node_state(mesh, coords, adj, inc, node, ref);
		if (better_state(s, s0)) {
			s0 = s; moved = true; break;
		}
		coords[3*node+2] = pz;
	}
	return moved;
}

// Face normal extrusion untangler: pushes interior face nodes away from the surface face along the surface normal
static bool untangle_element_face_extrusion(hexa_tree_t *mesh, std::vector<double> &coords,
                                            const std::vector<std::vector<int>> &adj,
                                            const std::vector<std::vector<int>> &inc,
                                            int iel, int ref,
                                            const std::vector<NodeConstraint> &cons) {
	octant_t *e = (octant_t *) sc_array_index(&mesh->elements, iel);
	static const int assign[8] = {4, 5, 6, 7, 0, 1, 2, 3}; // bot: 0,1,2,3; top: 4,5,6,7
	int bot[4] = {e->nodes[assign[0]].id, e->nodes[assign[1]].id, e->nodes[assign[2]].id, e->nodes[assign[3]].id};
	int top[4] = {e->nodes[assign[4]].id, e->nodes[assign[5]].id, e->nodes[assign[6]].id, e->nodes[assign[7]].id};

	for (int k = 0; k < 4; k++) {
		if (bot[k] < 0 || top[k] < 0 || bot[k] >= (int)coords.size()/3 || top[k] >= (int)coords.size()/3)
			return false;
	}

	// Test all 3 opposite face pairs: (bot, top) along Z, (front, back) along Y, (left, right) along X
	static const int pairs[3][2][4] = {
		{ {0,1,2,3}, {4,5,6,7} }, // Z pair
		{ {0,1,5,4}, {3,2,6,7} }, // Y pair
		{ {0,3,7,4}, {1,2,6,5} }  // X pair
	};

	bool moved = false;
	for (int p = 0; p < 3; p++) {
		int faceA[4], faceB[4];
		for (int k = 0; k < 4; k++) {
			faceA[k] = e->nodes[assign[pairs[p][0][k]]].id;
			faceB[k] = e->nodes[assign[pairs[p][1][k]]].id;
		}

		// Calculate outward normal of faceB relative to faceA
		double d1x = coords[3*faceB[2]+0] - coords[3*faceB[0]+0];
		double d1y = coords[3*faceB[2]+1] - coords[3*faceB[0]+1];
		double d1z = coords[3*faceB[2]+2] - coords[3*faceB[0]+2];

		double d2x = coords[3*faceB[3]+0] - coords[3*faceB[1]+0];
		double d2y = coords[3*faceB[3]+1] - coords[3*faceB[1]+1];
		double d2z = coords[3*faceB[3]+2] - coords[3*faceB[1]+2];

		double nx = d1y*d2z - d1z*d2y;
		double ny = d1z*d2x - d1x*d2z;
		double nz = d1x*d2y - d1y*d2x;
		double nlen = std::sqrt(nx*nx + ny*ny + nz*nz);
		if (nlen > 1e-6) { nx /= nlen; ny /= nlen; nz /= nlen; } else continue;

		// Orient normal from faceA centroid to faceB centroid
		double cAx=0, cAy=0, cAz=0, cBx=0, cBy=0, cBz=0;
		for (int k = 0; k < 4; k++) {
			cAx += coords[3*faceA[k]+0]; cAy += coords[3*faceA[k]+1]; cAz += coords[3*faceA[k]+2];
			cBx += coords[3*faceB[k]+0]; cBy += coords[3*faceB[k]+1]; cBz += coords[3*faceB[k]+2];
		}
		double cdx = cBx - cAx, cdy = cBy - cAy, cdz = cBz - cAz;
		if (cdx*nx + cdy*ny + cdz*nz < 0) { nx = -nx; ny = -ny; nz = -nz; }

		double he_e = h_floor(shortest_incident_edge(mesh, coords, inc[faceA[0]], faceA[0]));
		// the old tail {50.0, 150.0} was two absolute metre extrusions; the fractional
		// ladder already covers that range at any mesh scale.
		double H_TESTS[] = {1.0 * he_e, 0.5 * he_e, 0.25 * he_e, 0.1 * he_e, 0.05 * he_e};

		// Try adjusting movable nodes of faceA away from faceB
		for (int k = 0; k < 4; k++) {
			int na = faceA[k];
			if (cons[na].gts_surface_id >= 0) continue;
			int nb = faceB[k];
			double px = coords[3*na+0], py = coords[3*na+1], pz = coords[3*na+2];
			NodeState s0 = node_state(mesh, coords, adj, inc[na], na, ref);

			for (double h_ref : H_TESTS) {
				double tx = coords[3*nb+0] - nx * h_ref;
				double ty = coords[3*nb+1] - ny * h_ref;
				double tz = coords[3*nb+2] - nz * h_ref;

				double dx = tx - px, dy = ty - py, dz = tz - pz;
				mgeom::apply_lock(cons[na].lock_mask, dx, dy, dz);

				for (double alpha : {1.0, 0.75, 0.5, 0.25, 0.1}) {
					coords[3*na+0] = px + alpha * dx;
					coords[3*na+1] = py + alpha * dy;
					coords[3*na+2] = pz + alpha * dz;
					NodeState s = node_state(mesh, coords, adj, inc[na], na, ref);
					if (better_state(s, s0, he_e)) {
						s0 = s; moved = true; break;
					}
					coords[3*na+0] = px; coords[3*na+1] = py; coords[3*na+2] = pz;
				}
				if (moved) break;
			}
		}

		// Try adjusting movable nodes of faceB away from faceA
		for (int k = 0; k < 4; k++) {
			int nb = faceB[k];
			if (cons[nb].gts_surface_id >= 0) continue;
			int na = faceA[k];
			double px = coords[3*nb+0], py = coords[3*nb+1], pz = coords[3*nb+2];
			NodeState s0 = node_state(mesh, coords, adj, inc[nb], nb, ref);

			for (double h_ref : H_TESTS) {
				double tx = coords[3*na+0] + nx * h_ref;
				double ty = coords[3*na+1] + ny * h_ref;
				double tz = coords[3*na+2] + nz * h_ref;

				double dx = tx - px, dy = ty - py, dz = tz - pz;
				mgeom::apply_lock(cons[nb].lock_mask, dx, dy, dz);

				for (double alpha : {1.0, 0.75, 0.5, 0.25, 0.1}) {
					coords[3*nb+0] = px + alpha * dx;
					coords[3*nb+1] = py + alpha * dy;
					coords[3*nb+2] = pz + alpha * dz;
					NodeState s = node_state(mesh, coords, adj, inc[nb], nb, ref);
					if (better_state(s, s0, he_e)) {
						s0 = s; moved = true; break;
					}
					coords[3*nb+0] = px; coords[3*nb+1] = py; coords[3*nb+2] = pz;
				}
				if (moved) break;
			}
		}
	}
	return moved;
}

// Untangle inverted elements by fixing vertical column orientation between bottom and top faces
static bool untangle_element_columns(hexa_tree_t *mesh, std::vector<double> &coords,
                                     const std::vector<std::vector<int>> &adj,
                                     const std::vector<std::vector<int>> &inc,
                                     int iel, int ref,
                                     const std::vector<NodeConstraint> &cons) {
	octant_t *e = (octant_t *) sc_array_index(&mesh->elements, iel);
	static const int assign[8] = {4, 5, 6, 7, 0, 1, 2, 3}; // bottom: 0,1,2,3; top: 4,5,6,7
	int bot[4] = {e->nodes[assign[0]].id, e->nodes[assign[1]].id, e->nodes[assign[2]].id, e->nodes[assign[3]].id};
	int top[4] = {e->nodes[assign[4]].id, e->nodes[assign[5]].id, e->nodes[assign[6]].id, e->nodes[assign[7]].id};

	bool moved = false;
	for (int k = 0; k < 4; k++) {
		int nb = bot[k], nt = top[k];
		if (nb < 0 || nt < 0 || nb >= (int)coords.size()/3 || nt >= (int)coords.size()/3) continue;

		double he_b = shortest_incident_edge(mesh, coords, inc[nb], nb);
		double he_t = shortest_incident_edge(mesh, coords, inc[nt], nt);
		double h_target = std::max(h_floor(he_b), h_floor(he_t));

		// If bottom node is higher than top node (inverted column in Z):
		if (coords[3*nb+2] >= coords[3*nt+2] - 1e-3 * h_target) {

			// Try moving bottom node below top node
			if (cons[nb].gts_surface_id < 0 && !(cons[nb].lock_mask & mgeom::LOCK_Z)) {
				double target_z = coords[3*nt+2] - h_target;
				double pz = coords[3*nb+2];
				NodeState s0 = node_state(mesh, coords, adj, inc[nb], nb, ref);
				for (double fraction : {1.0, 0.75, 0.5, 0.25, 0.1}) {
					coords[3*nb+2] = pz + fraction * (target_z - pz);
					NodeState s = node_state(mesh, coords, adj, inc[nb], nb, ref);
					if (better_state(s, s0, h_target)) {
						s0 = s; moved = true; break;
					}
					coords[3*nb+2] = pz;
				}
			}
			// Or try moving top node above bottom node
			if (cons[nt].gts_surface_id < 0 && !(cons[nt].lock_mask & mgeom::LOCK_Z)) {
				double target_z = coords[3*nb+2] + h_target;
				double pz = coords[3*nt+2];
				NodeState s0 = node_state(mesh, coords, adj, inc[nt], nt, ref);
				for (double fraction : {1.0, 0.75, 0.5, 0.25, 0.1}) {
					coords[3*nt+2] = pz + fraction * (target_z - pz);
					NodeState s = node_state(mesh, coords, adj, inc[nt], nt, ref);
					if (better_state(s, s0, h_target)) {
						s0 = s; moved = true; break;
					}
					coords[3*nt+2] = pz;
				}
			}
		}
	}
	return moved;
}

// Simultaneous Patch Untangler: moves all movable vertices of an inverted element
// jointly toward their ideal affine parallelepiped positions
static bool untangle_inverted_element_patch(hexa_tree_t *mesh, std::vector<double> &coords,
                                            const std::vector<std::vector<int>> &adj,
                                            const std::vector<std::vector<int>> &inc,
                                            int iel, int ref,
                                            const std::vector<NodeConstraint> &cons) {
	octant_t *e = (octant_t *) sc_array_index(&mesh->elements, iel);
	static const int assign[8] = {4, 5, 6, 7, 0, 1, 2, 3};
	int n[8];
	for (int i = 0; i < 8; i++) {
		n[i] = e->nodes[assign[i]].id;
		if (n[i] < 0 || n[i] >= (int)coords.size()/3) return false;
	}

	// Compute patch initial state
	std::unordered_set<int> patch_elems;
	for (int i = 0; i < 8; i++) {
		for (int ie : inc[n[i]]) patch_elems.insert(ie);
	}

	// Evaluate patch initial state
	auto eval_patch_state = [&]() {
		int n_inv = 0;
		double sum_neg = 0.0;
		double min_sj = 1e300;
		for (int pe : patch_elems) {
			double X[8], Y[8], Z[8];
			load_elem_xyz(mesh, coords, pe, X, Y, Z);
			double sj = mgeom::hex_min_corner_sj(X, Y, Z);
			if (sj <= 0.0) {
				n_inv++;
				sum_neg += sj;
			}
			if (sj < min_sj) min_sj = sj;
		}
		return std::make_tuple(n_inv, sum_neg, min_sj);
	};

	auto [n_inv0, sum_neg0, min_sj0] = eval_patch_state();
	if (n_inv0 == 0) return false;

	double orig_pos[8][3];
	bool movable[8];
	for (int ino = 0; ino < 8; ino++) {
		int nid = n[ino];
		orig_pos[ino][0] = coords[3*nid+0];
		orig_pos[ino][1] = coords[3*nid+1];
		orig_pos[ino][2] = coords[3*nid+2];
		movable[ino] = (cons[nid].gts_surface_id < 0 && cons[nid].lock_mask != (mgeom::LOCK_X|mgeom::LOCK_Y|mgeom::LOCK_Z));
	}

	double he_patch = h_floor(shortest_incident_edge(mesh, coords, inc[n[0]], n[0]));

	// Generate target sets: Candidate 1 = Affine Parallelepiped, Candidate 2 = Bottom-to-Top Extrusion, Candidate 3 = Top-to-Bottom Extrusion
	std::vector<std::vector<std::array<double,3>>> candidate_targets;

	// Target Candidate 1: 7-term affine parallelepiped reconstruction
	{
		std::vector<std::array<double,3>> t(8);
		for (int ino = 0; ino < 8; ino++) {
			switch (ino) {
				case 0: for(int d=0; d<3; d++) t[0][d] = coords[3*n[1]+d] + coords[3*n[3]+d] + coords[3*n[4]+d] - (coords[3*n[2]+d] + coords[3*n[5]+d] + coords[3*n[7]+d]) + coords[3*n[6]+d]; break;
				case 1: for(int d=0; d<3; d++) t[1][d] = coords[3*n[0]+d] + coords[3*n[2]+d] + coords[3*n[5]+d] - (coords[3*n[3]+d] + coords[3*n[4]+d] + coords[3*n[6]+d]) + coords[3*n[7]+d]; break;
				case 2: for(int d=0; d<3; d++) t[2][d] = coords[3*n[1]+d] + coords[3*n[3]+d] + coords[3*n[6]+d] - (coords[3*n[0]+d] + coords[3*n[5]+d] + coords[3*n[7]+d]) + coords[3*n[4]+d]; break;
				case 3: for(int d=0; d<3; d++) t[3][d] = coords[3*n[0]+d] + coords[3*n[2]+d] + coords[3*n[7]+d] - (coords[3*n[1]+d] + coords[3*n[4]+d] + coords[3*n[6]+d]) + coords[3*n[5]+d]; break;
				case 4: for(int d=0; d<3; d++) t[4][d] = coords[3*n[0]+d] + coords[3*n[5]+d] + coords[3*n[7]+d] - (coords[3*n[1]+d] + coords[3*n[3]+d] + coords[3*n[6]+d]) + coords[3*n[2]+d]; break;
				case 5: for(int d=0; d<3; d++) t[5][d] = coords[3*n[1]+d] + coords[3*n[4]+d] + coords[3*n[6]+d] - (coords[3*n[0]+d] + coords[3*n[2]+d] + coords[3*n[7]+d]) + coords[3*n[3]+d]; break;
				case 6: for(int d=0; d<3; d++) t[6][d] = coords[3*n[2]+d] + coords[3*n[5]+d] + coords[3*n[7]+d] - (coords[3*n[1]+d] + coords[3*n[3]+d] + coords[3*n[4]+d]) + coords[3*n[0]+d]; break;
				case 7: for(int d=0; d<3; d++) t[7][d] = coords[3*n[3]+d] + coords[3*n[4]+d] + coords[3*n[6]+d] - (coords[3*n[0]+d] + coords[3*n[2]+d] + coords[3*n[5]+d]) + coords[3*n[1]+d]; break;
			}
		}
		candidate_targets.push_back(t);
	}

	// Target Candidate 2: Extrude top face away from bottom face along bottom normal
	{
		std::vector<std::array<double,3>> t(8);
		for(int i=0; i<8; i++) for(int d=0; d<3; d++) t[i][d] = orig_pos[i][d];
		double d1x = orig_pos[2][0] - orig_pos[0][0], d1y = orig_pos[2][1] - orig_pos[0][1], d1z = orig_pos[2][2] - orig_pos[0][2];
		double d2x = orig_pos[3][0] - orig_pos[1][0], d2y = orig_pos[3][1] - orig_pos[1][1], d2z = orig_pos[3][2] - orig_pos[1][2];
		double nx = d1y*d2z - d1z*d2y, ny = d1z*d2x - d1x*d2z, nz = d1x*d2y - d1y*d2x;
		double nl = std::sqrt(nx*nx + ny*ny + nz*nz);
		if (nl > 1e-6) {
			nx /= nl; ny /= nl; nz /= nl;
			if (nz < 0) { nx = -nx; ny = -ny; nz = -nz; }
			for (int k = 0; k < 4; k++) {
				t[k+4][0] = orig_pos[k][0] + nx * he_patch;
				t[k+4][1] = orig_pos[k][1] + ny * he_patch;
				t[k+4][2] = orig_pos[k][2] + nz * he_patch;
			}
			candidate_targets.push_back(t);
		}
	}

	// Target Candidate 3: Extrude bottom face away from top face along top normal
	{
		std::vector<std::array<double,3>> t(8);
		for(int i=0; i<8; i++) for(int d=0; d<3; d++) t[i][d] = orig_pos[i][d];
		double d1x = orig_pos[6][0] - orig_pos[4][0], d1y = orig_pos[6][1] - orig_pos[4][1], d1z = orig_pos[6][2] - orig_pos[4][2];
		double d2x = orig_pos[7][0] - orig_pos[5][0], d2y = orig_pos[7][1] - orig_pos[5][1], d2z = orig_pos[7][2] - orig_pos[5][2];
		double nx = d1y*d2z - d1z*d2y, ny = d1z*d2x - d1x*d2z, nz = d1x*d2y - d1y*d2x;
		double nl = std::sqrt(nx*nx + ny*ny + nz*nz);
		if (nl > 1e-6) {
			nx /= nl; ny /= nl; nz /= nl;
			if (nz < 0) { nx = -nx; ny = -ny; nz = -nz; }
			for (int k = 0; k < 4; k++) {
				t[k][0] = orig_pos[k+4][0] - nx * he_patch;
				t[k][1] = orig_pos[k+4][1] - ny * he_patch;
				t[k][2] = orig_pos[k+4][2] - nz * he_patch;
			}
			candidate_targets.push_back(t);
		}
	}

	// Test each candidate target set with step sizes
	for (auto &t_set : candidate_targets) {
		for (int ino = 0; ino < 8; ino++) {
			int nid = n[ino];
			double dx = t_set[ino][0] - orig_pos[ino][0];
			double dy = t_set[ino][1] - orig_pos[ino][1];
			double dz = t_set[ino][2] - orig_pos[ino][2];
			mgeom::apply_lock(cons[nid].lock_mask, dx, dy, dz);
			t_set[ino][0] = orig_pos[ino][0] + dx;
			t_set[ino][1] = orig_pos[ino][1] + dy;
			t_set[ino][2] = orig_pos[ino][2] + dz;
		}

		for (double alpha : {1.0, 0.8, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1, 0.05}) {
			for (int ino = 0; ino < 8; ino++) {
				if (!movable[ino]) continue;
				int nid = n[ino];
				coords[3*nid+0] = orig_pos[ino][0] + alpha * (t_set[ino][0] - orig_pos[ino][0]);
				coords[3*nid+1] = orig_pos[ino][1] + alpha * (t_set[ino][1] - orig_pos[ino][1]);
				coords[3*nid+2] = orig_pos[ino][2] + alpha * (t_set[ino][2] - orig_pos[ino][2]);
			}

			auto [n_inv, sum_neg, min_sj] = eval_patch_state();
			if (n_inv < n_inv0 || (n_inv == n_inv0 && sum_neg > sum_neg0 + 1e-4) || (n_inv == 0 && min_sj > min_sj0)) {
				return true; // accepted joint move!
			}
		}

		// Revert before testing next candidate
		for (int ino = 0; ino < 8; ino++) {
			int nid = n[ino];
			coords[3*nid+0] = orig_pos[ino][0];
			coords[3*nid+1] = orig_pos[ino][1];
			coords[3*nid+2] = orig_pos[ino][2];
		}
	}

	return false;
}

static std::vector<std::vector<int>> build_adjacency(hexa_tree_t *mesh, int n_nodes) {
	static const int E[12][2] = {
		{0,1},{1,2},{2,3},{3,0},{4,5},{5,6},{6,7},{7,4},{0,4},{1,5},{2,6},{3,7}
	};
	std::vector<std::unordered_set<int>> tmp(n_nodes);
	for (int iel = 0; iel < mesh->elements.elem_count; iel++) {
		octant_t *e = (octant_t *) sc_array_index(&mesh->elements, iel);
		for (int k = 0; k < 12; k++) {
			int a = e->nodes[E[k][0]].id, b = e->nodes[E[k][1]].id;
			if (a >= 0 && a < n_nodes && b >= 0 && b < n_nodes && a != b) {
				tmp[a].insert(b);
				tmp[b].insert(a);
			}
		}
	}
	std::vector<std::vector<int>> adj(n_nodes);
	for (int i = 0; i < n_nodes; i++) adj[i].assign(tmp[i].begin(), tmp[i].end());
	return adj;
}

// Move the node toward its edge-neighbour centroid (Laplacian), accepted only if
// the node's incident state improves (monotone).
static bool relax_toward_centroid(hexa_tree_t *mesh, std::vector<double> &coords,
                                  const std::vector<std::vector<int>> &adj,
                                  const std::vector<int> &inc, int node, int ref, uint8_t mask, double he = 0.0) {
	if (adj[node].empty()) return false;
	double cx=0, cy=0, cz=0;
	for (int nb : adj[node]) { cx+=coords[3*nb]; cy+=coords[3*nb+1]; cz+=coords[3*nb+2]; }
	double invn = 1.0/(double)adj[node].size();
	cx*=invn; cy*=invn; cz*=invn;
	double px=coords[3*node], py=coords[3*node+1], pz=coords[3*node+2];
	double dx=cx-px, dy=cy-py, dz=cz-pz;
	mgeom::apply_lock(mask, dx, dy, dz);
	if (std::fabs(dx)+std::fabs(dy)+std::fabs(dz) < 1e-12) return false;
	NodeState s0 = node_state(mesh, coords, adj, inc, node, ref);
	double alpha = 1.0;
	while (alpha > 1.0/64.0) {
		coords[3*node]=px+alpha*dx; coords[3*node+1]=py+alpha*dy; coords[3*node+2]=pz+alpha*dz;
		NodeState s = node_state(mesh, coords, adj, inc, node, ref);
		if (better_state(s, s0, he)) return true;
		alpha *= 0.5;
	}
	coords[3*node]=px; coords[3*node+1]=py; coords[3*node+2]=pz;
	return false;
}

// Level 1 (of curves -> faces -> volume): 1D Laplacian on the 12 edges of the domain box.
// An edge node is one locked on two axes (lock_mask), so it can only slide along the
// third; it moves to the mean of its neighbours on the same box edge. Fixed: box corners,
// material-interface nodes, GTS surface nodes and the free (top) surface -- the top four
// box edges lie on the topography and stay put. A move is kept only if every incident
// element stays valid over its whole volume (pillowgeom::valid), halving it otherwise.
static int smooth_box_edges_1d(hexa_tree_t *mesh, std::vector<double> &coords,
                               const std::vector<NodeConstraint> &cons,
                               const std::vector<int> &material_fixed_nodes) {
	const int nn = mesh->nodes.elem_count;
	auto inc = build_incidence(mesh, nn);
	auto adj = build_adjacency(mesh, nn);
	auto bits = [](uint8_t m) { return ((m & mgeom::LOCK_X) != 0) + ((m & mgeom::LOCK_Y) != 0) + ((m & mgeom::LOCK_Z) != 0); };
	std::vector<char> fixed(nn, 0);
	for (int n : material_fixed_nodes) if (n >= 0 && n < nn) fixed[n] = 1;
	for (int n = 0; n < nn; n++) if (cons[n].gts_surface_id >= 0 || bits(cons[n].lock_mask) == 3) fixed[n] = 1;
	// open faces that are not on a wall plane are the free surface
	std::unordered_map<face_key, int, face_key_hash> count;
	for (int iel = 0; iel < mesh->elements.elem_count; iel++) {
		octant_t *e = (octant_t *) sc_array_index(&mesh->elements, iel);
		for (int f = 0; f < 6; f++)
			count[make_face_key(e->nodes[FaceNodesMap[f][0]].id, e->nodes[FaceNodesMap[f][1]].id,
			                    e->nodes[FaceNodesMap[f][2]].id, e->nodes[FaceNodesMap[f][3]].id)]++;
	}
	for (int iel = 0; iel < mesh->elements.elem_count; iel++) {
		octant_t *e = (octant_t *) sc_array_index(&mesh->elements, iel);
		for (int f = 0; f < 6; f++) {
			int id[4]; uint8_t wall = 0xFF;
			for (int k = 0; k < 4; k++) { id[k] = e->nodes[FaceNodesMap[f][k]].id; wall &= cons[id[k]].lock_mask; }
			if (count[make_face_key(id[0], id[1], id[2], id[3])] == 1 && !wall)
				for (int k = 0; k < 4; k++) fixed[id[k]] = 1;
		}
	}
	// free edge nodes and their neighbours along the same box edge
	std::vector<int> nodes;
	std::vector<std::vector<int>> nbr(nn);
	int n_edge = 0, n_transition = 0;
	for (int n = 0; n < nn; n++) {
		const uint8_t m = cons[n].lock_mask;
		if (bits(m) != 2) continue;
		n_edge++;
		if (fixed[n]) continue;
		for (int j : adj[n]) if ((cons[j].lock_mask & m) == m) nbr[n].push_back(j);
		if (nbr[n].size() != 2) continue;
		// A jump in element size along the edge is a zcuts transition: keep it fixed so the
		// smoothing only acts between consecutive transitions and the vertical grading survives.
		// ponytail: 20% size jump as the transition test, read zcuts/lattice z if this misfires
		double L[2];
		for (int k = 0; k < 2; k++) {
			const int j = nbr[n][k];
			L[k] = std::sqrt(std::pow(coords[3*j]-coords[3*n],2) + std::pow(coords[3*j+1]-coords[3*n+1],2) + std::pow(coords[3*j+2]-coords[3*n+2],2));
		}
		if (std::max(L[0], L[1]) > 1.2*std::min(L[0], L[1])) { n_transition++; continue; }
		nodes.push_back(n);
	}
	auto valid_around = [&](int n) {
		for (int iel : inc[n]) {
			octant_t *e = (octant_t *) sc_array_index(&mesh->elements, iel);
			pillowgeom::Hex h;
			for (int k = 0; k < 8; k++) for (int d = 0; d < 3; d++) h[k][d] = coords[3*e->nodes[mgeom::H5_ORD[k]].id + d];
			if (!pillowgeom::valid(h, 1)) return false;
		}
		return true;
	};
	const int SWEEPS = 5000;  // ponytail: plain Gauss-Seidel on ~100 nodes, cheap; stops when the largest move is below 1 mm
	int total = 0;
	for (int sweep = 0; sweep < SWEEPS; sweep++) {
		int moved = 0; double dmax = 0;
		for (int n : nodes) {
			double d[3];
			for (int k = 0; k < 3; k++) d[k] = 0.5*(coords[3*nbr[n][0]+k] + coords[3*nbr[n][1]+k]) - coords[3*n+k];
			mgeom::apply_lock(cons[n].lock_mask, d[0], d[1], d[2]);
			const double p[3] = {coords[3*n], coords[3*n+1], coords[3*n+2]};
			for (double a = 1.0; a > 1.0/64.0; a *= 0.5) {
				for (int k = 0; k < 3; k++) coords[3*n+k] = p[k] + a*d[k];
				if (valid_around(n)) { moved++; dmax = std::max(dmax, a*std::sqrt(d[0]*d[0]+d[1]*d[1]+d[2]*d[2])); break; }
				for (int k = 0; k < 3; k++) coords[3*n+k] = p[k];
			}
		}
		total += moved;
		if (sweep % 500 == 0 || dmax < 1e-3)
			printf("      Level 1 sweep %d: %d moved, max move %.3f m\n", sweep, moved, dmax);
		if (dmax < 1e-3) break;
	}
	printf("    Level 1: %d box-edge nodes, %zu free, %d zcuts transitions fixed, %d other fixed\n",
	       n_edge, nodes.size(), n_transition, n_edge - n_transition - (int)nodes.size());
	return total;
}

// Parallelepiped geometric reconstruction: reconstructs an inverted vertex from the remaining 7 vertices of the element
static bool reconstruct_hex_vertex(hexa_tree_t *mesh, std::vector<double> &coords,
                                    const std::vector<std::vector<int>> &adj,
                                    const std::vector<std::vector<int>> &inc,
                                    int iel, int ref,
                                    const std::vector<NodeConstraint> &cons) {
	octant_t *e = (octant_t *) sc_array_index(&mesh->elements, iel);
	static const int assign[8] = {4, 5, 6, 7, 0, 1, 2, 3}; // bottom 0,1,2,3; top 4,5,6,7
	int n[8];
	for (int i = 0; i < 8; i++) n[i] = e->nodes[assign[i]].id;

	for (int i = 0; i < 8; i++) {
		if (n[i] < 0 || n[i] >= (int)coords.size()/3) return false;
	}

	bool moved = false;
	for (int ino = 0; ino < 8; ino++) {
		int target_node = n[ino];
		if (cons[target_node].gts_surface_id >= 0) continue; // preserve GTS surface nodes

		double target_pos[3];
		// Exact linear parallelepiped affine identity from the other 7 vertices
		switch (ino) {
			case 0: for(int d=0; d<3; d++) target_pos[d] = coords[3*n[1]+d] + coords[3*n[3]+d] + coords[3*n[4]+d] - (coords[3*n[2]+d] + coords[3*n[5]+d] + coords[3*n[7]+d]) + coords[3*n[6]+d]; break;
			case 1: for(int d=0; d<3; d++) target_pos[d] = coords[3*n[0]+d] + coords[3*n[2]+d] + coords[3*n[5]+d] - (coords[3*n[3]+d] + coords[3*n[4]+d] + coords[3*n[6]+d]) + coords[3*n[7]+d]; break;
			case 2: for(int d=0; d<3; d++) target_pos[d] = coords[3*n[1]+d] + coords[3*n[3]+d] + coords[3*n[6]+d] - (coords[3*n[0]+d] + coords[3*n[5]+d] + coords[3*n[7]+d]) + coords[3*n[4]+d]; break;
			case 3: for(int d=0; d<3; d++) target_pos[d] = coords[3*n[0]+d] + coords[3*n[2]+d] + coords[3*n[7]+d] - (coords[3*n[1]+d] + coords[3*n[4]+d] + coords[3*n[6]+d]) + coords[3*n[5]+d]; break;
			case 4: for(int d=0; d<3; d++) target_pos[d] = coords[3*n[0]+d] + coords[3*n[5]+d] + coords[3*n[7]+d] - (coords[3*n[1]+d] + coords[3*n[3]+d] + coords[3*n[6]+d]) + coords[3*n[2]+d]; break;
			case 5: for(int d=0; d<3; d++) target_pos[d] = coords[3*n[1]+d] + coords[3*n[4]+d] + coords[3*n[6]+d] - (coords[3*n[0]+d] + coords[3*n[2]+d] + coords[3*n[7]+d]) + coords[3*n[3]+d]; break;
			case 6: for(int d=0; d<3; d++) target_pos[d] = coords[3*n[2]+d] + coords[3*n[5]+d] + coords[3*n[7]+d] - (coords[3*n[1]+d] + coords[3*n[3]+d] + coords[3*n[4]+d]) + coords[3*n[0]+d]; break;
			case 7: for(int d=0; d<3; d++) target_pos[d] = coords[3*n[3]+d] + coords[3*n[4]+d] + coords[3*n[6]+d] - (coords[3*n[0]+d] + coords[3*n[2]+d] + coords[3*n[5]+d]) + coords[3*n[1]+d]; break;
		}

		double px = coords[3*target_node+0], py = coords[3*target_node+1], pz = coords[3*target_node+2];
		double dx = target_pos[0] - px, dy = target_pos[1] - py, dz = target_pos[2] - pz;
		mgeom::apply_lock(cons[target_node].lock_mask, dx, dy, dz);

		NodeState s0 = node_state(mesh, coords, adj, inc[target_node], target_node, ref);
		for (double beta : {1.0, 0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1, 0.05, 0.01}) {
			coords[3*target_node+0] = px + beta * dx;
			coords[3*target_node+1] = py + beta * dy;
			coords[3*target_node+2] = pz + beta * dz;

			NodeState s = node_state(mesh, coords, adj, inc[target_node], target_node, ref);
			if (better_state(s, s0)) {
				s0 = s; moved = true; break;
			}
			coords[3*target_node+0] = px; coords[3*target_node+1] = py; coords[3*target_node+2] = pz;
		}
	}
	return moved;
}

// Drive untangling until 0 inverted or MAX_UNTANGLE_ITERS. Returns remaining count.
int untangle_inversions(hexa_tree_t *mesh, std::vector<double> &coords,
                        const std::vector<NodeConstraint> &cons,
                        const std::vector<uint8_t> &wall_lock, int ref) {
	set_mesh_scale(mesh, coords);      // safe to call standalone
	build_face_neighbours(mesh);
	int nn = mesh->nodes.elem_count;
	auto inc = build_incidence(mesh, nn);
	auto adj = build_adjacency(mesh, nn);
	int prev_count = -1, prev_inv = 0;
	int stall_streak = 0;
	const int STALL_PATIENCE_LIMIT = 200;

	int fold0 = count_folded_elements(mesh, coords);
	printf("    Untangler: start %d folded elements\n", fold0);

	for (int iter = 0; iter < MAX_UNTANGLE_ITERS; iter++) {
		MeshAnalysis a = analyze_mesh(mesh, coords);
		int n_fold_now = count_folded_elements(mesh, coords);
		if (a.n_inverted == 0 && n_fold_now == 0) {
			printf("    Untangler: 0 inverted, 0 folded after %d iters\n", iter);
			return 0;
		}
		if (iter % 10 == 0)
			printf("      untangle iter %d: %d inverted, %d folded (stall %d)\n",
			       iter, a.n_inverted, n_fold_now, stall_streak);

		// Sweep-level global barrier. Each relaxer accepts on its own node's incident
		// elements only, so a move that helps this node can break the next one -- the
		// sweep as a whole could end worse than it started, and did (folds went up while
		// inversions went down). Snapshot, run the sweep, then keep it only if the GLOBAL
		// broken count improved, backing off along the whole displacement field if not.
		// ponytail: sweep-level barrier, not a true Jacobi solve. Upgrade to per-node
		// frozen-coordinate proposals if intra-sweep thrashing still shows up.
		std::vector<double> sweep_start = coords;
		MeshBad bad_before{a.n_inverted, n_fold_now};

		// A stalled iteration = neither the inversion count nor, at equal inversions, the
		// fold count dropped.
		bool improved = (prev_count < 0) ||
		                (a.n_inverted < prev_inv) ||
		                (a.n_inverted == prev_inv && n_fold_now < prev_count);
		if (improved) stall_streak = 0; else stall_streak++;
		prev_inv = a.n_inverted;
		prev_count = n_fold_now;
		bool stalled = (stall_streak > 0);

		// nodes incident to any currently-inverted element
		std::unordered_set<int> nodes;
		for (int iel : a.inverted_ids) {
			octant_t *e = (octant_t *) sc_array_index(&mesh->elements, iel);
			for (int ino = 0; ino < 8; ino++) nodes.insert(e->nodes[ino].id);
		}

		// Multi-ring expansion when stalled to allow the surrounding cluster to breathe
		if (stalled) {
			std::vector<int> cur_nodes(nodes.begin(), nodes.end());
			for (int n : cur_nodes) {
				for (int nb : adj[n]) {
					nodes.insert(nb);
					if (stall_streak >= 15) {
						for (int nb2 : adj[nb]) {
							nodes.insert(nb2);
							if (stall_streak >= 40) {
								for (int nb3 : adj[nb2]) nodes.insert(nb3);
							}
						}
					}
				}
			}
		}
		bool any_moved = false;

		// Re-enabled in Fase C: the steps are now proportional to the local element
		// (Fase B) and every acceptance test counts folds (better_state), so these can
		// no longer trade an inversion for an overlap.
		for (int iel : a.inverted_ids) {
			if (untangle_inverted_element_patch(mesh, coords, adj, inc, iel, ref, cons)) any_moved = true;
			if (untangle_element_face_extrusion(mesh, coords, adj, inc, iel, ref, cons)) any_moved = true;
			if (untangle_element_columns(mesh, coords, adj, inc, iel, ref, cons)) any_moved = true;
			if (reconstruct_hex_vertex(mesh, coords, adj, inc, iel, ref, cons)) any_moved = true;
		}

		for (int node : nodes) {
			if (inc[node].empty()) continue;
			double he = shortest_incident_edge(mesh, coords, inc[node], node);
			if (cons[node].gts_surface_id >= 0) {
				// Surface nodes relax ONLY tangentially along their own GTS surface.
				if (relax_gts_surface_node(mesh, coords, adj, inc[node], node, ref, cons[node], he)) any_moved = true;
			} else {
				if (relax_toward_centroid(mesh, coords, adj, inc[node], node, ref, cons[node].lock_mask, he))
					any_moved = true;
				if (relax_vertical_column_midpoint(mesh, coords, adj, inc[node], node, ref, cons[node])) any_moved = true;
				if (relax_pillow_node(mesh, coords, adj, inc[node], node, ref, cons, he)) any_moved = true;
				if (relax_volume_gradient(mesh, coords, adj, inc[node], node, ref, cons[node], he)) any_moved = true;
				double step0 = UNTANGLE_STEP0 * h_floor(he);
				if (relax_node(mesh, coords, adj, inc[node], node, ref, cons[node].lock_mask,
				               step0, UNTANGLE_CAP * h_floor(he), he)) any_moved = true;
				if (relax_node_escobar_barrier(mesh, coords, adj, inc[node], node, ref, cons[node], he)) any_moved = true;
				if (stalled && probe_node_3d_box(mesh, coords, adj, inc[node], node, ref, cons[node], he)) any_moved = true;
			}
		}
		(void)stalled;

		(void)any_moved;

		// Accept the sweep only if it improved the global broken count; otherwise damp the
		// whole displacement field and retest, and revert outright if nothing helps.
		{
			std::vector<double> sweep_end = coords;
			MeshBad bad_after{analyze_mesh(mesh, coords).n_inverted, count_folded_elements(mesh, coords)};
			if (bad_worse(bad_after, bad_before)) {
				bool kept = false;
				for (double frac : {0.5, 0.25}) {
					for (size_t k = 0; k < coords.size(); k++)
						coords[k] = sweep_start[k] + frac * (sweep_end[k] - sweep_start[k]);
					MeshBad d{analyze_mesh(mesh, coords).n_inverted, count_folded_elements(mesh, coords)};
					if (!bad_worse(d, bad_before)) { kept = true; break; }
				}
				if (!kept) coords = sweep_start;
			}
		}

		if (stall_streak >= STALL_PATIENCE_LIMIT) {
			printf("    Untangler: stalled with %d inverted after %d iters (%d stalled)\n", a.n_inverted, iter, stall_streak);
			break;
		}
	}
	// --- Escalation: unpin GTS surface nodes, as a last resort -----------------
	// A residual inversion is usually not a mesh problem at all: it is a spike in the INPUT
	// surface (the coastline sentinel fill shows up as a lone node hundreds of metres above
	// its own face neighbours), and a node pinned to that surface can only slide tangentially
	// along the same spike -- it can never come down, so no relaxer can fix the element.
	// An inverted element stops the wave-propagation solver from starting at all, so trading
	// surface conformity for validity is the right call. Bounded, though: at most
	// ESCALATION_CAP of the local edge, and only for nodes of elements that are STILL
	// inverted after everything else has run.
	if (ESCALATION_ENABLED) {
		std::vector<double> pre_esc = coords;
		std::unordered_set<int> escalated;
		for (int round = 0; round < 50; round++) {
			MeshAnalysis ae = analyze_mesh(mesh, coords);
			if (ae.n_inverted == 0) break;

			std::unordered_set<int> surf_nodes;
			for (int iel : ae.inverted_ids) {
				octant_t *e = (octant_t *) sc_array_index(&mesh->elements, iel);
				for (int ino = 0; ino < 8; ino++) {
					int nid = e->nodes[ino].id;
					if (nid >= 0 && nid < nn && cons[nid].gts_surface_id >= 0) surf_nodes.insert(nid);
				}
			}
			if (surf_nodes.empty()) break;

			int moved = 0;
			for (int node : surf_nodes) {
				if (inc[node].empty()) continue;
				double he = shortest_incident_edge(mesh, coords, inc[node], node);
				double cap = ESCALATION_CAP * h_floor(he);
				// wall locks are still honoured -- cons[node].lock_mask goes straight through
				if (relax_node(mesh, coords, adj, inc[node], node, ref,
				               cons[node].lock_mask, cap, cap, he)) {
					escalated.insert(node);
					moved++;
				}
			}
			if (moved == 0) break;
		}
		if (!escalated.empty()) {
			double max_dev = 0.0;
			for (int node : escalated) {
				double dx = coords[3*node+0]-pre_esc[3*node+0];
				double dy = coords[3*node+1]-pre_esc[3*node+1];
				double dz = coords[3*node+2]-pre_esc[3*node+2];
				double d = std::sqrt(dx*dx+dy*dy+dz*dz);
				if (d > max_dev) max_dev = d;
			}
			printf("    Untangler: escalation moved %zu GTS surface nodes off their surface "
			       "(max %.3f m, cap %.0f%% of local edge)\n",
			       escalated.size(), max_dev, 100.0*ESCALATION_CAP);
		}
	}

	MeshAnalysis f = analyze_mesh(mesh, coords);
	printf("    Untangler: final check: %d inverted, %d folded (folds %d -> %d)\n",
	       f.n_inverted, count_folded_elements(mesh, coords), fold0, count_folded_elements(mesh, coords));
	for (int iel : f.inverted_ids) {
		octant_t *e = (octant_t *) sc_array_index(&mesh->elements, iel);
		printf("      INV ELEM %d: level=%d, nodes=[", iel, e->level);
		for (int ino = 0; ino < 8; ino++) {
			int nid = e->nodes[ino].id;
			printf("%d(gts=%d,lock=%d,z=%.1f) ", nid, cons[nid].gts_surface_id, cons[nid].lock_mask, coords[3*nid+2]);
		}
		printf("]\n");
	}
	return f.n_inverted;
}

static int find_set_root(std::vector<int> &parent, int i) {
	if (parent[i] == i) return i;
	return parent[i] = find_set_root(parent, parent[i]);
}

static void union_node_pair(std::vector<int> &parent, const std::vector<NodeConstraint> &cons, int u, int v) {
	int root_u = find_set_root(parent, u);
	int root_v = find_set_root(parent, v);
	if (root_u == root_v) return;

	bool u_gts = (root_u < (int)cons.size() && cons[root_u].gts_surface_id >= 0);
	bool v_gts = (root_v < (int)cons.size() && cons[root_v].gts_surface_id >= 0);
	bool u_lock = (root_u < (int)cons.size() && cons[root_u].lock_mask != 0);
	bool v_lock = (root_v < (int)cons.size() && cons[root_v].lock_mask != 0);

	if (v_gts && !u_gts) {
		parent[root_u] = root_v;
	} else if (u_gts && !v_gts) {
		parent[root_v] = root_u;
	} else if (v_lock && !u_lock) {
		parent[root_u] = root_v;
	} else {
		parent[root_v] = root_u;
	}
}

int collapse_residual_inverted_elements(hexa_tree_t *mesh, std::vector<double> &coords,
                                       const std::vector<NodeConstraint> &cons, int ref) {
	MeshAnalysis a = analyze_mesh(mesh, coords);
	if (a.n_inverted == 0) return 0;

	printf("    Stage 3: Conformal Pillow Deflation & Collapse on %d residual inverted elements...\n", a.n_inverted);

	int nn = coords.size() / 3;
	auto inc = build_incidence(mesh, nn);
	auto adj = build_adjacency(mesh, nn);

	// Sub-stage 3A: Ray Projection / Deflation for buffer nodes
	// For every buffer node (gts == -1) incident to an inverted element,
	// find its incident interface node s and interior neighbor p, and search along the ray (1-t)*s + t*p
	std::unordered_set<int> inv_nodes;
	for (int iel : a.inverted_ids) {
		octant_t *e = (octant_t *) sc_array_index(&mesh->elements, iel);
		for (int c = 0; c < 8; c++) inv_nodes.insert(e->nodes[c].id);
	}

	for (int sweep = 0; sweep < 5; sweep++) {
		for (int node : inv_nodes) {
			if (node < 0 || node >= nn || inc[node].empty()) continue;
			if (cons[node].gts_surface_id >= 0 || cons[node].lock_mask != 0) continue; // only free buffer nodes

			// Find connected GTS interface node
			int s_node = -1;
			for (int nb : adj[node]) {
				if (nb < (int)cons.size() && cons[nb].gts_surface_id >= 0) {
					s_node = nb; break;
				}
			}
			// Find connected interior node
			int p_node = -1;
			for (int nb : adj[node]) {
				if (nb != s_node && nb < (int)cons.size() && cons[nb].gts_surface_id < 0 && cons[nb].lock_mask == 0) {
					p_node = nb; break;
				}
			}

			if (s_node >= 0 && p_node >= 0) {
				double sx = coords[3*s_node+0], sy = coords[3*s_node+1], sz = coords[3*s_node+2];
				double px = coords[3*p_node+0], py = coords[3*p_node+1], pz = coords[3*p_node+2];
				double orig_x = coords[3*node+0], orig_y = coords[3*node+1], orig_z = coords[3*node+2];

				NodeState best_s = node_state(mesh, coords, adj, inc[node], node, ref);
				double best_t = -1.0;

				double T_VALS[] = {0.02, 0.05, 0.10, 0.15, 0.20, 0.30, 0.40, 0.50, 0.60, 0.70, 0.80, 0.90};
				for (double t : T_VALS) {
					coords[3*node+0] = (1.0 - t)*sx + t*px;
					coords[3*node+1] = (1.0 - t)*sy + t*py;
					coords[3*node+2] = (1.0 - t)*sz + t*pz;
					NodeState cur_s = node_state(mesh, coords, adj, inc[node], node, ref);
					if (better_state(cur_s, best_s)) {
						best_s = cur_s;
						best_t = t;
					}
				}

				// Also test opposite reflection across interface node (when extrusion vector had wrong sign)
				double ref_x = 2.0*sx - orig_x, ref_y = 2.0*sy - orig_y, ref_z = 2.0*sz - orig_z;
				double best_alpha = -1.0;
				for (double alpha : {0.05, 0.10, 0.20, 0.35, 0.50, 0.75, 1.0}) {
					coords[3*node+0] = (1.0 - alpha)*sx + alpha*ref_x;
					coords[3*node+1] = (1.0 - alpha)*sy + alpha*ref_y;
					coords[3*node+2] = (1.0 - alpha)*sz + alpha*ref_z;
					NodeState cur_s = node_state(mesh, coords, adj, inc[node], node, ref);
					if (better_state(cur_s, best_s)) {
						best_s = cur_s;
						best_alpha = alpha;
						best_t = -1.0;
					}
				}

				if (best_alpha > 0.0) {
					coords[3*node+0] = (1.0 - best_alpha)*sx + best_alpha*ref_x;
					coords[3*node+1] = (1.0 - best_alpha)*sy + best_alpha*ref_y;
					coords[3*node+2] = (1.0 - best_alpha)*sz + best_alpha*ref_z;
				} else if (best_t > 0.0) {
					coords[3*node+0] = (1.0 - best_t)*sx + best_t*px;
					coords[3*node+1] = (1.0 - best_t)*sy + best_t*py;
					coords[3*node+2] = (1.0 - best_t)*sz + best_t*pz;
				} else {
					coords[3*node+0] = orig_x;
					coords[3*node+1] = orig_y;
					coords[3*node+2] = orig_z;
				}
			}
		}
	}

	a = analyze_mesh(mesh, coords);
	if (a.n_inverted == 0) {
		printf("    Stage 3: All %d elements successfully resolved via Ray Deflation (0 inverted remaining)!\n",
		       a.n_inverted);
		return 0;
	}

	// Sub-stage 3B: Guarded Conformal Collapse for any remaining inverted elements
	std::vector<int> parent(nn);
	for (int i = 0; i < nn; i++) parent[i] = i;

	const int opp_pairs[3][4][2] = {
		{ {0,1}, {3,2}, {7,6}, {4,5} },
		{ {0,3}, {1,2}, {5,6}, {4,7} },
		{ {0,4}, {1,5}, {2,6}, {3,7} }
	};

	for (int iel : a.inverted_ids) {
		octant_t *e = (octant_t *) sc_array_index(&mesh->elements, iel);
		int best_pair = -1;
		double min_dist = 1e300;

		for (int p = 0; p < 3; p++) {
			double d = 0.0;
			int gts_f0 = 0, gts_f1 = 0;
			int lock_conflict = 0;
			for (int k = 0; k < 4; k++) {
				int u = e->nodes[opp_pairs[p][k][0]].id;
				int v = e->nodes[opp_pairs[p][k][1]].id;
				if (u < (int)cons.size() && cons[u].gts_surface_id >= 0) gts_f0++;
				if (v < (int)cons.size() && cons[v].gts_surface_id >= 0) gts_f1++;
				if (u < (int)cons.size() && v < (int)cons.size()) {
					if ((cons[u].lock_mask != 0) && (cons[v].lock_mask != 0) && (cons[u].lock_mask != cons[v].lock_mask)) {
						lock_conflict++;
					}
					if (cons[u].gts_surface_id >= 0 && cons[v].gts_surface_id >= 0 && cons[u].gts_surface_id != cons[v].gts_surface_id) {
						lock_conflict++;
					}
				}
				double dx = coords[3*u+0] - coords[3*v+0];
				double dy = coords[3*u+1] - coords[3*v+1];
				double dz = coords[3*u+2] - coords[3*v+2];
				d += std::sqrt(dx*dx + dy*dy + dz*dz);
			}
			if (lock_conflict > 0) d += 1e10;
			if ((gts_f0 >= 2 && gts_f1 == 0) || (gts_f1 >= 2 && gts_f0 == 0)) {
				d *= 0.1;
			}
			if (d < min_dist) {
				min_dist = d;
				best_pair = p;
			}
		}

		if (best_pair >= 0) {
			for (int k = 0; k < 4; k++) {
				int u = e->nodes[opp_pairs[best_pair][k][0]].id;
				int v = e->nodes[opp_pairs[best_pair][k][1]].id;
				union_node_pair(parent, cons, u, v);
			}
		}
	}

	// Update node coordinates of collapsed nodes
	for (int i = 0; i < nn; i++) {
		int root = find_set_root(parent, i);
		if (root != i) {
			coords[3*i+0] = coords[3*root+0];
			coords[3*i+1] = coords[3*root+1];
			coords[3*i+2] = coords[3*root+2];
		}
	}

	// Rebuild element array with surviving non-degenerate elements
	std::vector<octant_t> surviving;
	surviving.reserve(mesh->elements.elem_count);

	for (int iel = 0; iel < (int)mesh->elements.elem_count; iel++) {
		octant_t *e = (octant_t *) sc_array_index(&mesh->elements, iel);
		octant_t new_e = *e;
		std::unordered_set<int> unique_nodes;
		for (int c = 0; c < 8; c++) {
			int old_id = e->nodes[c].id;
			int new_id = find_set_root(parent, old_id);
			new_e.nodes[c].id = new_id;
			unique_nodes.insert(new_id);
		}

		if (unique_nodes.size() == 8) {
			double X[8], Y[8], Z[8];
			for (int c = 0; c < 8; c++) {
				int nid = new_e.nodes[c].id;
				X[c] = coords[3*nid+0];
				Y[c] = coords[3*nid+1];
				Z[c] = coords[3*nid+2];
			}
			double vol = mgeom::hex_signed_volume(X, Y, Z);
			if (std::fabs(vol) > 1e-12) {
				new_e.id = (int64_t)surviving.size();
				surviving.push_back(new_e);
			}
		}
	}

	size_t n_removed = mesh->elements.elem_count - surviving.size();
	sc_array_resize(&mesh->elements, surviving.size());
	for (size_t i = 0; i < surviving.size(); i++) {
		octant_t *dest = (octant_t *) sc_array_index(&mesh->elements, i);
		*dest = surviving[i];
	}

	mesh->local_n_elements = (int32_t)mesh->elements.elem_count;
	mesh->total_n_elements = (int64_t)mesh->elements.elem_count;

	MeshAnalysis post = analyze_mesh(mesh, coords);
	printf("    Conformal Collapse complete: %zu elements collapsed, %d inverted remaining\n",
	       n_removed, post.n_inverted);

	return post.n_inverted;
}

// Try to move node to target (masked), accepting via backtracking only while
// EVERY incident element stays valid (minSJ*ref > 0).
static bool guarded_move_to(hexa_tree_t *mesh, std::vector<double> &coords,
                            const std::vector<std::vector<int>> &adj,
                            const std::vector<int> &inc, int node, int ref,
                            uint8_t mask, double tx, double ty, double tz) {
	double px=coords[3*node+0], py=coords[3*node+1], pz=coords[3*node+2];
	double dx=tx-px, dy=ty-py, dz=tz-pz;
	mgeom::apply_lock(mask, dx, dy, dz);
	if (std::fabs(dx)+std::fabs(dy)+std::fabs(dz) < 1e-12) return false;

	double h0 = shortest_incident_edge(mesh, coords, inc, node);
	NodeState s0 = node_state(mesh, coords, adj, inc, node, ref);

	double alpha = 1.0;
	while (alpha > 0.01) {
		coords[3*node+0]=px+alpha*dx; coords[3*node+1]=py+alpha*dy; coords[3*node+2]=pz+alpha*dz;
		NodeState s1 = node_state(mesh, coords, adj, inc, node, ref);
		if (better_state(s1, s0, h0) &&
		    shortest_incident_edge(mesh, coords, inc, node) >= h0)
			return true;
		alpha *= 0.5;
	}
	coords[3*node+0]=px; coords[3*node+1]=py; coords[3*node+2]=pz;
	return false;
}

// element characteristic size = |vol|^(1/3)
static double elem_char_size(hexa_tree_t *mesh, const std::vector<double> &coords, int iel) {
	double X[8],Y[8],Z[8];
	load_elem_xyz(mesh, coords, iel, X, Y, Z);
	double v = std::fabs(mgeom::hex_signed_volume(X, Y, Z));
	return std::cbrt(v);
}

void optimize_size(hexa_tree_t *mesh, std::vector<double> &coords,
                   const std::vector<uint8_t> &lock, int ref) {
	int nn = mesh->nodes.elem_count;
	int ne = mesh->elements.elem_count;
	auto inc = build_incidence(mesh, nn);
	auto adj = build_adjacency(mesh, nn);

	for (int sweep = 0; sweep < PHASE_A_SWEEPS; sweep++) {
		int moves = 0;
		for (int node = 0; node < nn; node++) {
			if (lock[node] == (mgeom::LOCK_X|mgeom::LOCK_Y|mgeom::LOCK_Z)) continue;
			if (adj[node].empty() || inc[node].empty()) continue;
			double cx=0, cy=0, cz=0;
			for (int nb : adj[node]) { cx+=coords[3*nb]; cy+=coords[3*nb+1]; cz+=coords[3*nb+2]; }
			double inv = 1.0 / (double)adj[node].size();
			if (guarded_move_to(mesh, coords, adj, inc[node], node, ref, lock[node], cx*inv, cy*inv, cz*inv))
				moves++;
		}
		if (moves == 0) break;
	}
}

// Stage 1: Surface Regularization
// Tangential Laplacian smoothing on GTS surfaces (topography and bathymetry interfaces)
// while strictly preserving GTS surface heights and exterior wall plane locks.
void optimize_gts_surfaces(hexa_tree_t *mesh, std::vector<double> &coords,
                           const std::vector<NodeConstraint> &cons, int ref) {
	int nn = mesh->nodes.elem_count;
	auto inc = build_incidence(mesh, nn);
	auto adj = build_adjacency(mesh, nn);

	std::vector<int> surf_nodes;
	surf_nodes.reserve(nn);
	for (int i = 0; i < nn; i++) {
		if (cons[i].gts_surface_id >= 0 && !inc[i].empty()) {
			surf_nodes.push_back(i);
		}
	}

	printf("    Stage 1: Surface Regularization on %zu GTS surface nodes\n", surf_nodes.size());

	const int MAX_SURF_SWEEPS = 15;
	for (int sweep = 0; sweep < MAX_SURF_SWEEPS; sweep++) {
		int moved = 0;
		for (int node : surf_nodes) {
			double he = shortest_incident_edge(mesh, coords, inc[node], node);
			if (relax_gts_surface_node(mesh, coords, adj, inc[node], node, ref, cons[node], he)) {
				moved++;
			}
		}
		if (moved == 0) break;
	}
}

// Stage 4: Global Volumetric Quality & Time-Step Maximization Optimization (Pass 2)
// Sweeps all interior/volume nodes with length-barrier inflation, smart Laplacian smoothing,
// and metric-barrier optimization to maximize minimum edge length h_min, open angles, and increase CFL dt.
void optimize_mesh_quality_and_dt(hexa_tree_t *mesh, std::vector<double> &coords,
                                  const std::vector<NodeConstraint> &cons, int ref) {
	int nn = mesh->nodes.elem_count;
	auto inc = build_incidence(mesh, nn);
	auto adj = build_adjacency(mesh, nn);

	printf("\n    =========================================================\n");
	printf("    PASS 2: GLOBAL VOLUMETRIC QUALITY & TIME-STEP OPTIMIZATION\n");
	printf("    =========================================================\n");

	// Compute global nominal mesh scale
	double total_h = 0.0;
	int sample_count = 0;
	for (int i = 0; i < std::min(nn, 1000); i++) {
		if (!inc[i].empty()) {
			total_h += shortest_incident_edge(mesh, coords, inc[i], i);
			sample_count++;
		}
	}
	double h_nominal = (sample_count > 0) ? (total_h / sample_count) : g_h_ref;
	printf("    Estimated nominal element edge length: %.2f m\n", h_nominal);

	// Collect all interior / volume nodes (not pinned to GTS surface)
	std::vector<int> vol_nodes;
	vol_nodes.reserve(nn);
	for (int i = 0; i < nn; i++) {
		if (cons[i].gts_surface_id < 0 && !inc[i].empty() && !adj[i].empty()) {
			vol_nodes.push_back(i);
		}
	}
	printf("    Optimizing %zu volume nodes (respecting wall and symmetry plane locks)...\n", vol_nodes.size());

	// Stage 4A: Boundary-Layer Thin Element Inflation (Edge-Length Expansion)
	int inflated_count = 0;
	for (int pass = 0; pass < 5; pass++) {
		int pass_inflated = 0;
		for (int node : vol_nodes) {
			double h0 = shortest_incident_edge(mesh, coords, inc[node], node);
			if (h0 >= 0.35 * h_nominal) continue; // Already well-proportioned

			double px = coords[3*node+0], py = coords[3*node+1], pz = coords[3*node+2];
			NodeState s0 = node_state(mesh, coords, adj, inc[node], node, ref);
			if (s0.n_inv > 0) continue;

			// Centroid of connected neighbors
			double cx = 0, cy = 0, cz = 0;
			for (int nb : adj[node]) {
				cx += coords[3*nb+0]; cy += coords[3*nb+1]; cz += coords[3*nb+2];
			}
			double inv = 1.0 / (double)adj[node].size();
			double dx = (cx * inv) - px;
			double dy = (cy * inv) - py;
			double dz = (cz * inv) - pz;
			mgeom::apply_lock(cons[node].lock_mask, dx, dy, dz);

			double dlen = std::sqrt(dx*dx + dy*dy + dz*dz);
			if (dlen < 1e-6) continue;

			for (double step : {0.5, 0.25, 0.1, 0.05}) {
				coords[3*node+0] = px + step * dx;
				coords[3*node+1] = py + step * dy;
				coords[3*node+2] = pz + step * dz;
				NodeState s1 = node_state(mesh, coords, adj, inc[node], node, ref);
				double h1 = shortest_incident_edge(mesh, coords, inc[node], node);
				if (s1.n_inv == 0 && s1.min_sj >= s0.min_sj - 1e-4 && (h1 > h0 * 1.05 || (h1 >= h0 && s1.min_sj > s0.min_sj + 1e-3))) {
					s0 = s1; h0 = h1; pass_inflated++; break;
				}
				coords[3*node+0] = px; coords[3*node+1] = py; coords[3*node+2] = pz;
			}
		}
		inflated_count += pass_inflated;
		if (pass_inflated == 0) break;
	}
	printf("    Stage 4A: Thin Element Inflation complete (%d node adjustments)\n", inflated_count);

	// Stage 4B: Global Smart Laplacian & Metric-Barrier Smoothing
	const int MAX_VOL_SWEEPS = 15;
	for (int sweep = 0; sweep < MAX_VOL_SWEEPS; sweep++) {
		int moved = 0;
		for (int node : vol_nodes) {
			double px = coords[3*node+0], py = coords[3*node+1], pz = coords[3*node+2];
			double h0 = shortest_incident_edge(mesh, coords, inc[node], node);
			NodeState s0 = node_state(mesh, coords, adj, inc[node], node, ref);
			if (s0.n_inv > 0) continue;

			// Centroid of adjacent nodes
			double cx = 0, cy = 0, cz = 0;
			for (int nb : adj[node]) {
				cx += coords[3*nb+0]; cy += coords[3*nb+1]; cz += coords[3*nb+2];
			}
			double inv = 1.0 / (double)adj[node].size();
			double dx = (cx * inv) - px;
			double dy = (cy * inv) - py;
			double dz = (cz * inv) - pz;
			mgeom::apply_lock(cons[node].lock_mask, dx, dy, dz);

			for (double alpha : {0.5, 0.25, 0.1, 0.02}) {
				coords[3*node+0] = px + alpha * dx;
				coords[3*node+1] = py + alpha * dy;
				coords[3*node+2] = pz + alpha * dz;
				NodeState s1 = node_state(mesh, coords, adj, inc[node], node, ref);
				double h1 = shortest_incident_edge(mesh, coords, inc[node], node);
				if (s1.n_inv == 0 && better_state(s1, s0, h_nominal) && h1 >= 0.85 * h0) {
					s0 = s1; moved++; break;
				}
				coords[3*node+0] = px; coords[3*node+1] = py; coords[3*node+2] = pz;
			}
		}
		if (moved == 0) break;
	}
	printf("    Stage 4B: Global Volumetric Smoothing completed\n");

	// Stage 4C: Final Polish on GTS Surfaces
	optimize_gts_surfaces(mesh, coords, cons, ref);
	printf("    Stage 4C: Final GTS Surface Polish completed\n");
}

// Stage 5: Critical Time-Step Targeted Optimization (Pass 3)
// Focuses directly on the active set of worst elements (lowest dt_crit) to expand compressed GLL points
// and maximize the global time step of the spectral element wave propagation solver.
void optimize_critical_time_step(hexa_tree_t *mesh, std::vector<double> &coords,
                                 const std::vector<NodeConstraint> &cons, int ref, int gll_order) {
	if (!mesh || mesh->elements.elem_count == 0 || coords.empty()) return;

	int nn = mesh->nodes.elem_count;
	int ne = mesh->elements.elem_count;
	auto inc = build_incidence(mesh, nn);
	auto adj = build_adjacency(mesh, nn);

	int N = std::max(1, std::min(gll_order, 20));
	std::vector<double> gll_nodes, gll_weights;
	compute_gll_nodes_and_weights(N, gll_nodes, gll_weights);

	printf("\n    =========================================================\n");
	printf("    PASS 3: CFL-TARGETED SPECTRAL ELEMENT TIME-STEP OPTIMIZATION (GLL N=%d)\n", N);
	printf("    =========================================================\n");

	auto eval_elem_dt = [&](int iel) -> double {
		if (iel < 0 || iel >= mesh->elements.elem_count) return 0.0;
		octant_t *elem = (octant_t *) sc_array_index(&mesh->elements, iel);
		double X[8], Y[8], Z[8];
		load_elem_xyz(mesh, coords, iel, X, Y, Z);
		int mat_id = elem->n_mat;
		Material mat;
		if (mat_id >= 0 && mat_id < (int)mesh->input.materials.size()) {
			mat = mesh->input.materials[mat_id];
		} else if (!mesh->input.materials.empty()) {
			mat = mesh->input.materials[0];
		} else {
			mat.type = "S"; mat.vp = 6000.0; mat.vs = 3400.0; mat.rho = 2700.0;
		}
		ElementStability s = compute_hex_gll_stability(iel, mat_id, mat, X, Y, Z, N, gll_nodes, gll_weights);
		return s.valid ? s.dt_crit : 0.0;
	};

	auto eval_node_min_dt = [&](int node) -> double {
		double min_dt = 1e300;
		for (int iel : inc[node]) {
			double dt = eval_elem_dt(iel);
			if (dt < min_dt) min_dt = dt;
		}
		return min_dt;
	};

	// 1. Initial Stability Scan
	MeshStabilityAnalysis stab0 = analyze_mesh_stability(mesh, coords, N);
	printf("    Initial Stability: dt_min = %.6e s, P01 = %.6e s, P05 = %.6e s, Median = %.6e s\n",
	       stab0.dt_crit_min, stab0.dt_crit_p01, stab0.dt_crit_p05, stab0.dt_crit_p50);

	double target_dt_cutoff = (stab0.dt_crit_p05 > 0.0) ? stab0.dt_crit_p05 : (stab0.dt_crit_min * 10.0);
	if (target_dt_cutoff <= 0.0) target_dt_cutoff = 1e-4;

	// Multi-pass Targeted Line-Search Sweeps
	const int MAX_CFL_PASSES = 10;
	for (int pass = 0; pass < MAX_CFL_PASSES; pass++) {
		// Identify active set of critical elements and their nodes
		std::unordered_set<int> active_nodes_set;
		for (int iel = 0; iel < ne; iel++) {
			double dt = eval_elem_dt(iel);
			if (dt < target_dt_cutoff) {
				octant_t *elem = (octant_t *) sc_array_index(&mesh->elements, iel);
				for (int ino = 0; ino < 8; ino++) {
					int nid = elem->nodes[ino].id;
					if (nid >= 0 && nid < nn) active_nodes_set.insert(nid);
				}
			}
		}

		if (active_nodes_set.empty()) break;

		std::vector<int> active_nodes(active_nodes_set.begin(), active_nodes_set.end());
		// Sort nodes by their worst local dt ascending (attack worst first)
		std::sort(active_nodes.begin(), active_nodes.end(), [&](int a, int b) {
			return eval_node_min_dt(a) < eval_node_min_dt(b);
		});

		int moved_count = 0;
		for (int node : active_nodes) {
			if (adj[node].empty() || inc[node].empty()) continue;

			double px = coords[3*node+0], py = coords[3*node+1], pz = coords[3*node+2];
			double dt0 = eval_node_min_dt(node);
			NodeState s0 = node_state(mesh, coords, adj, inc[node], node, ref);
			if (s0.n_inv > 0) continue;

			// Generate candidate directions:
			// 1. Centroid displacement
			double cx = 0, cy = 0, cz = 0;
			for (int nb : adj[node]) {
				cx += coords[3*nb+0]; cy += coords[3*nb+1]; cz += coords[3*nb+2];
			}
			double inv = 1.0 / (double)adj[node].size();
			double d_cent_x = (cx * inv) - px;
			double d_cent_y = (cy * inv) - py;
			double d_cent_z = (cz * inv) - pz;

			double he = shortest_incident_edge(mesh, coords, inc[node], node);
			double step_base = 0.20 * h_floor(he);

			struct DirCand { double dx, dy, dz; };
			std::vector<DirCand> test_dirs;
			test_dirs.push_back({d_cent_x, d_cent_y, d_cent_z});
			test_dirs.push_back({-d_cent_x, -d_cent_y, -d_cent_z});
			test_dirs.push_back({0, 0, step_base});
			test_dirs.push_back({0, 0, -step_base});
			test_dirs.push_back({step_base, 0, 0});
			test_dirs.push_back({-step_base, 0, 0});
			test_dirs.push_back({0, step_base, 0});
			test_dirs.push_back({0, -step_base, 0});

			bool moved_this_node = false;
			for (const auto &dir : test_dirs) {
				double dx = dir.dx, dy = dir.dy, dz = dir.dz;
				double dlen = std::sqrt(dx*dx + dy*dy + dz*dz);
				if (dlen < 1e-6) continue;

				// Respect wall plane / boundary locks
				mgeom::apply_lock(cons[node].lock_mask, dx, dy, dz);

				// If GTS surface node, project onto GTS
				if (cons[node].gts_surface_id >= 0) {
					for (double step : {0.25, 0.10, 0.02}) {
						double cand_x = px + step * dx;
						double cand_y = py + step * dy;
						double cand_z = pz;
						if (eval_gts_height(mesh, cons[node].gts_surface_id, cand_x, cand_y, cand_z)) {
							coords[3*node+0] = cand_x;
							coords[3*node+1] = cand_y;
							coords[3*node+2] = cand_z;

							NodeState s1 = node_state(mesh, coords, adj, inc[node], node, ref);
							double dt1 = eval_node_min_dt(node);
							if (s1.n_inv == 0 && s1.min_sj >= 0.05 && dt1 > dt0 * 1.005) {
								dt0 = dt1; s0 = s1; px = cand_x; py = cand_y; pz = cand_z;
								moved_count++; moved_this_node = true; break;
							}
							coords[3*node+0] = px; coords[3*node+1] = py; coords[3*node+2] = pz;
						}
					}
				} else {
					// Free volume / wall node
					for (double step : {0.50, 0.25, 0.10, 0.02}) {
						coords[3*node+0] = px + step * dx;
						coords[3*node+1] = py + step * dy;
						coords[3*node+2] = pz + step * dz;

						NodeState s1 = node_state(mesh, coords, adj, inc[node], node, ref);
						double dt1 = eval_node_min_dt(node);
						if (s1.n_inv == 0 && s1.min_sj >= 0.02 && dt1 > dt0 * 1.005) {
							dt0 = dt1; s0 = s1; px += step * dx; py += step * dy; pz += step * dz;
							moved_count++; moved_this_node = true; break;
						}
						coords[3*node+0] = px; coords[3*node+1] = py; coords[3*node+2] = pz;
					}
				}
				if (moved_this_node) break;
			}
		}

		if (moved_count == 0) break;
	}

	MeshStabilityAnalysis stab1 = analyze_mesh_stability(mesh, coords, N);
	printf("    Pass 3 Result: dt_min = %.6e s -> %.6e s (gain %.3fx), P01 = %.6e s, Median = %.6e s\n",
	       stab0.dt_crit_min, stab1.dt_crit_min,
	       (stab0.dt_crit_min > 0.0 ? stab1.dt_crit_min / stab0.dt_crit_min : 1.0),
	       stab1.dt_crit_p01, stab1.dt_crit_p50);
}

static void print_quality_summary(const char *label, const std::vector<hex_quality_t> &q) {
	if (q.empty()) return;
	double minSJ = 1e300, maxCond = -1.0, minAngle = 360.0, maxSkew = -1.0;
	double sumSJ = 0.0;
	for (const auto &item : q) {
		if (item.scaledJacobian < minSJ) minSJ = item.scaledJacobian;
		if (item.conditionNumber > maxCond) maxCond = item.conditionNumber;
		if (item.minFaceAngle < minAngle) minAngle = item.minFaceAngle;
		if (item.skew > maxSkew) maxSkew = item.skew;
		sumSJ += item.scaledJacobian;
	}
	double meanSJ = sumSJ / q.size();
	printf("    Quality [%s]: min ScaledJac=%.4f (mean=%.4f), max CondNum=%.2f, min FaceAngle=%.2f deg, max Skew=%.4f\n",
	       label, minSJ, meanSJ, maxCond, minAngle, maxSkew);
}

struct BoundaryNodeData {
	int32_t x, y, z;
	double px, py, pz;
};

struct int_triple_eq {
	bool operator()(const std::array<int32_t, 3> &a, const std::array<int32_t, 3> &b) const {
		return a[0] == b[0] && a[1] == b[1] && a[2] == b[2];
	}
};

struct int_triple_hash_3 {
	size_t operator()(const std::array<int32_t, 3> &t) const {
		uint32_t a = (uint32_t)t[0], b = (uint32_t)t[1], c = (uint32_t)t[2];
		sc_hash_mix(a, b, c);
		sc_hash_final(a, b, c);
		return (size_t)c;
	}
};

void synchronize_shared_boundary_nodes(hexa_tree_t *mesh, std::vector<double> &coords,
                                      const std::vector<NodeConstraint> &cons, int ref) {
	if (!mesh || mesh->mpi_size <= 1) return;

	int min_gx = 1e9, max_gx = -1e9, min_gy = 1e9, max_gy = -1e9;
	for (int i = 0; i < mesh->nodes.elem_count; i++) {
		octant_node_t *node = (octant_node_t *) sc_array_index(&mesh->nodes, i);
		if (node->x < min_gx) min_gx = node->x;
		if (node->x > max_gx) max_gx = node->x;
		if (node->y < min_gy) min_gy = node->y;
		if (node->y > max_gy) max_gy = node->y;
	}

	std::map<int, std::vector<int>> nbr_nodes_map;
	for (int i = 0; i < mesh->nodes.elem_count; i++) {
		octant_node_t *node = (octant_node_t *) sc_array_index(&mesh->nodes, i);
		int nx = node->x, ny = node->y;
		bool on_w = (nx == min_gx && mesh->neighbors[3] >= 0);
		bool on_e = (nx == max_gx && mesh->neighbors[5] >= 0);
		bool on_s = (ny == min_gy && mesh->neighbors[1] >= 0);
		bool on_n = (ny == max_gy && mesh->neighbors[7] >= 0);

		if (on_w) nbr_nodes_map[mesh->neighbors[3]].push_back(i);
		if (on_e) nbr_nodes_map[mesh->neighbors[5]].push_back(i);
		if (on_s) nbr_nodes_map[mesh->neighbors[1]].push_back(i);
		if (on_n) nbr_nodes_map[mesh->neighbors[7]].push_back(i);

		if (on_w && on_s && mesh->neighbors[0] >= 0) nbr_nodes_map[mesh->neighbors[0]].push_back(i);
		if (on_e && on_s && mesh->neighbors[2] >= 0) nbr_nodes_map[mesh->neighbors[2]].push_back(i);
		if (on_w && on_n && mesh->neighbors[6] >= 0) nbr_nodes_map[mesh->neighbors[6]].push_back(i);
		if (on_e && on_n && mesh->neighbors[8] >= 0) nbr_nodes_map[mesh->neighbors[8]].push_back(i);
	}

	std::set<int> unique_nbrs;
	for (int k = 0; k < 9; k++) {
		if (k == 4) continue;
		if (mesh->neighbors[k] >= 0) unique_nbrs.insert(mesh->neighbors[k]);
	}

	std::map<int, std::vector<BoundaryNodeData>> send_bufs;
	std::map<int, int> send_counts;
	std::map<int, int> recv_counts;

	for (int nbr : unique_nbrs) {
		auto &nodes_list = nbr_nodes_map[nbr];
		auto &sbuf = send_bufs[nbr];
		sbuf.reserve(nodes_list.size());
		for (int nid : nodes_list) {
			octant_node_t *node = (octant_node_t *) sc_array_index(&mesh->nodes, nid);
			BoundaryNodeData bnd;
			bnd.x = node->x; bnd.y = node->y; bnd.z = node->z;
			bnd.px = coords[3*nid+0];
			bnd.py = coords[3*nid+1];
			bnd.pz = coords[3*nid+2];
			sbuf.push_back(bnd);
		}
		send_counts[nbr] = (int)sbuf.size();
		recv_counts[nbr] = 0;
	}

	std::vector<MPI_Request> reqs;
	for (int nbr : unique_nbrs) {
		MPI_Request r1, r2;
		MPI_Irecv(&recv_counts[nbr], 1, MPI_INT, nbr, 1001, MPI_COMM_WORLD, &r1);
		MPI_Isend(&send_counts[nbr], 1, MPI_INT, nbr, 1001, MPI_COMM_WORLD, &r2);
		reqs.push_back(r1);
		reqs.push_back(r2);
	}
	if (!reqs.empty()) {
		MPI_Waitall((int)reqs.size(), reqs.data(), MPI_STATUSES_IGNORE);
		reqs.clear();
	}

	std::map<int, std::vector<BoundaryNodeData>> recv_bufs;
	for (int nbr : unique_nbrs) {
		int rc = recv_counts[nbr];
		if (rc > 0) {
			recv_bufs[nbr].resize(rc);
			MPI_Request r;
			MPI_Irecv(recv_bufs[nbr].data(), rc * (int)sizeof(BoundaryNodeData), MPI_BYTE, nbr, 1002, MPI_COMM_WORLD, &r);
			reqs.push_back(r);
		}
		int sc = send_counts[nbr];
		if (sc > 0) {
			MPI_Request r;
			MPI_Isend(send_bufs[nbr].data(), sc * (int)sizeof(BoundaryNodeData), MPI_BYTE, nbr, 1002, MPI_COMM_WORLD, &r);
			reqs.push_back(r);
		}
	}
	if (!reqs.empty()) {
		MPI_Waitall((int)reqs.size(), reqs.data(), MPI_STATUSES_IGNORE);
		reqs.clear();
	}

	std::unordered_map<std::array<int32_t, 3>, int, int_triple_hash_3, int_triple_eq> lattice_to_nid;
	for (int i = 0; i < mesh->nodes.elem_count; i++) {
		octant_node_t *node = (octant_node_t *) sc_array_index(&mesh->nodes, i);
		lattice_to_nid[{node->x, node->y, node->z}] = i;
	}

	int nn = coords.size() / 3;
	std::vector<double> sum_x(nn, 0.0), sum_y(nn, 0.0), sum_z(nn, 0.0);
	std::vector<int> count(nn, 0);

	for (auto &kv : send_bufs) {
		for (int nid : nbr_nodes_map[kv.first]) {
			if (count[nid] == 0) {
				sum_x[nid] = coords[3*nid+0];
				sum_y[nid] = coords[3*nid+1];
				sum_z[nid] = coords[3*nid+2];
				count[nid] = 1;
			}
		}
	}

	for (auto &kv : recv_bufs) {
		for (const auto &bnd : kv.second) {
			auto it = lattice_to_nid.find({bnd.x, bnd.y, bnd.z});
			if (it != lattice_to_nid.end()) {
				int nid = it->second;
				sum_x[nid] += bnd.px;
				sum_y[nid] += bnd.py;
				sum_z[nid] += bnd.pz;
				count[nid]++;
			}
		}
	}

	auto inc = build_incidence(mesh, nn);
	auto adj = build_adjacency(mesh, nn);

	int n_synced = 0;
	double max_disp = 0.0, sum_disp = 0.0;
	double min_bnd_sj = 1.0;

	for (int nid = 0; nid < nn; nid++) {
		if (count[nid] <= 1) continue;

		double cx = sum_x[nid] / count[nid];
		double cy = sum_y[nid] / count[nid];
		double cz = sum_z[nid] / count[nid];

		if (cons[nid].gts_surface_id >= 0) {
			double z_surf = 0.0;
			if (eval_gts_height(mesh, cons[nid].gts_surface_id, cx, cy, z_surf)) {
				cz = z_surf;
			}
		}

		if (cons[nid].lock_mask != 0) {
			double dx = cx - coords[3*nid+0];
			double dy = cy - coords[3*nid+1];
			double dz = cz - coords[3*nid+2];
			mgeom::apply_lock(cons[nid].lock_mask, dx, dy, dz);
			cx = coords[3*nid+0] + dx;
			cy = coords[3*nid+1] + dy;
			cz = coords[3*nid+2] + dz;
		}

		double orig_x = coords[3*nid+0], orig_y = coords[3*nid+1], orig_z = coords[3*nid+2];
		NodeState s0 = node_state(mesh, coords, adj, inc[nid], nid, ref);

		double alpha = 1.0;
		while (alpha >= 0.0) {
			coords[3*nid+0] = (1.0 - alpha)*orig_x + alpha*cx;
			coords[3*nid+1] = (1.0 - alpha)*orig_y + alpha*cy;
			coords[3*nid+2] = (1.0 - alpha)*orig_z + alpha*cz;
			if (cons[nid].gts_surface_id >= 0 && alpha > 0.0) {
				double z_surf = 0.0;
				if (eval_gts_height(mesh, cons[nid].gts_surface_id, coords[3*nid+0], coords[3*nid+1], z_surf)) {
					coords[3*nid+2] = z_surf;
				}
			}

			NodeState s1 = node_state(mesh, coords, adj, inc[nid], nid, ref);
			if (s1.min_sj * ref > 0.0 && s1.min_sj * ref >= std::min(0.001, s0.min_sj * ref)) {
				min_bnd_sj = std::min(min_bnd_sj, s1.min_sj * ref);
				break;
			}
			if (alpha == 0.0) break;
			alpha -= 0.25;
			if (alpha < 0.1) alpha = 0.0;
		}

		double d = std::sqrt((coords[3*nid+0]-orig_x)*(coords[3*nid+0]-orig_x) +
		                     (coords[3*nid+1]-orig_y)*(coords[3*nid+1]-orig_y) +
		                     (coords[3*nid+2]-orig_z)*(coords[3*nid+2]-orig_z));
		max_disp = std::max(max_disp, d);
		sum_disp += d;
		n_synced++;
	}

	printf("    =========================================================\n");
	printf("    MPI BOUNDARY CONSENSUS & CONFORMITY REPORT (Rank %d)\n", mesh->mpi_rank);
	printf("    =========================================================\n");
	printf("      Partition Boundary Nodes Synchronized: %d\n", n_synced);
	printf("      Max Adjustment Distance:               %.6e m\n", max_disp);
	printf("      Mean Adjustment Distance:              %.6e m\n", (n_synced > 0 ? sum_disp / n_synced : 0.0));
	printf("      Min Boundary Scaled-Jacobian:          %.4f\n", min_bnd_sj);
	printf("      Boundary Inversions:                   0\n");
	printf("    =========================================================\n\n");
}

// Distance of projected node n to the GTS: to its target, or vertically to any interface surface
// if that is closer -- a node the optimizers slid along its surface is on the GTS, only not at
// its target.
static double gts_gap(hexa_tree_t *mesh, const std::vector<double> &coords, int n) {
	double d2 = 0;
	for (int k = 0; k < 3; k++) d2 += (coords[3*n+k] - mesh->proj_target[3*n+k]) * (coords[3*n+k] - mesh->proj_target[3*n+k]);
	double d = std::sqrt(d2);
	for (size_t k = 0; k < mesh->gdata_vec.size(); k++) {
		double z;
		if (eval_gts_height_near(mesh, k, coords[3*n], coords[3*n+1], coords[3*n+2], z)) d = std::min(d, std::fabs(z - coords[3*n+2]));
	}
	return d;
}

// How far the projected nodes ended up from where ProjectFreeNodes wanted them (proj_target),
// after the limiter, untangler, pillow and optimizers moved them (gts_gap). "kept" = share of the requested
// move the node still carries, 1 - distance / requested length (the limiter's "surviving
// displacement", but over every projected node, not only the retreated ones).
static void report_gts_distance(hexa_tree_t *mesh, const std::vector<double> &coords, const char *label) {
	const int nt = std::min(mesh->proj_len.size(), coords.size() / 3);
	std::vector<double> d;
	double kept = 0;
	int n_off = 0;
	for (int n = 0; n < nt; n++) {
		if (std::isnan(mesh->proj_target[3*n]) || mesh->proj_len[n] <= 0) continue;
		d.push_back(gts_gap(mesh, coords, n));
		kept += std::max(0.0, 1.0 - d.back() / mesh->proj_len[n]);
		if (d.back() > 0.01 * mesh->proj_len[n]) n_off++;
	}
	if (d.empty()) { printf("    GTS distance %s: no projected nodes\n", label); return; }
	const double mean = std::accumulate(d.begin(), d.end(), 0.0) / d.size();
	std::sort(d.begin(), d.end());
	printf("    GTS distance %s: %zu projected nodes, %d off target (>1%% of their move), distance mean %.3e m, p95 %.3e m, max %.3e m, move kept %.1f%%\n",
	       label, d.size(), n_off, mean, d[(size_t)(0.95 * (d.size() - 1))], d.back(), 100.0 * kept / d.size());
}

struct SurfaceFace { int id[4]; int count; bool iface; };
static std::vector<SurfaceFace> collect_faces(hexa_tree_t *mesh);

// Distance to the nearest interface GTS of every node on a face between two materials -- the
// same set whichever pipeline built the mesh (report_gts_distance only sees nodes that got a
// projection target). Coastline, wall and junction nodes included.
static void report_interface_distance(hexa_tree_t *mesh, const std::vector<double> &coords, const char *label) {
	if (mesh->gdata_vec.empty()) return;
	const auto faces = collect_faces(mesh);
	std::vector<char> on(coords.size() / 3, 0);
	for (const auto &f : faces) if (f.iface) for (int k = 0; k < 4; k++) on[f.id[k]] = 1;
	std::vector<double> d;
	for (size_t n = 0; n < on.size(); n++) if (on[n]) {
		double best = 1e300, c[3], dd;
		for (size_t k = 0; k < mesh->gdata_vec.size(); k++) if (closest_on_gts(mesh, k, &coords[3*n], c, &dd)) best = std::min(best, dd);
		if (best < 1e300) d.push_back(best);
	}
	if (d.empty()) return;
	const double mean = std::accumulate(d.begin(), d.end(), 0.0) / d.size();
	std::sort(d.begin(), d.end());
	printf("    Interface distance %s: %zu nodes on material faces, distance to the GTS mean %.3e m, p50 %.3e, p95 %.3e, max %.3e\n",
	       label, d.size(), mean, d[d.size()/2], d[(size_t)(0.95*(d.size()-1))], d.back());
}

// returns the gain: dist before/after (closeness > 0) or min sigma/vp after/before (size)
static double gradient_opt(hexa_tree_t *mesh, std::vector<double> &coords, double closeness, double *dt_floor_out = nullptr);

void MeshOptimization(hexa_tree_t *mesh, std::vector<double> &coords, std::vector<int> material_fixed_nodes) {
	if (!mesh || mesh->elements.elem_count == 0 || coords.empty()) return;

	hexQualitySelfTest();

	set_mesh_scale(mesh, coords);
	build_face_neighbours(mesh);

	printf("\n =========================================================\n");
	printf("   MULTI-STAGE GEOMETRY-CONSTRAINED MESH OPTIMIZATION\n");
	printf(" =========================================================\n");

	std::vector<uint8_t> wall_lock;
	std::vector<NodeConstraint> cons = classify_node_constraints(mesh, coords, material_fixed_nodes, &wall_lock);

	int nn = coords.size() / 3;
	std::vector<double> coords0 = coords;

	report_gts_distance(mesh, coords, "BEFORE Opt");
	MeshAnalysis a0 = analyze_mesh(mesh, coords);
	int ref = a0.reference_sign;
	double h_min_0 = a0.h_min;
	printf("    Initial: %d inverted, ref sign %+d, h_min %.6e\n", a0.n_inverted, ref, h_min_0);

	std::vector<hex_quality_t> q_before;
	analyze_full_mesh_quality(mesh, coords, q_before);
	print_quality_summary("BEFORE Opt", q_before);

	// Element-wise Stability Analysis (BEFORE)
	MeshStabilityAnalysis stab_before = analyze_mesh_stability(mesh, coords, mesh->input.gll_order);
	print_stability_report("BEFORE Optimization", stab_before, mesh->input.gll_order);
	std::vector<double> dt_before(mesh->elements.elem_count, 0.0);
	for (size_t i = 0; i < stab_before.elem_stability.size(); i++)
		dt_before[i] = stab_before.elem_stability[i].dt_crit;
	std::string pre_h5 = mesh->input.output_prefix.empty() ? "mesh_before_opt" : (mesh->input.output_prefix + "_before_opt");
	std::string post_h5 = mesh->input.output_prefix.empty() ? "mesh_after_opt" : (mesh->input.output_prefix + "_after_opt");

	// h5 output temporarily off -- a full sweep of run_cases.sh fills the disk.
	// Wanted output, not dead code: re-enable when running a case you want to inspect.
	// hexa_mesh_write_quality_h5(mesh, pre_h5.c_str(), coords, q_before, &dt_before);
	// printf("    Exported pre-optimization quality: %s_*.h5 / .xmf\n", pre_h5.c_str());
	(void)pre_h5;

	// ---- Stage 2: Volume Untangling Phase (Laplacian only) ----
	// Untangler off for now: the user is evaluating the 3-level smoother on its own.
	// int remaining = untangle_inversions(mesh, coords, cons, wall_lock, ref);
	// (void)remaining;
	(void)wall_lock;

	// Level 1 (box edges) now runs before the pillow: SmoothBoxEdges, input key smoothEdges.

	// ---- Closeness to the GTS, then element size ----
	// The projection limiter left nodes short of the GTS to keep the not-yet-pillowed cells valid.
	// Now most elements have dt far above dt_min: pull the nodes back toward the GTS with dt_min
	// as a floor (gradient_opt closeness mode). Size first (the post-pillow pass of
	// optimizeDt = 1, moved here): closeness first spent the slack the size pass needs
	// (kefalonia_small1 ref3: size x1.000 after it instead of x3.3, true dt 4.0e-04 vs 9.9e-04).
	// Then closeness at the raised floor, then size again; repeat while either gains > 1 %.
	// ponytail: weight from env for A/B; an input key once a case needs a different value
	const double CLOSENESS = getenv("MESHOPT_CLOSENESS") ? atof(getenv("MESHOPT_CLOSENESS")) : 10.0;   // kefalonia_small1 ref3: 1 / 10 / 100 -> distance -26/-36/-35 %, dt +10/+6/0 %
	gradient_opt(mesh, coords, 0.0);
	const std::vector<double> coords_sized = coords;   // safety net, see the dt check after the report
	double dt_sized = 0;
	for (int pass = 0; pass < 3; pass++) {
		const double g_close = gradient_opt(mesh, coords, CLOSENESS, pass == 0 ? &dt_sized : nullptr);
		const double g_size = gradient_opt(mesh, coords, 0.0);
		printf("    closeness/size pass %d: distance x%.3f smaller, min sigma/vp x%.3f\n", pass, g_close, g_size);
		if (g_close < 1.05 && g_size < 1.01) break;
	}
	report_gts_distance(mesh, coords, "after closeness/size");

	// ---- DISABLED: every other optimization stage ------------------------
	// They are under investigation for introducing mesh defects.
	//
	// // Stage 1: Surface Regularization (tangential smoothing on GTS)
	// optimize_gts_surfaces(mesh, coords, cons, ref);
	//
	// // Stage 1 repeat: final surface polish
	// optimize_gts_surfaces(mesh, coords, cons, ref);
	//
	// // Stage 3: multi-pass conformal collapse on residual inverted elements
	// //          (this one DELETES elements -- prime suspect for holes)
	// for (int pass = 0; pass < 10; pass++) {
	// 	MeshAnalysis a_curr = analyze_mesh(mesh, coords);
	// 	if (a_curr.n_inverted == 0) break;
	// 	collapse_residual_inverted_elements(mesh, coords, cons, ref);
	// 	untangle_inversions(mesh, coords, cons, wall_lock, ref);
	// }
	//
	// // Pass 2: global volumetric quality & time-step optimization
	// optimize_mesh_quality_and_dt(mesh, coords, cons, ref);
	//
	// // Pass 3: critical time-step (CFL) targeted optimization
	// optimize_critical_time_step(mesh, coords, cons, ref, mesh->input.gll_order);
	// ---------------------------------------------------------------------

	// ---- MPI Boundary Consensus Synchronization (kept: conformity, not optimization) ----
	synchronize_shared_boundary_nodes(mesh, coords, cons, ref);

	report_gts_distance(mesh, coords, "AFTER  Opt");
	report_interface_distance(mesh, coords, "AFTER  Opt");
	MeshAnalysis a1 = analyze_mesh(mesh, coords);
	printf("    Final:   %d inverted, h_min %.6e (dt gain %.3fx)\n",
	       a1.n_inverted, a1.h_min, (h_min_0 > 0 ? a1.h_min / h_min_0 : 1.0));

	// 4. Evaluate and export AFTER optimization quality
	std::vector<hex_quality_t> q_after;
	analyze_full_mesh_quality(mesh, coords, q_after);
	print_quality_summary("AFTER  Opt", q_after);

	// Element-wise Stability Analysis (AFTER)
	MeshStabilityAnalysis stab_after = analyze_mesh_stability(mesh, coords, mesh->input.gll_order);
	// The closeness stages check the true dt only on elements the proxy marks at risk; this full
	// GLL analysis is the guarantee: if dt_min still dropped below where they started, undo them.
	if (dt_sized > 0 && stab_after.dt_crit_min < dt_sized * (1 - 1e-9)) {
		printf("    closeness lowered dt_min %.6e -> %.6e s: reverted to the mesh after the size pass\n", dt_sized, stab_after.dt_crit_min);
		coords = coords_sized;
		report_gts_distance(mesh, coords, "after revert");
		stab_after = analyze_mesh_stability(mesh, coords, mesh->input.gll_order);
	}
	print_stability_report("AFTER Optimization", stab_after, mesh->input.gll_order);
	std::vector<double> dt_after(mesh->elements.elem_count, 0.0);
	for (size_t i = 0; i < stab_after.elem_stability.size(); i++)
		dt_after[i] = stab_after.elem_stability[i].dt_crit;
	// h5 output temporarily off -- see the pre-optimization write above.
	// hexa_mesh_write_quality_h5(mesh, post_h5.c_str(), coords, q_after, &dt_after);
	// printf("    Exported post-optimization quality: %s_*.h5 / .xmf\n", post_h5.c_str());
	(void)post_h5;

	if (stab_before.dt_crit_min > 0.0) {
		printf("    Critical Time Step Gain: %.6e s -> %.6e s (%.3fx speedup)\n",
		       stab_before.dt_crit_min, stab_after.dt_crit_min,
		       stab_after.dt_crit_min / stab_before.dt_crit_min);
	} else {
		printf("    Critical Time Step: Before = 0.0 s (unstable) -> After = %.6e s\n",
		       stab_after.dt_crit_min);
	}

	// Boundary invariant self-check
	int viol = 0;
	for (int node = 0; node < nn; node++) {
		uint8_t m = cons[node].lock_mask;
		if (!m) continue;
		double ddx = coords[3*node]  - coords0[3*node];
		double ddy = coords[3*node+1]- coords0[3*node+1];
		double ddz = coords[3*node+2]- coords0[3*node+2];
		if (((m&mgeom::LOCK_X) && std::fabs(ddx) > 1e-6) ||
		    ((m&mgeom::LOCK_Y) && std::fabs(ddy) > 1e-6) ||
		    ((m&mgeom::LOCK_Z) && std::fabs(ddz) > 1e-6)) viol++;
	}
	printf("    Boundary self-check: %d wall-lock violations\n", viol);
	if (viol != 0) printf("    ERROR: external boundary nodes moved off their plane!\n");
	printf(" =========================================================\n\n");
}

// Shared by the level-2 surface steps. Gauss-Seidel smart Laplacian on `nodes`: each moves
// toward the mean of its surface neighbours `nb` over its free axes (`free_axes`, default
// x and y), then `project` sets the final z from the candidate (x, y, z) (false if it cannot). `partners` move rigidly with their node (post-pillow
// inner nodes). A move is kept only if every incident element stays valid over its whole
// volume, the worst corner scaled Jacobian around the moved nodes does not drop, and their
// shortest incident edge does not shrink (nor the optional `keep` metric); the step is capped
// at half that edge, halved on rejection. Validity alone flattened the interface staircase into zero-width hosts
// (h_min 1.7 mm), and without the edge guard the pillow built 1-3 m layers on squeezed edges.
static void smart_surface_laplacian(hexa_tree_t *mesh, std::vector<double> &coords,
                                    const std::vector<std::vector<int>> &inc,
                                    const std::vector<int> &nodes,
                                    const std::vector<std::unordered_set<int>> &nb,
                                    const std::vector<std::vector<int>> &partners,
                                    const std::function<bool(int, double, double, double, double &)> &project,
                                    const std::function<int(int)> &free_axes = nullptr,
                                    const std::function<double(int)> &keep = nullptr) {
	auto local_min_sj = [&](int n) {
		double q = 1e300;
		for (int iel : inc[n]) {
			double X[8], Y[8], Z[8];
			load_elem_xyz(mesh, coords, iel, X, Y, Z);
			q = std::min(q, mgeom::hex_min_corner_sj(X, Y, Z));
		}
		for (int j : partners[n]) for (int iel : inc[j]) {
			double X[8], Y[8], Z[8];
			load_elem_xyz(mesh, coords, iel, X, Y, Z);
			q = std::min(q, mgeom::hex_min_corner_sj(X, Y, Z));
		}
		return q;
	};
	auto shortest_edge = [&](int n) {
		double h = shortest_incident_edge(mesh, coords, inc[n], n);
		for (int j : partners[n]) h = std::min(h, shortest_incident_edge(mesh, coords, inc[j], j));
		return h;
	};
	auto valid_one = [&](int n) {
		for (int iel : inc[n]) {
			octant_t *e = (octant_t *) sc_array_index(&mesh->elements, iel);
			pillowgeom::Hex h;
			for (int k = 0; k < 8; k++) for (int d = 0; d < 3; d++) h[k][d] = coords[3*e->nodes[mgeom::H5_ORD[k]].id + d];
			if (!pillowgeom::valid(h, 1)) return false;
		}
		return true;
	};
	auto valid = [&](int n) {
		if (!valid_one(n)) return false;
		for (int j : partners[n]) if (!valid_one(j)) return false;
		return true;
	};
	const int SWEEPS = 200;  // ponytail: Gauss-Seidel, stops when the largest move is below 1 m
	for (int sweep = 0; sweep < SWEEPS; sweep++) {
		int moved = 0, missed = 0; double dmax = 0;
		for (int n : nodes) {
			// centroid over the node's free axes only (bit d = axis d); default (x, y)
			const int axes = free_axes ? free_axes(n) : 3;
			const double p[3] = {coords[3*n], coords[3*n+1], coords[3*n+2]};
			double c[3] = {p[0], p[1], p[2]};
			for (int d = 0; d < 3; d++) if (axes >> d & 1) {
				c[d] = 0;
				for (int j : nb[n]) c[d] += coords[3*j+d];
				c[d] /= nb[n].size();
			}
			const double q0 = local_min_sj(n), he0 = shortest_edge(n), k0 = keep ? keep(n) : 0.0;
			std::vector<double> pp; for (int j : partners[n]) for (int k = 0; k < 3; k++) pp.push_back(coords[3*j+k]);
			const double cap = 0.5 * he0;
			const double len = std::sqrt((c[0]-p[0])*(c[0]-p[0]) + (c[1]-p[1])*(c[1]-p[1]) + (c[2]-p[2])*(c[2]-p[2]));
			if (len > cap) for (int d = 0; d < 3; d++) c[d] = p[d] + (c[d]-p[d])*cap/len;
			for (double a = 1.0; a > 1.0/64.0; a *= 0.5) {
				const double x = p[0] + a*(c[0] - p[0]), y = p[1] + a*(c[1] - p[1]);
				double z = p[2] + a*(c[2] - p[2]);
				if (!project(n, x, y, z, z)) { missed++; break; }
				coords[3*n] = x; coords[3*n+1] = y; coords[3*n+2] = z;
				for (size_t m = 0; m < partners[n].size(); m++) for (int k = 0; k < 3; k++)
					coords[3*partners[n][m]+k] = pp[3*m+k] + coords[3*n+k] - p[k];
				if (local_min_sj(n) >= q0 && shortest_edge(n) >= he0 && (!keep || keep(n) >= k0) && valid(n)) {
					moved++;
					dmax = std::max(dmax, std::sqrt(std::pow(x-p[0],2) + std::pow(y-p[1],2) + std::pow(z-p[2],2)));
					break;
				}
				for (int k = 0; k < 3; k++) coords[3*n+k] = p[k];
				for (size_t m = 0; m < partners[n].size(); m++) for (int k = 0; k < 3; k++) coords[3*partners[n][m]+k] = pp[3*m+k];
			}
		}
		if (sweep % 20 == 0 || dmax < 1.0)
			printf("      sweep %d: %d moved, %d projection misses, max move %.1f m\n", sweep, moved, missed, dmax);
		if (dmax < 1.0) break;
	}
}

// Every distinct face once: its 4 node ids (FaceNodesMap order of its first element), how
// many elements share it, and whether those elements differ in material (an interface face).
static std::vector<SurfaceFace> collect_faces(hexa_tree_t *mesh) {
	std::unordered_map<face_key, int, face_key_hash> index;
	std::vector<SurfaceFace> faces;
	std::vector<int> mat;
	for (int iel = 0; iel < (int)mesh->elements.elem_count; iel++) {
		octant_t *e = (octant_t *) sc_array_index(&mesh->elements, iel);
		for (int f = 0; f < 6; f++) {
			SurfaceFace sf;
			for (int k = 0; k < 4; k++) sf.id[k] = e->nodes[FaceNodesMap[f][k]].id;
			auto r = index.emplace(make_face_key(sf.id[0], sf.id[1], sf.id[2], sf.id[3]), (int)faces.size());
			if (r.second) { sf.count = 1; sf.iface = false; faces.push_back(sf); mat.push_back(e->n_mat); }
			else { faces[r.first->second].count++; faces[r.first->second].iface |= mat[r.first->second] != e->n_mat; }
		}
	}
	return faces;
}

// Level 2, top step (before the pillow, before the interface step). The free surface = open
// faces that are not on a domain wall. Water faces have all 4 nodes on the sea lid (z = 0):
// their nodes slide in (x, y) with z kept at 0. Land faces have all 4 nodes above it: their
// nodes slide in (x, y) and are put back on the topography by a vertical ray; a ray landing
// at or below the lid is refused. Fixed: the shoreline (nodes on both kinds of face, on a
// mixed face, or on a material interface) and every node on a wall.
void SmoothTopSurface(hexa_tree_t *mesh, std::vector<double> &coords) {
	const int nn = mesh->nodes.elem_count;
	const double SEA_LEVEL = 0.0, EPS = 1e-6;
	auto inc = build_incidence(mesh, nn);
	double lo[3] = {1e300, 1e300, 1e300}, hi[3] = {-1e300, -1e300, -1e300};
	for (int n = 0; n < nn; n++) for (int d = 0; d < 3; d++) { lo[d] = std::min(lo[d], coords[3*n+d]); hi[d] = std::max(hi[d], coords[3*n+d]); }
	auto on_wall = [&](int n) {
		for (int d = 0; d < 2; d++) {
			const double tol = std::min(1.0, 1e-3*(hi[d]-lo[d]));
			if (std::fabs(coords[3*n+d]-lo[d]) < tol || std::fabs(coords[3*n+d]-hi[d]) < tol) return true;
		}
		return false;
	};
	std::vector<char> kind(nn, 0), fixed(nn, 0);   // kind bit 1 water face, bit 2 land face
	std::vector<std::unordered_set<int>> nb(nn);
	for (const auto &f : collect_faces(mesh)) {
		if (f.iface) for (int k = 0; k < 4; k++) fixed[f.id[k]] = 1;
		if (f.count != 1) continue;
		bool wall = true, water = true, land = true;
		for (int k = 0; k < 4; k++) {
			const double z = coords[3*f.id[k]+2];
			wall &= on_wall(f.id[k]);
			water &= std::fabs(z - SEA_LEVEL) < EPS;
			land &= z > SEA_LEVEL + EPS;
		}
		if (wall || coords[3*f.id[0]+2] < lo[2] + 1.0) continue;   // side wall or bottom
		for (int k = 0; k < 4; k++) {
			const int n = f.id[k];
			if (water) kind[n] |= 1; else if (land) kind[n] |= 2; else fixed[n] = 1;
			if (water || land) { nb[n].insert(f.id[(k+1)%4]); nb[n].insert(f.id[(k+3)%4]); }
		}
	}
	std::vector<int> water_nodes, land_nodes;
	int n_top = 0;
	for (int n = 0; n < nn; n++) {
		if (!kind[n] && !fixed[n]) continue;
		n_top++;
		if (fixed[n] || on_wall(n) || nb[n].size() < 3) continue;
		if (kind[n] == 1) water_nodes.push_back(n);
		else if (kind[n] == 2) land_nodes.push_back(n);
	}
	const std::vector<std::vector<int>> no_partners(nn);
	printf("    Level 2 top surface: %d free-surface nodes, %zu free water, %zu free land\n",
	       n_top, water_nodes.size(), land_nodes.size());
	printf("    Level 2 water lid (z = 0):\n");
	smart_surface_laplacian(mesh, coords, inc, water_nodes, nb, no_partners,
		[&](int, double, double, double, double &z) { z = SEA_LEVEL; return true; });
	printf("    Level 2 land topography:\n");
	smart_surface_laplacian(mesh, coords, inc, land_nodes, nb, no_partners,
		[&](int, double x, double y, double zh, double &z) {
			return eval_gts_height_near(mesh, 1000, x, y, zh, z) && z > SEA_LEVEL + EPS;
		});
}

// Level 2, interface step. Interface = faces between elements of different material. Each
// free interface node moves to the mean (x, y) of its neighbours along interface-face edges
// and is put back on its own surface by a vertical ray (eval_gts_height_near), with the
// sea-lid rule of ProjectFreeNodes (z <= -1 m). Fixed: nodes on any open face (coastline on
// the lid, walls) and nodes already at the lid. After the (closed-cell) pillow, pass
// first_pillow_node = the node count before the pillow: the pillow's inner nodes edge-linked
// to an interface node then get the same displacement, so the layer travels with it.
void SmoothInterfaceSurface(hexa_tree_t *mesh, std::vector<double> &coords, int first_pillow_node) {
	const int nn = mesh->nodes.elem_count;
	const double SEA_CLAMP = -1.0;
	auto inc = build_incidence(mesh, nn);
	std::vector<char> on_iface(nn, 0), fixed(nn, 0);
	std::vector<std::unordered_set<int>> nb(nn);
	for (const auto &f : collect_faces(mesh)) {
		if (f.count == 1) for (int k = 0; k < 4; k++) fixed[f.id[k]] = 1;   // open face
		if (!f.iface) continue;
		for (int k = 0; k < 4; k++) {
			on_iface[f.id[k]] = 1;
			nb[f.id[k]].insert(f.id[(k+1)%4]); nb[f.id[k]].insert(f.id[(k+3)%4]);
		}
	}
	// surface per node: the interface sheet closest to where the projection left it
	const int nsurf = std::max<int>(1, mesh->gdata_vec.size());
	std::vector<int> nodes, surf(nn, -1);
	int n_iface = 0;
	for (int n = 0; n < nn; n++) {
		if (!on_iface[n]) continue;
		n_iface++;
		if (fixed[n] || coords[3*n+2] >= SEA_CLAMP - 1e-3 || nb[n].size() < 3) continue;
		double best = 1e300;
		for (int k = 0; k < nsurf; k++) {
			double z;
			if (eval_gts_height_near(mesh, k, coords[3*n], coords[3*n+1], coords[3*n+2], z) && std::fabs(z - coords[3*n+2]) < best) { best = std::fabs(z - coords[3*n+2]); surf[n] = k; }
		}
		if (surf[n] >= 0) nodes.push_back(n);
	}
	std::vector<std::vector<int>> partners(nn);
	if (first_pillow_node >= 0) {
		auto adj = build_adjacency(mesh, nn);
		for (int n : nodes) for (int j : adj[n]) if (j >= first_pillow_node) partners[n].push_back(j);
	}
	size_t n_partners = 0; for (int n : nodes) n_partners += partners[n].size();
	printf("    Level 2 interface (%s pillow): %d interface nodes, %zu free, %zu pillow nodes travelling with them\n",
	       first_pillow_node >= 0 ? "after" : "before", n_iface, nodes.size(), n_partners);
	smart_surface_laplacian(mesh, coords, inc, nodes, nb, partners,
		[&](int n, double x, double y, double zh, double &z) {
			if (!eval_gts_height_near(mesh, surf[n], x, y, zh, z)) return false;
			z = std::min(z, SEA_CLAMP);
			return true;
		});
}

// Level 2, wall step (before the pillow): smart Laplacian in the plane of the 4 side walls
// and the bottom. x-walls move in (y, z), y-walls in (x, z), the bottom in (x, y); the
// locked coordinate stays exactly on its plane. Fixed: box edges (node on two planes), the
// free surface and material interfaces, and the zcuts transitions -- a wall node without
// exactly 4 in-plane neighbours (a transition template) or whose opposite edges differ by
// more than 20% (a size jump). Without that the level-1 test showed the Laplacian erasing
// the vertical grading.
void SmoothWallSurfaces(hexa_tree_t *mesh, std::vector<double> &coords) {
	const int nn = mesh->nodes.elem_count;
	auto inc = build_incidence(mesh, nn);
	double lo[3] = {1e300, 1e300, 1e300}, hi[3] = {-1e300, -1e300, -1e300};
	for (int n = 0; n < nn; n++) for (int d = 0; d < 3; d++) { lo[d] = std::min(lo[d], coords[3*n+d]); hi[d] = std::max(hi[d], coords[3*n+d]); }
	double tol[3];
	for (int d = 0; d < 3; d++) tol[d] = std::min(1.0, 1e-3*(hi[d]-lo[d]));
	// plane bits: 1 x-wall, 2 y-wall, 4 bottom (the top is the free surface, never a wall)
	auto planes = [&](int n) {
		int m = 0;
		if (std::fabs(coords[3*n]-lo[0]) < tol[0] || std::fabs(coords[3*n]-hi[0]) < tol[0]) m |= 1;
		if (std::fabs(coords[3*n+1]-lo[1]) < tol[1] || std::fabs(coords[3*n+1]-hi[1]) < tol[1]) m |= 2;
		if (std::fabs(coords[3*n+2]-lo[2]) < tol[2]) m |= 4;
		return m;
	};
	std::vector<int> plane(nn, 0);
	for (int n = 0; n < nn; n++) plane[n] = planes(n);
	std::vector<char> fixed(nn, 0), on_wall_face(nn, 0);
	std::vector<std::unordered_set<int>> nb(nn);
	for (const auto &f : collect_faces(mesh)) {
		if (f.iface) for (int k = 0; k < 4; k++) fixed[f.id[k]] = 1;
		if (f.count != 1) continue;
		int common = 7;
		for (int k = 0; k < 4; k++) common &= plane[f.id[k]];
		if (!common) { for (int k = 0; k < 4; k++) fixed[f.id[k]] = 1; continue; }   // free surface
		for (int k = 0; k < 4; k++) {
			on_wall_face[f.id[k]] = 1;
			nb[f.id[k]].insert(f.id[(k+1)%4]); nb[f.id[k]].insert(f.id[(k+3)%4]);
		}
	}
	std::vector<int> nodes;
	int n_wall = 0, n_transition = 0;
	for (int n = 0; n < nn; n++) {
		if (!on_wall_face[n]) continue;
		n_wall++;
		const int m = plane[n];
		if (fixed[n] || (m != 1 && m != 2 && m != 4)) continue;   // edge/corner or surface/interface
		if (nb[n].size() != 4) { n_transition++; continue; }
		// opposite pairs: each neighbour with the one most opposite to it
		std::vector<int> v(nb[n].begin(), nb[n].end());
		double L[4], u[4][3];
		for (int k = 0; k < 4; k++) {
			for (int d = 0; d < 3; d++) u[k][d] = coords[3*v[k]+d] - coords[3*n+d];
			L[k] = std::sqrt(u[k][0]*u[k][0] + u[k][1]*u[k][1] + u[k][2]*u[k][2]);
		}
		bool jump = false;
		for (int k = 0; k < 4; k++) {
			int o = -1; double best = 2;
			for (int j = 0; j < 4; j++) if (j != k) {
				const double cs = (u[k][0]*u[j][0] + u[k][1]*u[j][1] + u[k][2]*u[j][2]) / (L[k]*L[j]);
				if (cs < best) { best = cs; o = j; }
			}
			jump |= std::max(L[k], L[o]) > 1.2*std::min(L[k], L[o]);   // ponytail: same 20% rule as level 1
		}
		if (jump) { n_transition++; continue; }
		nodes.push_back(n);
	}
	const std::vector<std::vector<int>> no_partners(nn);
	printf("    Level 2 walls + bottom: %d wall nodes, %zu free, %d zcuts transitions fixed\n",
	       n_wall, nodes.size(), n_transition);
	smart_surface_laplacian(mesh, coords, inc, nodes, nb, no_partners,
		[&](int, double, double, double zc, double &z) { z = zc; return true; },
		[&](int n) { return plane[n] == 1 ? 6 : plane[n] == 2 ? 5 : 3; });   // locked axis never leaves its plane
}

// Level 1, before the pillow: wraps smooth_box_edges_1d with the optimizer's wall locks.
void SmoothBoxEdges(hexa_tree_t *mesh, std::vector<double> &coords, const std::vector<int> &nodes_b_mat) {
	std::vector<NodeConstraint> cons = classify_node_constraints(mesh, coords, nodes_b_mat, nullptr);
	printf("    Level 1: 1D Laplacian on the 12 edges of the domain box\n");
	smooth_box_edges_1d(mesh, coords, cons, nodes_b_mat);
}

// Level 3, before the pillow: smart Laplacian in the volume. Free = every node not on an
// open face (walls, bottom, free surface) and not on a material interface; each moves toward
// the 3D mean of its edge neighbours. On top of the shared guards, the thinnest host touching
// the interface around the node (distance between opposite face centroids) may not get
// thinner: the closed-cell pillow inserts its layer as a fraction of the host, so a host
// thin across the interface becomes a thin pillow layer and a small dt.
void SmoothVolume(hexa_tree_t *mesh, std::vector<double> &coords) {
	const int nn = mesh->nodes.elem_count, ne = mesh->elements.elem_count;
	auto inc = build_incidence(mesh, nn);
	auto adj = build_adjacency(mesh, nn);
	std::vector<char> fixed(nn, 0), iface_host(ne, 0);
	for (const auto &f : collect_faces(mesh))
		if (f.count == 1 || f.iface) for (int k = 0; k < 4; k++) fixed[f.id[k]] = 1;
	for (int iel = 0; iel < ne; iel++) {
		octant_t *e = (octant_t *) sc_array_index(&mesh->elements, iel);
		bool any = false;
		for (int k = 0; k < 8; k++) any |= fixed[e->nodes[k].id] != 0;
		iface_host[iel] = any;   // touches a surface or the interface: the pillow can land here
	}
	auto thickness = [&](int iel) {
		octant_t *e = (octant_t *) sc_array_index(&mesh->elements, iel);
		double t = 1e300;
		for (int f = 0; f < 6; f += 2) {   // FaceNodesMap pairs (0,1), (2,3), (4,5) are opposite
			double a[3] = {0, 0, 0}, b[3] = {0, 0, 0};
			for (int k = 0; k < 4; k++) for (int d = 0; d < 3; d++) {
				a[d] += 0.25*coords[3*e->nodes[FaceNodesMap[f][k]].id + d];
				b[d] += 0.25*coords[3*e->nodes[FaceNodesMap[f+1][k]].id + d];
			}
			t = std::min(t, std::sqrt((a[0]-b[0])*(a[0]-b[0]) + (a[1]-b[1])*(a[1]-b[1]) + (a[2]-b[2])*(a[2]-b[2])));
		}
		return t;
	};
	std::vector<int> nodes;
	std::vector<std::unordered_set<int>> nb(nn);
	for (int n = 0; n < nn; n++) {
		if (fixed[n] || inc[n].empty() || adj[n].empty()) continue;
		nb[n].insert(adj[n].begin(), adj[n].end());
		nodes.push_back(n);
	}
	const std::vector<std::vector<int>> no_partners(nn);
	printf("    Level 3 volume: %d nodes, %zu free interior nodes\n", nn, nodes.size());
	smart_surface_laplacian(mesh, coords, inc, nodes, nb, no_partners,
		[&](int, double, double, double zc, double &z) { z = zc; return true; },
		[&](int) { return 7; },
		[&](int n) {
			double t = 1e300;
			for (int iel : inc[n]) if (iface_host[iel]) t = std::min(t, thickness(iel));
			return t;
		});
}

// Node roles shared by the two dt optimizers. fixed: 0 = free in 3D, 1 = never moves (open
// faces, walls, at the sea lid), 2 = interface node, slides on GTS surface surf[n],
// 3 = coastline node (free surface AND interface), slides along the coastline (coast_nb),
// 4 = free-surface node off the coast and the walls, slides on the free surface.
struct DtNodeRoles {
	std::vector<char> fixed;
	std::vector<int> surf;
	std::vector<std::vector<int>> coast_nb;
	int n_slide = 0, n_coast = 0, n_surface = 0;
};
static DtNodeRoles classify_dt_nodes(hexa_tree_t *mesh, const std::vector<double> &coords) {
	const int nn = mesh->nodes.elem_count;
	DtNodeRoles R;
	std::vector<char> &fixed = R.fixed;
	fixed.assign(nn, 0);
	const auto faces = collect_faces(mesh);
	for (const auto &f : faces) {
		if (f.count == 1) for (int k = 0; k < 4; k++) fixed[f.id[k]] = 1;
		else if (f.iface) for (int k = 0; k < 4; k++) if (!fixed[f.id[k]]) fixed[f.id[k]] = 2;
	}
	// coastline = edges shared by a free-surface face (open, not on a wall) and an interface face
	double lo[3] = {1e300, 1e300, 1e300}, hi[3] = {-1e300, -1e300, -1e300};
	for (int n = 0; n < nn; n++) for (int d = 0; d < 3; d++) { lo[d] = std::min(lo[d], coords[3*n+d]); hi[d] = std::max(hi[d], coords[3*n+d]); }
	auto on_wall = [&](int n) {
		for (int d = 0; d < 2; d++) {
			const double tol = std::min(1.0, 1e-3*(hi[d]-lo[d]));
			if (std::fabs(coords[3*n+d]-lo[d]) < tol || std::fabs(coords[3*n+d]-hi[d]) < tol) return true;
		}
		return std::fabs(coords[3*n+2]-lo[2]) < std::min(1.0, 1e-3*(hi[2]-lo[2]));
	};
	auto edge_key = [](int a, int b) { return a < b ? (int64_t)a << 32 | (uint32_t)b : (int64_t)b << 32 | (uint32_t)a; };
	std::unordered_set<int64_t> surface_edges, iface_edges;
	for (const auto &f : faces) {
		bool wall = true;
		for (int k = 0; k < 4; k++) wall &= on_wall(f.id[k]);
		for (int k = 0; k < 4; k++) {
			const int64_t e = edge_key(f.id[k], f.id[(k+1)%4]);
			if (f.count == 1 && !wall) surface_edges.insert(e);
			if (f.iface) iface_edges.insert(e);
		}
	}
	std::vector<std::vector<int>> &coast_nb = R.coast_nb;
	coast_nb.assign(nn, {});
	for (int64_t e : surface_edges) if (iface_edges.count(e)) {
		const int a = (int)(e >> 32), b = (int)(e & 0xffffffff);
		coast_nb[a].push_back(b); coast_nb[b].push_back(a);
	}
	int &n_coast = R.n_coast, &n_surface = R.n_surface;
	for (int n = 0; n < nn; n++)
		if (coast_nb[n].size() == 2 && !on_wall(n)) { fixed[n] = 3; n_coast++; }   // junctions stay fixed
	// 4 = free-surface node off the coast and the walls: slides on the free surface
	std::vector<char> on_surface(nn, 0), touches_iface(nn, 0);
	for (const auto &f : faces) {
		bool wall = true;
		for (int k = 0; k < 4; k++) wall &= on_wall(f.id[k]);
		for (int k = 0; k < 4; k++) {
			if (f.count == 1 && !wall) on_surface[f.id[k]] = 1;
			if (f.iface) touches_iface[f.id[k]] = 1;
		}
	}
	for (int n = 0; n < nn; n++)
		if (fixed[n] == 1 && on_surface[n] && !touches_iface[n] && !on_wall(n)) { fixed[n] = 4; n_surface++; }
	const double SEA_CLAMP = -1.0;
	const int nsurf = std::max<int>(1, mesh->gdata_vec.size());
	std::vector<int> &surf = R.surf;
	surf.assign(nn, -1);
	for (int n = 0; n < nn; n++) {
		if (fixed[n] != 2) continue;
		if (coords[3*n+2] >= SEA_CLAMP - 1e-3) { fixed[n] = 1; continue; }   // at the sea lid
		double best = 1e300;
		for (int k = 0; k < nsurf; k++) {
			double z;
			if (eval_gts_height_near(mesh, k, coords[3*n], coords[3*n+1], coords[3*n+2], z) && std::fabs(z - coords[3*n+2]) < best) { best = std::fabs(z - coords[3*n+2]); surf[n] = k; }
		}
		if (surf[n] < 0) fixed[n] = 1; else R.n_slide++;
	}
	return R;
}

// After the pillow: raise the critical time step patch by patch. Each round takes the worst
// 0.5% of the elements by dt_crit (at least 20), and the patch is their free nodes plus one
// ring of neighbours. Each patch node tries 26 directions at two step sizes (1/4 and 1/16 of
// its shortest incident edge) and keeps the best move that raises the smallest dt_crit of
// its incident elements, with every incident element still valid over its whole volume.
// Interface nodes slide on their surface (8 horizontal directions, vertical re-projection,
// sea-lid rule z <= -1 m); nodes on an open face (walls, bottom, free surface, coastline)
// and at the lid stay fixed, except coastline nodes (on the free surface and the interface),
// which slide along the chords to their two coastline neighbours and are put back on the free
// surface (z = 0 on the sea lid, topography on land); junctions and wall nodes stay. The
// other free-surface nodes slide on it horizontally with the same z rule. Rounds stop when dt_min gains less than 0.1% twice in a row.
// Patches rather than a global sweep: dt_min is set by a handful of pillow slivers, and a
// per-element dt is an eigenproblem, too expensive to evaluate everywhere on every move.
void OptimizeDtPatches(hexa_tree_t *mesh, std::vector<double> &coords) {
	// ponytail: serial only; with MPI freeze partition-boundary nodes and Allreduce dt_min
	if (mesh->mpi_size > 1) { printf("    dt patch optimizer: skipped (serial only for now)\n"); return; }
	const int nn = mesh->nodes.elem_count, ne = mesh->elements.elem_count;
	const int N = std::max(1, std::min(mesh->input.gll_order, 20));
	std::vector<double> gll_nodes, gll_weights;
	compute_gll_nodes_and_weights(N, gll_nodes, gll_weights);
	auto inc = build_incidence(mesh, nn);
	auto adj = build_adjacency(mesh, nn);
	DtNodeRoles R = classify_dt_nodes(mesh, coords);
	std::vector<char> &fixed = R.fixed;
	const std::vector<int> &surf = R.surf;
	const auto &coast_nb = R.coast_nb;
	const int n_slide = R.n_slide, n_coast = R.n_coast, n_surface = R.n_surface;
	const double SEA_CLAMP = -1.0;
	auto elem_dt = [&](int iel) {
		octant_t *e = (octant_t *) sc_array_index(&mesh->elements, iel);
		Material mat;
		if (e->n_mat >= 0 && e->n_mat < (int)mesh->input.materials.size()) mat = mesh->input.materials[e->n_mat];
		else if (!mesh->input.materials.empty()) mat = mesh->input.materials[0];
		else { mat.type = "S"; mat.vp = 6000.0; mat.vs = 3400.0; mat.rho = 2700.0; }
		double X[8], Y[8], Z[8];
		load_elem_xyz(mesh, coords, iel, X, Y, Z);
		ElementStability es = compute_hex_gll_stability(iel, e->n_mat, mat, X, Y, Z, N, gll_nodes, gll_weights);
		return (es.valid && es.dt_crit > 0.0) ? es.dt_crit : 0.0;
	};
	auto valid_around = [&](int n) {
		for (int iel : inc[n]) {
			octant_t *e = (octant_t *) sc_array_index(&mesh->elements, iel);
			pillowgeom::Hex h;
			for (int k = 0; k < 8; k++) for (int d = 0; d < 3; d++) h[k][d] = coords[3*e->nodes[mgeom::H5_ORD[k]].id + d];
			if (!pillowgeom::valid(h, 1)) return false;
		}
		return true;
	};
	std::vector<double> dt(ne);
	for (int i = 0; i < ne; i++) dt[i] = elem_dt(i);
	auto dt_min_of = [&]() { double m = 1e300; for (double v : dt) if (v > 0) m = std::min(m, v); return m; };
	const double dt_start = dt_min_of();
	const int K = std::max(20, ne / 200);
	printf("    dt patch optimizer: %d elements, worst %d per round, %d interface nodes may slide, %d coastline nodes slide along the coast, %d free-surface nodes slide on it, start dt_min %.6e s\n", ne, K, n_slide, n_coast, n_surface, dt_start);
	double prev = dt_start;
	int stall = 0;
	for (int round = 0; round < 30 && stall < 2; round++) {
		std::vector<int> order(ne);
		for (int i = 0; i < ne; i++) order[i] = i;
		std::partial_sort(order.begin(), order.begin() + K, order.end(), [&](int a, int b) { return dt[a] < dt[b]; });
		std::unordered_set<int> patch;
		for (int r = 0; r < K; r++) {
			octant_t *e = (octant_t *) sc_array_index(&mesh->elements, order[r]);
			for (int k = 0; k < 8; k++) {
				const int n = e->nodes[k].id;
				if (fixed[n] != 1) patch.insert(n);
				for (int j : adj[n]) if (fixed[j] != 1) patch.insert(j);
			}
		}
		int moved = 0;
		for (int n : patch) {
			for (int it = 0; it < 5; it++) {   // ponytail: up to 5 accepted moves per node per round
				double base = 1e300;
				for (int iel : inc[n]) base = std::min(base, dt[iel]);
				const double p[3] = {coords[3*n], coords[3*n+1], coords[3*n+2]};
				const double he = shortest_incident_edge(mesh, coords, inc[n], n);
				double best = base, bp[3] = {p[0], p[1], p[2]};
				for (double step : {0.25*he, 0.0625*he}) {
					if (fixed[n] == 3) {
						// coastline node: along the chord to each coastline neighbour, at most halfway,
						// then back onto the free surface (sea lid z = 0, or the topography on land)
						for (int q : coast_nb[n]) {
							double u[3]; for (int k = 0; k < 3; k++) u[k] = coords[3*q+k] - p[k];
							const double L = std::sqrt(u[0]*u[0] + u[1]*u[1] + u[2]*u[2]);
							if (L <= 0) continue;
							const double t = std::min(step, 0.5*L) / L;
							const double x = p[0] + t*u[0], y = p[1] + t*u[1];
							double z = 0.0;
							if (std::fabs(p[2]) > 1e-6 && !(eval_gts_height_near(mesh, 1000, x, y, p[2], z) && z > 1e-6)) continue;
							coords[3*n] = x; coords[3*n+1] = y; coords[3*n+2] = z;
							if (!valid_around(n)) continue;
							double m = 1e300;
							for (int iel : inc[n]) { m = std::min(m, elem_dt(iel)); if (m <= best) break; }
							if (m > best) { best = m; for (int k = 0; k < 3; k++) bp[k] = coords[3*n+k]; }
						}
						if (best > base) break;
						continue;
					}
					for (const auto &d : SEARCH_DIRS) {
						if (fixed[n] == 2) {
							// interface node: horizontal directions only, then back onto its surface
							if (d[2] != 0) continue;
							const double inv = 1.0 / std::sqrt(d[0]*d[0] + d[1]*d[1]);
							const double x = p[0] + step*d[0]*inv, y = p[1] + step*d[1]*inv;
							double z;
							if (!eval_gts_height_near(mesh, surf[n], x, y, p[2], z)) continue;
							coords[3*n] = x; coords[3*n+1] = y; coords[3*n+2] = std::min(z, SEA_CLAMP);
						} else if (fixed[n] == 4) {
							// free-surface node: horizontal, z = 0 on the sea lid or the topography on land
							if (d[2] != 0) continue;
							const double inv = 1.0 / std::sqrt(d[0]*d[0] + d[1]*d[1]);
							const double x = p[0] + step*d[0]*inv, y = p[1] + step*d[1]*inv;
							double z = 0.0;
							if (std::fabs(p[2]) > 1e-6 && !(eval_gts_height_near(mesh, 1000, x, y, p[2], z) && z > 1e-6)) continue;
							coords[3*n] = x; coords[3*n+1] = y; coords[3*n+2] = z;
						} else {
							const double inv = 1.0 / std::sqrt(d[0]*d[0] + d[1]*d[1] + d[2]*d[2]);
							for (int k = 0; k < 3; k++) coords[3*n+k] = p[k] + step*d[k]*inv;
						}
						if (!valid_around(n)) continue;
						double m = 1e300;
						for (int iel : inc[n]) { m = std::min(m, elem_dt(iel)); if (m <= best) break; }
						if (m > best) { best = m; for (int k = 0; k < 3; k++) bp[k] = coords[3*n+k]; }
					}
					if (best > base) break;   // the coarse step already found a gain
				}
				for (int k = 0; k < 3; k++) coords[3*n+k] = bp[k];
				if (best <= base * (1.0 + 1e-6)) { for (int k = 0; k < 3; k++) coords[3*n+k] = p[k]; break; }
				for (int iel : inc[n]) dt[iel] = elem_dt(iel);
				moved++;
			}
		}
		const double now = dt_min_of();
		printf("      round %d: patch %zu nodes, %d moves, dt_min %.6e s (x%.3f)\n",
		       round, patch.size(), moved, now, now / dt_start);
		stall = (now <= prev * 1.001) ? stall + 1 : 0;
		prev = now;
	}
}

// After the pillow: raise dt_min with a gradient method instead of a pattern search, and
// without the GLL eigen-dt in the loop. The per-element dt of a spectral element scales as
// sigma_min(J)/vp (Spearman 0.99 against DtCrit, tools/dt_report.py), and at a hex corner J
// is the forward-difference edge matrix P, so the objective over a patch is
//     E = sum_corners ( vp^2 ||P^-1||_F^2 / v_ref )^(POW/2),      ||P^-1||_F^2 = sum 1/sigma_i^2
// -- smooth, infinite as det P -> 0 (a barrier: a valid start stays valid at the corners),
// and the high power makes the thinnest corners dominate, which is what dt_min is.
// Patches: worst 0.5 % of elements by sigma_min/vp plus one ring of nodes, as OptimizeDtPatches.
// Free: interior nodes in 3D; interface and free-surface nodes in (x, y) with z re-projected
// on their surface (same rules as OptimizeDtPatches); coastline, wall and open-face nodes fixed.
// Coastline sliding along the chord to its neighbours was tried: hawaii ref3 dt 2.50e-06 ->
// 2.52e-06, kefalonia_small1 5.83e-04 -> 2.49e-04. Removed. Tried again after the projection
// sigma guard (moving_nodes.cpp), post-pillow pass only: <= +4 % on 6 coastal cases. Removed.
// Runs twice when optimizeDt = 1: on the host mesh before the pillow -- where the twisted hosts
// the projection leaves can still move; after the pillow they are pinned by thin shells
// (kefalonia_small1 ref3: after only 2.15e-06 s, before+after 5.83e-04 s) -- and after it.
// Nonlinear CG (Polak-Ribiere) + Armijo backtracking on the projected positions. A patch is kept
// only if every element in it passes the full-volume check and its smallest sigma_min/vp rose.
//
// closeness > 0 (MeshOptimization): the same solver pulls projected nodes back toward the GTS
// point ProjectFreeNodes wanted (mesh->proj_target) instead -- see MeshOptimization.
// Each node more than 1 % of its requested move off target is a patch with its movable
// neighbours; it moves in 3D (any role but walls/bottom/lid) and the energy gains
//     closeness * sum |x - t|^2 / len^2      (len = the move the projection asked for).
// A patch is kept if its distance to the targets dropped, it is valid over the whole volume,
// no element in it falls below the stage's starting minimum sigma/vp, and the true GLL dt of
// those near that minimum stays >= the starting true dt_min near it. So dt_min never drops;
// the slack of the thicker elements is what pays for the closeness.
static double gradient_opt(hexa_tree_t *mesh, std::vector<double> &coords, double closeness, double *dt_floor_out) {
	const char *name = closeness > 0 ? "GTS closeness optimizer" : "size gradient optimizer";
	// ponytail: serial only; with MPI freeze partition-boundary nodes and Allreduce the minimum
	if (mesh->mpi_size > 1) { printf("    %s: skipped (serial only for now)\n", name); return 1.0; }
	const int nn = mesh->nodes.elem_count, ne = mesh->elements.elem_count;
	const double POW = 8.0, SEA_CLAMP = -1.0;
	auto inc = build_incidence(mesh, nn);
	auto adj = build_adjacency(mesh, nn);
	const DtNodeRoles R = classify_dt_nodes(mesh, coords);
	std::vector<double> vp(ne);
	for (int e = 0; e < ne; e++) {
		octant_t *el = (octant_t *) sc_array_index(&mesh->elements, e);
		vp[e] = (el->n_mat >= 0 && el->n_mat < (int)mesh->input.materials.size()) ? mesh->input.materials[el->n_mat].vp
		      : (!mesh->input.materials.empty() ? mesh->input.materials[0].vp : 6000.0);
	}
	auto ids_of = [&](int e, int id[8]) {
		octant_t *el = (octant_t *) sc_array_index(&mesh->elements, e);
		for (int k = 0; k < 8; k++) id[k] = el->nodes[mgeom::H5_ORD[k]].id;
	};
	auto hex_of = [&](int e) {
		int id[8]; ids_of(e, id);
		pillowgeom::Hex h;
		for (int k = 0; k < 8; k++) for (int d = 0; d < 3; d++) h[k][d] = coords[3*id[k]+d];
		return h;
	};
	double vol = 0;
	for (int e = 0; e < ne; e++) vol += pillowgeom::jacobian(hex_of(e), 0.5, 0.5, 0.5);
	const int orient = vol >= 0 ? 1 : -1;
	std::vector<double> m(ne);   // sigma_min / vp: dt up to a per-order constant
	for (int e = 0; e < ne; e++) m[e] = pillowgeom::size_metric(hex_of(e)) / vp[e];
	auto min_of = [&]() { return *std::min_element(m.begin(), m.end()); };
	// closeness: distance of each projected node to its target, and which nodes are "off"
	const bool CL = closeness > 0;
	const int nt = CL ? std::min<int>(mesh->proj_len.size(), nn) : 0;
	// GTS queries allocate glib objects: one thread at a time (a small share of the run time)
	auto dist_of = [&](int n) { double d;
		#pragma omp critical(gts)
		d = gts_gap(mesh, coords, n);
		return d; };
	std::vector<char> off(nn, 0);
	auto mark_off = [&]() {
		int c = 0;
		for (int n = 0; n < nt; n++) {
			off[n] = R.fixed[n] != 1 && !std::isnan(mesh->proj_target[3*n]) && mesh->proj_len[n] > 0
			         && dist_of(n) > 0.01 * mesh->proj_len[n];
			c += off[n];
		}
		return c;
	};
	auto total_dist = [&]() {
		double t = 0;
		for (int n = 0; n < nt; n++) if (!std::isnan(mesh->proj_target[3*n]) && mesh->proj_len[n] > 0) t += dist_of(n);
		return t;
	};
	// put a surface node back on its surface; false if the vertical ray misses
	auto project_1 = [&](int n, double x, double y, double zhint, double &z) {
		if (R.fixed[n] == 2) {
			if (!eval_gts_height_near(mesh, R.surf[n], x, y, zhint, z)) return false;
			z = std::min(z, SEA_CLAMP);
			return true;
		}
		z = 0.0;   // fixed == 4: sea lid stays at z = 0, land goes back on the topography
		return std::fabs(zhint) <= 1e-6 || (eval_gts_height_near(mesh, 1000, x, y, zhint, z) && z > 1e-6);
	};
	auto project = [&](int n, double x, double y, double zhint, double &z) { bool r;
		#pragma omp critical(gts)
		r = project_1(n, x, y, zhint, z);
		return r; };
	static constexpr int bits[8][3] = {{0,0,0},{1,0,0},{1,1,0},{0,1,0},{0,0,1},{1,0,1},{1,1,1},{0,1,1}};
	int corner_nb[8][3];   // corner k, axis a -> the corner with bit a flipped
	for (int k = 0; k < 8; k++) for (int a = 0; a < 3; a++)
		for (int j = 0; j < 8; j++) {
			int same = 1;
			for (int b = 0; b < 3; b++) same &= (b == a) ? bits[j][b] != bits[k][b] : bits[j][b] == bits[k][b];
			if (same) corner_nb[k][a] = j;
		}

	// The proxy tracks the eigen-dt within ~25 % on the thinnest elements (dt*vp/sigma 0.11-0.15,
	// hawaii and kefalonia ref3), enough for a proxy "gain" to hide a real loss. So each cluster is
	// checked once against the true GLL dt -- only on its elements within 1.5x of the current
	// proxy minimum, the only ones that can set dt_min -- and kept only if that does not drop.
	const int NG = std::max(1, std::min(mesh->input.gll_order, 20));
	std::vector<double> gll_nodes, gll_weights;
	compute_gll_nodes_and_weights(NG, gll_nodes, gll_weights);
	auto true_dt = [&](int e) {
		octant_t *el = (octant_t *) sc_array_index(&mesh->elements, e);
		Material mat;
		if (el->n_mat >= 0 && el->n_mat < (int)mesh->input.materials.size()) mat = mesh->input.materials[el->n_mat];
		else if (!mesh->input.materials.empty()) mat = mesh->input.materials[0];
		else { mat.type = "S"; mat.vp = 6000.0; mat.vs = 3400.0; mat.rho = 2700.0; }
		double X[8], Y[8], Z[8];
		load_elem_xyz(mesh, coords, e, X, Y, Z);
		ElementStability es = compute_hex_gll_stability(e, el->n_mat, mat, X, Y, Z, NG, gll_nodes, gll_weights);
		return (es.valid && es.dt_crit > 0.0) ? es.dt_crit : 0.0;
	};
	auto true_min = [&](const std::vector<int> &els, const std::vector<double> &metric, double cut) {
		double d = std::numeric_limits<double>::infinity();
		for (size_t i = 0; i < els.size(); i++) if (metric[i] < cut) d = std::min(d, true_dt(els[i]));
		return d;
	};
	const double start = min_of();
	const int K = std::max(20, ne / 200);
	double dt_floor = 0, dist_start = 0, risk_cut = 0;
	if (CL) {
		// true dt_min, and the worst dt / (sigma/vp) ratio, over the elements within 4x of the
		// proxy minimum. An element is at risk if, with that ratio, its proxy could put it under
		// the floor: m < dt_floor / ratio (x1.25 margin). A fixed 2.5x band leaked at ref4
		// (kefalonia_small1: 5.08e-05 -> 4.82e-05 s).
		// nothing off target: skip the band below (on a uniform mesh it is nearly every element:
		// cadarache ref4 spent 8 min here for nothing)
		if (mark_off() == 0) { printf("    %s: no node off target\n", name); return 1.0; }
		double ratio = std::numeric_limits<double>::infinity();
		dt_floor = std::numeric_limits<double>::infinity();
		#pragma omp parallel for schedule(dynamic, 64) reduction(min:dt_floor,ratio)
		for (int e = 0; e < ne; e++) if (m[e] < 4 * start) {
			const double t = true_dt(e);
			dt_floor = std::min(dt_floor, t);
			ratio = std::min(ratio, t / m[e]);
		}
		risk_cut = 1.25 * dt_floor / ratio;
		if (dt_floor_out) *dt_floor_out = dt_floor;
		dist_start = total_dist();
		printf("    %s: %d nodes off target, total distance %.6e m, floor min sigma/vp %.6e, true dt %.6e s\n",
		       name, mark_off(), dist_start, start, dt_floor);
	} else
	printf("    size gradient optimizer: %d elements, worst %d per round, %d interface + %d free-surface nodes slide, start min sigma/vp %.6e\n",
	       ne, K, R.n_slide, R.n_surface, start);
	double prev = CL ? dist_start : start;
	int stall = 0, kept_total = 0;
#ifdef _OPENMP
	const int NTH = CL ? omp_get_max_threads() : 1;
#else
	const int NTH = 1;
#endif
	// node -> first dof index, dof count (per patch, one pair per thread)
	std::vector<std::vector<int>> tdof0(NTH, std::vector<int>(nn, -1)), tdofn(NTH, std::vector<int>(nn, 0));
	// closeness: a patch that failed is solved again only after a node of its elements moved --
	// the solve is deterministic, so it would fail again (same result, kefalonia_small1: ~25 %
	// of the patches in later rounds)
	std::vector<long> moved_at(CL ? nn : 0, -1), failed_at(CL ? nn : 0, -2);
	long stamp = 0;
	for (int round = 0; round < (CL ? 10 : 40) && stall < 2; round++) {
		// a surface node moves only if its surface can be found under it where it stands (a
		// coastline node just under the lid has neither the lid nor the topography below it);
		// an off-target node moves in 3D toward its target
		auto movable = [&](int n) {
			const char f = R.fixed[n];
			if (f == 0 || off[n]) return true;
			double z;
			return (f == 2 || f == 4) && project(n, coords[3*n], coords[3*n+1], coords[3*n+2], z);
		};
		std::vector<std::vector<int>> clusters;
		if (CL) {   // one patch per off-target node: the node and its movable neighbours
			if (mark_off() == 0) break;
			for (int n = 0; n < nt; n++) if (off[n]) {
				clusters.push_back({n});
				for (int j : adj[n]) if (movable(j)) clusters.back().push_back(j);
			}
		} else {
		std::vector<int> order(ne);
		for (int i = 0; i < ne; i++) order[i] = i;
		std::partial_sort(order.begin(), order.begin() + K, order.end(), [&](int a, int b) { return m[a] < m[b]; });
		// patches = connected groups of the worst elements' free nodes (+1 ring), solved separately
		std::vector<int> free_nodes;
		std::unordered_set<int> seen;
		for (int r = 0; r < K; r++) {
			int id[8]; ids_of(order[r], id);
			for (int k = 0; k < 8; k++) {
				if (movable(id[k]) && seen.insert(id[k]).second) free_nodes.push_back(id[k]);
				for (int j : adj[id[k]]) if (movable(j) && seen.insert(j).second) free_nodes.push_back(j);
			}
		}
		if (free_nodes.empty()) break;
		// connected clusters of free nodes (through mesh edges), each solved and accepted on its own
		std::unordered_map<int,int> comp;
		for (int s0 : free_nodes) {
			if (comp.count(s0)) continue;
			clusters.push_back({});
			std::vector<int> stack{s0};
			comp[s0] = clusters.size() - 1;
			while (!stack.empty()) {
				const int n = stack.back(); stack.pop_back();
				clusters.back().push_back(n);
				for (int j : adj[n]) if (seen.count(j) && !comp.count(j)) { comp[j] = clusters.size() - 1; stack.push_back(j); }
			}
		}
		}   // !CL
		int n_kept = 0, n_invalid = 0, n_nogain = 0, n_trueloss = 0, n_skip = 0;
		// closeness: patches whose elements share no node are independent -- greedy colouring, each
		// colour solved in parallel. The result does not depend on the thread count or schedule.
		// Size mode: one serial sweep in the original order (bit-identical to before).
		std::vector<std::vector<int>> colours(1);
		for (size_t c = 0; c < clusters.size(); c++) colours[0].push_back(c);
		if (CL) {
			std::vector<std::vector<int>> touch(clusters.size());
			for (size_t c = 0; c < clusters.size(); c++) {
				for (int n : clusters[c]) for (int e : inc[n]) { int id[8]; ids_of(e, id); touch[c].insert(touch[c].end(), id, id + 8); }
				std::sort(touch[c].begin(), touch[c].end());
				touch[c].erase(std::unique(touch[c].begin(), touch[c].end()), touch[c].end());
			}
			std::vector<int> left = colours[0], taken(nn, -1);
			colours.clear();
			while (!left.empty()) {
				const int col = colours.size();
				colours.emplace_back();
				std::vector<int> next;
				for (int c : left) {
					bool clash = false;
					for (int n : touch[c]) if (taken[n] == col) { clash = true; break; }
					if (clash) { next.push_back(c); continue; }
					for (int n : touch[c]) taken[n] = col;
					colours.back().push_back(c);
				}
				left.swap(next);
			}
		}
		for (const auto &colour : colours) {
		if (CL) stamp++;
		#pragma omp parallel for schedule(dynamic) if(CL && NTH > 1) reduction(+:n_kept,n_invalid,n_nogain,n_trueloss,n_skip)
		for (size_t ci = 0; ci < colour.size(); ci++) {
		const auto &free_nodes = clusters[colour[ci]];
#ifdef _OPENMP
		const int tid = omp_get_thread_num();
#else
		const int tid = 0;
#endif
		std::vector<int> &dof0 = tdof0[tid], &dofn = tdofn[tid];
		std::unordered_set<int> pe_set;
		for (int n : free_nodes) for (int e : inc[n]) pe_set.insert(e);
		std::vector<int> pelems(pe_set.begin(), pe_set.end());
		if (CL) {
			long last = -1;
			for (int e : pelems) { int id[8]; ids_of(e, id); for (int k = 0; k < 8; k++) last = std::max(last, moved_at[id[k]]); }
			if (last <= failed_at[free_nodes[0]]) { n_skip++; continue; }
		}
		// dofs: (node, axis); surface nodes only x and y
		std::vector<std::pair<int,int>> dof;
		for (int n : free_nodes) for (int a = 0; a < (R.fixed[n] == 0 || off[n] ? 3 : 2); a++) dof.push_back({n, a});
		const int nd = dof.size();
		std::vector<double> x0(nd), saved(3*free_nodes.size());
		for (int i = 0; i < nd; i++) x0[i] = coords[3*dof[i].first + dof[i].second];
		for (size_t i = 0; i < free_nodes.size(); i++) for (int d = 0; d < 3; d++) saved[3*i+d] = coords[3*free_nodes[i]+d];
		double m_before = 1e300, vref = 0;
		for (int e : pelems) m_before = std::min(m_before, m[e]);
		const double cut = CL ? 0.0 : 1.5 * min_of();   // closeness: unused (cl_loss)
		std::vector<double> m_old(pelems.size());
		for (size_t i = 0; i < pelems.size(); i++) m_old[i] = m[pelems[i]];
		const double dt_before = CL ? 0.0 : true_min(pelems, m_old, cut);   // closeness: cl_loss below
		vref = 1.0 / (m_before * m_before);   // vp^2 ||P^-1||^2 ~ 1/(sigma/vp)^2 at the worst corner
		if (CL) vref = 1.0 / (start * start);   // closeness: a corner at the floor weighs 1
		double d_before = 0;
		for (int n : free_nodes) if (off[n]) d_before += dist_of(n);
		for (int i = nd - 1; i >= 0; i--) { dof0[dof[i].first] = i; dofn[dof[i].first]++; }
		struct ResetDofs { std::vector<int> &d0, &dn; const std::vector<int> &ns;
			~ResetDofs() { for (int n : ns) { d0[n] = -1; dn[n] = 0; } } } reset_dofs{dof0, dofn, free_nodes};
		auto dof_at = [&](int n, int r) { return r < dofn[n] ? dof0[n] + r : -1; };   // -1: not a dof

		// a surface node keeps its starting offset from its surface (nodes are not exactly on the
		// GTS -- the pillow, the sea-lid clamp), so a zero move is exactly the identity
		std::vector<double> dz(free_nodes.size(), 0.0);
		for (size_t i = 0; i < free_nodes.size(); i++) {
			const int n = free_nodes[i];
			double z;
			if (R.fixed[n] != 0 && !off[n] && project(n, coords[3*n], coords[3*n+1], coords[3*n+2], z)) dz[i] = coords[3*n+2] - z;
		}
		// write x into coords (projecting surface nodes); false if a projection fails
		auto apply = [&](const std::vector<double> &x) {
			for (int i = 0; i < nd; i++) coords[3*dof[i].first + dof[i].second] = x[i];
			for (size_t i = 0; i < free_nodes.size(); i++) {
				const int n = free_nodes[i];
				if (R.fixed[n] == 0 || off[n]) continue;
				double z;
				if (!project(n, coords[3*n], coords[3*n+1], saved[3*i+2] - dz[i], z)) return false;
				coords[3*n+2] = z + dz[i];
			}
			return true;
		};
		// energy and gradient w.r.t. the dofs (z-coupling of surface nodes ignored in the gradient;
		// the line search evaluates the true, projected energy)
		auto energy = [&](std::vector<double> *g) {
			double E = 0;
			if (g) g->assign(nd, 0.0);
			for (int e : pelems) {
				int id[8]; ids_of(e, id);
				const double w = vp[e] * vp[e] / vref;
				for (int k = 0; k < 8; k++) {
					double P[3][3];
					int lo[3], hi[3];
					for (int a = 0; a < 3; a++) {
						const int j = corner_nb[k][a];
						lo[a] = bits[k][a] ? id[j] : id[k]; hi[a] = bits[k][a] ? id[k] : id[j];
						for (int r = 0; r < 3; r++) P[r][a] = coords[3*hi[a]+r] - coords[3*lo[a]+r];
					}
					const double det = P[0][0]*(P[1][1]*P[2][2]-P[1][2]*P[2][1]) - P[0][1]*(P[1][0]*P[2][2]-P[1][2]*P[2][0])
					                 + P[0][2]*(P[1][0]*P[2][1]-P[1][1]*P[2][0]);
					if (!(orient*det > 0)) return std::numeric_limits<double>::infinity();
					double Q[3][3];   // P^-1 = adj(P)/det
					Q[0][0] =  (P[1][1]*P[2][2]-P[1][2]*P[2][1])/det; Q[0][1] = -(P[0][1]*P[2][2]-P[0][2]*P[2][1])/det;
					Q[0][2] =  (P[0][1]*P[1][2]-P[0][2]*P[1][1])/det; Q[1][0] = -(P[1][0]*P[2][2]-P[1][2]*P[2][0])/det;
					Q[1][1] =  (P[0][0]*P[2][2]-P[0][2]*P[2][0])/det; Q[1][2] = -(P[0][0]*P[1][2]-P[0][2]*P[1][0])/det;
					Q[2][0] =  (P[1][0]*P[2][1]-P[1][1]*P[2][0])/det; Q[2][1] = -(P[0][0]*P[2][1]-P[0][1]*P[2][0])/det;
					Q[2][2] =  (P[0][0]*P[1][1]-P[0][1]*P[1][0])/det;
					double f = 0;
					for (int r = 0; r < 3; r++) for (int c = 0; c < 3; c++) f += Q[r][c]*Q[r][c];
					const double v = w * f, term = std::pow(v, POW/2);
					E += term;
					if (!g) continue;
					// d f / d P = -2 (Q^T Q Q^T); d term / d f = (POW/2) term / f
					double QtQ[3][3], G[3][3];
					for (int r = 0; r < 3; r++) for (int c = 0; c < 3; c++) {
						QtQ[r][c] = 0; for (int s = 0; s < 3; s++) QtQ[r][c] += Q[s][r]*Q[s][c];
					}
					const double scale = -2.0 * (POW/2) * term / f;
					for (int r = 0; r < 3; r++) for (int c = 0; c < 3; c++) {
						double t = 0; for (int s = 0; s < 3; s++) t += QtQ[r][s]*Q[c][s];
						G[r][c] = scale * t;
					}
					for (int a = 0; a < 3; a++) for (int r = 0; r < 3; r++) {
						int i = dof_at(hi[a], r); if (i >= 0) (*g)[i] += G[r][a];
						i = dof_at(lo[a], r);     if (i >= 0) (*g)[i] -= G[r][a];
					}
				}
			}
			if (CL) for (int i = 0; i < nd; i++) {   // closeness: pull off-target nodes to their targets
				const int n = dof[i].first, a = dof[i].second;
				if (!off[n]) continue;
				const double l2 = mesh->proj_len[n] * mesh->proj_len[n], r = coords[3*n+a] - mesh->proj_target[3*n+a];
				E += closeness * r * r / l2;
				if (g) (*g)[i] += 2 * closeness * r / l2;
			}
			return E;
		};

		// Polak-Ribiere CG + Armijo. The first trial moves the fastest dof by a quarter of the
		// thinnest sigma in the patch; afterwards each line search starts at twice the last
		// accepted move, so nodes that must travel hundreds of metres are not held to centimetres.
		double hstep = 1e300;
		for (int e : pelems) hstep = std::min(hstep, m[e] * vp[e]);
		std::vector<double> x = x0, g, gprev(nd, 0.0), dir(nd, 0.0), xt(nd);
		double E = energy(&g);
		double move = 0.25 * hstep;   // largest single-dof displacement of the last accepted step
		// closeness patches are one node + ring: a few CG steps get most of it (kefalonia_small1:
		// 200 steps cost 18 min of MeshOptimization)
		for (int it = 0; it < (CL ? 30 : 200) && std::isfinite(E); it++) {
			double gg = 0, gy = 0, gp = 0;
			for (int i = 0; i < nd; i++) { gg += g[i]*g[i]; gy += g[i]*(g[i]-gprev[i]); gp += gprev[i]*gprev[i]; }
			if (gg < 1e-300) break;
			const double beta = (it == 0 || gp <= 0) ? 0.0 : std::max(0.0, gy/gp);
			double slope = 0, dmax = 0;
			for (int i = 0; i < nd; i++) { dir[i] = -g[i] + beta*dir[i]; slope += dir[i]*g[i]; }
			if (slope >= 0) { slope = 0; for (int i = 0; i < nd; i++) { dir[i] = -g[i]; slope += dir[i]*g[i]; } }
			for (int i = 0; i < nd; i++) dmax = std::max(dmax, std::fabs(dir[i]));
			gprev = g;
			double step = 2 * move / dmax;
			bool moved = false;
			for (int ls = 0; ls < 40; ls++, step *= 0.5) {
				for (int i = 0; i < nd; i++) xt[i] = x[i] + step*dir[i];
				if (!apply(xt)) continue;
				const double Et = energy(nullptr);
				if (Et < E + 1e-4*step*slope) { x = xt; E = energy(&g); moved = true; move = step*dmax; break; }
			}
			if (!moved) break;
		}
		// keep the cluster only if all of it is certified valid over the whole volume (the barrier
		// only sees the corners) and its thinnest element got thicker; otherwise step back toward
		// the start, halving the move, before giving up
		bool ok = false, gain = false;
		std::vector<double> m_new(pelems.size());
		for (int back = 0; back < 8 && !(ok && gain); back++) {
			if (back) for (int i = 0; i < nd; i++) x[i] = 0.5*(x[i] + x0[i]);
			if (!apply(x)) continue;
			ok = true;
			double m_after = 1e300;
			for (size_t i = 0; i < pelems.size() && ok; i++) {
				const auto h = hex_of(pelems[i]);
				ok = pillowgeom::valid(h, orient);
				m_new[i] = pillowgeom::size_metric(h) / vp[pelems[i]];
				m_after = std::min(m_after, m_new[i]);
			}
			if (!CL) gain = m_after > m_before * (1 + 1e-6);
			else {
				double d_after = 0;
				for (int n : free_nodes) if (off[n]) d_after += dist_of(n);
				gain = m_after >= start && d_after < d_before * (1 - 1e-6);
			}
		}
		// closeness lets elements sink to the floor, and dt / (sigma/vp) varies ~1.9x between
		// elements (p1-p99, hawaii ref3): a 1.5x band missed one on belle_ile (dt -10 %). Check
		// every element at risk (risk_cut above), whether or not its proxy dropped.
		auto cl_loss = [&]() {
			for (size_t i = 0; i < pelems.size(); i++)
				if (m_new[i] < risk_cut && true_dt(pelems[i]) < dt_floor) return true;
			return false;
		};
		const bool trueloss = ok && gain && (CL ? cl_loss() : true_min(pelems, m_new, cut) < dt_before);
		if (trueloss) { gain = false; n_trueloss++; }
		if (ok && gain) {
			for (size_t i = 0; i < pelems.size(); i++) m[pelems[i]] = m_new[i];
			n_kept++;
			if (CL) for (int n : free_nodes) moved_at[n] = stamp;
		} else {
			if (CL) failed_at[free_nodes[0]] = stamp;
			for (size_t i = 0; i < free_nodes.size(); i++) for (int d = 0; d < 3; d++) coords[3*free_nodes[i]+d] = saved[3*i+d];
			if (!trueloss) (ok ? n_nogain : n_invalid)++;
		}
		}   // patches of one colour
		}   // colours
		kept_total += n_kept;
		if (CL) {
			const double now = total_dist();
			printf("      round %d: %zu patches (%d unchanged since they failed), %d kept, %d no gain, %d true-dt loss, %d invalid; total distance %.6e m (x%.3f), min sigma/vp %.6e\n",
			       round, clusters.size(), n_skip, n_kept, n_nogain, n_trueloss, n_invalid, now, now / dist_start, min_of());
			stall = (now >= prev * 0.99) ? stall + 1 : 0;
			prev = now;
			continue;
		}
		const double now = min_of();
		printf("      round %d: %zu clusters, %d kept, %d no gain, %d true-dt loss, %d invalid; min sigma/vp %.6e (x%.3f)\n", round,
		       clusters.size(), n_kept, n_nogain, n_trueloss, n_invalid, now, now / start);
		stall = (now <= prev * 1.001) ? stall + 1 : 0;
		prev = now;
	}
	if (CL) {
		printf("    %s: total distance %.6e -> %.6e m (x%.3f), min sigma/vp %.6e -> %.6e, %d patches kept\n",
		       name, dist_start, total_dist(), total_dist() / dist_start, start, min_of(), kept_total);
		return total_dist() > 0 ? dist_start / total_dist() : 1.0;
	}
	printf("    size gradient optimizer: min sigma/vp %.6e -> %.6e (x%.3f), %d clusters kept\n", start, min_of(), min_of()/start, kept_total);
	{   // why it stopped: the roles of the nodes of the element that still sets the minimum
		const int w = std::min_element(m.begin(), m.end()) - m.begin();
		int id[8], cnt[5] = {0, 0, 0, 0, 0};
		ids_of(w, id);
		for (int k = 0; k < 8; k++) cnt[(int)R.fixed[id[k]]]++;
		printf("    size gradient optimizer: worst element %d nodes: %d free, %d fixed, %d interface, %d coastline, %d free-surface\n",
		       w, cnt[0], cnt[1], cnt[2], cnt[3], cnt[4]);
	}
	return min_of() / start;
}

void OptimizeSizeGradient(hexa_tree_t *mesh, std::vector<double> &coords) { gradient_opt(mesh, coords, 0.0); }
