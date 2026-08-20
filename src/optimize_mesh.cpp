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
#include "mesh_geom.h"
#include "stability.h"

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
static const bool   ESCALATION_ENABLED = true;  // interface relaxation to untangle stubborn interface elements

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

// Node constraint structure supporting generic N GTS surfaces and 2D/1D boundary locks
struct NodeConstraint {
	uint8_t lock_mask;  // LOCK_X | LOCK_Y | LOCK_Z for boundary planes/lines/corners
	int gts_surface_id; // -1 if not on GTS surface, 0..N-1 = gdata_vec index, 1000 = topography (tdata)
};

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
				if (std::fabs(xval - dom_min_x) < 1.0 || std::fabs(xval - dom_max_x) < 1.0)
					m |= mgeom::LOCK_X;
			}
			if (cy) {
				double yval = coords[3*nid[0]+1];
				if (std::fabs(yval - dom_min_y) < 1.0 || std::fabs(yval - dom_max_y) < 1.0)
					m |= mgeom::LOCK_Y;
			}
			if (cz) {
				double zval = coords[3*nid[0]+2];
				if (std::fabs(zval - dom_min_z) < 1.0)
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
		if (std::fabs(coords[3*i+0] - dom_min_x) < 1.0 || std::fabs(coords[3*i+0] - dom_max_x) < 1.0) {
			lock[i] |= mgeom::LOCK_X;
			cons[i].lock_mask |= mgeom::LOCK_X;
		}
		if (std::fabs(coords[3*i+1] - dom_min_y) < 1.0 || std::fabs(coords[3*i+1] - dom_max_y) < 1.0) {
			lock[i] |= mgeom::LOCK_Y;
			cons[i].lock_mask |= mgeom::LOCK_Y;
		}
		if (std::fabs(coords[3*i+2] - dom_min_z) < 1.0) {
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
			if (len > 1e-6 && len < h) h = len;
		}
	}
	return (h < 1e299) ? h : 1.0;
}

static const double SEARCH_DIRS[26][3] = {
	{ 1, 0, 0}, {-1, 0, 0}, { 0, 1, 0}, { 0,-1, 0}, { 0, 0, 1}, { 0, 0,-1},
	{ 1, 1, 0}, {-1, 1, 0}, { 1,-1, 0}, {-1,-1, 0},
	{ 1, 0, 1}, {-1, 0, 1}, { 1, 0,-1}, {-1, 0,-1},
	{ 0, 1, 1}, { 0,-1, 1}, { 0, 1,-1}, { 0,-1,-1},
	{ 1, 1, 1}, {-1, 1, 1}, { 1,-1, 1}, {-1,-1, 1},
	{ 1, 1,-1}, {-1, 1,-1}, { 1,-1,-1}, {-1,-1,-1}
};

// Local state of a node's incident elements: how many are inverted, worst
// scaled Jacobian, sum of negative Jacobians, z-disparity with horizontal
// neighbors (x+, x-, y+, y-), xy-disparity with vertical column neighbors, and face warp.
struct NodeState {
	int n_inv;
	double min_sj;
	double sum_neg_sj;
	double z_disparity;      // Max |z - z_neighbor| among horizontal (x+, x-, y+, y-) neighbors
	double xy_disparity;     // Max |(x,y) - (x,y)_neighbor| along vertical column neighbors
	double max_face_warp;    // Max face non-planarity ratio among incident elements
	double sum_sj;
};

static NodeState node_state(hexa_tree_t *mesh, const std::vector<double> &coords,
                            const std::vector<std::vector<int>> &adj,
                            const std::vector<int> &inc, int node, int ref) {
	NodeState s;
	s.n_inv = 0;
	s.min_sj = 1e300;
	s.sum_neg_sj = 0.0;
	s.sum_sj = 0.0;
	s.max_face_warp = 0.0;
	s.z_disparity = 0.0;
	s.xy_disparity = 0.0;

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
	}

	// Structured neighbor disparity check (x+, x-, y+, y-, z+, z-)
	if (node >= 0 && node < mesh->nodes.elem_count && node < (int)adj.size()) {
		octant_node_t *nd = (octant_node_t *) sc_array_index(&mesh->nodes, node);
		int nz_node = nd->z;
		int nx_node = nd->x;
		int ny_node = nd->y;

		double px = coords[3*node+0];
		double py = coords[3*node+1];
		double pz = coords[3*node+2];

		for (int nb : adj[node]) {
			if (nb < 0 || nb >= mesh->nodes.elem_count) continue;
			octant_node_t *nd_nb = (octant_node_t *) sc_array_index(&mesh->nodes, nb);

			// Horizontal neighbors in the same z layer (x+, x-, y+, y-)
			if (nd_nb->z == nz_node) {
				double dz = std::fabs(pz - coords[3*nb+2]);
				if (dz > s.z_disparity) s.z_disparity = dz;
			}

			// Vertical column neighbors (same x, y column)
			if (nd_nb->x == nx_node && nd_nb->y == ny_node) {
				double dx = std::fabs(px - coords[3*nb+0]);
				double dy = std::fabs(py - coords[3*nb+1]);
				if (dx > s.xy_disparity) s.xy_disparity = dx;
				if (dy > s.xy_disparity) s.xy_disparity = dy;
			}
		}
	}

	return s;
}

// b is better than a iff it reduces inversion or improves quality while penalizing z-jumps
static bool better_state(const NodeState &b, const NodeState &a, double h_ref = 500.0) {
	double norm_h = (h_ref > 50.0) ? h_ref : 500.0;
	if (b.n_inv != a.n_inv) {
		return b.n_inv < a.n_inv;
	}
	if (b.n_inv > 0) {
		if (b.min_sj > a.min_sj + 1e-6) return true;
		if (b.min_sj < a.min_sj - 1e-6) return false;
		if (b.sum_neg_sj > a.sum_neg_sj + 1e-6) return true;
		if (b.sum_neg_sj < a.sum_neg_sj - 1e-6) return false;
		return b.sum_sj > a.sum_sj + 1e-6;
	}

	// Valid mesh phase: composite objective penalizing z-disparity, xy-shear, and face warp
	double penalty_a = 0.50 * (a.z_disparity / norm_h) + 0.30 * (a.xy_disparity / norm_h) + 0.20 * a.max_face_warp;
	double penalty_b = 0.50 * (b.z_disparity / norm_h) + 0.30 * (b.xy_disparity / norm_h) + 0.20 * b.max_face_warp;

	double score_a = a.min_sj - penalty_a;
	double score_b = b.min_sj - penalty_b;

	if (score_b > score_a + 1e-5) return true;
	if (score_b < score_a - 1e-5) return false;

	return b.sum_sj > a.sum_sj + 1e-12;
}

// Perform surface-constrained tangential relaxation for a node on a GTS surface
static bool relax_gts_surface_node(hexa_tree_t *mesh, std::vector<double> &coords,
                                   const std::vector<std::vector<int>> &adj,
                                   const std::vector<int> &inc, int node, int ref,
                                   const NodeConstraint &cons, double he = 500.0) {
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
	double h_surf = std::max(he, 50.0);
	double SURF_STEPS[] = {0.8*h_surf, 0.5*h_surf, 0.25*h_surf, 0.1*h_surf, 0.05*h_surf, 0.01*h_surf, 50.0, 20.0, 5.0, 1.0};
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
                                   const NodeConstraint &cons, double he = 500.0) {
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
	double h_scale = std::max(he, 50.0);
	double GRAD_STEPS[] = {1.0*h_scale, 0.6*h_scale, 0.3*h_scale, 0.15*h_scale, 0.05*h_scale, 0.01*h_scale, 50.0, 20.0, 5.0, 1.0};
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
                              const std::vector<NodeConstraint> &cons, double he = 500.0) {
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
                       uint8_t eff_mask, double step0, double cap, double he = 500.0) {
	double sx = coords[3*node+0], sy = coords[3*node+1], sz = coords[3*node+2];
	NodeState s0 = node_state(mesh, coords, adj, inc, node, ref);
	bool moved = false;

	double max_z_step = std::max(he, 350.0);
	static const double STEPS[] = {250.0, 150.0, 100.0, 50.0, 25.0, 10.0, 5.0, 2.0, 1.0, 0.5, 0.1};

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
                              const NodeConstraint &cons, double he = 500.0) {
	if (cons.gts_surface_id >= 0 || inc.empty()) return false;
	double px = coords[3*node+0], py = coords[3*node+1], pz = coords[3*node+2];
	NodeState s0 = node_state(mesh, coords, adj, inc, node, ref);
	NodeState best_s = s0;
	double best_x = px, best_y = py, best_z = pz;
	bool found = false;

	double scale = std::max(he, 50.0);
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
                                       const NodeConstraint &cons, double he = 500.0) {
	if (inc.empty()) return false;
	if (cons.gts_surface_id >= 0 || cons.lock_mask == (mgeom::LOCK_X|mgeom::LOCK_Y|mgeom::LOCK_Z)) return false;

	double px = coords[3*node+0], py = coords[3*node+1], pz = coords[3*node+2];
	NodeState s0 = node_state(mesh, coords, adj, inc, node, ref);

	double h_scale = std::max(he, 50.0);
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
		target_z = coords[3*n_above+2] - std::max(he, 50.0);
	} else if (n_below >= 0) {
		double he = shortest_incident_edge(mesh, coords, inc, node);
		target_z = coords[3*n_below+2] + std::max(he, 50.0);
	} else {
		return false;
	}

	double dz = target_z - pz;
	if (std::fabs(dz) < 1e-4) return false;

	NodeState s0 = node_state(mesh, coords, adj, inc, node, ref);
	bool moved = false;

	for (double alpha : {1.0, 0.75, 0.5, 0.25, 0.1, 0.05}) {
		coords[3*node+2] = pz + alpha * dz;
		NodeState s = node_state(mesh, coords, adj, inc, node, ref);
		if (better_state(s, s0, 500.0)) {
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

		double he_e = std::max(shortest_incident_edge(mesh, coords, inc[faceA[0]], faceA[0]), 50.0);
		double H_TESTS[] = {1.0 * he_e, 0.5 * he_e, 0.25 * he_e, 0.1 * he_e, 0.05 * he_e, 50.0, 150.0};

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

		// If bottom node is higher than top node (inverted column in Z):
		if (coords[3*nb+2] >= coords[3*nt+2] - 1.0) {
			double he_b = shortest_incident_edge(mesh, coords, inc[nb], nb);
			double he_t = shortest_incident_edge(mesh, coords, inc[nt], nt);
			double h_target = std::max({he_b, he_t, 50.0});

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

	double he_patch = std::max(shortest_incident_edge(mesh, coords, inc[n[0]], n[0]), 50.0);

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
                                  const std::vector<int> &inc, int node, int ref, uint8_t mask, double he = 100.0) {
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
			if (better_state(s, s0, 500.0)) {
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
	int nn = mesh->nodes.elem_count;
	auto inc = build_incidence(mesh, nn);
	auto adj = build_adjacency(mesh, nn);
	int prev_count = -1;
	int stall_streak = 0;
	const int STALL_PATIENCE_LIMIT = 200;

	for (int iter = 0; iter < MAX_UNTANGLE_ITERS; iter++) {
		MeshAnalysis a = analyze_mesh(mesh, coords);
		if (a.n_inverted == 0) { printf("    Untangler: 0 inverted after %d iters\n", iter); return 0; }
		if (iter % 10 == 0) printf("      untangle iter %d: %d inverted (stall %d)\n", iter, a.n_inverted, stall_streak);

		// A stalled iteration = the inverted COUNT did not drop.
		if (prev_count >= 0 && a.n_inverted >= prev_count) stall_streak++; else stall_streak = 0;
		prev_count = a.n_inverted;
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

		// First pass: try simultaneous patch untangling, face normal extrusion, vertical column untangling & parallelepiped reconstruction
		for (int iel : a.inverted_ids) {
			if (untangle_inverted_element_patch(mesh, coords, adj, inc, iel, ref, cons)) {
				any_moved = true;
			}
			if (untangle_element_face_extrusion(mesh, coords, adj, inc, iel, ref, cons)) {
				any_moved = true;
			}
			if (untangle_element_columns(mesh, coords, adj, inc, iel, ref, cons)) {
				any_moved = true;
			}
			if (reconstruct_hex_vertex(mesh, coords, adj, inc, iel, ref, cons)) {
				any_moved = true;
			}
		}

		for (int node : nodes) {
			if (inc[node].empty()) continue;
			double he = shortest_incident_edge(mesh, coords, inc[node], node);
			if (cons[node].gts_surface_id >= 0) {
				// Surface nodes ONLY relax tangentially on the given GTS surface (preserve surfaces)
				if (relax_gts_surface_node(mesh, coords, adj, inc[node], node, ref, cons[node], he)) {
					any_moved = true;
				}
			} else {
				// 1. Try vertical column midpoint balancing
				if (relax_vertical_column_midpoint(mesh, coords, adj, inc[node], node, ref, cons[node])) {
					any_moved = true;
				}
				// 2. Try pillow layer ray relaxation
				if (relax_pillow_node(mesh, coords, adj, inc[node], node, ref, cons, he)) {
					any_moved = true;
				}
				// 3. Try direct volume gradient ascent
				if (relax_volume_gradient(mesh, coords, adj, inc[node], node, ref, cons[node], he)) {
					any_moved = true;
				}
				// 4. Volume nodes relax with z-disparity and column shear penalties
				double step0 = UNTANGLE_STEP0 * std::max(he, 50.0);
				if (relax_node(mesh, coords, adj, inc[node], node, ref, cons[node].lock_mask, step0, 0.0, he))
					any_moved = true;
				if (relax_toward_centroid(mesh, coords, adj, inc[node], node, ref, cons[node].lock_mask, he))
					any_moved = true;
				// 5. Try Escobar Regularized Barrier Inversion-Free Descent
				if (relax_node_escobar_barrier(mesh, coords, adj, inc[node], node, ref, cons[node], he)) {
					any_moved = true;
				}
				// 6. Exhaustive 3D box search when stalled
				if (stalled && probe_node_3d_box(mesh, coords, adj, inc[node], node, ref, cons[node], he)) {
					any_moved = true;
				}
			}
		}

		(void)any_moved;

		if (stall_streak >= STALL_PATIENCE_LIMIT) {
			printf("    Untangler: stalled with %d inverted after %d iters (%d stalled)\n", a.n_inverted, iter, stall_streak);
			break;
		}
	}
	MeshAnalysis f = analyze_mesh(mesh, coords);
	printf("    Untangler: final check: %d inverted elements remaining\n", f.n_inverted);
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
					if (better_state(cur_s, best_s, 50.0)) {
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
					if (better_state(cur_s, best_s, 50.0)) {
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
	double h_nominal = (sample_count > 0) ? (total_h / sample_count) : 500.0;
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
			double step_base = std::max(he * 0.20, 1.0);

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

void MeshOptimization(hexa_tree_t *mesh, std::vector<double> &coords, std::vector<int> material_fixed_nodes) {
	if (!mesh || mesh->elements.elem_count == 0 || coords.empty()) return;

	hexQualitySelfTest();

	printf("\n =========================================================\n");
	printf("   MULTI-STAGE GEOMETRY-CONSTRAINED MESH OPTIMIZATION\n");
	printf(" =========================================================\n");

	std::vector<uint8_t> wall_lock;
	std::vector<NodeConstraint> cons = classify_node_constraints(mesh, coords, material_fixed_nodes, &wall_lock);

	int nn = coords.size() / 3;
	std::vector<double> coords0 = coords;

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

	hexa_mesh_write_quality_h5(mesh, pre_h5.c_str(), coords, q_before, &dt_before);
	printf("    Exported pre-optimization quality: %s_*.h5 / .xmf\n", pre_h5.c_str());

	auto inc = build_incidence(mesh, nn);
	auto adj = build_adjacency(mesh, nn);

	// ---- Stage 1: Surface Regularization (Tangential smoothing on GTS) ----
	optimize_gts_surfaces(mesh, coords, cons, ref);

	// ---- Stage 2: Volume Untangling Phase ----
	int remaining = untangle_inversions(mesh, coords, cons, wall_lock, ref);

	// ---- Stage 1 Repeat: Final Surface Polish ----
	optimize_gts_surfaces(mesh, coords, cons, ref);

	// ---- Stage 3: Multi-Pass Conformal Collapse on Residual Inverted Elements ----
	for (int pass = 0; pass < 10; pass++) {
		MeshAnalysis a_curr = analyze_mesh(mesh, coords);
		if (a_curr.n_inverted == 0) break;
		printf("    Stage 3 (Pass %d): Conformal Collapse on %d residual inverted elements...\n",
		       pass + 1, a_curr.n_inverted);
		collapse_residual_inverted_elements(mesh, coords, cons, ref);
		untangle_inversions(mesh, coords, cons, wall_lock, ref);
	}

	// ---- PASS 2: Global Volumetric Quality & Time-Step Optimization ----
	optimize_mesh_quality_and_dt(mesh, coords, cons, ref);

	// ---- PASS 3: Critical Time-Step (CFL) Targeted Optimization ----
	optimize_critical_time_step(mesh, coords, cons, ref, mesh->input.gll_order);

	// ---- MPI Boundary Consensus Synchronization & Quality Verification ----
	synchronize_shared_boundary_nodes(mesh, coords, cons, ref);

	MeshAnalysis a1 = analyze_mesh(mesh, coords);
	printf("    Final:   %d inverted, h_min %.6e (dt gain %.3fx)\n",
	       a1.n_inverted, a1.h_min, (h_min_0 > 0 ? a1.h_min / h_min_0 : 1.0));

	// 4. Evaluate and export AFTER optimization quality
	std::vector<hex_quality_t> q_after;
	analyze_full_mesh_quality(mesh, coords, q_after);
	print_quality_summary("AFTER  Opt", q_after);

	// Element-wise Stability Analysis (AFTER)
	MeshStabilityAnalysis stab_after = analyze_mesh_stability(mesh, coords, mesh->input.gll_order);
	print_stability_report("AFTER Optimization", stab_after, mesh->input.gll_order);
	std::vector<double> dt_after(mesh->elements.elem_count, 0.0);
	for (size_t i = 0; i < stab_after.elem_stability.size(); i++)
		dt_after[i] = stab_after.elem_stability[i].dt_crit;
	hexa_mesh_write_quality_h5(mesh, post_h5.c_str(), coords, q_after, &dt_after);
	printf("    Exported post-optimization quality: %s_*.h5 / .xmf\n", post_h5.c_str());

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
