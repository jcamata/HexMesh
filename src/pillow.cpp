// Previous implementation preserved verbatim for reference. The active implementation
// is in pillow_rebuild.cpp; it validates a complete insertion before changing the mesh.
#if 0
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <vector>
#include <iostream>
#include <unordered_map>
#include <array>
#include <unordered_set>
#include <algorithm>
#include <sc.h>
#include <sc_containers.h>
#include "hexa.h"
#include "mesh_geom.h"
#include "verify_mesh.h"
#include "pillow_untangle.h"

/*
 * Double-Layer Pillowing for Spectral Hexahedral Elements (Conforming Mesh for Continuous Galerkin).
 * 
 * 1. Factor N = 2 Grid Scaling:
 *    Doubling the integer grid coordinates (x,y,z) -> (2x,2y,2z) places original
 *    nodes at even integers (0, 2, 4...), allowing newly inserted pillow layer
 *    buffer nodes to be allocated at exact integer positions (1, 3, 5...)
 *    without any fractional coordinates.
 * 
 * 2. Axis-Aligned Pure Face Extrusion (Zero Cross-Axis Overlap):
 *    Pillowing buffer nodes are extruded strictly along the orthogonal face normal
 *    direction (x, y, or z) of each interface quad face.
 *    This prevents diagonal offsets at step corners, ensuring no element extrudes
 *    into neighboring elements along the wrong axis.
 */

struct int_triple {
	int x, y, z;
	bool operator==(const int_triple &o) const { return x == o.x && y == o.y && z == o.z; }
};

struct int_triple_hash {
	size_t operator()(const int_triple &t) const {
		uint32_t a = (uint32_t)t.x, b = (uint32_t)t.y, c = (uint32_t)t.z;
		sc_hash_mix(a, b, c);
		sc_hash_final(a, b, c);
		return (size_t)c;
	}
};

static const int FaceNormal[6][3] = {
	{-1,  0,  0}, // Face 0 (x-)
	{+1,  0,  0}, // Face 1 (x+)
	{ 0, -1,  0}, // Face 2 (y-)
	{ 0, +1,  0}, // Face 3 (y+)
	{ 0,  0, -1}, // Face 4 (un-flipped 2026-08-19: identical 836/836 inverted-element ids on the
	              // real bathymetry case with or without this flip -- proven irrelevant there;
	              // fixes the flat-interface case 836->0, so keeping the natural sign)
	{ 0,  0, +1}  // Face 5 (see Face 4 note)
};

struct face_key {
	int nodes[4];
	bool operator==(const face_key &o) const {
		return nodes[0] == o.nodes[0] && nodes[1] == o.nodes[1] &&
		       nodes[2] == o.nodes[2] && nodes[3] == o.nodes[3];
	}
};

struct face_key_hash {
	size_t operator()(const face_key &fk) const {
		uint32_t a = (uint32_t)fk.nodes[0], b = (uint32_t)fk.nodes[1], c = (uint32_t)fk.nodes[2];
		sc_hash_mix(a, b, c);
		c ^= (uint32_t)fk.nodes[3];
		sc_hash_final(a, b, c);
		return (size_t)c;
	}
};

static face_key make_face_key(const int n[4]) {
	face_key fk;
	fk.nodes[0] = n[0]; fk.nodes[1] = n[1]; fk.nodes[2] = n[2]; fk.nodes[3] = n[3];
	std::sort(fk.nodes, fk.nodes + 4);
	return fk;
}

static inline int64_t pillow_key(int node_id, int mat_id) {
	return ((int64_t)node_id << 32) | (uint32_t)mat_id;
}

// Structure to store interface quad face info. A and B are the two materials
// meeting at this quad, oriented deterministically by n_mat (A = lower value)
// so run-to-run output doesn't depend on element scan order.
typedef struct {
	int elemA_id, elemB_id; // the two elements on either side
	int matA, matB;         // their actual materials (matA < matB)
	int faceA, faceB;       // face index in elemA / elemB
	int quad_nodes[4];      // global node IDs of the quad
	int normal[3]; // face normal vector (nx, ny, nz) pointing from the A side to the B side
} interface_quad_t;

// Helper to evaluate minimum normalized corner cross product (sine of angle) across all 4 corners of a quad.
// Positive (>0) indicates convex; <=0 indicates reflex (>=180 deg) or self-intersecting bowtie.
static double quad_min_corner_s(const interface_quad_t &q, const std::vector<double> &coords) {
	double p[4][3];
	for (int k = 0; k < 4; k++) {
		int nid = q.quad_nodes[k];
		for (int d = 0; d < 3; d++) p[k][d] = coords[3*nid+d];
	}
	double d02[3] = { p[2][0]-p[0][0], p[2][1]-p[0][1], p[2][2]-p[0][2] };
	double d13[3] = { p[3][0]-p[1][0], p[3][1]-p[1][1], p[3][2]-p[1][2] };
	double N[3] = {
		d02[1]*d13[2] - d02[2]*d13[1],
		d02[2]*d13[0] - d02[0]*d13[2],
		d02[0]*d13[1] - d02[1]*d13[0]
	};
	double n_len = std::sqrt(N[0]*N[0] + N[1]*N[1] + N[2]*N[2]);
	if (n_len < 1e-12) return -1.0;

	double min_s = 1e300;
	for (int k = 0; k < 4; k++) {
		double v_prev[3] = { p[(k+3)%4][0]-p[k][0], p[(k+3)%4][1]-p[k][1], p[(k+3)%4][2]-p[k][2] };
		double v_next[3] = { p[(k+1)%4][0]-p[k][0], p[(k+1)%4][1]-p[k][1], p[(k+1)%4][2]-p[k][2] };
		double C[3] = {
			v_next[1]*v_prev[2] - v_next[2]*v_prev[1],
			v_next[2]*v_prev[0] - v_next[0]*v_prev[2],
			v_next[0]*v_prev[1] - v_next[1]*v_prev[0]
		};
		double l_prev = std::sqrt(v_prev[0]*v_prev[0] + v_prev[1]*v_prev[1] + v_prev[2]*v_prev[2]);
		double l_next = std::sqrt(v_next[0]*v_next[0] + v_next[1]*v_next[1] + v_next[2]*v_next[2]);
		double denom = n_len * l_prev * l_next;
		double s = (denom > 1e-12) ? (C[0]*N[0] + C[1]*N[1] + C[2]*N[2]) / denom : -1.0;
		if (s < min_s) min_s = s;
	}
	return min_s;
}

// The feasible cone of a shared pillow buffer node, in one place because two stages need the
// same answer: the extrusion field (which has to aim inside it) and the node-pinch despeckle
// (which has to detect when it is empty).
//
// At an interface corner the pillow hex's corner Jacobian is (e1 x e2) . d = cn . d, so a single
// buffer node serving every interface quad that meets this node on this material side needs one
// direction d with cn_q . d > 0 for all of them. The deepest such direction maximises
// min_q (cn_q . d) -- the Chebyshev centre of the spherical cone. Its optimum is supported by
// one, two or three of the constraints: a single normal, the bisector normalize(n_i + n_j) of a
// pair, or for a triple the direction equidistant from all three, which solves
// n_i.d = n_j.d = n_l.d, i.e. d is along (n_i - n_j) x (n_j - n_l). NOT the normalised sum --
// that is the average, and PILLOW_CONE_SELFCHECK caught the sphere beating it by 0.588.
//
// Returns the best achievable margin (<= 0 means the cone is EMPTY: no direction is valid, so at
// least one incident hex must invert whatever the mesher does) and, if asked, the direction.
static double cone_max_min(const std::vector<std::array<double,3>> &N,
                           std::array<double,3> *dir_out = NULL)
{
	const std::vector<std::array<double,3>> &U = N;
	if (U.empty()) return 1.0;                      // no constraint: anything is valid

	double bm = -1e300;
	std::array<double,3> best = {0.0, 0.0, 1.0};
	auto margin = [&U](const std::array<double,3> &d) {
		double m = 1e300;
		for (size_t i = 0; i < U.size(); i++) {
			double t = U[i][0]*d[0] + U[i][1]*d[1] + U[i][2]*d[2];
			if (t < m) m = t;
		}
		return m;
	};
	auto consider = [&](std::array<double,3> v) {
		double l = std::sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);
		if (l < 1e-12) return;
		for (int k = 0; k < 3; k++) v[k] /= l;
		double m = margin(v);
		if (m > bm) { bm = m; best = v; }
	};
	const size_t nk = U.size();
	std::array<double,3> tot = {0.0, 0.0, 0.0};
	for (size_t i = 0; i < nk; i++) {
		consider(U[i]);
		for (int k = 0; k < 3; k++) tot[k] += U[i][k];
		for (size_t j = i+1; j < nk; j++) {
			consider({U[i][0]+U[j][0], U[i][1]+U[j][1], U[i][2]+U[j][2]});
			for (size_t l = j+1; l < nk; l++) {
				double u[3], w[3], cr[3];
				for (int k = 0; k < 3; k++) { u[k] = U[i][k]-U[j][k]; w[k] = U[j][k]-U[l][k]; }
				cr[0] = u[1]*w[2]-u[2]*w[1];
				cr[1] = u[2]*w[0]-u[0]*w[2];
				cr[2] = u[0]*w[1]-u[1]*w[0];
				consider({ cr[0],  cr[1],  cr[2]});
				consider({-cr[0], -cr[1], -cr[2]});
			}
		}
	}
	consider(tot);
	if (dir_out) *dir_out = best;
	return bm;
}

// Which way is "into each material side" for this interface quad: the quad's own geometric normal
// (from the diagonals, the robust choice on a warped quad) plus the sign that turns it into
// material A and into material B, each decided by that host element's own interior node.
//
// Shared because two stages must agree exactly on it: the extrusion field, which aims the buffer
// node along it, and the cone-opening convexification, which moves interface nodes to make room
// for that aim. If they drifted, one would be opening a cone the other does not measure.
// Returns false on a degenerate quad, which contributes nothing to either.
static bool quad_side_signs(hexa_tree_t *mesh, const std::vector<double> &coords,
                            const interface_quad_t &q, double ng[3], int &sA, int &sB)
{
	const double *p[4];
	for (int n = 0; n < 4; n++) p[n] = &coords[3*q.quad_nodes[n]];
	double d02[3], d13[3];
	for (int k = 0; k < 3; k++) { d02[k] = p[2][k]-p[0][k]; d13[k] = p[3][k]-p[1][k]; }
	ng[0] = d02[1]*d13[2]-d02[2]*d13[1];
	ng[1] = d02[2]*d13[0]-d02[0]*d13[2];
	ng[2] = d02[0]*d13[1]-d02[1]*d13[0];
	double nl = std::sqrt(ng[0]*ng[0]+ng[1]*ng[1]+ng[2]*ng[2]);
	if (nl < 1e-12) return false;
	for (int k = 0; k < 3; k++) ng[k] /= nl;

	octant_t *eA = (octant_t *) sc_array_index(&mesh->elements, q.elemA_id);
	octant_t *eB = (octant_t *) sc_array_index(&mesh->elements, q.elemB_id);
	int iA = eA->nodes[FaceNodesMap_inv[q.faceA][0]].id;
	int iB = eB->nodes[FaceNodesMap_inv[q.faceB][0]].id;
	double dotA = 0.0, dotB = 0.0;
	for (int k = 0; k < 3; k++) {
		dotA += ng[k] * (coords[3*iA+k] - p[0][k]);
		dotB += ng[k] * (coords[3*iB+k] - p[0][k]);
	}
	sA = (dotA < 0.0 ? -1 : 1);
	sB = (dotB < 0.0 ? -1 : 1);
	return true;
}

// The corner normal of an interface quad at its corner n: e1 x e2 of the two in-plane edges,
// taken raw. On a reflex corner it genuinely opposes the quad's diagonal normal, and that is the
// information the cone needs, not something to normalise away. Multiply by sA / sB from
// quad_side_signs to get it oriented into a material side; the pillow hex's corner Jacobian there
// is then exactly cn . d.
static void quad_corner_normal(const std::vector<double> &coords,
                               const interface_quad_t &q, int n, double cn[3])
{
	int v = q.quad_nodes[n], a = q.quad_nodes[(n+1)%4], b = q.quad_nodes[(n+3)%4];
	double e1[3], e2[3];
	for (int k = 0; k < 3; k++) {
		e1[k] = coords[3*a+k] - coords[3*v+k];
		e2[k] = coords[3*b+k] - coords[3*v+k];
	}
	cn[0] = e1[1]*e2[2]-e1[2]*e2[1];
	cn[1] = e1[2]*e2[0]-e1[0]*e2[2];
	cn[2] = e1[0]*e2[1]-e1[1]*e2[0];
}

// Convexify interface quads prior to pillow buffer extrusion (Metodologia 1).
// Eliminates reflex corners (theta >= 180 deg) and bowtie crossings by moving interface nodes
// along the GTS surface, strictly ensuring:
// 1. Coupled star monotonicity: no neighboring interface quad in the node's star is degraded.
// 2. Zero host element inversions.
// 3. Exact GTS surface reprojection.
// 4. Zero hanging nodes (C0 conforming continuous nodal shift).
static int ConvexifyInterfaceQuads(hexa_tree_t *mesh, std::vector<double> &coords,
                                   const std::vector<interface_quad_t> &interface_quads,
                                   int ref)
{
	const size_t ne = mesh->elements.elem_count;
	std::unordered_map<int, std::vector<int>> n2e;
	for (size_t iel = 0; iel < ne; iel++) {
		octant_t *e = (octant_t *) sc_array_index(&mesh->elements, iel);
		for (int k = 0; k < 8; k++) {
			n2e[e->nodes[k].id].push_back((int)iel);
		}
	}

	// Bounding box for domain boundary walls
	double dom_min_x = 1e300, dom_max_x = -1e300;
	double dom_min_y = 1e300, dom_max_y = -1e300;
	for (size_t i = 0; i < coords.size()/3; i++) {
		double x = coords[3*i+0], y = coords[3*i+1];
		if (x < dom_min_x) dom_min_x = x;
		if (x > dom_max_x) dom_max_x = x;
		if (y < dom_min_y) dom_min_y = y;
		if (y > dom_max_y) dom_max_y = y;
	}
	const double wall_x = std::min(1.0, 1e-3 * (dom_max_x - dom_min_x));
	const double wall_y = std::min(1.0, 1e-3 * (dom_max_y - dom_min_y));

	// Precompute h_min, node-to-quad star mapping, and interface neighbor graph
	std::unordered_map<int, double> iface_hmin;
	std::unordered_map<int, std::vector<int>> n2q;
	std::unordered_map<int, std::vector<int>> n2nbr;
	for (size_t iq = 0; iq < interface_quads.size(); iq++) {
		const auto &q = interface_quads[iq];
		for (int n = 0; n < 4; n++) {
			int n1 = q.quad_nodes[n];
			int n2 = q.quad_nodes[(n+1)%4];
			n2q[n1].push_back((int)iq);
			n2nbr[n1].push_back(n2);
			n2nbr[n2].push_back(n1);

			double dx = coords[3*n1+0] - coords[3*n2+0];
			double dy = coords[3*n1+1] - coords[3*n2+1];
			double dz = coords[3*n1+2] - coords[3*n2+2];
			double elen = std::sqrt(dx*dx + dy*dy + dz*dz);
			auto it1 = iface_hmin.find(n1);
			if (it1 == iface_hmin.end() || elen < it1->second) iface_hmin[n1] = elen;
			auto it2 = iface_hmin.find(n2);
			if (it2 == iface_hmin.end() || elen < it2->second) iface_hmin[n2] = elen;
		}
	}

	// Unique neighbors in n2nbr
	for (auto &kv : n2nbr) {
		std::sort(kv.second.begin(), kv.second.end());
		kv.second.erase(std::unique(kv.second.begin(), kv.second.end()), kv.second.end());
	}

	// Associate interface nodes with closest GTS surface
	std::unordered_map<int, int> node_gts;
	for (auto &kv : iface_hmin) {
		int nid = kv.first;
		double px = coords[3*nid+0], py = coords[3*nid+1], pz = coords[3*nid+2];
		int best_g = -1;
		double best_dz = 1e300;
		for (size_t g = 0; g < mesh->gdata_vec.size(); g++) {
			double z_eval = pz;
			if (eval_gts_height_near(mesh, (int)g, px, py, pz, z_eval)) {
				double dz = std::fabs(pz - z_eval);
				if (dz < best_dz) { best_dz = dz; best_g = (int)g; }
			}
		}
		node_gts[nid] = (best_dz < 2.0) ? best_g : -1;
	}

	auto check_elem_inverted = [&](int iel) -> bool {
		double X[8], Y[8], Z[8];
		load_elem_xyz(mesh, coords, iel, X, Y, Z);
		return mgeom::is_inverted(mgeom::hex_signed_volume(X, Y, Z),
		                          mgeom::hex_min_corner_sj(X, Y, Z), ref);
	};

	// The feasible cone at a node, which is the quantity the extrusion field will later have to
	// aim inside: one shared buffer node per (node, material side) must satisfy cn . d > 0 for
	// every interface quad meeting it there. Section 14 measured the residual's cones as empty by
	// less than 0.05 -- within ~3 degrees of feasible -- so they are marginal geometry, and node
	// motion on the GTS surface is the lever that fits. Material reassignment was tried and made
	// things worse; the numbers are at the despeckle.
	auto eval_node_cone = [&](int nid) -> double {
		auto it = n2q.find(nid);
		if (it == n2q.end() || it->second.empty()) return 1.0;
		std::unordered_map<int, std::vector<std::array<double,3>>> by_mat;
		for (int iq : it->second) {
			const interface_quad_t &q = interface_quads[iq];
			double ng[3]; int sA, sB;
			if (!quad_side_signs(mesh, coords, q, ng, sA, sB)) continue;
			int n = -1;
			for (int k = 0; k < 4; k++) if (q.quad_nodes[k] == nid) { n = k; break; }
			if (n < 0) continue;
			double cn[3];
			quad_corner_normal(coords, q, n, cn);
			double l = std::sqrt(cn[0]*cn[0] + cn[1]*cn[1] + cn[2]*cn[2]);
			if (l < 1e-12) continue;
			for (int k = 0; k < 3; k++) cn[k] /= l;
			by_mat[q.matA].push_back({sA*cn[0], sA*cn[1], sA*cn[2]});
			by_mat[q.matB].push_back({sB*cn[0], sB*cn[1], sB*cn[2]});
		}
		double worst = 1e300;
		for (auto &kv : by_mat) {
			double m = cone_max_min(kv.second);
			if (m < worst) worst = m;
		}
		return (worst > 1e299) ? 1.0 : worst;
	};

	// Moving a node changes the corner normals of its ring too -- the corner normal at a
	// neighbour uses the edge back to this node -- so the monotonicity test has to cover the ring,
	// exactly as the convexity test covers every quad in the star rather than one corner.
	auto eval_cone_star = [&](int nid) -> double {
		double worst = eval_node_cone(nid);
		auto it = n2nbr.find(nid);
		if (it != n2nbr.end())
			for (int nb2 : it->second) {
				double m = eval_node_cone(nb2);
				if (m < worst) worst = m;
			}
		return worst;
	};

	// OFF by default, and the reason is a negative result worth keeping (section 15). Turned on at
	// PILLOW_CONE_OPEN=0 it does exactly what it promises -- every empty cone in the set
	// disappears, on every mesh, 6 -> 0 at ref3 -- and it changes the inverted count by ZERO:
	// ref3 stays at 8. That retires the attribution in sections 13.5 and 14 that the residual was
	// "forced by" those cones. It was correlation, not cause.
	//
	// What it does move, ref3, against the default: neighbour folds 25 -> 22, and dt case by case
	// -- hyeres 8.6x BETTER, hawaii 1.3x better, seven meshes unchanged, mauna_loa_small 2x worse
	// and Argostoli 14x worse (1.44e-06 -> 1.01e-07). That last one is why this is not the
	// default. At 0.05, where it also moves merely-tight cones, ref3 goes 8 -> 11 inverted.
	//
	// Kept behind the knob rather than deleted: on a mesh set where hyeres-like cases dominate and
	// Argostoli-like ones do not, PILLOW_CONE_OPEN=0 is a real option.
	double CONE_OPEN = -2.0;      // <= -1 never fires: cone margins are bounded below by -1
	if (const char *ev = getenv("PILLOW_CONE_OPEN")) CONE_OPEN = atof(ev);
	// Off has to mean off: the extra candidates and the extra acceptance term below change which
	// move the ORIGINAL reflex path picks even when no cone is ever out of tolerance (measured:
	// ref3 went 25 -> 28 neighbour folds with the threshold already disabled). So gate them.
	const bool cone_obj = (CONE_OPEN > -1.0);

	auto eval_star_worst_s = [&](int nid) -> double {
		double worst = 1e300;
		auto it = n2q.find(nid);
		if (it == n2q.end()) return 1.0;
		for (int iq : it->second) {
			double s = quad_min_corner_s(interface_quads[iq], coords);
			if (s < worst) worst = s;
		}
		return worst;
	};

	int n_convexified = 0;

	// Sweeps over interface nodes
	for (int sweep = 0; sweep < 6; sweep++) {
		int sweep_moves = 0;

		for (auto &kv : n2q) {
			int nid = kv.first;
			double px = coords[3*nid+0], py = coords[3*nid+1], pz = coords[3*nid+2];

			// Domain boundary check
			if (std::fabs(px - dom_min_x) < wall_x || std::fabs(px - dom_max_x) < wall_x) continue;
			if (std::fabs(py - dom_min_y) < wall_y || std::fabs(py - dom_max_y) < wall_y) continue;

			double worst_s_before = eval_star_worst_s(nid);
			double cone_before    = cone_obj ? eval_cone_star(nid) : 1.0;
			// Adjust nodes with a reflex corner (the original objective) or with a cone too tight
			// for a shared buffer node to fit (the new one). Either alone is reason to move.
			if (worst_s_before >= 0.05 && cone_before >= CONE_OPEN) continue;

			// Find worst quad in star to determine chord target
			int worst_qid = -1;
			double worst_q_s = 1e300;
			for (int iq : kv.second) {
				double s = quad_min_corner_s(interface_quads[iq], coords);
				if (s < worst_q_s) { worst_q_s = s; worst_qid = iq; }
			}
			if (worst_qid < 0) continue;

			const auto &wq = interface_quads[worst_qid];
			int rk = -1;
			for (int k = 0; k < 4; k++) {
				if (wq.quad_nodes[k] == nid) { rk = k; break; }
			}
			if (rk < 0) continue;

			double wp[4][3];
			for (int k = 0; k < 4; k++) {
				int qnid = wq.quad_nodes[k];
				for (int d = 0; d < 3; d++) wp[k][d] = coords[3*qnid+d];
			}

			// Candidate 1: Chord convexification target
			double m_chord[3] = {
				0.5 * (wp[(rk+3)%4][0] + wp[(rk+1)%4][0]),
				0.5 * (wp[(rk+3)%4][1] + wp[(rk+1)%4][1]),
				0.5 * (wp[(rk+3)%4][2] + wp[(rk+1)%4][2])
			};
			double opp[3] = { wp[(rk+2)%4][0], wp[(rk+2)%4][1], wp[(rk+2)%4][2] };
			double target1[3] = {
				m_chord[0] + 0.15 * (m_chord[0] - opp[0]),
				m_chord[1] + 0.15 * (m_chord[1] - opp[1]),
				m_chord[2] + 0.15 * (m_chord[2] - opp[2])
			};

			// Candidate 2: Star Laplacian centroid
			double target2[3] = { 0.0, 0.0, 0.0 };
			const auto &nbrs = n2nbr[nid];
			if (!nbrs.empty()) {
				for (int nbid : nbrs) {
					target2[0] += coords[3*nbid+0];
					target2[1] += coords[3*nbid+1];
					target2[2] += coords[3*nbid+2];
				}
				target2[0] /= nbrs.size();
				target2[1] /= nbrs.size();
				target2[2] /= nbrs.size();
			} else {
				target2[0] = target1[0]; target2[1] = target1[1]; target2[2] = target1[2];
			}

			// Candidate 3: Blend
			double target3[3] = {
				0.5 * (target1[0] + target2[0]),
				0.5 * (target1[1] + target2[1]),
				0.5 * (target1[2] + target2[2])
			};

			// Candidates 4-7: plain in-plane steps. The three above are aimed at convexity and
			// say nothing useful when the corner is already convex but the cone is tight, which
			// is the whole new case. A four-point pattern search is the cheap robust answer --
			// the cone's dependence on the node position has no closed form worth deriving for a
			// handful of nodes per mesh.
			const double step = 0.25 * iface_hmin[nid];
			double target4[3] = { px + step, py, pz };
			double target5[3] = { px - step, py, pz };
			double target6[3] = { px, py + step, pz };
			double target7[3] = { px, py - step, pz };

			const double (*cands[7])[3] = {
				(const double(*)[3])target3, (const double(*)[3])target1, (const double(*)[3])target2,
				(const double(*)[3])target4, (const double(*)[3])target5,
				(const double(*)[3])target6, (const double(*)[3])target7 };

			double best_tx = px, best_ty = py, best_tz = pz;
			double best_s_after = worst_s_before;
			double best_cone_after = cone_before;
			bool move_found = false;

			int gid = node_gts[nid];
			double max_d = 0.25 * iface_hmin[nid];

			for (int ic = 0; ic < (cone_obj ? 7 : 3); ic++) {
				double dx = (*cands[ic])[0] - px;
				double dy = (*cands[ic])[1] - py;
				double dist = std::sqrt(dx*dx + dy*dy);
				if (dist > max_d && dist > 1e-12) {
					dx *= max_d / dist;
					dy *= max_d / dist;
				}
				if (std::sqrt(dx*dx + dy*dy) < 1e-4 * iface_hmin[nid]) continue;

				for (double alpha : {1.0, 0.5, 0.25}) {
					double cur_tx = px + alpha * dx;
					double cur_ty = py + alpha * dy;
					double cur_tz = pz;
					bool proj_ok = false;

					if (gid >= 0) {
						if (eval_gts_height_near(mesh, gid, cur_tx, cur_ty, pz, cur_tz)) {
							proj_ok = true;
						}
					} else {
						proj_ok = true;
					}
					if (!proj_ok) continue;

					// Temporarily apply candidate
					coords[3*nid+0] = cur_tx;
					coords[3*nid+1] = cur_ty;
					coords[3*nid+2] = cur_tz;

					double cand_s    = eval_star_worst_s(nid);
					double cand_cone = cone_obj ? eval_cone_star(nid) : 1.0;

					// Two objectives, one monotonicity rule: a move must strictly improve one of
					// them and must not degrade the other. Opening a cone by creating a reflex
					// corner, or convexifying a quad by closing a neighbour's cone, both just move
					// the problem, and this is what stops it.
					bool better_s    = (cand_s > best_s_after + 1e-4)
					                && (!cone_obj || cand_cone >= best_cone_after - 1e-6);
					bool better_cone = cone_obj && (cand_cone > best_cone_after + 1e-6)
					                && (cand_s >= best_s_after - 1e-6);
					if (better_s || better_cone) {
						// Verify incident host elements validity
						bool host_ok = true;
						for (int ie : n2e[nid]) {
							if (check_elem_inverted(ie)) { host_ok = false; break; }
						}
						if (host_ok) {
							best_s_after = cand_s;
							best_cone_after = cand_cone;
							best_tx = cur_tx; best_ty = cur_ty; best_tz = cur_tz;
							move_found = true;
						}
					}

					// Revert for next evaluation
					coords[3*nid+0] = px;
					coords[3*nid+1] = py;
					coords[3*nid+2] = pz;

					if (move_found && best_s_after > 0.05) break; // Reached clean convexity
				}
				if (move_found && best_s_after > 0.05) break;
			}

			if (move_found) {
				coords[3*nid+0] = best_tx;
				coords[3*nid+1] = best_ty;
				coords[3*nid+2] = best_tz;
				sweep_moves++;
				n_convexified++;
			}
		}

		if (sweep_moves == 0) break;
	}

	return n_convexified;
}

void ApplyDoublePillowing(hexa_tree_t *mesh, std::vector<double> &coords, std::vector<int> &nodes_b_mat) {
	if (!mesh || mesh->elements.elem_count == 0) return;

	printf("    Applying Double-Layer Pillowing (Axis-Aligned Pure Face Extrusion)...\n");

	const int N_FACTOR = 2;

	// Buffer-layer thickness as a FRACTION of the physical edge from the interface
	// node toward the material-side interior node (found via the reference-mesh
	// topology, FaceNodesMap_inv). Placing the buffer along the real (deformed)
	// edge keeps every created hex valid on steep bathymetry, unlike an
	// axis-aligned fixed delta which inverts. 0<f<1; 0.25 = thin layer.
	const double PILLOW_FRACTION = 0.15;

	// Interpolation fraction for the NEW face/edge/vertex buffer placement (try_create_node_v2):
	// how far from the interface node toward its real full-step neighbour. 0.5 = exact midpoint;
	// user requested testing 0.45 (slightly toward the interface) after a visual check.
	const double V2_FRACTION = 0.45;

	// Step 7 pull-back: smallest fraction of the originally-placed buffer offset we accept
	// before giving up on a buffer node. NOT a numerical floor -- a visual one: below ~0.4 the
	// pillow sheet at the coastline reads as a crack/hole in the rendered mesh (user rejected
	// the 0.05 floor for exactly that). An inverted fat element is repairable downstream and
	// visible as geometry; a collapsed valid one is neither.
	const double ALPHA_MIN = 0.4;

	// 1. Scale integer lattice coordinates of all nodes and elements by N=2
	for (int ino = 0; ino < mesh->nodes.elem_count; ino++) {
		octant_node_t *node = (octant_node_t *) sc_array_index(&mesh->nodes, ino);
		node->x *= N_FACTOR;
		node->y *= N_FACTOR;
		node->z *= N_FACTOR;
	}
	for (int iel = 0; iel < mesh->elements.elem_count; iel++) {
		octant_t *elem = (octant_t *) sc_array_index(&mesh->elements, iel);
		elem->x *= N_FACTOR;
		elem->y *= N_FACTOR;
		elem->z *= N_FACTOR;
		for (int ino = 0; ino < 8; ino++) {
			elem->nodes[ino].x *= N_FACTOR;
			elem->nodes[ino].y *= N_FACTOR;
			elem->nodes[ino].z *= N_FACTOR;
		}
	}

	size_t initial_node_count = mesh->nodes.elem_count;
	size_t initial_elem_count = mesh->elements.elem_count;

	int min_gx = mesh->x_start * N_FACTOR;
	int max_gx = mesh->x_end * N_FACTOR;
	int min_gy = mesh->y_start * N_FACTOR;
	int max_gy = mesh->y_end * N_FACTOR;
	int min_gz = 0;
	int max_gz = 3 * mesh->max_z * N_FACTOR;

	// Position -> node id lookup over ORIGINAL (pre-pillow) nodes only, built once, never
	// updated. Used to find the real mesh neighbour a buffer node's physical position should be
	// interpolated toward (see get_or_create_node below). Deliberately NOT extended to buffer
	// nodes as they're created (unlike the "NO position dedup" note on buffer nodes further
	// down, which is about a different, already-diagnosed failure mode): confirmed by direct
	// measurement on 2026-08-19 that original nodes have zero position collisions in this mesh
	// (549917/549917 unique), so this lookup is unambiguous by construction, and since it never
	// contains buffer nodes there's nothing to keep in sync as new ones are created.
	std::unordered_map<int_triple, int, int_triple_hash> orig_pos_hash;
	orig_pos_hash.reserve(initial_node_count);
	for (size_t i = 0; i < initial_node_count; i++) {
		octant_node_t *nd = (octant_node_t *) sc_array_index(&mesh->nodes, i);
		orig_pos_hash[{nd->x, nd->y, nd->z}] = (int)i;
	}

	std::vector<interface_quad_t> interface_quads;

	struct face_ref {
		int elem_id;
		int face_idx;
	};

	// 1b. Material despeckle. A one-cell-wide spike or notch in the material
	// field gives its interface nodes incident quad normals on OPPOSING axes
	// (+x and -x), so the feasible cone for the single shared buffer node is
	// empty and the incident pillow hexes cannot all be valid -- no amount of
	// untangling fixes them. Flip any element that disagrees with 5+ of its 6
	// face neighbours; two rounds clear nested spikes.
	{
		std::unordered_map<face_key, face_ref, face_key_hash> nmap;
		nmap.reserve(mesh->elements.elem_count * 4);
		std::vector<std::array<int,6>> nb(mesh->elements.elem_count);
		for (auto &a : nb) a.fill(-1);
		for (int iel = 0; iel < mesh->elements.elem_count; iel++) {
			octant_t *elem = (octant_t *) sc_array_index(&mesh->elements, iel);
			for (int iface = 0; iface < 6; iface++) {
				int fn[4];
				for (int n = 0; n < 4; n++) fn[n] = elem->nodes[FaceNodesMap[iface][n]].id;
				face_key key = make_face_key(fn);
				auto it = nmap.find(key);
				if (it == nmap.end()) nmap[key] = { iel, iface };
				else {
					nb[iel][iface] = it->second.elem_id;
					nb[it->second.elem_id][it->second.face_idx] = iel;
					nmap.erase(it);
				}
			}
		}
		int nflip_total = 0;
		for (int round = 0; round < 2; round++) {
			std::vector<int> flip;
			for (int iel = 0; iel < mesh->elements.elem_count; iel++) {
				octant_t *elem = (octant_t *) sc_array_index(&mesh->elements, iel);
				int ndiff = 0;
				for (int f = 0; f < 6; f++) {
					int j = nb[iel][f];
					if (j < 0) continue;
					octant_t *o = (octant_t *) sc_array_index(&mesh->elements, j);
					if (o->n_mat != elem->n_mat) ndiff++;
				}
				if (ndiff >= 4) flip.push_back(iel);
			}
			if (flip.empty()) break;
			for (int iel : flip) {
				octant_t *elem = (octant_t *) sc_array_index(&mesh->elements, iel);
				// Reassign to whichever material is the majority among the (up to 6)
				// face neighbours, ties broken by lowest n_mat -- a fixed 0<->1 toggle
				// only makes sense with exactly two materials.
				std::unordered_map<int,int> count;
				for (int f = 0; f < 6; f++) {
					int j = nb[iel][f];
					if (j < 0) continue;
					octant_t *o = (octant_t *) sc_array_index(&mesh->elements, j);
					count[o->n_mat]++;
				}
				int best_mat = elem->n_mat, best_count = -1;
				for (auto &kv : count) {
					if (kv.second > best_count || (kv.second == best_count && kv.first < best_mat)) {
						best_count = kv.second; best_mat = kv.first;
					}
				}
				elem->n_mat = best_mat;
			}
			nflip_total += (int)flip.size();
		}
		printf("    Material despeckle: %d one-cell spikes/notches flipped\n", nflip_total);

		// 1c. Node-level pinch despeckle. The element rule above only catches a cell that
		// disagrees with 4+ of its 6 face neighbours, and on real bathymetry it fires zero
		// times. The configuration that actually leaves an EMPTY FEASIBLE CONE is thinner than
		// that: a material reaching a node from two directions that are not face-connected
		// around that node -- a step corner meeting diagonally, a one-cell notch at the
		// coastline. Such a node's incident quads carry normals on both ends of an axis, so the
		// single shared buffer node cannot sit inside every incident hex, at any thickness or
		// direction. Measured after the extrusion field and the local repair: 20 of 20 residual
		// inverted elements on Argostoli ref3 and on hyeres ref3 touch one of these nodes. It is
		// the last family and it is topological -- the fix is to remove the pinch from the
		// material field, one cell at a time, before the interface quads are built.
		// Tried and rejected: replacing this face-connectivity proxy with the exact cone test
		// (cone_max_min over the incident interface quads' corner normals at each node, which is
		// literally "can one shared buffer node serve them all"). At a threshold of 0 it fires
		// 122 times on Argostoli ref3 against this rule's 40 and makes things WORSE -- 4 -> 13
		// inverted -- because reassigning a whole cell's material to close a cone shreds the
		// coastline. At any threshold below -0.05 it never fires at all. Both facts say the same
		// thing: the cones that matter are empty by less than 0.05, i.e. within ~3 degrees of
		// feasible. They are marginal geometry, not material pinches, and material is the wrong
		// lever for them. See section 14.
		int npinch_total = 0;
		// A cell is reassigned at most once: flipping the small component can create a pinch on
		// the other side, and without this the same two cells swap material every round forever.
		std::unordered_set<int> flipped_once;
		for (int round = 0; round < 8; round++) {
			// nodes sitting on a material interface
			std::unordered_set<int> iface_nodes;
			for (int iel = 0; iel < (int)mesh->elements.elem_count; iel++) {
				octant_t *elem = (octant_t *) sc_array_index(&mesh->elements, iel);
				for (int f = 0; f < 6; f++) {
					int j = nb[iel][f];
					if (j < 0) continue;
					octant_t *o = (octant_t *) sc_array_index(&mesh->elements, j);
					if (o->n_mat == elem->n_mat) continue;
					for (int n = 0; n < 4; n++) iface_nodes.insert(elem->nodes[FaceNodesMap[f][n]].id);
				}
			}
			// every cell incident to those nodes (not just the interface-adjacent ones: a cell
			// deeper in the material is what connects two apparent components)
			std::unordered_map<int, std::vector<int>> n2e;
			for (int iel = 0; iel < (int)mesh->elements.elem_count; iel++) {
				octant_t *elem = (octant_t *) sc_array_index(&mesh->elements, iel);
				for (int k = 0; k < 8; k++) {
					int nid = elem->nodes[k].id;
					if (iface_nodes.count(nid)) n2e[nid].push_back(iel);
				}
			}
			std::unordered_set<int> flip;      // cells to reassign this round
			for (auto &kv : n2e) {
				const std::vector<int> &cells = kv.second;
				if (cells.size() < 3) continue;
				std::unordered_map<int, std::vector<int>> by_mat;
				for (int ie : cells) {
					octant_t *e = (octant_t *) sc_array_index(&mesh->elements, ie);
					by_mat[e->n_mat].push_back(ie);
				}
				if (by_mat.size() < 2) continue;
				for (auto &mk : by_mat) {
					const std::vector<int> &S = mk.second;
					if (S.size() < 2) continue;          // a single cell cannot be pinched
					// connected components of S under face adjacency, restricted to S
					std::unordered_set<int> inS(S.begin(), S.end());
					std::unordered_set<int> seen;
					std::vector<std::vector<int>> comp;
					for (int seed : S) {
						if (seen.count(seed)) continue;
						std::vector<int> stack(1, seed), cur;
						seen.insert(seed);
						while (!stack.empty()) {
							int ie = stack.back(); stack.pop_back();
							cur.push_back(ie);
							for (int f = 0; f < 6; f++) {
								int j = nb[ie][f];
								if (j < 0 || !inS.count(j) || seen.count(j)) continue;
								seen.insert(j); stack.push_back(j);
							}
						}
						comp.push_back(cur);
					}
					if (comp.size() < 2) continue;       // face-connected around the node: fine
					// pinch: keep the largest component, flip the others
					size_t big = 0;
					for (size_t i = 1; i < comp.size(); i++)
						if (comp[i].size() > comp[big].size()) big = i;
					for (size_t i = 0; i < comp.size(); i++)
						if (i != big) for (int ie : comp[i])
							if (!flipped_once.count(ie)) flip.insert(ie);
				}
			}

			if (flip.empty()) break;
			for (int iel : flip) {
				octant_t *elem = (octant_t *) sc_array_index(&mesh->elements, iel);
				std::unordered_map<int,int> count;    // majority of the face neighbours, as above
				for (int f = 0; f < 6; f++) {
					int j = nb[iel][f];
					if (j < 0) continue;
					octant_t *o = (octant_t *) sc_array_index(&mesh->elements, j);
					if (o->n_mat != elem->n_mat) count[o->n_mat]++;
				}
				int best_mat = elem->n_mat, best_count = -1;
				for (auto &kv : count) {
					if (kv.second > best_count || (kv.second == best_count && kv.first < best_mat)) {
						best_count = kv.second; best_mat = kv.first;
					}
				}
				elem->n_mat = best_mat;
				flipped_once.insert(iel);
			}
			npinch_total += (int)flip.size();
		}
		printf("    Node-pinch despeckle: %d cells reassigned (empty-cone nodes)\n", npinch_total);
	}

	// 2. Closed Manifold Scan: Pair ALL faces across mesh->elements

	std::unordered_map<face_key, face_ref, face_key_hash> fmap;
	fmap.reserve(mesh->elements.elem_count * 4);

	for (int iel = 0; iel < mesh->elements.elem_count; iel++) {
		octant_t *elem = (octant_t *) sc_array_index(&mesh->elements, iel);
		for (int iface = 0; iface < 6; iface++) {
			int fn[4];
			for (int n = 0; n < 4; n++) {
				fn[n] = elem->nodes[FaceNodesMap[iface][n]].id;
			}
			face_key key = make_face_key(fn);
			auto it = fmap.find(key);
			if (it == fmap.end()) {
				fmap[key] = { iel, iface };
			} else {
				octant_t *prev_elem = (octant_t *) sc_array_index(&mesh->elements, it->second.elem_id);
				if (elem->n_mat != prev_elem->n_mat) {
					interface_quad_t quad;
					// Lower n_mat is always "A", regardless of scan order -- keeps
					// interface_quads (and everything derived from it) reproducible
					// across runs of the same input.
					bool curA = (elem->n_mat < prev_elem->n_mat);
					if (curA) {
						quad.elemA_id = elem->id;      quad.matA = elem->n_mat;
						quad.elemB_id = prev_elem->id; quad.matB = prev_elem->n_mat;
						quad.faceA = iface;
						quad.faceB = it->second.face_idx;
					} else {
						quad.elemA_id = prev_elem->id; quad.matA = prev_elem->n_mat;
						quad.elemB_id = elem->id;      quad.matB = elem->n_mat;
						quad.faceA = it->second.face_idx;
						quad.faceB = iface;
					}
					quad.normal[0] = FaceNormal[quad.faceA][0];
					quad.normal[1] = FaceNormal[quad.faceA][1];
					quad.normal[2] = FaceNormal[quad.faceA][2];

					octant_t *eA = (octant_t *) sc_array_index(&mesh->elements, quad.elemA_id);
					for (int n = 0; n < 4; n++) {
						quad.quad_nodes[n] = eA->nodes[FaceNodesMap[quad.faceA][n]].id;
					}
					interface_quads.push_back(quad);
				}
				fmap.erase(it);
			}
		}
	}

	printf("    Pillowing scan: scanned %d elements, %zu complete interface quads found\n",
	       (int)mesh->elements.elem_count, interface_quads.size());

	if (interface_quads.empty()) return;

	// Pre-extrusion convexification of interface quads on GTS surface
	if (getenv("NO_CONVEXIFY") == NULL) {
		double vol_sum = 0.0;
		double X[8], Y[8], Z[8];
		for (size_t ie = 0; ie < mesh->elements.elem_count; ie++) {
			load_elem_xyz(mesh, coords, (int)ie, X, Y, Z);
			vol_sum += mgeom::hex_signed_volume(X, Y, Z);
		}
		int ref_host = mgeom::reference_sign(vol_sum);
		int n_conv = ConvexifyInterfaceQuads(mesh, coords, interface_quads, ref_host);
		printf("    Pre-extrusion quad convexification: %d node adjustments on surface\n", n_conv);
	}

	// 3. Helper to create unique buffer nodes with pure axis-aligned offset in `coords`

	// Per-buffer-node bookkeeping for the adaptive shrink pass (step 7, after element
	// creation): the source interface node, the UNSCALED averaged offset direction*magnitude
	// (ox,oy,oz -- what PILLOW_FRACTION multiplies), and a per-node scale in (0,1] applied on
	// top of PILLOW_FRACTION, so a buffer node whose resulting pillow hex turns out invalid or
	// folded onto its parent can be pulled back toward the interface without touching
	// connectivity (safe/conforming by construction -- only ever moves a position, same
	// property the untangler relies on).
	std::unordered_map<int, int> buffer_orig;
	std::unordered_map<int, std::array<double,3>> buffer_dir;
	std::unordered_map<int, double> buffer_scale;

	// NO position dedup. Pass B calls the two functions below exactly once per (interface node,
	// material side), which is what keeps the layer conforming. Merging by integer target
	// instead fused buffer nodes of DISTINCT interface nodes wherever a bathymetry step or a
	// domain-edge clamp put two targets in the same cell, giving the incident original hex two
	// identical corner ids -- degenerate, and unfixable by any node movement. That was every one
	// of the surviving inverted originals.
	auto push_buffer_node = [&](int orig_node_id, int target_x, int target_y, int target_z,
	                             double px_new, double py_new, double pz_new,
	                             double ox, double oy, double oz) -> int {
		octant_node_t *orig_node = (octant_node_t *) sc_array_index(&mesh->nodes, orig_node_id);
		int orig_color = orig_node->color;

		octant_node_t *new_node = (octant_node_t *) sc_array_push(&mesh->nodes);
		int new_id = (int)(mesh->nodes.elem_count - 1);
		new_node->id = new_id;
		new_node->x = target_x;
		new_node->y = target_y;
		new_node->z = target_z;
		new_node->fixed = 0;
		new_node->color = orig_color;

		coords.push_back(px_new);
		coords.push_back(py_new);
		coords.push_back(pz_new);

		buffer_orig[new_id] = orig_node_id;
		buffer_dir[new_id] = {ox, oy, oz};
		buffer_scale[new_id] = 1.0;

		return new_id;
	};

	// ORIGINAL placement: integer target = one lattice step along the face normal (unchanged,
	// used only for dedup/topology bookkeeping); physical = interface node + fraction f of the
	// AVERAGED interior offset (mean of interior-minus-interface over all incident interface
	// faces on this material side). Averaging aligns neighbouring buffers into a congruent slab,
	// so the pillow hex does not twist on steep bathymetry -- kept as the fallback for the cases
	// the new face/edge/vertex placement below can't handle (genuine same-axis opposition, or a
	// missing real neighbour to interpolate toward).
	auto get_or_create_node = [&](int orig_node_id, int offset_dx, int offset_dy, int offset_dz,
	                              double ox, double oy, double oz) -> int {
		octant_node_t *orig_node = (octant_node_t *) sc_array_index(&mesh->nodes, orig_node_id);
		int target_x = std::clamp(orig_node->x + offset_dx, min_gx, max_gx);
		int target_y = std::clamp(orig_node->y + offset_dy, min_gy, max_gy);
		int target_z = std::clamp(orig_node->z + offset_dz, min_gz, max_gz);

		double px = coords[3 * orig_node_id + 0];
		double py = coords[3 * orig_node_id + 1];
		double pz = coords[3 * orig_node_id + 2];
		double px_new = px + PILLOW_FRACTION * ox;
		double py_new = py + PILLOW_FRACTION * oy;
		double pz_new = pz + PILLOW_FRACTION * oz;

		return push_buffer_node(orig_node_id, target_x, target_y, target_z, px_new, py_new, pz_new, ox, oy, oz);
	};

	// NEW placement: classify the (interface node, material side) by how many distinct AXES its
	// incident quads' normals span (1 = face, 2 = edge, 3 = vertex -- true same-axis opposition,
	// e.g. a 1-cell spike the despeckle pass missed, is a separate, already-diagnosed case, not
	// handled here). Integer target = orig + one +-1 step per involved axis (the odd position
	// between two even original nodes). Physical target = the REAL mesh node a full 2-step away
	// in that same direction (found via orig_pos_hash, confirmed collision-free on this mesh),
	// interpolated at t=0.5 -- not the old thin PILLOW_FRACTION=0.15 sliver: since this now
	// anchors to an actual neighbour instead of an averaged, possibly-inconsistent direction,
	// a well-proportioned half-size element is more robust (more slack before a small direction
	// error flips its Jacobian sign) than a thin one, and leaves less repair work for the
	// optimizer afterward. Returns -1 if this (node, side) can't be placed this way (caller
	// falls back to get_or_create_node above): true opposition, or the full-step neighbour
	// doesn't exist (e.g. near a refinement/domain edge).
	auto try_create_node_v2 = [&](int orig_node_id, uint8_t mask) -> int {
		if (((mask & 3) == 3) || ((mask & 12) == 12) || ((mask & 48) == 48)) return -1; // true opposition
		int dx = (mask & 1) ? 1 : (mask & 2) ? -1 : 0;
		int dy = (mask & 4) ? 1 : (mask & 8) ? -1 : 0;
		int dz = (mask & 16) ? 1 : (mask & 32) ? -1 : 0;
		if (dx == 0 && dy == 0 && dz == 0) return -1; // no incident quad recorded (shouldn't happen)

		octant_node_t *orig_node = (octant_node_t *) sc_array_index(&mesh->nodes, orig_node_id);
		int target_x = std::clamp(orig_node->x + dx, min_gx, max_gx);
		int target_y = std::clamp(orig_node->y + dy, min_gy, max_gy);
		int target_z = std::clamp(orig_node->z + dz, min_gz, max_gz);

		int_triple neighbor_key = { orig_node->x + 2*dx, orig_node->y + 2*dy, orig_node->z + 2*dz };
		auto it = orig_pos_hash.find(neighbor_key);
		if (it == orig_pos_hash.end()) return -1; // no real neighbour there, fall back

		int neighbor_id = it->second;
		double px = coords[3*orig_node_id+0], py = coords[3*orig_node_id+1], pz = coords[3*orig_node_id+2];
		double nx_ = coords[3*neighbor_id+0], ny_ = coords[3*neighbor_id+1], nz_ = coords[3*neighbor_id+2];
		double px_new = px + V2_FRACTION*(nx_-px), py_new = py + V2_FRACTION*(ny_-py), pz_new = pz + V2_FRACTION*(nz_-pz);
		double ox = nx_-px, oy = ny_-py, oz = nz_-pz; // unscaled, for buffer_dir bookkeeping only

		return push_buffer_node(orig_node_id, target_x, target_y, target_z, px_new, py_new, pz_new, ox, oy, oz);
	};

	// 4. Map every interface node, per material side it borders, to its buffer node.
	// Keyed by (node id, material id) instead of a fixed mat0/mat1 pair, so a node
	// sitting at a junction of 3+ materials gets one buffer per material side, not
	// just two.
	auto pillow_key = [](int node_id, int mat_id) -> int64_t {
		return ((int64_t)node_id << 32) | (uint32_t)mat_id;
	};
	std::unordered_map<int64_t, int> pillow_map;

	// Pass A: accumulate, per (interface node, material side), the interior
	// offset (interior-minus-interface) over ALL incident interface faces that
	// border that material. {sum_x, sum_y, sum_z, count, sum_magnitude}.
	std::unordered_map<int64_t, std::array<double,5>> acc;
	auto accumulate = [&](std::array<double,5>& a, int inode, int nid) {
		double dx=coords[3*inode+0]-coords[3*nid+0];
		double dy=coords[3*inode+1]-coords[3*nid+1];
		double dz=coords[3*inode+2]-coords[3*nid+2];
		a[0]+=dx; a[1]+=dy; a[2]+=dz; a[3]+=1.0; a[4]+=std::sqrt(dx*dx+dy*dy+dz*dz);
	};
	for (auto &q : interface_quads) {
		octant_t *eA = (octant_t *) sc_array_index(&mesh->elements, q.elemA_id);
		octant_t *eB = (octant_t *) sc_array_index(&mesh->elements, q.elemB_id);
		for (int n = 0; n < 4; n++) {
			int nid = q.quad_nodes[n];
			accumulate(acc[pillow_key(nid, q.matA)], eA->nodes[FaceNodesMap_inv[q.faceA][n]].id, nid);
			accumulate(acc[pillow_key(nid, q.matB)], eB->nodes[FaceNodesMap_inv[q.faceB][n]].id, nid);
		}
	}

	// Per-(node, material side) bitmask of which of the 6 axis directions (+x,-x,+y,-y,+z,-z,
	// bits 1/2/4/8/16/32) are represented among incident quads -- SIGNED per side (matA gets
	// -normal's bit, matB gets +normal's bit, matching get_or_create_node's existing -nx.../
	// +nx... convention), so it directly says which way THIS side's buffer should extrude.
	// Feeds two things: (a) try_create_node_v2's face/edge/vertex placement below, (b) the
	// opposing-normals warning (a (node,mat) side whose bits include BOTH signs on the same
	// axis has an EMPTY feasible cone: no single buffer node can keep every incident pillow hex
	// valid, same-axis 180-degree opposition -- e.g. a 1-cell spike the despeckle pass missed --
	// is not fixable by the face/edge/vertex scheme either, since that's not a face/edge/vertex
	// pattern at all; those cases fall back to the old averaged-direction placement).
	std::unordered_map<int64_t, uint8_t> side_mask;
	for (auto &q : interface_quads) {
		uint8_t bA = 0, bB = 0;
		int nx = q.normal[0], ny = q.normal[1], nz = q.normal[2];
		if (-nx > 0) bA |= 1; if (-nx < 0) bA |= 2;
		if (-ny > 0) bA |= 4; if (-ny < 0) bA |= 8;
		if (-nz > 0) bA |= 16; if (-nz < 0) bA |= 32;
		if (nx > 0) bB |= 1; if (nx < 0) bB |= 2;
		if (ny > 0) bB |= 4; if (ny < 0) bB |= 8;
		if (nz > 0) bB |= 16; if (nz < 0) bB |= 32;
		for (int n = 0; n < 4; n++) {
			side_mask[pillow_key(q.quad_nodes[n], q.matA)] |= bA;
			side_mask[pillow_key(q.quad_nodes[n], q.matB)] |= bB;
		}
	}
	{
		int nopp = 0;
		for (auto &kv : side_mask) {
			uint8_t m = kv.second;
			if (((m & 3) == 3) || ((m & 12) == 12) || ((m & 48) == 48)) nopp++;
		}
		if (nopp > 0) {
			printf("    WARNING: %d / %zu (node, material) sides have opposing quad normals "
			       "(empty feasible cone) -- their pillow hexes cannot all be valid\n",
			       nopp, side_mask.size());
			// which nodes they are, so the residual inversions can be attributed to them
			// (same role as the invmap_*.csv dumps): node id, material side, axis mask, position
			FILE *fc = fopen("pillow_cone.csv", "w");
			if (fc) {
				fprintf(fc, "node,mat,mask,x,y,z\n");
				for (auto &kv : side_mask) {
					uint8_t m = kv.second;
					if (!(((m & 3) == 3) || ((m & 12) == 12) || ((m & 48) == 48))) continue;
					int nid = (int)(kv.first >> 32), mat = (int)(uint32_t)kv.first;
					fprintf(fc, "%d,%d,%d,%.3f,%.3f,%.3f\n", nid, mat, (int)m,
					        coords[3*nid+0], coords[3*nid+1], coords[3*nid+2]);
				}
				fclose(fc);
			}
		}
	}

	// --- Geometric extrusion field: direction and length of every buffer node ---------------
	// Direction = average of the incident quads' GEOMETRIC normals, oriented into the material
	// side, then Laplacian-smoothed over the interface graph; length = ETA times the smaller of
	// the shortest incident quad edge and the distance to the host element's interior node, so
	// the thickness follows the local element size instead of the distance to some lattice
	// neighbour. Both limits are needed: the in-plane one keeps the layer from being wider than
	// the quad, the interior one keeps it from punching through the host element, which is thin
	// in z near the surface while the interface quad can be kilometres wide (hawaii, mauna_loa). Both reference implementations do exactly this
	// (HexGen_Hex2Spline's Pillow(), Marechal's boundary layers, IMR 2016), and it targets the
	// two signatures measured on the inverted elements of this mesh (tools/pillow_report.py):
	// a pillow hex is never inverted while its four corner offsets agree within 60 degrees, and
	// the inversion rate climbs with the thickness spread inside a single hex. The old
	// lattice-based placement stays behind PILLOW_LEGACY=1 as an A/B control.
	const bool legacy_placement = (getenv("PILLOW_LEGACY") != NULL);
	double PILLOW_ETA = 0.45;                 // thickness as a fraction of the shortest quad edge
	                                          // (swept 0.45-0.75 in section 16: 0.45 wins, raising
	                                          // it fails the inverted-count constraint everywhere)
	if (const char *ev = getenv("PILLOW_ETA")) PILLOW_ETA = atof(ev);
	const int    DIR_SMOOTH_ITERS = 3;
	const double DIR_SMOOTH_W     = 0.5;
	// Floor for the host-room cap, as a fraction of the size-based nominal thickness. 0.2, not the
	// original 0.4: with the interior limit actually working (the sqrt below), a 16-point sweep of
	// ETA x LEN_FLOOR over the 11 ref3 meshes puts the only passing points at ETA 0.45 with a floor
	// of 0.2 or 0.1, and 0.2 dominates 0.1 on every mesh. Raising the floor is monotonically bad
	// -- 0.4/0.6/0.8 give 10/21/30 inverted -- which makes sense: the floor forces thickness back
	// exactly where the cap has just concluded there is no room, so it pushes host elements out.
	// Raising ETA is bad too: 0.55, 0.65 and 0.75 fail at every floor. Section 16 has the grid.
	double LEN_FLOOR              = 0.2;
	if (const char *ev = getenv("PILLOW_LEN_FLOOR")) LEN_FLOOR = atof(ev);

	std::unordered_map<int64_t, std::array<double,3>> ext_dir;
	std::unordered_map<int64_t, double> ext_len;
	std::unordered_map<int64_t, std::vector<int64_t>> ext_nbr;
	// The validity constraint the direction has to satisfy, collected alongside the direction
	// itself: for each interface quad incident to this (node, material side), the CORNER normal
	// at this node, e1 x e2 of the two in-plane edges, oriented into the side. The pillow hex's
	// corner Jacobian here is exactly (e1 x e2) . d, so the layer is valid at this corner only
	// while cn . d > 0. Per corner, not per quad: on a warped or reflex quad the four corner
	// normals do not agree with the quad's diagonal normal, and it is the corner ones that decide.
	std::unordered_map<int64_t, std::vector<std::array<double,3>>> ext_cn;
	{
		auto add_side = [&](int64_t key, const double d[3], const double cn[3], double len,
		                    int64_t n1, int64_t n2) {
			auto &e = ext_dir[key];
			e[0] += d[0]; e[1] += d[1]; e[2] += d[2];
			double cl = std::sqrt(cn[0]*cn[0] + cn[1]*cn[1] + cn[2]*cn[2]);
			if (cl > 1e-12) ext_cn[key].push_back({cn[0]/cl, cn[1]/cl, cn[2]/cl});
			// Tried and rejected: also pushing `d` (the quad's diagonal normal, oriented into
			// this side) as a second constraint, to couple the four corners of a quad against
			// each other rather than only against their own corner normals. It tightens the
			// cone without buying anything -- same inverted count over the 11 ref3 meshes, and
			// it costs dt: kefalonia 4.2e-5 -> 2.3e-6 s, hyeres 7.2e-7 -> 2.8e-7, and Argostoli
			// goes from 5 to 14 neighbour folds. The per-corner constraint is the one that pays.
			auto it = ext_len.find(key);
			if (it == ext_len.end() || len < it->second) ext_len[key] = len;
			std::vector<int64_t> &nb = ext_nbr[key];
			nb.push_back(n1); nb.push_back(n2);
		};
		for (auto &q : interface_quads) {
			const double *p[4];
			for (int n = 0; n < 4; n++) p[n] = &coords[3*q.quad_nodes[n]];
			double emin = 1e300;
			for (int n = 0; n < 4; n++) {
				const double *a = p[n], *b = p[(n+1)%4];
				double e = std::sqrt((a[0]-b[0])*(a[0]-b[0]) + (a[1]-b[1])*(a[1]-b[1]) +
				                     (a[2]-b[2])*(a[2]-b[2]));
				if (e < emin) emin = e;
			}
			double ng[3]; int sA, sB;
			if (!quad_side_signs(mesh, coords, q, ng, sA, sB) || emin >= 1e299) continue;
			octant_t *eA = (octant_t *) sc_array_index(&mesh->elements, q.elemA_id);
			octant_t *eB = (octant_t *) sc_array_index(&mesh->elements, q.elemB_id);
			double dA[3], dB[3];
			for (int k = 0; k < 3; k++) { dA[k] = sA*ng[k]; dB[k] = sB*ng[k]; }
			for (int n = 0; n < 4; n++) {
				int nid = q.quad_nodes[n];
				int nb1 = q.quad_nodes[(n+1)%4], nb3 = q.quad_nodes[(n+3)%4];
				// this corner's own normal, e1 x e2 of the two in-plane edges, flipped into each
				// material side by the same sign the side direction got. Taken raw: on a reflex
				// corner it genuinely opposes the quad's diagonal normal, and that is the
				// information the cone repair below needs, not something to normalise away.
				double cn[3], cnA[3], cnB[3];
				quad_corner_normal(coords, q, n, cn);
				for (int k = 0; k < 3; k++) { cnA[k] = sA*cn[k]; cnB[k] = sB*cn[k]; }
				// how much room this corner actually has on each side: the edge from the
				// interface node to the host element's own interior node
				int inA = eA->nodes[FaceNodesMap_inv[q.faceA][n]].id;
				int inB = eB->nodes[FaceNodesMap_inv[q.faceB][n]].id;
				// DISTANCE to the host interior node, not its square. This accumulated a*a and
				// was then compared against emin, a length -- a dimensional error that silently
				// disabled the whole interior limit on every metre-scale mesh (hA ~ 1e4 against
				// emin ~ 1e2, so the min always picked emin) while inverting the meaning on the
				// sub-metre `plane` case (hA ~ 1e-4 < emin, layer ~100x thinner than designed).
				double hA = 0.0, hB = 0.0;
				for (int k = 0; k < 3; k++) {
					double a = coords[3*inA+k] - p[n][k], b = coords[3*inB+k] - p[n][k];
					hA += a*a; hB += b*b;
				}
				hA = std::sqrt(hA); hB = std::sqrt(hB);
				add_side(pillow_key(nid, q.matA), dA, cnA, PILLOW_ETA*std::min(emin, hA),
				         pillow_key(nb1, q.matA), pillow_key(nb3, q.matA));
				add_side(pillow_key(nid, q.matB), dB, cnB, PILLOW_ETA*std::min(emin, hB),
				         pillow_key(nb1, q.matB), pillow_key(nb3, q.matB));
			}
		}
		auto renorm = [](std::array<double,3> &v) -> bool {
			double l = std::sqrt(v[0]*v[0]+v[1]*v[1]+v[2]*v[2]);
			if (l < 1e-12) return false;
			v[0] /= l; v[1] /= l; v[2] /= l;
			return true;
		};
		for (auto &kv : ext_dir) renorm(kv.second);
		// Laplacian smoothing over the interface graph, per material side -- the step that
		// removes the >60 deg disagreement between the four corners of a quad.
		for (int it = 0; it < DIR_SMOOTH_ITERS; it++) {
			std::unordered_map<int64_t, std::array<double,3>> next = ext_dir;
			for (auto &kv : ext_dir) {
				auto nb = ext_nbr.find(kv.first);
				if (nb == ext_nbr.end() || nb->second.empty()) continue;
				std::array<double,3> avg = {0.0, 0.0, 0.0};
				int cnt = 0;
				for (size_t j = 0; j < nb->second.size(); j++) {
					auto it2 = ext_dir.find(nb->second[j]);
					if (it2 == ext_dir.end()) continue;
					avg[0] += it2->second[0]; avg[1] += it2->second[1]; avg[2] += it2->second[2];
					cnt++;
				}
				if (cnt == 0) continue;
				std::array<double,3> v;
				for (int k = 0; k < 3; k++)
					v[k] = (1.0-DIR_SMOOTH_W)*kv.second[k] + DIR_SMOOTH_W*avg[k]/(double)cnt;
				if (renorm(v)) next[kv.first] = v;
			}
			ext_dir.swap(next);
		}

		// --- cone repair: put the direction back inside the feasible set ---------------------
		// Neither of the two steps above knows about signs. The accumulator sums raw normals, so
		// at a coastline step where a node's incident corner normals on the SAME material side
		// oppose each other the sum cancels and renorm leaves a near-zero vector in the map --
		// which the Laplacian then spreads to the neighbours at weight 0.5 for three rounds. And
		// the Laplacian itself averages across sharp features and across the two shores of a
		// one-cell channel, which are edge-neighbours in the interface graph while carrying
		// nearly opposite normals. Either way the result can leave the feasible cone and put the
		// buffer node on the WRONG SIDE of its own interface quad.
		//
		// This was the whole residual (tools/pillow_residual.py on the previous run): 100 % of the
		// inverted pillow hexes had an offset more than 45 degrees -- median 110 to 135, up to 165
		// -- from a direction that would have worked, while a valid direction existed for every
		// single residual quad with a margin of ~0.99. Inversion rate on the hexes whose four
		// offsets straddle their own quad: 26-74 %, against 0.000-0.076 % on all the others; the
		// two meshes with no such hex had no inverted element at all. The legacy path already
		// guards this (avg_offset zeroes opposing axes and falls back); the field path lost it.
		//
		// The fix is the constrained problem stated honestly: maximise the worst corner margin
		// min_q (cn_q . d) over unit d -- the Chebyshev centre of the spherical cone. With a
		// handful of quads meeting at a node the optimum is supported by one, two or three of the
		// constraints, so enumerating the normals, their pairwise and their triple sums reaches
		// it. Applied ONLY where the field has actually left the cone, so the smoothing of step
		// 3 -- which is what removed the >60 degree disagreement in the first place -- is kept
		// everywhere it is doing its job. Where the cone is genuinely empty (a one-cell notch the
		// node-pinch despeckle of section 5 did not reach) no direction scores positive and the
		// best available is still the least bad one.
		double CONE_MARGIN = 0.05;            // ~3 deg of slack, calibration knob
		if (const char *ev = getenv("PILLOW_CONE_MARGIN")) CONE_MARGIN = atof(ev);
		// PILLOW_CONE_SELFCHECK=1 brute-forces the same max-min over a 2000-point Fibonacci
		// sphere and reports the worst amount by which cone_max_min above was beaten. That
		// enumeration claims the optimum is supported by one, two or three constraints; if that
		// were wrong the sphere would find a better direction and the gap would be large. The
		// sphere itself only resolves to ~0.05 rad, so a few 1e-3 either way is the sampling,
		// not an error. Same role as UNTANGLE_SELFCHECK in pillow_untangle.cpp.
		const bool cone_selfcheck = (getenv("PILLOW_CONE_SELFCHECK") != NULL);
		double worst_gap = 0.0;
		int n_reaimed = 0, n_empty_cone = 0;
		for (auto &kv : ext_dir) {
			auto itc = ext_cn.find(kv.first);
			if (itc == ext_cn.end() || itc->second.empty()) continue;
			const std::vector<std::array<double,3>> &N = itc->second;
			auto margin = [&N](const std::array<double,3> &d) {
				double m = 1e300;
				for (size_t i = 0; i < N.size(); i++) {
					double t = N[i][0]*d[0] + N[i][1]*d[1] + N[i][2]*d[2];
					if (t < m) m = t;
				}
				return m;
			};
			// A direction the two steps above failed to normalise is still sitting in the map
			// as a near-zero vector, and its margin is ~0 -- which BEATS every real candidate
			// whenever the cone is empty, so the garbage would survive and create_from_field
			// would then reject the key (|d| < 1e-12) and fall back to avg_offset, the legacy
			// lattice guess with none of the guarantees here. Treat a non-unit direction as
			// having no margin at all, so a real direction always wins.
			double dl0 = std::sqrt(kv.second[0]*kv.second[0] + kv.second[1]*kv.second[1]
			                     + kv.second[2]*kv.second[2]);
			double m0 = (dl0 > 0.5) ? margin(kv.second) : -1e300;
			if (m0 > CONE_MARGIN) continue;               // already comfortably inside

			// The deepest direction the cone allows, from the shared helper above.
			std::array<double,3> best;
			double bm = cone_max_min(N, &best);
			if (bm <= m0) { best = kv.second; bm = m0; }
			const double enum_bm = bm;   // the enumeration's own depth, before the projection
			                             // trades some of it away -- this is what the selfcheck
			                             // has to be compared against, not the final bm
			// Aim the smoothed
			// direction at HALF that depth rather than at the bare acceptance threshold: a
			// projection that stops on the cone boundary is grazing, and the first perturbation
			// downstream pushes it back out (measured: targeting CONE_MARGIN gave 38 neighbour
			// folds against 19 for a deeper aim, at 9 inverted against 8). Half the achievable
			// depth is self-scaling -- generous in a wide cone, modest in a tight one -- and
			// needs no second constant. If the projection cannot reach it, the enumeration's
			// direction stands, which is also what happens when the cone is empty.
			if (dl0 > 0.5 && bm > 0.0) {
				const double target = 0.5 * bm;
				std::array<double,3> pd = {kv.second[0]/dl0, kv.second[1]/dl0, kv.second[2]/dl0};
				for (int it2 = 0; it2 < 32; it2++) {
					size_t wi = 0; double wm = 1e300;
					for (size_t i = 0; i < N.size(); i++) {
						double t = N[i][0]*pd[0] + N[i][1]*pd[1] + N[i][2]*pd[2];
						if (t < wm) { wm = t; wi = i; }
					}
					if (wm >= target) break;
					for (int k = 0; k < 3; k++) pd[k] += (target - wm) * N[wi][k];
					double pl = std::sqrt(pd[0]*pd[0] + pd[1]*pd[1] + pd[2]*pd[2]);
					if (pl < 1e-12) break;
					for (int k = 0; k < 3; k++) pd[k] /= pl;
				}
				if (margin(pd) >= target) { best = pd; bm = margin(pd); }
			}
			if (cone_selfcheck) {
				double bf = -1e300;
				for (int sp = 0; sp < 2000; sp++) {
					double z = 1.0 - 2.0*(sp + 0.5)/2000.0;
					double rr = std::sqrt(std::max(0.0, 1.0 - z*z));
					double phi = M_PI*(3.0 - std::sqrt(5.0))*sp;
					std::array<double,3> v = {rr*std::cos(phi), rr*std::sin(phi), z};
					double m = margin(v);
					if (m > bf) bf = m;
				}
				if (bf - enum_bm > worst_gap) worst_gap = bf - enum_bm;
			}
			if (bm > m0) { kv.second = best; n_reaimed++; }
			if (bm <= 0.0) n_empty_cone++;
		}
		printf("    Extrusion cone repair: %d (node, side) directions re-aimed into the feasible "
		       "cone, %d with an empty cone (no direction valid)\n", n_reaimed, n_empty_cone);
		if (cone_selfcheck)
			printf("    Cone repair selfcheck: worst margin the sphere search beat the "
			       "enumeration by: %.3e\n", worst_gap);

		std::unordered_map<int64_t, double> nominal = ext_len;   // size-based thickness, pre-cap

		// Second pass, with the smoothed directions known: cap the length so the offset also
		// keeps the HOST element valid. This is Marechal's "adjust the vector size until the
		// mother hex stays valid" (IMR 2016, fig. 8), in closed form. Writing the offset in the
		// basis of the host's three edges at that corner, d = c1 e1 + c2 e2 + c3 e3, the corner
		// Jacobians of the host scale as (1 - L*sum(ci)) and (1 - L*ci), so L < ETA/max(...)
		// keeps every one of them positive. Without this a smoothed interface normal on a steep
		// flank points outside the host's own edge cone and slides the corner out of its cell:
		// the layer comes out perfect and the host hexes turn inside out (hawaii, mauna_loa).
		for (auto &q : interface_quads) {
			octant_t *eh[2] = { (octant_t *) sc_array_index(&mesh->elements, q.elemA_id),
			                    (octant_t *) sc_array_index(&mesh->elements, q.elemB_id) };
			const int face[2] = { q.faceA, q.faceB };
			for (int n = 0; n < 4; n++) {
				int nid = q.quad_nodes[n];
				const int64_t key[2] = { pillow_key(nid, q.matA), pillow_key(nid, q.matB) };
				for (int sd = 0; sd < 2; sd++) {
					auto itd = ext_dir.find(key[sd]);
					auto itl = ext_len.find(key[sd]);
					if (itd == ext_dir.end() || itl == ext_len.end()) continue;
					int li = FaceNodesMap[face[sd]][n];          // corner in elem->nodes order
					int k  = (li + 4) % 8;                       // same corner in CORNER_NB order
					const double *p0 = &coords[3*eh[sd]->nodes[li].id];
					double E[3][3];                              // columns = the three edges
					for (int j = 0; j < 3; j++) {
						int nb = eh[sd]->nodes[mgeom::H5_ORD[mgeom::CORNER_NB[k][j]]].id;
						for (int r = 0; r < 3; r++) E[r][j] = coords[3*nb+r] - p0[r];
					}
					double det = E[0][0]*(E[1][1]*E[2][2]-E[1][2]*E[2][1])
					           - E[0][1]*(E[1][0]*E[2][2]-E[1][2]*E[2][0])
					           + E[0][2]*(E[1][0]*E[2][1]-E[1][1]*E[2][0]);
					double en[3];                                // relative conditioning test:
					for (int j = 0; j < 3; j++)                  // an ill-conditioned corner
						en[j] = std::sqrt(E[0][j]*E[0][j] + E[1][j]*E[1][j] + E[2][j]*E[2][j]);
					if (std::fabs(det) < 1e-9 * en[0]*en[1]*en[2]) continue;   // gives no usable bound
					const std::array<double,3> &d = itd->second;
					double c[3];                                 // Cramer: E c = d
					for (int j = 0; j < 3; j++) {
						double M[3][3];
						for (int r = 0; r < 3; r++) for (int cc = 0; cc < 3; cc++)
							M[r][cc] = (cc == j) ? d[r] : E[r][cc];
						c[j] = (M[0][0]*(M[1][1]*M[2][2]-M[1][2]*M[2][1])
						      - M[0][1]*(M[1][0]*M[2][2]-M[1][2]*M[2][0])
						      + M[0][2]*(M[1][0]*M[2][1]-M[1][1]*M[2][0])) / det;
					}
					double worst = c[0] + c[1] + c[2];
					for (int j = 0; j < 3; j++) if (c[j] > worst) worst = c[j];
					if (worst <= 1e-12) continue;                // host only grows: no bound
					// Never collapse the layer to a sliver to save a host element: below ~0.4 of
					// the nominal thickness the pillow sheet reads as a crack in the rendered
					// mesh (the ALPHA_MIN=0.05 episode). An inverted fat element is repairable
					// downstream and visible; a collapsed valid one is neither.
					double lim = std::max(PILLOW_ETA / worst, LEN_FLOOR * nominal[key[sd]]);
					if (lim < itl->second) itl->second = lim;
				}
			}
		}
	}

	// Buffer node from the extrusion field. The integer lattice target is the same one the old
	// placement used -- topology bookkeeping and every downstream consumer of node->x/y/z depend
	// on it -- only the physical position changes. Returns -1 (caller falls back) when this
	// (node, side) has no usable field entry.
	auto create_from_field = [&](int orig_nid, int64_t key, uint8_t mask,
	                             int sgn, int nx, int ny, int nz) -> int {
		int dx = ((mask & 3)  == 3)  ? 0 : (mask & 1)  ? 1 : (mask & 2)  ? -1 : 0;
		int dy = ((mask & 12) == 12) ? 0 : (mask & 4)  ? 1 : (mask & 8)  ? -1 : 0;
		int dz = ((mask & 48) == 48) ? 0 : (mask & 16) ? 1 : (mask & 32) ? -1 : 0;
		if (dx == 0 && dy == 0 && dz == 0) { dx = sgn*nx; dy = sgn*ny; dz = sgn*nz; }

		auto itd = ext_dir.find(key);
		auto itl = ext_len.find(key);
		if (itd == ext_dir.end() || itl == ext_len.end() || itl->second <= 0.0) return -1;
		std::array<double,3> d = itd->second;
		double dl = std::sqrt(d[0]*d[0] + d[1]*d[1] + d[2]*d[2]);
		if (dl < 1e-12) return -1;

		octant_node_t *on = (octant_node_t *) sc_array_index(&mesh->nodes, orig_nid);
		int tx = std::clamp(on->x + dx, min_gx, max_gx);
		int ty = std::clamp(on->y + dy, min_gy, max_gy);
		int tz = std::clamp(on->z + dz, min_gz, max_gz);

		double L = itl->second;
		double ox = L*d[0]/dl, oy = L*d[1]/dl, oz = L*d[2]/dl;
		double px = coords[3*orig_nid+0], py = coords[3*orig_nid+1], pz = coords[3*orig_nid+2];
		return push_buffer_node(orig_nid, tx, ty, tz, px+ox, py+oy, pz+oz, ox, oy, oz);
	};

	// Averaged offset with a magnitude FLOOR: use the mean interior DIRECTION
	// (congruent slab, no twist) but keep the mean interior DISTANCE as the
	// length, so cancelling directions at ridges/valleys do not collapse the
	// layer to zero thickness. Fallback to the non-opposing axis normal if the
	// averaged direction vanishes or has opposing quad normals.
	auto avg_offset = [&](const std::array<double,5>& a, uint8_t mask, int sgn, int nx, int ny, int nz) -> std::array<double,3> {
		double c = a[3] > 0 ? a[3] : 1.0;
		double dx = a[0]/c, dy = a[1]/c, dz = a[2]/c;
		double mag = a[4]/c;                          // mean interior distance (~cell)
		if (mag < 1.0) mag = 50.0;

		// Zero out components along axes with true opposition (+ and - on same axis)
		if ((mask & 3) == 3) dx = 0.0;
		if ((mask & 12) == 12) dy = 0.0;
		if ((mask & 48) == 48) dz = 0.0;

		double dl = std::sqrt(dx*dx + dy*dy + dz*dz);
		if (dl > 1e-6) { return { dx/dl*mag, dy/dl*mag, dz/dl*mag }; }

		// Fallback to non-opposed axis normal
		double fx = ((mask & 3) == 3) ? 0.0 : (double)sgn * nx;
		double fy = ((mask & 12) == 12) ? 0.0 : (double)sgn * ny;
		double fz = ((mask & 48) == 48) ? 0.0 : (double)sgn * nz;
		double fl = std::sqrt(fx*fx + fy*fy + fz*fz);
		if (fl > 1e-6) { return { fx/fl*mag, fy/fl*mag, fz/fl*mag }; }

		// If all components were opposed or zero, extrude along the quad's specific normal
		double qx = (double)sgn * nx, qy = (double)sgn * ny, qz = (double)sgn * nz;
		double ql = std::sqrt(qx*qx + qy*qy + qz*qz);
		if (ql > 1e-6) { return { qx/ql * (0.3*mag), qy/ql * (0.3*mag), qz/ql * (0.3*mag) }; }

		return { 0.0, 0.0, (double)sgn * 0.3 * mag };
	};

	// Pass B: create one shared buffer node per interface node per material side it borders
	// (keeps the mesh conforming: no hanging nodes). Try the face/edge/vertex placement first;
	// fall back to the averaged-direction placement for the cases it can't handle (true
	// same-axis opposition, or no real neighbour to interpolate toward).
	int n_v2 = 0, n_fallback = 0;
	for (auto &q : interface_quads) {
		int nx = q.normal[0];
		int ny = q.normal[1];
		int nz = q.normal[2];

		for (int n = 0; n < 4; n++) {
			int orig_nid = q.quad_nodes[n];
			int64_t keyA = pillow_key(orig_nid, q.matA);
			if (pillow_map.find(keyA) == pillow_map.end()) {
				int idA = legacy_placement ? try_create_node_v2(orig_nid, side_mask[keyA])
				                           : create_from_field(orig_nid, keyA, side_mask[keyA], -1, nx, ny, nz);
				if (idA >= 0) { n_v2++; }
				else {
					n_fallback++;
					std::array<double,3> oA = avg_offset(acc[keyA], side_mask[keyA], -1, nx, ny, nz);
					int off_x = ((side_mask[keyA] & 3) == 3) ? 0 : -nx;
					int off_y = ((side_mask[keyA] & 12) == 12) ? 0 : -ny;
					int off_z = ((side_mask[keyA] & 48) == 48) ? 0 : -nz;
					idA = get_or_create_node(orig_nid, off_x, off_y, off_z, oA[0], oA[1], oA[2]);
				}
				pillow_map[keyA] = idA;
			}
			int64_t keyB = pillow_key(orig_nid, q.matB);
			if (pillow_map.find(keyB) == pillow_map.end()) {
				int idB = legacy_placement ? try_create_node_v2(orig_nid, side_mask[keyB])
				                           : create_from_field(orig_nid, keyB, side_mask[keyB], +1, nx, ny, nz);
				if (idB >= 0) { n_v2++; }
				else {
					n_fallback++;
					std::array<double,3> oB = avg_offset(acc[keyB], side_mask[keyB], +1, nx, ny, nz);
					int off_x = ((side_mask[keyB] & 3) == 3) ? 0 : +nx;
					int off_y = ((side_mask[keyB] & 12) == 12) ? 0 : +ny;
					int off_z = ((side_mask[keyB] & 48) == 48) ? 0 : +nz;
					idB = get_or_create_node(orig_nid, off_x, off_y, off_z, oB[0], oB[1], oB[2]);
				}
				pillow_map[keyB] = idB;
			}
		}
	}
	printf("    Buffer node placement: %d %s, %d fell back to averaged-direction\n",
	       n_v2, legacy_placement ? "face/edge/vertex (legacy)" : "extrusion field",
	       n_fallback);

	// 5. Global Element Remapping for Conformity (NO HANGING NODES)
	// Remap every original element to its (node, n_mat) buffer, whatever n_mat is.
	size_t n_orig_elems = mesh->elements.elem_count;
	for (size_t iel = 0; iel < n_orig_elems; iel++) {
		octant_t *elem = (octant_t *) sc_array_index(&mesh->elements, iel);
		for (int ino = 0; ino < 8; ino++) {
			auto it = pillow_map.find(pillow_key(elem->nodes[ino].id, elem->n_mat));
			if (it != pillow_map.end()) {
				elem->nodes[ino].id = it->second;
				octant_node_t *pn = (octant_node_t *) sc_array_index(&mesh->nodes, it->second);
				elem->nodes[ino].x = pn->x;
				elem->nodes[ino].y = pn->y;
				elem->nodes[ino].z = pn->z;
			}
		}
	}

	// 6. Create Double Pillow Layer Hexahedra using Exact Topological Face-Nodes Mapping
	int n_created_pillow = 0;
	for (size_t iq = 0; iq < interface_quads.size(); iq++) {
		interface_quad_t q = interface_quads[iq];

		int matA_nodes[4], matB_nodes[4];
		for (int n = 0; n < 4; n++) {
			int orig_nid = q.quad_nodes[n];
			matA_nodes[n] = pillow_map[pillow_key(orig_nid, q.matA)];
			matB_nodes[n] = pillow_map[pillow_key(orig_nid, q.matB)];
		}

		// Create Pillow Element A (matA side)
		octant_t *elemA = (octant_t *) sc_array_index(&mesh->elements, q.elemA_id);
		int8_t lvlA = elemA->level;
		int32_t exA = elemA->x, eyA = elemA->y, ezA = elemA->z;

		octant_t *pelemA = (octant_t *) sc_array_push(&mesh->elements);
		memset(pelemA, 0, sizeof(octant_t));
		pelemA->id = (int64_t)(mesh->elements.elem_count - 1);
		pelemA->n_mat = q.matA;
		pelemA->level = lvlA;
		pelemA->x = exA; pelemA->y = eyA; pelemA->z = ezA;

		// Exact 3D Topological Assignment using FaceNodesMap and FaceNodesMap_inv
		for (int k = 0; k < 4; k++) {
			int lo = FaceNodesMap[q.faceA][k];     // face nodes (original interface)
			int li = FaceNodesMap_inv[q.faceA][k]; // opposite face nodes (matA buffer)

			pelemA->nodes[lo].id = q.quad_nodes[k];
			octant_node_t *n_lo = (octant_node_t *) sc_array_index(&mesh->nodes, q.quad_nodes[k]);
			pelemA->nodes[lo].x = n_lo->x; pelemA->nodes[lo].y = n_lo->y; pelemA->nodes[lo].z = n_lo->z;
			pelemA->nodes[lo].fixed = n_lo->fixed; pelemA->nodes[lo].color = n_lo->color;

			pelemA->nodes[li].id = matA_nodes[k];
			octant_node_t *n_li = (octant_node_t *) sc_array_index(&mesh->nodes, matA_nodes[k]);
			pelemA->nodes[li].x = n_li->x; pelemA->nodes[li].y = n_li->y; pelemA->nodes[li].z = n_li->z;
			pelemA->nodes[li].fixed = n_li->fixed; pelemA->nodes[li].color = n_li->color;
		}
		n_created_pillow++;

		// Create Pillow Element B (matB side)
		octant_t *elemB = (octant_t *) sc_array_index(&mesh->elements, q.elemB_id);
		int8_t lvlB = elemB->level;
		int32_t exB = elemB->x, eyB = elemB->y, ezB = elemB->z;

		octant_t *pelemB = (octant_t *) sc_array_push(&mesh->elements);
		memset(pelemB, 0, sizeof(octant_t));
		pelemB->id = (int64_t)(mesh->elements.elem_count - 1);
		pelemB->n_mat = q.matB;
		pelemB->level = lvlB;
		pelemB->x = exB; pelemB->y = eyB; pelemB->z = ezB;

		// Exact 3D Topological Assignment using FaceNodesMap and FaceNodesMap_inv
		for (int k = 0; k < 4; k++) {
			int lo = FaceNodesMap[q.faceB][k];     // face nodes (original interface)
			int li = FaceNodesMap_inv[q.faceB][k]; // opposite face nodes (matB buffer)

			pelemB->nodes[lo].id = q.quad_nodes[k];
			octant_node_t *n_lo = (octant_node_t *) sc_array_index(&mesh->nodes, q.quad_nodes[k]);
			pelemB->nodes[lo].x = n_lo->x; pelemB->nodes[lo].y = n_lo->y; pelemB->nodes[lo].z = n_lo->z;
			pelemB->nodes[lo].fixed = n_lo->fixed; pelemB->nodes[lo].color = n_lo->color;

			pelemB->nodes[li].id = matB_nodes[k];
			octant_node_t *n_li = (octant_node_t *) sc_array_index(&mesh->nodes, matB_nodes[k]);
			pelemB->nodes[li].x = n_li->x; pelemB->nodes[li].y = n_li->y; pelemB->nodes[li].z = n_li->z;
			pelemB->nodes[li].fixed = n_li->fixed; pelemB->nodes[li].color = n_li->color;
		}
		n_created_pillow++;
	}


	// 7. Validity-driven pull-back of buffer nodes (conforming by construction: only positions
	// move, never connectivity). A buffer node whose incident hexes came out inverted is slid
	// back along its own placement segment toward its interface node until every incident hex is
	// valid, or ALPHA_MIN is reached. This is the local equivalent of "collapsing" the bad
	// element -- the layer thins where it does not fit instead of the mesh losing an element.
	if (getenv("NO_PULLBACK") == NULL) {
		double vol_sum = 0.0;
		double X[8], Y[8], Z[8];
		for (size_t ie = 0; ie < mesh->elements.elem_count; ie++) {
			load_elem_xyz(mesh, coords, (int)ie, X, Y, Z);
			vol_sum += mgeom::hex_signed_volume(X, Y, Z);
		}
		int ref = mgeom::reference_sign(vol_sum);

		// buffer node -> incident elements
		std::unordered_map<int, std::vector<int>> incid;
		for (size_t ie = 0; ie < mesh->elements.elem_count; ie++) {
			octant_t *e = (octant_t *) sc_array_index(&mesh->elements, ie);
			for (int k = 0; k < 8; k++) {
				int id = e->nodes[k].id;
				if (id >= (int)initial_node_count) incid[id].push_back((int)ie);
			}
		}

		// Worst corner Jacobian (sign-corrected) over a node's incident elements, and how many
		// of them are inverted. This is the whole objective: raise the worst corner above zero.
		auto score = [&](const std::vector<int> &els, int &bad) -> double {
			double worst = 1e300;
			bad = 0;
			for (int ie : els) {
				load_elem_xyz(mesh, coords, ie, X, Y, Z);
				double v = mgeom::hex_signed_volume(X, Y, Z);
				double sj = mgeom::hex_min_corner_sj(X, Y, Z) * ref;
				if (mgeom::is_inverted(v, mgeom::hex_min_corner_sj(X, Y, Z), ref)) bad++;
				if (sj < worst) worst = sj;
			}
			return worst;
		};

		// Feasible set of one buffer node: on the segment side of its interface node, thickness
		// within [ALPHA_MIN, 1] of the placed one, direction within ANGLE_MAX of the extrusion
		// direction. Any candidate is projected back into it, so every accepted position keeps
		// the layer visible and shaped -- the constraint the unconstrained node smoothing in
		// PillowingInterface.cpp lacked when it crumpled the coastline.
		auto project = [&](const double p0[3], const double u0[3], double t_nom,
		                   double a_min, double a_cos, double c[3]) {
			double w[3] = { c[0]-p0[0], c[1]-p0[1], c[2]-p0[2] };
			double l = std::sqrt(w[0]*w[0] + w[1]*w[1] + w[2]*w[2]);
			if (l < 1e-12) { for (int k = 0; k < 3; k++) c[k] = p0[k] + a_min*t_nom*u0[k]; return; }
			double wh[3] = { w[0]/l, w[1]/l, w[2]/l };
			double cs = wh[0]*u0[0] + wh[1]*u0[1] + wh[2]*u0[2];
			if (cs < a_cos) {                           // rotate back onto the cone boundary
				double t[3] = { wh[0]-cs*u0[0], wh[1]-cs*u0[1], wh[2]-cs*u0[2] };
				double tl = std::sqrt(t[0]*t[0] + t[1]*t[1] + t[2]*t[2]);
				double sn = std::sqrt(std::max(0.0, 1.0 - a_cos*a_cos));
				for (int k = 0; k < 3; k++)
					wh[k] = a_cos*u0[k] + (tl > 1e-12 ? sn*t[k]/tl : 0.0);
			}
			if (l > t_nom)        l = t_nom;            // never thicker than placed
			if (l < a_min*t_nom)  l = a_min*t_nom;      // never below this level's floor
			for (int k = 0; k < 3; k++) c[k] = p0[k] + l*wh[k];
		};

		int n_moved = 0, n_fixed = 0, n_relaxed = 0;
		double min_alpha = 1.0;

		for (int round = 0; round < 3; round++) {
			int moved_this_round = 0;
			for (auto &kv : incid) {
				int bid = kv.first;
				auto ob = buffer_orig.find(bid);
				if (ob == buffer_orig.end()) continue;
				int bad0; double s0 = score(kv.second, bad0);
				if (bad0 == 0) continue;

				int oid = ob->second;
				const double p0[3] = { coords[3*oid+0], coords[3*oid+1], coords[3*oid+2] };
				double w0[3] = { coords[3*bid+0]-p0[0], coords[3*bid+1]-p0[1], coords[3*bid+2]-p0[2] };
				double t_cur = std::sqrt(w0[0]*w0[0] + w0[1]*w0[1] + w0[2]*w0[2]);
				if (t_cur < 1e-12) continue;
				const double u0[3] = { w0[0]/t_cur, w0[1]/t_cur, w0[2]/t_cur };
				// thickness as placed: the current one divided by whatever scale is already applied
				double t_nom = t_cur / std::max(buffer_scale[bid], 1e-9);

				// Pattern search: 6 axis moves plus grow/shrink along the extrusion direction,
				// step halved whenever no candidate improves. Purely local -- this node's
				// position against this node's incident elements, nothing else is read or written.
				const double dirs[8][3] = {
					{ 1,0,0},{-1,0,0},{0, 1,0},{0,-1,0},{0,0, 1},{0,0,-1},
					{ u0[0], u0[1], u0[2]}, {-u0[0],-u0[1],-u0[2]}
				};
				double best[3] = { coords[3*bid+0], coords[3*bid+1], coords[3*bid+2] };
				double best_s = s0; int best_bad = bad0;
				bool improved_any = false;

				// Escalation ladder: search inside the visual constraints first (thickness >=
				// ALPHA_MIN of the placed one, direction within 45 deg). Only for the nodes that
				// would otherwise stay inverted, retry with a lower floor and a wider cone --
				// and there accept a move ONLY if it strictly removes an inverted element, never
				// for a Jacobian gain. A dozen thinner cells at isolated spots is a different
				// thing from thinning the whole coastline, which is what ALPHA_MIN = 0.05 did.
				const double LADDER[2][2] = { { ALPHA_MIN, 0.7071 },    // 0.40, 45 deg
				                              { 0.25,      0.5    } };  // 0.25, 60 deg
				for (int lvl = 0; lvl < 2 && best_bad > 0; lvl++) {
					const double a_min = LADDER[lvl][0], a_cos = LADDER[lvl][1];
					double step = 0.5 * t_cur;
					for (int it = 0; it < 8 && best_bad > 0; it++) {
						bool improved = false;
						for (int d = 0; d < 8; d++) {
							double c[3] = { best[0] + step*dirs[d][0],
							                best[1] + step*dirs[d][1],
							                best[2] + step*dirs[d][2] };
							project(p0, u0, t_nom, a_min, a_cos, c);
							for (int k = 0; k < 3; k++) coords[3*bid+k] = c[k];
							int bad; double sc = score(kv.second, bad);
							bool take = (lvl == 0) ? (bad < best_bad || (bad == best_bad && sc > best_s))
							                       : (bad < best_bad);
							if (take) {
								best_bad = bad; best_s = sc;
								for (int k = 0; k < 3; k++) best[k] = c[k];
								improved = true; improved_any = true;
								if (lvl > 0) n_relaxed++;
							}
						}
						if (!improved) step *= 0.5;
						if (step < 1e-3 * t_nom) break;
					}
				}

				for (int k = 0; k < 3; k++) coords[3*bid+k] = best[k];
				if (improved_any) {
					double tb = std::sqrt((best[0]-p0[0])*(best[0]-p0[0]) +
					                      (best[1]-p0[1])*(best[1]-p0[1]) +
					                      (best[2]-p0[2])*(best[2]-p0[2]));
					buffer_scale[bid] = tb / t_nom;
					if (buffer_scale[bid] < min_alpha) min_alpha = buffer_scale[bid];
					moved_this_round++;
					if (best_bad == 0) n_fixed++;
				}
			}
			n_moved += moved_this_round;
			if (moved_this_round == 0) break;
		}

		// Coupled leftovers: a buffer node shared by two hexes that spoil each other cannot be
		// fixed by moving it alone, whatever the stencil. Hand those clusters to the patch-local
		// elliptic solve, which moves the whole cluster at once (see src/pillow_untangle.cpp).
		if (getenv("NO_ELLIPTIC") == NULL) {
			std::unordered_map<int, BufferAnchor> anchors;
			for (auto &kv : buffer_orig) {
				int bid = kv.first, oid = kv.second;
				BufferAnchor a;
				a.orig_id = oid;
				double w[3];
				for (int k = 0; k < 3; k++) {
					a.orig_xyz[k] = coords[3*oid+k];
					w[k] = coords[3*bid+k] - a.orig_xyz[k];
				}
				double l = std::sqrt(w[0]*w[0] + w[1]*w[1] + w[2]*w[2]);
				if (l < 1e-12) continue;
				for (int k = 0; k < 3; k++) a.u[k] = w[k]/l;
				a.t_nom = l / std::max(buffer_scale[bid], 1e-9);
				anchors[bid] = a;
			}

			// Precompute interface node anchors for tangential sliding on the GTS surface
			double dom_min_x = 1e300, dom_max_x = -1e300;
			double dom_min_y = 1e300, dom_max_y = -1e300;
			for (size_t i = 0; i < coords.size()/3; i++) {
				double x = coords[3*i+0], y = coords[3*i+1];
				if (x < dom_min_x) dom_min_x = x;
				if (x > dom_max_x) dom_max_x = x;
				if (y < dom_min_y) dom_min_y = y;
				if (y > dom_max_y) dom_max_y = y;
			}
			const double wall_x = std::min(1.0, 1e-3 * (dom_max_x - dom_min_x));
			const double wall_y = std::min(1.0, 1e-3 * (dom_max_y - dom_min_y));

			std::unordered_map<int, double> iface_hmin;
			for (const auto &q : interface_quads) {
				for (int n = 0; n < 4; n++) {
					int n1 = q.quad_nodes[n];
					int n2 = q.quad_nodes[(n+1)%4];
					double dx = coords[3*n1+0] - coords[3*n2+0];
					double dy = coords[3*n1+1] - coords[3*n2+1];
					double dz = coords[3*n1+2] - coords[3*n2+2];
					double elen = std::sqrt(dx*dx + dy*dy + dz*dz);
					auto it1 = iface_hmin.find(n1);
					if (it1 == iface_hmin.end() || elen < it1->second) iface_hmin[n1] = elen;
					auto it2 = iface_hmin.find(n2);
					if (it2 == iface_hmin.end() || elen < it2->second) iface_hmin[n2] = elen;
				}
			}

			std::unordered_map<int, InterfaceAnchor> iface_anchors;
			for (auto &kv : iface_hmin) {
				int oid = kv.first;
				InterfaceAnchor ia;
				ia.orig_xyz[0] = coords[3*oid+0];
				ia.orig_xyz[1] = coords[3*oid+1];
				ia.orig_xyz[2] = coords[3*oid+2];
				ia.max_disp = 0.20 * kv.second;
				ia.lock_mask = 0;
				ia.buf_count = 0;
				if (std::fabs(ia.orig_xyz[0] - dom_min_x) < wall_x || std::fabs(ia.orig_xyz[0] - dom_max_x) < wall_x)
					ia.lock_mask |= mgeom::LOCK_X;
				if (std::fabs(ia.orig_xyz[1] - dom_min_y) < wall_y || std::fabs(ia.orig_xyz[1] - dom_max_y) < wall_y)
					ia.lock_mask |= mgeom::LOCK_Y;

				int best_k = -1;
				double best_dz = 1e300;
				for (size_t k = 0; k < mesh->gdata_vec.size(); k++) {
					double z_eval = ia.orig_xyz[2];
					if (eval_gts_height_near(mesh, (int)k, ia.orig_xyz[0], ia.orig_xyz[1], ia.orig_xyz[2], z_eval)) {
						double dz = std::fabs(ia.orig_xyz[2] - z_eval);
						if (dz < best_dz) {
							best_dz = dz;
							best_k = (int)k;
						}
					}
				}
				// ONLY associate with GTS surface if the node actually lies on it (dz < 2.0 m)
				if (best_dz < 2.0) {
					ia.gts_id = best_k;
				} else {
					ia.gts_id = -1; // flat interface / not on GTS surface: keeps z = orig_xyz[2]
				}
				iface_anchors[oid] = ia;
			}
			for (auto &kv : buffer_orig) {
				int bid = kv.first, oid = kv.second;
				auto it = iface_anchors.find(oid);
				if (it != iface_anchors.end() && it->second.buf_count < 4) {
					it->second.buf_ids[it->second.buf_count++] = bid;
				}
			}

			int npatch = 0;
			int nfix = EllipticPatchUntangle(mesh, coords, ref, anchors, iface_anchors, 0.25, 0.5, &npatch);
			printf("    Elliptic patch untangle: %d patches, %d elements recovered\n", npatch, nfix);
		}

		int still_bad = 0;
		for (size_t ie = 0; ie < mesh->elements.elem_count; ie++) {
			load_elem_xyz(mesh, coords, (int)ie, X, Y, Z);
			if (mgeom::is_inverted(mgeom::hex_signed_volume(X, Y, Z),
			                       mgeom::hex_min_corner_sj(X, Y, Z), ref)) still_bad++;
		}
		printf("    Pillow repair: %d buffer nodes moved (%d cleared all incidents, %d needed the "
		       "relaxed level, min alpha %.2f), %d elements still inverted\n",
		       n_moved, n_fixed, n_relaxed, min_alpha, still_bad);
	}

	// Update local and total mesh counts
	mesh->local_n_nodes = (int32_t)mesh->nodes.elem_count;
	mesh->total_n_nodes = (int64_t)mesh->nodes.elem_count;
	mesh->local_n_elements = (int32_t)mesh->elements.elem_count;
	mesh->total_n_elements = (int64_t)mesh->elements.elem_count;

	// Register all interface nodes in nodes_b_mat for downstream output
	std::unordered_set<int> bset(nodes_b_mat.begin(), nodes_b_mat.end());
	for (auto &q : interface_quads) {
		for (int n = 0; n < 4; n++) {
			int nid = q.quad_nodes[n];
			if (bset.insert(nid).second) {
				nodes_b_mat.push_back(nid);
			}
		}
	}

	// part_nodes was sized for the pre-pillow node count in hexa_mesh.cpp; every consumer
	// (the VTK writer, Apply_material, the optimizers) indexes it by current node id, so it
	// must grow with the node array -- otherwise they read past the allocation, which only
	// segfaults when the following page happens to be unmapped. Same re-allocation the other
	// node-creating passes already do (PillowingInterface.cpp, hexa_pml.cpp).
	free(mesh->part_nodes);
	mesh->part_nodes = (int32_t *) malloc(mesh->local_n_nodes * sizeof(int32_t));
	for (int ino = 0; ino < mesh->local_n_nodes; ino++)
		mesh->part_nodes[ino] = mesh->mpi_rank;

	printf("    Pillowing summary: %zu new nodes created, %d new pillow elements added\n",
	       mesh->nodes.elem_count - initial_node_count,
	       n_created_pillow);
	printf("    Total mesh count: %d elements, %d nodes (%d interface nodes registered)\n",
	       mesh->local_n_elements, mesh->local_n_nodes, (int)nodes_b_mat.size());
}
#endif // archived double-layer extrusion
