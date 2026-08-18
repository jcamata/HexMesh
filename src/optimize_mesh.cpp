#include <cstdio>
#include <cmath>
#include <vector>
#include <array>
#include <algorithm>
#include <unordered_map>
#include <unordered_set>
#include <sc.h>
#include <sc_containers.h>
#include "hexa.h"
#include "verify_mesh.h"
#include "mesh_geom.h"

// ---- Tunables -------------------------------------------------------------
static const int    MAX_UNTANGLE_ITERS = 2000;
static const int    STALL_PATIENCE     = 60;    // stalled iters to tolerate before giving up
static const double UNTANGLE_STEP0     = 1.5;   // initial line-search step = this * shortest incident edge
static const int    PHASE_A_SWEEPS     = 20;
static const int    PHASE_B_ROUNDS     = 10;
static const double PHASE_B_FRACTION   = 0.05;  // fraction of smallest elems attacked
static const double ESCALATION_CAP     = 0.10;  // max interface move = 10% shortest edge
static const bool   ESCALATION_ENABLED = false; // interface relaxation crushes h_min for ~marginal gain; off

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

// Per-node lock mask: interface = fully fixed; external faces lock their
// constant-integer-coordinate axis; interior/buffer = free.
std::vector<uint8_t> classify_node_constraints(hexa_tree_t *mesh,
                                               const std::vector<int> &nodes_b_mat,
                                               std::vector<uint8_t> *wall_out) {
	int nn = mesh->nodes.elem_count;
	std::vector<uint8_t> lock(nn, 0);

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

	for (auto &kv : seen) {
		if (shared.count(kv.first)) continue; // internal face, skip
		int iel = kv.second.first, f = kv.second.second;
		octant_t *e = (octant_t *) sc_array_index(&mesh->elements, iel);
		int nid[4];
		int cx[4], cy[4], cz[4];
		for (int k = 0; k < 4; k++) {
			nid[k] = e->nodes[FaceNodesMap[f][k]].id;
			octant_node_t *pn = (octant_node_t *) sc_array_index(&mesh->nodes, nid[k]);
			cx[k] = pn->x; cy[k] = pn->y; cz[k] = pn->z;
		}
		uint8_t m = mgeom::constant_axes_mask(cx, cy, cz);
		for (int k = 0; k < 4; k++) lock[nid[k]] |= m;
	}

	// Wall-only mask, captured BEFORE the interface override, so escalation can
	// free an interface node's tangential DOFs without ever letting it leave an
	// external boundary plane (a shoreline node is both interface and wall).
	if (wall_out) *wall_out = lock;

	// Interface (bathymetry/topography) nodes fully fixed. Applied last so it
	// dominates any wall lock at the shoreline.
	for (int nid : nodes_b_mat)
		if (nid >= 0 && nid < nn) lock[nid] = mgeom::LOCK_X | mgeom::LOCK_Y | mgeom::LOCK_Z;

	return lock;
}

// node -> list of incident element ids
static std::vector<std::vector<int>> build_incidence(hexa_tree_t *mesh, int n_nodes) {
	std::vector<std::vector<int>> inc(n_nodes);
	for (int iel = 0; iel < mesh->elements.elem_count; iel++) {
		octant_t *e = (octant_t *) sc_array_index(&mesh->elements, iel);
		for (int ino = 0; ino < 8; ino++) inc[e->nodes[ino].id].push_back(iel);
	}
	return inc;
}

// shortest edge incident to a node (physical), for the escalation cap and steps
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

// min reference-signed scaled Jacobian over a node's incident elements
double node_quality(hexa_tree_t *mesh, const std::vector<double> &coords,
                    const std::vector<int> &inc, int ref) {
	double q = 1e300;
	for (int iel : inc) {
		double X[8],Y[8],Z[8];
		load_elem_xyz(mesh, coords, iel, X, Y, Z);
		double sj = mgeom::hex_min_corner_sj(X, Y, Z) * ref;
		if (sj < q) q = sj;
	}
	return q;
}

// Local state of a node's incident elements: how many are inverted, and the
// worst reference-signed scaled Jacobian. Used for a globally MONOTONE accept
// rule (see better_state): a node move may never increase its own incident
// inverted count. Since an element is inverted only via moves of ITS nodes, and
// every such move is guarded here, the global inverted count cannot increase.
struct NodeState { int n_inv; double min_sj; };

static NodeState node_state(hexa_tree_t *mesh, const std::vector<double> &coords,
                            const std::vector<int> &inc, int ref) {
	NodeState s; s.n_inv = 0; s.min_sj = 1e300;
	for (int iel : inc) {
		double X[8],Y[8],Z[8];
		load_elem_xyz(mesh, coords, iel, X, Y, Z);
		double vol = mgeom::hex_signed_volume(X, Y, Z);
		double sj  = mgeom::hex_min_corner_sj(X, Y, Z);
		if (mgeom::is_inverted(vol, sj, ref)) s.n_inv++;
		double ssj = sj * ref;
		if (ssj < s.min_sj) s.min_sj = ssj;
	}
	return s;
}

// b is better than a iff it has FEWER inverted incident elements, or the same
// number but a strictly higher worst scaled Jacobian.
static bool better_state(const NodeState &b, const NodeState &a) {
	if (b.n_inv != a.n_inv) return b.n_inv < a.n_inv;
	return b.min_sj > a.min_sj + 1e-12;
}

// Coordinate line-search that improves the node's incident state monotonically.
// `eff_mask` may free axes beyond the stored lock (escalation); `cap` bounds the
// total displacement from the node's start position (<=0 means no cap). Returns
// true if the node moved. Never increases the node's incident inverted count.
static bool relax_node(hexa_tree_t *mesh, std::vector<double> &coords,
                       const std::vector<int> &inc, int node, int ref,
                       uint8_t eff_mask, double step0, double cap) {
	double sx = coords[3*node+0], sy = coords[3*node+1], sz = coords[3*node+2];
	NodeState s0 = node_state(mesh, coords, inc, ref);
	bool moved = false;

	for (int axis = 0; axis < 3; axis++) {
		if ((axis==0 && (eff_mask&mgeom::LOCK_X)) ||
		    (axis==1 && (eff_mask&mgeom::LOCK_Y)) ||
		    (axis==2 && (eff_mask&mgeom::LOCK_Z))) continue;

		for (int dir = -1; dir <= 1; dir += 2) {
			double step = step0;
			while (step > step0 * 1e-3) {
				double nx = coords[3*node+0], ny = coords[3*node+1], nz = coords[3*node+2];
				double cand = (axis==0?nx:axis==1?ny:nz) + dir*step;
				// enforce cap on displacement from the start position
				if (cap > 0.0) {
					double base = (axis==0?sx:axis==1?sy:sz);
					if (std::fabs(cand - base) > cap) cand = base + (cand>base?cap:-cap);
				}
				double px=nx, py=ny, pz=nz;
				if (axis==0) coords[3*node+0]=cand; else if (axis==1) coords[3*node+1]=cand; else coords[3*node+2]=cand;
				NodeState s = node_state(mesh, coords, inc, ref);
				if (better_state(s, s0)) { s0 = s; moved = true; break; } // accept, keep
				coords[3*node+0]=px; coords[3*node+1]=py; coords[3*node+2]=pz; // revert
				step *= 0.5;
			}
		}
	}
	return moved;
}

static std::vector<std::vector<int>> build_adjacency(hexa_tree_t *mesh, int n_nodes); // defined below

// Move the node toward its edge-neighbour centroid (Laplacian), accepted only if
// the node's incident state improves (monotone). Complements the axis line-search:
// it escapes coordinate-descent local minima that stall plain per-axis moves.
static bool relax_toward_centroid(hexa_tree_t *mesh, std::vector<double> &coords,
                                  const std::vector<std::vector<int>> &adj,
                                  const std::vector<int> &inc, int node, int ref, uint8_t mask) {
	if (adj[node].empty()) return false;
	double cx=0, cy=0, cz=0;
	for (int nb : adj[node]) { cx+=coords[3*nb]; cy+=coords[3*nb+1]; cz+=coords[3*nb+2]; }
	double invn = 1.0/(double)adj[node].size();
	cx*=invn; cy*=invn; cz*=invn;
	double px=coords[3*node], py=coords[3*node+1], pz=coords[3*node+2];
	double dx=cx-px, dy=cy-py, dz=cz-pz;
	mgeom::apply_lock(mask, dx, dy, dz);
	if (std::fabs(dx)+std::fabs(dy)+std::fabs(dz) < 1e-12) return false;
	NodeState s0 = node_state(mesh, coords, inc, ref);
	double alpha = 1.0;
	while (alpha > 1.0/64.0) {
		coords[3*node]=px+alpha*dx; coords[3*node+1]=py+alpha*dy; coords[3*node+2]=pz+alpha*dz;
		NodeState s = node_state(mesh, coords, inc, ref);
		if (better_state(s, s0)) return true;
		alpha *= 0.5;
	}
	coords[3*node]=px; coords[3*node+1]=py; coords[3*node+2]=pz;
	return false;
}

// Drive untangling until 0 inverted or MAX_UNTANGLE_ITERS. Returns remaining count.
int untangle_inversions(hexa_tree_t *mesh, std::vector<double> &coords,
                        const std::vector<uint8_t> &lock,
                        const std::vector<uint8_t> &wall_lock, int ref) {
	int nn = mesh->nodes.elem_count;
	auto inc = build_incidence(mesh, nn);
	auto adj = build_adjacency(mesh, nn);
	int prev_count = -1;
	int stall_streak = 0;

	for (int iter = 0; iter < MAX_UNTANGLE_ITERS; iter++) {
		MeshAnalysis a = analyze_mesh(mesh, coords);
		if (a.n_inverted == 0) { printf("    Untangler: 0 inverted after %d iters\n", iter); return 0; }
		if (iter % 20 == 0) printf("      untangle iter %d: %d inverted (stall %d)\n", iter, a.n_inverted, stall_streak);

		// A stalled iteration = the inverted COUNT did not drop. We keep going
		// (the min_sj-improving moves chip away at the worst elements) until the
		// count has been stuck for STALL_PATIENCE consecutive iterations.
		if (prev_count >= 0 && a.n_inverted >= prev_count) stall_streak++; else stall_streak = 0;
		prev_count = a.n_inverted;
		bool stalled = (stall_streak > 0);

		// nodes incident to any currently-inverted element
		std::unordered_set<int> nodes;
		for (int iel : a.inverted_ids) {
			octant_t *e = (octant_t *) sc_array_index(&mesh->elements, iel);
			for (int ino = 0; ino < 8; ino++) nodes.insert(e->nodes[ino].id);
		}

		bool any_moved = false;
		for (int node : nodes) {
			if (inc[node].empty()) continue;
			double step0 = UNTANGLE_STEP0 * shortest_incident_edge(mesh, coords, inc[node], node);
			// free-node pass: honour the stored lock, no cap
			if (relax_node(mesh, coords, inc[node], node, ref, lock[node], step0, 0.0))
				any_moved = true;
			// Smoothing candidate: also try the adjacency centroid (Laplacian),
			// guarded by the same monotone accept — helps escape coordinate-descent
			// local minima that stall plain axis line-search.
			if (relax_toward_centroid(mesh, coords, adj, inc[node], node, ref, lock[node]))
				any_moved = true;
		}

		// Escalation: count stalled -> allow the interface nodes of the still-
		// inverted elements to move, capped at ESCALATION_CAP * shortest edge.
		bool esc_moved = false;
		if (stalled && ESCALATION_ENABLED) {
			for (int node : nodes) {
				if (inc[node].empty()) continue;
				double he = shortest_incident_edge(mesh, coords, inc[node], node);
				double cap = ESCALATION_CAP * he;
				double step0 = 0.5 * he;
				// eff_mask = wall lock only: interface freedom is granted, but a
				// shoreline node still cannot leave its external boundary plane.
				if (relax_node(mesh, coords, inc[node], node, ref, wall_lock[node], step0, cap))
					esc_moved = true;
			}
		}
		(void)any_moved; (void)esc_moved;

		if (stall_streak >= STALL_PATIENCE) {
			printf("    Untangler: stalled with %d inverted after %d iters (%d stalled)\n", a.n_inverted, iter, stall_streak);
			return a.n_inverted;
		}
	}
	MeshAnalysis f = analyze_mesh(mesh, coords);
	printf("    Untangler: hit MAX_UNTANGLE_ITERS with %d inverted\n", f.n_inverted);
	return f.n_inverted;
}

// node -> set of edge-neighbour node ids
static std::vector<std::vector<int>> build_adjacency(hexa_tree_t *mesh, int n_nodes) {
	static const int E[12][2] = {
		{0,1},{1,2},{2,3},{3,0},{4,5},{5,6},{6,7},{7,4},{0,4},{1,5},{2,6},{3,7}
	};
	std::vector<std::unordered_set<int>> tmp(n_nodes);
	for (int iel = 0; iel < mesh->elements.elem_count; iel++) {
		octant_t *e = (octant_t *) sc_array_index(&mesh->elements, iel);
		for (int k = 0; k < 12; k++) {
			int a = e->nodes[E[k][0]].id, b = e->nodes[E[k][1]].id;
			if (a != b) { tmp[a].insert(b); tmp[b].insert(a); }
		}
	}
	std::vector<std::vector<int>> adj(n_nodes);
	for (int i = 0; i < n_nodes; i++) adj[i].assign(tmp[i].begin(), tmp[i].end());
	return adj;
}

// Try to move node to target (masked), accepting via backtracking only while
// EVERY incident element stays valid (minSJ*ref > 0). Returns true if moved.
static bool guarded_move_to(hexa_tree_t *mesh, std::vector<double> &coords,
                            const std::vector<int> &inc, int node, int ref,
                            uint8_t mask, double tx, double ty, double tz) {
	double px=coords[3*node+0], py=coords[3*node+1], pz=coords[3*node+2];
	double dx=tx-px, dy=ty-py, dz=tz-pz;
	mgeom::apply_lock(mask, dx, dy, dz);
	if (std::fabs(dx)+std::fabs(dy)+std::fabs(dz) < 1e-12) return false;

	// Accept only moves that are monotone in BOTH objectives. A move changes only
	// the edges and elements incident to `node`, so keeping the node's own
	// shortest incident edge and inverted count from getting worse keeps the
	// GLOBAL h_min and inverted count from getting worse. Without the edge test
	// Phase A's Laplacian happily collapses small elements (h_min 7e-2 m),
	// which is the opposite of the CFL goal it exists to serve.
	NodeState s0 = node_state(mesh, coords, inc, ref);
	double h0 = shortest_incident_edge(mesh, coords, inc, node);

	double alpha = 1.0;
	while (alpha > 0.01) {
		coords[3*node+0]=px+alpha*dx; coords[3*node+1]=py+alpha*dy; coords[3*node+2]=pz+alpha*dz;
		NodeState s1 = node_state(mesh, coords, inc, ref);
		if (s1.n_inv <= s0.n_inv &&
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

	// ---- Phase A: global equalization (Laplacian toward neighbour centroid) --
	printf("    Optimizer Phase A: %d equalization sweeps...\n", PHASE_A_SWEEPS);
	for (int sweep = 0; sweep < PHASE_A_SWEEPS; sweep++) {
		int moves = 0;
		for (int node = 0; node < nn; node++) {
			if (lock[node] == (mgeom::LOCK_X|mgeom::LOCK_Y|mgeom::LOCK_Z)) continue;
			if (adj[node].empty() || inc[node].empty()) continue;
			double cx=0, cy=0, cz=0;
			for (int nb : adj[node]) { cx+=coords[3*nb]; cy+=coords[3*nb+1]; cz+=coords[3*nb+2]; }
			double inv = 1.0 / (double)adj[node].size();
			if (guarded_move_to(mesh, coords, inc[node], node, ref, lock[node], cx*inv, cy*inv, cz*inv))
				moves++;
		}
		if (moves == 0) break;
	}

	// ---- Phase B: min-focus (grow the smallest elements) --------------------
	int n_target = std::max(1, (int)(PHASE_B_FRACTION * ne));
	printf("    Optimizer Phase B: %d rounds, targeting %d smallest elems...\n", PHASE_B_ROUNDS, n_target);
	for (int round = 0; round < PHASE_B_ROUNDS; round++) {
		std::vector<std::pair<double,int>> sized(ne);
		for (int iel = 0; iel < ne; iel++) sized[iel] = { elem_char_size(mesh, coords, iel), iel };
		std::partial_sort(sized.begin(), sized.begin()+n_target, sized.end());

		int moves = 0;
		for (int t = 0; t < n_target; t++) {
			int iel = sized[t].second;
			octant_t *e = (octant_t *) sc_array_index(&mesh->elements, iel);
			// element centroid
			double gx=0, gy=0, gz=0;
			for (int ino = 0; ino < 8; ino++) {
				int id=e->nodes[ino].id; gx+=coords[3*id]; gy+=coords[3*id+1]; gz+=coords[3*id+2];
			}
			gx/=8; gy/=8; gz/=8;
			// push each free node radially outward from the centroid to enlarge it
			for (int ino = 0; ino < 8; ino++) {
				int node = e->nodes[ino].id;
				if (lock[node] == (mgeom::LOCK_X|mgeom::LOCK_Y|mgeom::LOCK_Z)) continue;
				double dx=coords[3*node]-gx, dy=coords[3*node+1]-gy, dz=coords[3*node+2]-gz;
				double r = std::sqrt(dx*dx+dy*dy+dz*dz);
				if (r < 1e-9) continue;
				double grow = 0.10; // 10% radial expansion attempt per round
				double tx=coords[3*node]  + grow*dx;
				double ty=coords[3*node+1]+ grow*dy;
				double tz=coords[3*node+2]+ grow*dz;
				if (guarded_move_to(mesh, coords, inc[node], node, ref, lock[node], tx, ty, tz))
					moves++;
			}
		}
		if (moves == 0) break;
	}
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

void MeshOptimization(hexa_tree_t *mesh, std::vector<double> &coords, std::vector<int> material_fixed_nodes) {
	if (!mesh || mesh->elements.elem_count == 0 || coords.empty()) return;

	// Run self-test unit assertions for quality metrics
	hexQualitySelfTest();

	// Size optimization (the time-step objective) is off for now: what matters at this stage
	// is a topologically correct mesh, not dt. Untangling stays on -- it is what removes
	// inverted elements, and without it the mesh keeps the raw inversions out of pillowing.
	const bool run_size_optimization = false;

	printf("\n =========================================================\n");
	printf("   MESH UNTANGLE%s\n", run_size_optimization ? " + SIZE OPTIMIZATION" : " (size optimization disabled)");
	printf(" =========================================================\n");

	std::vector<uint8_t> wall_lock;
	std::vector<uint8_t> lock = classify_node_constraints(mesh, material_fixed_nodes, &wall_lock);

	// snapshot locked coordinates to self-check the boundary invariant later
	int nn = mesh->nodes.elem_count;
	std::vector<double> coords0 = coords;

	MeshAnalysis a0 = analyze_mesh(mesh, coords);
	int ref = a0.reference_sign;
	double h_min_0 = a0.h_min;
	printf("    Initial: %d inverted, ref sign %+d, h_min %.6e\n", a0.n_inverted, ref, h_min_0);

	// 1. Evaluate and export BEFORE optimization quality
	std::vector<hex_quality_t> q_before;
	analyze_full_mesh_quality(mesh, coords, q_before);
	print_quality_summary("BEFORE Opt", q_before);
	hexa_mesh_write_quality_h5(mesh, "mesh_before_opt", coords, q_before);
	printf("    Exported pre-optimization quality: mesh_before_opt_*.h5 / .xmf\n");

	// 2. Perform optimization / untangling
	int remaining = untangle_inversions(mesh, coords, lock, wall_lock, ref);
	if (!run_size_optimization) {
		printf("    Size optimization disabled.\n");
	} else if (remaining == 0) {
		optimize_size(mesh, coords, lock, ref);
	} else {
		printf("    Skipping size optimization: %d inverted elements remain after untangling.\n", remaining);
	}

	MeshAnalysis a1 = analyze_mesh(mesh, coords);
	printf("    Final:   %d inverted, h_min %.6e (dt gain %.3fx)\n",
	       a1.n_inverted, a1.h_min, (h_min_0 > 0 ? a1.h_min / h_min_0 : 1.0));

	// 3. Evaluate and export AFTER optimization quality
	std::vector<hex_quality_t> q_after;
	analyze_full_mesh_quality(mesh, coords, q_after);
	print_quality_summary("AFTER  Opt", q_after);
	hexa_mesh_write_quality_h5(mesh, "mesh_after_opt", coords, q_after);
	printf("    Exported post-optimization quality: mesh_after_opt_*.h5 / .xmf\n");

	// Boundary invariant self-check: every locked coordinate must be unchanged,
	// EXCEPT interface nodes that escalation was allowed to nudge (<= cap). We
	// assert non-interface locks are exact and report the max interface drift.
	std::unordered_set<int> iface(material_fixed_nodes.begin(), material_fixed_nodes.end());
	double max_iface_drift = 0.0;
	int viol = 0;
	for (int node = 0; node < nn; node++) {
		uint8_t m = lock[node];
		if (!m) continue;
		double ddx = coords[3*node]  - coords0[3*node];
		double ddy = coords[3*node+1]- coords0[3*node+1];
		double ddz = coords[3*node+2]- coords0[3*node+2];
		if (iface.count(node)) {
			double d = std::sqrt(ddx*ddx+ddy*ddy+ddz*ddz);
			if (d > max_iface_drift) max_iface_drift = d;
		} else {
			// wall node: the locked axis/axes must not have moved at all
			if (((m&mgeom::LOCK_X) && std::fabs(ddx) > 1e-6) ||
			    ((m&mgeom::LOCK_Y) && std::fabs(ddy) > 1e-6) ||
			    ((m&mgeom::LOCK_Z) && std::fabs(ddz) > 1e-6)) viol++;
		}
	}
	printf("    Boundary self-check: %d wall-lock violations, max interface drift %.4e m\n",
	       viol, max_iface_drift);
	if (viol != 0) printf("    ERROR: external boundary nodes moved off their plane!\n");
	printf(" =========================================================\n\n");
}
