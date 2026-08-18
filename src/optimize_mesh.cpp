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

	double bbox_x1 = mesh->gdata.bbox ? mesh->gdata.bbox->x1 : 0.0;
	double bbox_x2 = mesh->gdata.bbox ? mesh->gdata.bbox->x2 : 0.0;
	double bbox_y1 = mesh->gdata.bbox ? mesh->gdata.bbox->y1 : 0.0;
	double bbox_y2 = mesh->gdata.bbox ? mesh->gdata.bbox->y2 : 0.0;
	double bbox_z1 = mesh->gdata.bbox ? mesh->gdata.bbox->z1 : 0.0;
	double bbox_z2 = mesh->gdata.bbox ? mesh->gdata.bbox->z2 : 0.0;

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
				if (std::fabs(xval - bbox_x1) < 1.0 || std::fabs(xval - bbox_x2) < 1.0)
					m |= mgeom::LOCK_X;
			}
			if (cy) {
				double yval = coords[3*nid[0]+1];
				if (std::fabs(yval - bbox_y1) < 1.0 || std::fabs(yval - bbox_y2) < 1.0)
					m |= mgeom::LOCK_Y;
			}
			if (cz) {
				double zval = coords[3*nid[0]+2];
				if (std::fabs(zval - bbox_z1) < 1.0 || std::fabs(zval - bbox_z2) < 1.0)
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

// Local state of a node's incident elements: how many are inverted, and the
// worst reference-signed scaled Jacobian.
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

// Perform surface-constrained tangential relaxation for a node on a GTS surface
static bool relax_gts_surface_node(hexa_tree_t *mesh, std::vector<double> &coords,
                                   const std::vector<std::vector<int>> &adj,
                                   const std::vector<int> &inc, int node, int ref,
                                   const NodeConstraint &cons) {
	if (inc.empty() || adj[node].empty() || cons.gts_surface_id < 0) return false;
	// Z always comes from the surface height below, never masked by apply_lock like
	// x/y are -- a Z-locked node cannot be legitimately relaxed onto the surface.
	if (cons.lock_mask & mgeom::LOCK_Z) return false;

	double cx = 0, cy = 0;
	int cnt = 0;
	for (int nb : adj[node]) {
		if (nb >= 0 && nb < (int)coords.size()/3) {
			cx += coords[3*nb+0];
			cy += coords[3*nb+1];
			cnt++;
		}
	}
	if (cnt == 0) return false;
	cx /= cnt; cy /= cnt;

	double px = coords[3*node+0], py = coords[3*node+1], pz = coords[3*node+2];
	double dx = cx - px, dy = cy - py, dz = 0.0;

	mgeom::apply_lock(cons.lock_mask, dx, dy, dz);
	if (std::fabs(dx) + std::fabs(dy) < 1e-12) return false;

	NodeState s0 = node_state(mesh, coords, inc, ref);
	double alpha = 1.0;
	while (alpha > 0.01) {
		double cand_x = px + alpha * dx;
		double cand_y = py + alpha * dy;
		double cand_z = pz;

		if (eval_gts_height(mesh, cons.gts_surface_id, cand_x, cand_y, cand_z)) {
			coords[3*node+0] = cand_x;
			coords[3*node+1] = cand_y;
			coords[3*node+2] = cand_z;

			NodeState s = node_state(mesh, coords, inc, ref);
			if (better_state(s, s0)) return true;
		}

		coords[3*node+0] = px;
		coords[3*node+1] = py;
		coords[3*node+2] = pz;
		alpha *= 0.5;
	}

	return false;
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


// Coordinate line-search that improves the node's incident state monotonically.
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
                        const std::vector<NodeConstraint> &cons,
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

		bool any_moved = false;
		for (int node : nodes) {
			if (inc[node].empty()) continue;
			if (cons[node].gts_surface_id >= 0) {
				if (relax_gts_surface_node(mesh, coords, adj, inc[node], node, ref, cons[node]))
					any_moved = true;
			} else {
				double step0 = UNTANGLE_STEP0 * shortest_incident_edge(mesh, coords, inc[node], node);
				if (relax_node(mesh, coords, inc[node], node, ref, cons[node].lock_mask, step0, 0.0))
					any_moved = true;
				if (relax_toward_centroid(mesh, coords, adj, inc[node], node, ref, cons[node].lock_mask))
					any_moved = true;
			}
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


// Try to move node to target (masked), accepting via backtracking only while
// EVERY incident element stays valid (minSJ*ref > 0).
static bool guarded_move_to(hexa_tree_t *mesh, std::vector<double> &coords,
                            const std::vector<int> &inc, int node, int ref,
                            uint8_t mask, double tx, double ty, double tz) {
	double px=coords[3*node+0], py=coords[3*node+1], pz=coords[3*node+2];
	double dx=tx-px, dy=ty-py, dz=tz-pz;
	mgeom::apply_lock(mask, dx, dy, dz);
	if (std::fabs(dx)+std::fabs(dy)+std::fabs(dz) < 1e-12) return false;

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
	hexa_mesh_write_quality_h5(mesh, "mesh_before_opt", coords, q_before);
	printf("    Exported pre-optimization quality: mesh_before_opt_*.h5 / .xmf\n");

	auto inc = build_incidence(mesh, nn);
	auto adj = build_adjacency(mesh, nn);

	// ---- Stage 1: GTS Surface & Boundary Regularization Phase ----
	int surf_moves = 0;
	printf("    Stage 1: Regularizing GTS surfaces & boundary walls (15 sweeps)...\n");
	for (int sweep = 0; sweep < 15; sweep++) {
		int moves = 0;
		for (int node = 0; node < nn; node++) {
			if (cons[node].gts_surface_id >= 0) {
				if (relax_gts_surface_node(mesh, coords, adj, inc[node], node, ref, cons[node]))
					moves++;
			} else if (cons[node].lock_mask != (mgeom::LOCK_X|mgeom::LOCK_Y|mgeom::LOCK_Z)) {
				if (relax_toward_centroid(mesh, coords, adj, inc[node], node, ref, cons[node].lock_mask))
					moves++;
			}
		}
		surf_moves += moves;
		if (moves == 0) break;
	}
	printf("    Stage 1 finished: %d node moves across GTS surfaces/walls.\n", surf_moves);

	// ---- Stage 2: Volume Untangling Phase ----
	int remaining = untangle_inversions(mesh, coords, cons, wall_lock, ref);

	// ---- Stage 3: Volume & Boundary Equalization Phase ----
	printf("    Stage 3: Volume & boundary equalization (10 sweeps)...\n");
	for (int sweep = 0; sweep < 10; sweep++) {
		int moves = 0;
		for (int node = 0; node < nn; node++) {
			if (cons[node].lock_mask == (mgeom::LOCK_X|mgeom::LOCK_Y|mgeom::LOCK_Z)) continue;
			if (cons[node].gts_surface_id >= 0) {
				if (relax_gts_surface_node(mesh, coords, adj, inc[node], node, ref, cons[node]))
					moves++;
			} else {
				if (relax_toward_centroid(mesh, coords, adj, inc[node], node, ref, cons[node].lock_mask))
					moves++;
			}
		}
		if (moves == 0) break;
	}

	MeshAnalysis a1 = analyze_mesh(mesh, coords);
	printf("    Final:   %d inverted, h_min %.6e (dt gain %.3fx)\n",
	       a1.n_inverted, a1.h_min, (h_min_0 > 0 ? a1.h_min / h_min_0 : 1.0));

	// 4. Evaluate and export AFTER optimization quality
	std::vector<hex_quality_t> q_after;
	analyze_full_mesh_quality(mesh, coords, q_after);
	print_quality_summary("AFTER  Opt", q_after);
	hexa_mesh_write_quality_h5(mesh, "mesh_after_opt", coords, q_after);
	printf("    Exported post-optimization quality: mesh_after_opt_*.h5 / .xmf\n");

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
