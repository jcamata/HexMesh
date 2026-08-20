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
		if (nopp > 0)
			printf("    WARNING: %d / %zu (node, material) sides have opposing quad normals "
			       "(empty feasible cone) -- their pillow hexes cannot all be valid\n",
			       nopp, side_mask.size());
	}

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
				int idA = try_create_node_v2(orig_nid, side_mask[keyA]);
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
				int idB = try_create_node_v2(orig_nid, side_mask[keyB]);
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
	printf("    Buffer node placement: %d face/edge/vertex, %d fell back to averaged-direction\n",
	       n_v2, n_fallback);

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
