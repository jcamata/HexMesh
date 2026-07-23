#include <gts.h>
#include <glib.h>
#include <vector>
#include <iostream>
using namespace std;
#include <set>
#include <map>
#include <tuple>
#include <algorithm>
#include <unordered_map>
#include <sc.h>
#include <sc_io.h>
#include <sc_containers.h>
#include <mpi.h>
#include <cassert>

#include "hexa.h"
#include "hilbert.h"

#include <ctime>
#include <chrono>

unsigned pillow_hash_fn(const void *v, const void *u)
{
	const pillow_t *q = (const pillow_t *)v;
	uint32_t a, b, c;

	a = (uint32_t)q->x;
	b = (uint32_t)q->y;
	c = (uint32_t)q->z;
	sc_hash_mix(a, b, c);
	sc_hash_final(a, b, c);
	return (unsigned)c;
}

int pillow_equal_fn(const void *v1, const void *v2, const void *u)
{
	const pillow_t *q1 = (const pillow_t *)v1;
	const pillow_t *q2 = (const pillow_t *)v2;
	return (q1->x == q2->x && q1->y == q2->y && q1->z == q2->z);
}

int AddPoint(hexa_tree_t *mesh, sc_hash_array_t *hash_nodes, GtsPoint *p, std::vector<double> &coords, int x, int y, int z)
{
	size_t position;
	octant_node_t *r;
	octant_node_t key;

	key.x = x;
	key.y = y;
	key.z = z;

	r = (octant_node_t *)sc_hash_array_insert_unique(hash_nodes, &key, &position);

	if (r != NULL)
	{
		r->x = x;
		r->y = y;
		r->z = z;
		r->id = mesh->nodes.elem_count;
		octant_node_t *n = (octant_node_t *)sc_array_push(&mesh->nodes);
		n->id = r->id;
		n->x = x;
		n->y = y;
		n->z = z;
		n->color = -1;
		n->fixed = 0;

		const double xx = p->x;
		const double yy = p->y;
		const double zz = p->z;

		coords.push_back(xx);
		coords.push_back(yy);
		coords.push_back(zz);
		return r->id;
	}
	else
	{
		r = (octant_node_t *)sc_array_index(&hash_nodes->a, position);
		return r->id;
	}
}

GtsPoint *LinearMapHex(const double *cord_in_ref, const double *cord_in_x, const double *cord_in_y, const double *cord_in_z)
{

	double N[8];
	GtsPoint *point;
	double out[3];

	N[0] = (1 - cord_in_ref[0]) * (1 - cord_in_ref[1]) * (1 - cord_in_ref[2]) / double(8);
	N[1] = (1 + cord_in_ref[0]) * (1 - cord_in_ref[1]) * (1 - cord_in_ref[2]) / double(8);
	N[2] = (1 + cord_in_ref[0]) * (1 + cord_in_ref[1]) * (1 - cord_in_ref[2]) / double(8);
	N[3] = (1 - cord_in_ref[0]) * (1 + cord_in_ref[1]) * (1 - cord_in_ref[2]) / double(8);

	N[4] = (1 - cord_in_ref[0]) * (1 - cord_in_ref[1]) * (1 + cord_in_ref[2]) / double(8);
	N[5] = (1 + cord_in_ref[0]) * (1 - cord_in_ref[1]) * (1 + cord_in_ref[2]) / double(8);
	N[6] = (1 + cord_in_ref[0]) * (1 + cord_in_ref[1]) * (1 + cord_in_ref[2]) / double(8);
	N[7] = (1 - cord_in_ref[0]) * (1 + cord_in_ref[1]) * (1 + cord_in_ref[2]) / double(8);

	out[0] = 0;
	out[1] = 0;
	out[2] = 0;

	for (int i = 0; i < 8; i++)
	{
		out[0] = N[i] * cord_in_x[i] + out[0];
		out[1] = N[i] * cord_in_y[i] + out[1];
		out[2] = N[i] * cord_in_z[i] + out[2];
	}

	point = gts_point_new(gts_point_class(), out[0], out[1], out[2]);

	return point;
}

void RedoNodeMapping(hexa_tree_t *mesh)
{

	int factor = 12;
	// just the multiplication
	//  it allow us add int points in the mesh
	//  keeping a structured mesh
	for (int iel = 0; iel < mesh->elements.elem_count; iel++)
	{
		octant_t *elem = (octant_t *)sc_array_index(&mesh->elements, iel);
		for (int ino = 0; ino < 8; ino++)
		{
			elem->nodes[ino].x = factor * elem->nodes[ino].x;
			elem->nodes[ino].y = factor * elem->nodes[ino].y;
			elem->nodes[ino].z = factor * elem->nodes[ino].z;
		}
		elem->x = 4 * elem->x;
		elem->y = 4 * elem->y;
		elem->z = 4 * elem->z;
	}

	for (int ino = 0; ino < mesh->nodes.elem_count; ino++)
	{
		octant_node_t *node = (octant_node_t *)sc_array_index(&mesh->nodes, ino);
		node->x = factor * node->x;
		node->y = factor * node->y;
		node->z = factor * node->z;
	}
	mesh->x_start = factor * mesh->x_start;
	mesh->y_start = factor * mesh->y_start;
	mesh->x_end = factor * mesh->x_end;
	mesh->y_end = factor * mesh->y_end;

	mesh->ncellx = 4 * mesh->ncellx;
	mesh->ncelly = 4 * mesh->ncelly;
	mesh->max_z = 4 * mesh->max_z;
}

void CopyPropEl(hexa_tree_t *mesh, int id, octant_t *elem1)
{

	octant_t *elem = (octant_t *)sc_array_index(&mesh->elements, id);

	elem1->level = -1;
	elem1->tem = elem->tem;
	elem1->pad = elem->pad;
	elem1->n_mat = elem->n_mat;
	elem1->pml_id = elem->pml_id;
	elem1->father = id;
	elem1->boundary = elem->boundary;

	for (int ino = 0; ino < 8; ino++)
	{
		elem1->nodes[ino].color = elem->nodes[ino].color;
		elem1->nodes[ino].fixed = elem->nodes[ino].fixed;
		elem1->nodes[ino].x = elem->nodes[ino].x;
		elem1->nodes[ino].y = elem->nodes[ino].y;
		elem1->nodes[ino].z = elem->nodes[ino].z;
	}

	for (int iedge = 0; iedge < 12; iedge++)
	{
		elem1->edge[iedge].coord[0] = elem->edge[iedge].coord[0];
		elem1->edge[iedge].coord[1] = elem->edge[iedge].coord[1];
		elem1->edge[iedge].id = elem->edge[iedge].id;
		elem1->edge[iedge].ref = elem->edge[iedge].ref;
	}

	for (int isurf = 0; isurf < 6; isurf++)
	{
		elem1->surf[isurf].ext = elem->surf[isurf].ext;
	}

	elem1->x = elem->x;
	elem1->y = elem->y;
	elem1->z = elem->z;
}

// Topology-based conformity check, keyed by exact sorted node IDs.
// Defect classes:
//   - degenerate elements: repeated node IDs inside one element
//   - faces shared by more than 2 elements
//   - interior orphan faces: multiplicity 1 but not on a domain boundary plane
// The previous lattice-coordinate check could not see torn faces: once a node
// is remapped, the two sides no longer share lattice keys and never collide
// in the hash, so tears were reported as 0.
struct TopoFaceKey {
	int n[4];
	bool operator==(const TopoFaceKey &o) const {
		return n[0]==o.n[0] && n[1]==o.n[1] && n[2]==o.n[2] && n[3]==o.n[3];
	}
};
struct TopoFaceKeyHash {
	size_t operator()(const TopoFaceKey &k) const {
		uint64_t h = 1469598103934665603ULL;
		for (int i = 0; i < 4; i++) {
			h ^= (uint64_t)(uint32_t)k.n[i];
			h *= 1099511628211ULL;
		}
		return (size_t)h;
	}
};

static TopoFaceKey MakeTopoFaceKey(const octant_t *el, int iface)
{
	TopoFaceKey k;
	for (int i = 0; i < 4; i++) k.n[i] = el->nodes[FaceNodesMap[iface][i]].id;
	std::sort(k.n, k.n + 4);
	return k;
}

static bool FaceOnDomainBoundary(hexa_tree_t *mesh, const octant_t *el, int iface)
{
	// A face is on the domain boundary when its 4 nodes share one lattice
	// boundary plane (same planes used by SurfaceIdentification).
	bool x0=true, x1=true, y0=true, y1=true, z0=true, z1=true;
	for (int i = 0; i < 4; i++) {
		const octant_node_t &nd = el->nodes[FaceNodesMap[iface][i]];
		x0 = x0 && (nd.x == mesh->x_start);
		x1 = x1 && (nd.x == mesh->x_end);
		y0 = y0 && (nd.y == mesh->y_start);
		y1 = y1 && (nd.y == mesh->y_end);
		z0 = z0 && (nd.z == 0);
		z1 = z1 && (nd.z == 3 * mesh->max_z);
	}
	return x0 || x1 || y0 || y1 || z0 || z1;
}

static int RunTopologyCheck(hexa_tree_t *mesh, const char *filename, const char *label)
{
	struct FaceInfo { int count, iel, iface; };
	std::unordered_map<TopoFaceKey, FaceInfo, TopoFaceKeyHash> fmap;
	fmap.reserve(mesh->elements.elem_count * 4);

	FILE *cf = fopen(filename, "w");
	int n_degenerate = 0;

	for (int iel = 0; iel < mesh->elements.elem_count; iel++) {
		octant_t *elem = (octant_t *)sc_array_index(&mesh->elements, iel);

		bool repeated = false;
		for (int a = 0; a < 8 && !repeated; a++)
			for (int b = a + 1; b < 8; b++)
				if (elem->nodes[a].id == elem->nodes[b].id) { repeated = true; break; }
		if (repeated) {
			n_degenerate++;
			fprintf(cf, "degenerate element iel=%d ids=", iel);
			for (int a = 0; a < 8; a++) fprintf(cf, " %d", elem->nodes[a].id);
			fprintf(cf, "\n");
		}

		for (int iface = 0; iface < 6; iface++) {
			TopoFaceKey key = MakeTopoFaceKey(elem, iface);
			auto it = fmap.find(key);
			if (it == fmap.end()) fmap[key] = {1, iel, iface};
			else it->second.count++;
		}
	}

	int n_orphan = 0, n_multi = 0, n_boundary = 0;
	for (auto &kv : fmap) {
		const FaceInfo &fi = kv.second;
		if (fi.count == 2) continue;
		octant_t *elem = (octant_t *)sc_array_index(&mesh->elements, fi.iel);
		if (fi.count == 1) {
			if (FaceOnDomainBoundary(mesh, elem, fi.iface)) { n_boundary++; continue; }
			n_orphan++;
			fprintf(cf, "orphan interior face iel=%d iface=%d n_mat=%d ids=%d %d %d %d\n",
					fi.iel, fi.iface, elem->n_mat,
					kv.first.n[0], kv.first.n[1], kv.first.n[2], kv.first.n[3]);
		} else {
			n_multi++;
			fprintf(cf, "face with multiplicity %d iel=%d iface=%d ids=%d %d %d %d\n",
					fi.count, fi.iel, fi.iface,
					kv.first.n[0], kv.first.n[1], kv.first.n[2], kv.first.n[3]);
		}
	}

	fprintf(cf, "%s: %d orphan interior, %d multiplicity>2, %d degenerate elements, %d boundary faces\n",
			label, n_orphan, n_multi, n_degenerate, n_boundary);
	fclose(cf);
	printf("     [%s] orphan interior faces: %d | mult>2: %d | degenerate elements: %d\n",
			label, n_orphan, n_multi, n_degenerate);
	return n_orphan + n_multi + n_degenerate;
}

// Checks that every element referencing a given node ID stores the same (x,y,z).
// A mismatch means MovingNodes updated coordinates in some elements but not others.
static int RunNodeConsistencyCheck(hexa_tree_t *mesh, const char *filename)
{
	struct NodePos { int x, y, z, iel_first; };
	std::unordered_map<int, NodePos> seen;
	seen.reserve(mesh->elements.elem_count * 4);

	FILE *cf = fopen(filename, "w");
	int n_mismatch = 0;

	for (int iel = 0; iel < mesh->elements.elem_count; iel++) {
		octant_t *elem = (octant_t *)sc_array_index(&mesh->elements, iel);
		for (int ino = 0; ino < 8; ino++) {
			int nid = elem->nodes[ino].id;
			int x = elem->nodes[ino].x, y = elem->nodes[ino].y, z = elem->nodes[ino].z;
			auto it = seen.find(nid);
			if (it == seen.end()) {
				seen[nid] = {x, y, z, iel};
			} else if (it->second.x != x || it->second.y != y || it->second.z != z) {
				n_mismatch++;
				fprintf(cf,
				        "mismatch node=%d: first seen iel=%d (%d,%d,%d) vs iel=%d ino=%d (%d,%d,%d)\n",
				        nid, it->second.iel_first,
				        it->second.x, it->second.y, it->second.z,
				        iel, ino, x, y, z);
			}
		}
	}
	fprintf(cf, "Node consistency: %d mismatch(es)\n", n_mismatch);
	fclose(cf);
	return n_mismatch;
}

static bool IsNodeInBoundaryHash(sc_hash_array_t *hash_b_mat, const octant_node_t &node)
{
	size_t position;
	octant_node_t key;
	key.x = node.x;
	key.y = node.y;
	key.z = node.z;
	return sc_hash_array_lookup(hash_b_mat, &key, &position);
}

static octant_node_t *LookupNodeInHash(sc_hash_array_t *hash_nodes, int x, int y, int z)
{
	size_t position;
	octant_node_t key;
	key.x = x;
	key.y = y;
	key.z = z;

	if (!sc_hash_array_lookup(hash_nodes, &key, &position)) {
		return NULL;
	}
	return (octant_node_t *)sc_array_index(&hash_nodes->a, position);
}

static int CountPillowableFacesInOctree(hexa_tree_t *mesh, octree_t *oct, sc_hash_array_t *hash_b_mat)
{
	int n_faces = 0;

	for (int iel = 0; iel < 8; iel++) {
		octant_t *elem = (octant_t *)sc_array_index(&mesh->elements, oct->id[iel]);
		for (int isurf = 0; isurf < 6; isurf++) {
			bool surf = true;
			for (int ive = 0; ive < 4; ive++) {
				int inode = FaceNodesMap[isurf][ive];
				if (!IsNodeInBoundaryHash(hash_b_mat, elem->nodes[inode])) {
					surf = false;
					break;
				}
			}
			if (surf) {
				n_faces++;
			}
		}
	}

	return n_faces;
}

void Pillowing(hexa_tree_t *mesh, std::vector<double> &coords, std::vector<int> &nodes_b_mat)
{
	bool deb = true;
	bool clamped = true;
	int n_orig = mesh->elements.elem_count;

	// -- hash_nodes: all original nodes by integer (x,y,z) -------------------
	sc_hash_array_t *hash_nodes = (sc_hash_array_t *)sc_hash_array_new(
		sizeof(octant_node_t), node_hash_fn, node_equal_fn, &clamped);
	for (int iel = 0; iel < n_orig; iel++) {
		octant_t *elem = (octant_t *)sc_array_index(&mesh->elements, iel);
		for (int ino = 0; ino < 8; ino++) {
			size_t pos;
			octant_node_t key;
			key.x = elem->nodes[ino].x; key.y = elem->nodes[ino].y; key.z = elem->nodes[ino].z;
			octant_node_t *r = (octant_node_t *)sc_hash_array_insert_unique(hash_nodes, &key, &pos);
			if (r) { r->x = key.x; r->y = key.y; r->z = key.z; r->id = elem->nodes[ino].id; }
		}
	}

	// -- hash_b_mat: material-interface nodes by (x,y,z) ---------------------
	sc_hash_array_t *hash_b_mat = (sc_hash_array_t *)sc_hash_array_new(
		sizeof(octant_node_t), node_hash_fn, node_equal_fn, &clamped);
	for (int i = 0; i < (int)nodes_b_mat.size(); i++) {
		size_t pos;
		octant_node_t *nd = (octant_node_t *)sc_array_index(&mesh->nodes, nodes_b_mat[i]);
		octant_node_t key; key.x = nd->x; key.y = nd->y; key.z = nd->z;
		octant_node_t *r = (octant_node_t *)sc_hash_array_insert_unique(hash_b_mat, &key, &pos);
		if (r) { r->x = key.x; r->y = key.y; r->z = key.z; r->id = nd->id; }
	}

	// -- pre-pillow diagnostics -----------------------------------------------
	if (deb) {
		char fn[80];
		sprintf(fn, "premesh_conformity_%04d_%04d.txt", mesh->mpi_size, mesh->mpi_rank);
		RunTopologyCheck(mesh, fn, "pre-pillow");
		sprintf(fn, "node_consistency_%04d_%04d.txt", mesh->mpi_size, mesh->mpi_rank);
		printf("     Node consistency: %d mismatch(es)\n", RunNodeConsistencyCheck(mesh, fn));
	}

	// -- STEP 1+2: interface faces, inward displacements, pinch resolution ----
	struct FaceRef { int iel, iface; };
	struct InterfacePair { int iel_a, iface_a; };
	std::vector<InterfacePair> ifaces;

	// face_push[f] = inward integer offset (factor-12 space, 1/3 element = 4).
	// 4 instead of 6 so that two pillow nodes pushed toward each other across
	// a one-element gap land on different lattice slots instead of merging.
	static const int face_push[6][3] = {
		{ 4,  0,  0}, {-4,  0,  0},
		{ 0,  4,  0}, { 0, -4,  0},
		{ 0,  0,  4}, { 0,  0, -4}
	};

	struct NodeDisp {
		int dx, dy, dz;
		bool face_used[6];
		int ref_iel, ref_ino;
	};
	std::unordered_map<int, NodeDisp> ndisp;

	// Real-world position of a pillow node, computed WITHOUT LinearMapHex's
	// trilinear blend across ref_elem's full 8 corners. That blend is only
	// accurate when ref_elem is close to a rectangular box; near steep
	// bathymetry ref_elem is frequently warped by MovingNodes' surface
	// projection, and trilinear-interpolating a "pure z push" through a warped
	// hex leaks it into x/y too (the hex's shear couples all three parametric
	// axes together) -- this is what actually produced the "should be same
	// x,y, different z" twist. It also let cr run past +-1 (extrapolating
	// beyond ref_elem into the next element) whenever ref_ino was already the
	// extremal corner for the push's own axis.
	// Fix: move along ref_elem's own REAL local edges instead. For each active
	// axis, take the ACTUAL coordinates of ref_ino's neighboring corner along
	// that axis (nodes[8] is in x/y/z-bit order: pairs (0,1)(3,2)(4,5)(7,6) on
	// x, (0,3)(1,2)(4,7)(5,6) on y, (0,4)(1,5)(2,6)(3,7) on z) and move the same
	// FRACTION of that real edge vector that dx/dy/dz is of the lattice edge
	// (+-12). Every axis only ever touches two real corners at a time, so a
	// warped opposite face of ref_elem can never leak into this corner's push.
	static const int xnb[8] = {1,0,3,2,5,4,7,6};
	static const int ynb[8] = {3,2,1,0,7,6,5,4};
	static const int znb[8] = {4,5,6,7,0,1,2,3};
	auto pillowPos = [&coords](const NodeDisp &nd, octant_t *ref)->std::array<double,3> {
		int ino = nd.ref_ino;
		int oid = ref->nodes[ino].id;
		double px = coords[3*oid+0], py = coords[3*oid+1], pz = coords[3*oid+2];
		auto axisPush = [&](int comp, int nbIdx) {
			if (comp == 0) return;
			int nbid = ref->nodes[nbIdx].id;
			// lattice delta along this axis between ino and its neighbor (+-12)
			double latdx = ref->nodes[nbIdx].x - ref->nodes[ino].x;
			double latdy = ref->nodes[nbIdx].y - ref->nodes[ino].y;
			double latdz = ref->nodes[nbIdx].z - ref->nodes[ino].z;
			double lat = latdx != 0 ? latdx : (latdy != 0 ? latdy : latdz);
			double t = (double)comp / lat;
			px += t * (coords[3*nbid+0] - coords[3*oid+0]);
			py += t * (coords[3*nbid+1] - coords[3*oid+1]);
			pz += t * (coords[3*nbid+2] - coords[3*oid+2]);
		};
		axisPush(nd.dx, xnb[ino]);
		axisPush(nd.dy, ynb[ino]);
		axisPush(nd.dz, znb[ino]);
		return { px, py, pz };
	};

	// Pinch resolution: a boundary node with zero net displacement (mat-0
	// sliver one element thick, pushed from both sides) or whose pillow
	// target slot is occupied makes the shrink-set boundary non-manifold
	// there — pillowing is topologically impossible. Resolve by dissolving
	// the mat-0 elements whose interface faces touch the pinched node into
	// mat-1 (the discrete interface moves one cell; those cells were
	// geometrically ambiguous anyway) and re-deriving the interface.
	// Iterate, since dissolving can expose new pinches.
	const int max_pinch_iters = 10;
	int dissolved_total = 0;
	// signed volume of a hex from its 8 corner coords (standard node order)
	auto hexvol8 = [](const double *X, const double *Y, const double *Z) -> double {
		static const int T[6][4]={{0,1,2,6},{0,2,3,6},{0,3,7,6},{0,7,4,6},{0,4,5,6},{0,5,1,6}};
		double v=0.0;
		for(auto &t:T){
			double ax=X[t[1]]-X[t[0]],ay=Y[t[1]]-Y[t[0]],az=Z[t[1]]-Z[t[0]];
			double bx=X[t[2]]-X[t[0]],by=Y[t[2]]-Y[t[0]],bz=Z[t[2]]-Z[t[0]];
			double cx=X[t[3]]-X[t[0]],cy=Y[t[3]]-Y[t[0]],cz=Z[t[3]]-Z[t[0]];
			v += (ax*(by*cz-bz*cy)-ay*(bx*cz-bz*cx)+az*(bx*cy-by*cx))/6.0;
		}
		return v;
	};
	for (int pinch_iter = 0; ; pinch_iter++) {
		ifaces.clear();
		ndisp.clear();

		// STEP 1: pair ALL faces by exact sorted node IDs; interface =
		// shared by one n_mat==0 and one n_mat!=0 element. No hash_b_mat
		// filter: the pillowed face set is then the complete boundary
		// between the material regions (closure invariant), so the global
		// node remap in STEP 5 can never tear an unpillowed mat-0/mat-1
		// face. (The old version required all 4 face nodes in hash_b_mat;
		// wherever the node projection skipped an octree, faces were
		// silently dropped and the mesh was torn at the patch perimeter.)
		{
			std::unordered_map<TopoFaceKey, FaceRef, TopoFaceKeyHash> fmap;
			fmap.reserve(n_orig * 4);

			for (int iel = 0; iel < n_orig; iel++) {
				octant_t *elem = (octant_t *)sc_array_index(&mesh->elements, iel);
				for (int iface = 0; iface < 6; iface++) {
					TopoFaceKey key = MakeTopoFaceKey(elem, iface);
					auto it = fmap.find(key);
					if (it == fmap.end()) {
						fmap[key] = {iel, iface};
						continue;
					}
					octant_t *prev = (octant_t *)sc_array_index(&mesh->elements, it->second.iel);
					bool cur0 = (elem->n_mat == 0);
					bool prv0 = (prev->n_mat == 0);
					if (cur0 != prv0) {
						if (cur0)
							ifaces.push_back({iel, iface});
						else
							ifaces.push_back({it->second.iel, it->second.iface});
					}
					fmap.erase(it);
				}
			}
		}

		// STEP 2: accumulate inward displacement per boundary node
		// (one push per face orientation).
		for (auto &ip : ifaces) {
			octant_t *elem = (octant_t *)sc_array_index(&mesh->elements, ip.iel_a);
			int iface = ip.iface_a;
			for (int k = 0; k < 4; k++) {
				int ino = FaceNodesMap[iface][k];
				int nid = elem->nodes[ino].id;
				auto it = ndisp.find(nid);
				if (it == ndisp.end()) {
					NodeDisp nd = {};
					nd.dx = face_push[iface][0]; nd.dy = face_push[iface][1]; nd.dz = face_push[iface][2];
					nd.face_used[iface] = true;
					nd.ref_iel = ip.iel_a; nd.ref_ino = ino;
					ndisp[nid] = nd;
				} else if (!it->second.face_used[iface]) {
					it->second.dx += face_push[iface][0];
					it->second.dy += face_push[iface][1];
					it->second.dz += face_push[iface][2];
					it->second.face_used[iface] = true;
				}
			}
		}

		// Clamp corner node displacement components so corner nodes don't overshoot
		for (auto &kv : ndisp) {
			NodeDisp &nd = kv.second;
			if (nd.dx > 4)  nd.dx = 4;   else if (nd.dx < -4)  nd.dx = -4;
			if (nd.dy > 4)  nd.dy = 4;   else if (nd.dy < -4)  nd.dy = -4;
			if (nd.dz > 4)  nd.dz = 4;   else if (nd.dz < -4)  nd.dz = -4;
		}

		// D5: predictive pillow-degeneracy. Compute the volume of the pillow
		// element that WOULD be built on each interface face (same geometry as
		// STEP 3/4: outer face at the original node positions, inner face at the
		// face_push offset via LinearMapHex). Flag mat-0 cells whose pillow would
		// be inverted or near-zero so they are dissolved here, BEFORE STEP 3 ever
		// creates the bad pillow. Node-moving cannot repair these afterwards
		// (shared-node oscillation), so the only fix is to not create them.
		std::set<int> degen_iel;
		{
			std::vector<double> pv(ifaces.size(), 0.0);
			std::vector<char> pcomp(ifaces.size(), 0);
			std::vector<double> apv; double ssum = 0.0;
			for (size_t fi = 0; fi < ifaces.size(); fi++) {
				octant_t *ea = (octant_t *)sc_array_index(&mesh->elements, ifaces[fi].iel_a);
				int iface = ifaces[fi].iface_a;
				double X[8], Y[8], Z[8]; bool full = true;
				for (int k = 0; k < 4; k++) {
					int lo = FaceNodesMap[iface][k];
					int li = FaceNodesMap_inv[iface][k];
					int nid = ea->nodes[lo].id;
					X[lo] = coords[3*nid]; Y[lo] = coords[3*nid+1]; Z[lo] = coords[3*nid+2];
					auto it = ndisp.find(nid);
					if (it == ndisp.end()) { full = false; break; }
					NodeDisp &nd = it->second;
					octant_t *ref = (octant_t *)sc_array_index(&mesh->elements, nd.ref_iel);
					std::array<double,3> pp = pillowPos(nd, ref);
					X[li] = pp[0]; Y[li] = pp[1]; Z[li] = pp[2];
				}
				if (!full) continue;                 // pinched node: left to the pinch logic
				pv[fi] = hexvol8(X, Y, Z); pcomp[fi] = 1;
				ssum += (pv[fi] > 0 ? 1.0 : -1.0);
				apv.push_back(pv[fi] < 0 ? -pv[fi] : pv[fi]);
			}
			if (!apv.empty()) {
				double sgnp = ssum >= 0 ? 1.0 : -1.0;
				std::nth_element(apv.begin(), apv.begin()+apv.size()/2, apv.end());
				double thrp = 0.01 * apv[apv.size()/2];   // < 1% of median pillow volume => degenerate
				for (size_t fi = 0; fi < ifaces.size(); fi++)
					if (pcomp[fi] && pv[fi]*sgnp < thrp) degen_iel.insert(ifaces[fi].iel_a);
			}
		}

		// Detect pinched nodes: zero net displacement, target slot occupied
		// by an existing node, or two boundary nodes sharing one target slot.
		std::set<int> pinched;

		// Non-manifold interface edges: an edge shared by more than 2
		// interface faces (diagonal mat-0 contact) would make the 4 pillow
		// elements built on those faces share one identical side quad
		// (multiplicity-4 face). Resolution is edge-targeted: dissolve only
		// the mat-0 elements whose interface face contains the edge itself
		// (node-based dissolution proved too broad and oscillated).
		std::set<std::pair<int, int> > pinched_edges;
		{
			std::map<std::pair<int, int>, int> edge_count;
			for (auto &ip : ifaces) {
				octant_t *ea = (octant_t *)sc_array_index(&mesh->elements, ip.iel_a);
				for (int k = 0; k < 4; k++) {
					int a = ea->nodes[FaceNodesMap[ip.iface_a][k]].id;
					int b = ea->nodes[FaceNodesMap[ip.iface_a][(k + 1) % 4]].id;
					if (a > b) std::swap(a, b);
					edge_count[std::make_pair(a, b)]++;
				}
			}
			for (auto &ec : edge_count) {
				if (ec.second > 2)
					pinched_edges.insert(ec.first);
			}
		}

		std::map<std::tuple<int, int, int>, int> target_slot;
		for (auto &kv : ndisp) {
			int nid = kv.first;
			NodeDisp &nd = kv.second;
			if (nd.dx == 0 && nd.dy == 0 && nd.dz == 0) {
				pinched.insert(nid);
				continue;
			}
			octant_t *ref_elem = (octant_t *)sc_array_index(&mesh->elements, nd.ref_iel);
			int nx = ref_elem->nodes[nd.ref_ino].x + nd.dx;
			int ny = ref_elem->nodes[nd.ref_ino].y + nd.dy;
			int nz = ref_elem->nodes[nd.ref_ino].z + nd.dz;
			if (LookupNodeInHash(hash_nodes, nx, ny, nz) != NULL) {
				pinched.insert(nid);
				continue;
			}
			auto slot = std::make_tuple(nx, ny, nz);
			auto sit = target_slot.find(slot);
			if (sit != target_slot.end()) {
				pinched.insert(nid);
				pinched.insert(sit->second);
			} else {
				target_slot[slot] = nid;
			}
		}

		if (false) { // Disabled pinch element dissolution per user request
			if (pinched.empty() && pinched_edges.empty() && degen_iel.empty()) break;
			if (pinch_iter >= max_pinch_iters) {
				printf("     WARNING: %d pinched nodes / %d non-manifold edges remain after %d pinch "
						"iterations; falling back to fresh-ID pillow nodes (locally degenerate geometry)\n",
						(int)pinched.size(), (int)pinched_edges.size(), max_pinch_iters);
				break;
			}

			int ndiss = 0;
			for (auto &ip : ifaces) {
				octant_t *ea = (octant_t *)sc_array_index(&mesh->elements, ip.iel_a);
				if (ea->n_mat != 0) continue;

				bool dissolve = false;
				for (int k = 0; k < 4 && !dissolve; k++) {
					int a = ea->nodes[FaceNodesMap[ip.iface_a][k]].id;
					if (pinched.count(a)) { dissolve = true; break; }
					int b = ea->nodes[FaceNodesMap[ip.iface_a][(k + 1) % 4]].id;
					int lo = a < b ? a : b, hi = a < b ? b : a;
					if (pinched_edges.count(std::make_pair(lo, hi))) dissolve = true;
				}
				if (degen_iel.count(ip.iel_a)) dissolve = true;   // D5: would-be-degenerate pillow
				if (dissolve) {
					ea->n_mat = 1;
					ndiss++;
				}
			}
			dissolved_total += ndiss;
			printf("     Pinch resolution iter %d: %d pinched nodes, %d non-manifold edges, "
					"dissolved %d sliver elements into mat-1\n",
					pinch_iter, (int)pinched.size(), (int)pinched_edges.size(), ndiss);
			if (ndiss == 0) break; // nothing left to dissolve; fallback handles the rest
		} else {
			break; // break out after building ifaces & ndisp without dissolving elements
		}
	}
	if (dissolved_total > 0)
		printf("     Pinch resolution: dissolved %d sliver elements in total\n", dissolved_total);
	printf("     Found %d interface faces\n", (int)ifaces.size());

	// -- STEP 3: create one pillow node per unique boundary node --------------
	// Pillow nodes ALWAYS get a fresh node ID. Reusing an existing node when
	// the lattice slot was occupied (old behaviour on zero net displacement)
	// collapsed pillow elements onto their base nodes, producing degenerate
	// elements and faces shared by 4 elements. If the slot is still taken
	// (unresolved pinch after the resolution loop), only the lattice
	// bookkeeping key is nudged; the real geometry comes from coords[].
	std::unordered_map<int, int> pillow_map;
	int n_nudged = 0;

	for (auto &kv : ndisp) {
		int nid = kv.first;
		NodeDisp &nd = kv.second;
		octant_t *ref_elem = (octant_t *)sc_array_index(&mesh->elements, nd.ref_iel);

		int nx = ref_elem->nodes[nd.ref_ino].x + nd.dx;
		int ny = ref_elem->nodes[nd.ref_ino].y + nd.dy;
		int nz = ref_elem->nodes[nd.ref_ino].z + nd.dz;

		std::array<double,3> pp = pillowPos(nd, ref_elem);

		size_t pos;
		octant_node_t key; key.x = nx; key.y = ny; key.z = nz;
		octant_node_t *ra = (octant_node_t *)sc_hash_array_insert_unique(hash_nodes, &key, &pos);
		if (ra == NULL) {
			n_nudged++;
			while (ra == NULL) {
				key.z += 1;
				ra = (octant_node_t *)sc_hash_array_insert_unique(hash_nodes, &key, &pos);
			}
		}
		ra->x = key.x; ra->y = key.y; ra->z = key.z;
		ra->id = (int)(hash_nodes->a.elem_count - 1);
		ra->fixed = 0; ra->color = 0;
		coords.push_back(pp[0]); coords.push_back(pp[1]); coords.push_back(pp[2]);
		pillow_map[nid] = ra->id;
	}
	if (n_nudged > 0)
		printf("     WARNING: %d pillow nodes created on nudged lattice slots (unresolved pinches)\n", n_nudged);
	printf("     Created %d pillow nodes\n", (int)pillow_map.size());

	// Smooth initial pillow node positions along the interface topology to untwist corner shear
	// Disabled: pure topological Laplacian smoothing with no anchor to the actual
	// bathymetry surface -- warps the interface layer away from the input geometry
	// (crumpled patches at the coastline in gmsh).
	if (false) {
		std::unordered_map<int, std::vector<int>> pillow_adj;
		pillow_adj.reserve(pillow_map.size());

		for (auto &ip : ifaces) {
			octant_t *ea = (octant_t *)sc_array_index(&mesh->elements, ip.iel_a);
			int iface = ip.iface_a;
			for (int k = 0; k < 4; k++) {
				int nid1 = ea->nodes[FaceNodesMap[iface][k]].id;
				int nid2 = ea->nodes[FaceNodesMap[iface][(k + 1) % 4]].id;
				auto it1 = pillow_map.find(nid1);
				auto it2 = pillow_map.find(nid2);
				if (it1 != pillow_map.end() && it2 != pillow_map.end()) {
					int p1 = it1->second, p2 = it2->second;
					pillow_adj[p1].push_back(p2);
					pillow_adj[p2].push_back(p1);
				}
			}
		}

		for (int sm_iter = 0; sm_iter < 3; sm_iter++) {
			std::unordered_map<int, std::array<double, 3>> sm_coords;
			for (auto &kv : pillow_map) {
				int p = kv.second;
				auto ait = pillow_adj.find(p);
				if (ait == pillow_adj.end() || ait->second.empty()) continue;
				double sx = 0, sy = 0, sz = 0;
				for (int nbr : ait->second) {
					sx += coords[3*nbr+0];
					sy += coords[3*nbr+1];
					sz += coords[3*nbr+2];
				}
				size_t deg = ait->second.size();
				sm_coords[p] = { sx / deg, sy / deg, sz / deg };
			}
			for (auto &kv : sm_coords) {
				int p = kv.first;
				coords[3*p+0] = 0.5 * (coords[3*p+0] + kv.second[0]);
				coords[3*p+1] = 0.5 * (coords[3*p+1] + kv.second[1]);
				coords[3*p+2] = 0.5 * (coords[3*p+2] + kv.second[2]);
			}
		}
	}

	// -- STEP 4: collect pillow element data BEFORE any element modification --
	// outer face (FaceNodesMap[iface])     = original positions, shared with mat-1
	// inner face (FaceNodesMap_inv[iface]) = pillow positions,   shared with mat-0
	struct PillemData {
		int n_mat, pad, tem, x, y, z, father;
		int8_t level, pml_id;
		bool boundary;
		octant_node_t nodes[8];
	};
	std::vector<PillemData> pillow_elems;
	pillow_elems.reserve(ifaces.size());

	for (auto &ip : ifaces) {
		octant_t *ea = (octant_t *)sc_array_index(&mesh->elements, ip.iel_a);
		int iface = ip.iface_a;
		PillemData pd;
		pd.n_mat = 0; pd.pad = ea->pad; pd.tem = ea->tem;
		pd.x = ea->x; pd.y = ea->y; pd.z = ea->z;
		pd.level = -1; pd.pml_id = ea->pml_id;
		pd.father = ip.iel_a; pd.boundary = ea->boundary;
		memset(pd.nodes, 0, sizeof(pd.nodes));

		for (int k = 0; k < 4; k++) {
			int lo = FaceNodesMap[iface][k];
			int li = FaceNodesMap_inv[iface][k];

			// outer face: original position (touches mat-1 neighbour)
			pd.nodes[lo] = ea->nodes[lo];

			// inner face: pillow position (will touch updated mat-0 element)
			int orig_nid = ea->nodes[lo].id;
			auto pit = pillow_map.find(orig_nid);
			if (pit == pillow_map.end()) {
				printf("ERROR: pillow_map missing nid=%d\n", orig_nid);
				pd.nodes[li] = ea->nodes[lo];
			} else {
				int pnid = pit->second;
				octant_node_t *pn = (octant_node_t *)sc_array_index(&hash_nodes->a, pnid);
				pd.nodes[li].id = pnid;
				pd.nodes[li].x = pn->x; pd.nodes[li].y = pn->y; pd.nodes[li].z = pn->z;
				pd.nodes[li].fixed = 0; pd.nodes[li].color = 0;
			}
		}
		pillow_elems.push_back(pd);
	}

	// -- STEP 5: update ALL original mat-0 elements ---------------------------
	// Global pillow_map: every mat-0 element sharing boundary node N gets the
	// same N_a => no cross-octree inconsistency, no post-pass needed.
	for (int iel = 0; iel < n_orig; iel++) {
		octant_t *elem = (octant_t *)sc_array_index(&mesh->elements, iel);
		if (elem->n_mat != 0) continue;
		for (int ino = 0; ino < 8; ino++) {
			auto it = pillow_map.find(elem->nodes[ino].id);
			if (it == pillow_map.end()) continue;
			octant_node_t *pn = (octant_node_t *)sc_array_index(&hash_nodes->a, it->second);
			elem->nodes[ino].id = it->second;
			elem->nodes[ino].x = pn->x; elem->nodes[ino].y = pn->y; elem->nodes[ino].z = pn->z;
		}
	}

	// -- STEP 6: push pillow elements into mesh -------------------------------
	for (auto &pd : pillow_elems) {
		octant_t *pelem = (octant_t *)sc_array_push(&mesh->elements);
		pelem->id      = (int64_t)(mesh->elements.elem_count - 1);
		pelem->n_mat   = pd.n_mat;  pelem->pad    = pd.pad;
		pelem->tem     = pd.tem;    pelem->x      = pd.x;
		pelem->y       = pd.y;      pelem->z      = pd.z;
		pelem->level   = pd.level;  pelem->pml_id = pd.pml_id;
		pelem->father  = pd.father; pelem->boundary = pd.boundary;
		memcpy(pelem->nodes, pd.nodes, sizeof(pelem->nodes));
		for (int ie = 0; ie < 12; ie++) { pelem->edge[ie].ref = false; pelem->edge[ie].id = 0; }
		for (int is = 0; is < 6;  is++) { pelem->surf[is].ext = false; }
	}

	// -- STEP 7: validity repair on the interface band ------------------------
	// Where the bathy interface is steep, the thin pillow layer folds and
	// inverts (negative volume), dragging the remapped mat-0 element with it
	// (these render as holes). Repair by pulling each offending pillow node back
	// toward its interface-origin position. The limit of pulling every pillow
	// node onto its origin is the pre-pillow mesh, which is valid by
	// construction (projection limiter), so halving toward it is monotone and
	// never inverts a currently-valid element. A pillow node is shared by the
	// pillow element AND the remapped mat-0 element, so moving it repairs both.
	// Only pillow nodes move; deeper mat-0/mat-1 nodes are untouched.
	{
		std::unordered_map<int,int> pillow_origin;        // pillow node id -> interface origin id
		pillow_origin.reserve(pillow_map.size()*2);
		for (auto &kv : pillow_map) pillow_origin[kv.second] = kv.first;

		auto hexvol = [&](octant_t *e)->double {
			// IMPORTANT: use the SAME node order the h5 writer emits
			// (assign_elem_nodes = {4,5,6,7,0,1,2,3} in hexa_h5.cpp). e->nodes is
			// NOT in the standard hex order assumed by the tet decomposition T, and
			// the reorder is not a clean z-flip, so computing T directly over
			// e->nodes gives a geometrically WRONG volume — the validity repair
			// then "sees" no inversions while the written mesh has many. Reordering
			// here makes this volume identical to what ParaView/the solver compute.
			static const int ord[8] = {4,5,6,7,0,1,2,3};
			double X[8],Y[8],Z[8];
			for(int i=0;i<8;i++){int id=e->nodes[ord[i]].id; X[i]=coords[3*id];Y[i]=coords[3*id+1];Z[i]=coords[3*id+2];}
			static const int T[6][4]={{0,1,2,6},{0,2,3,6},{0,3,7,6},{0,7,4,6},{0,4,5,6},{0,5,1,6}};
			double v=0.0;
			for(auto &t:T){
				double ax=X[t[1]]-X[t[0]],ay=Y[t[1]]-Y[t[0]],az=Z[t[1]]-Z[t[0]];
				double bx=X[t[2]]-X[t[0]],by=Y[t[2]]-Y[t[0]],bz=Z[t[2]]-Z[t[0]];
				double cx=X[t[3]]-X[t[0]],cy=Y[t[3]]-Y[t[0]],cz=Z[t[3]]-Z[t[0]];
				v += (ax*(by*cz-bz*cy) - ay*(bx*cz-bz*cx) + az*(bx*cy-by*cx))/6.0;
			}
			return v;
		};

		int ne = mesh->elements.elem_count;
		// global orientation sign + median |volume| (one pass, for the threshold)
		std::vector<double> av; av.reserve(ne);
		double ssum = 0.0;
		for(int iel=0; iel<ne; iel++){
			octant_t *e=(octant_t*)sc_array_index(&mesh->elements,iel);
			double v=hexvol(e); ssum += (v>0?1.0:-1.0); av.push_back(v<0?-v:v);
		}
		double sgn = ssum>=0 ? 1.0 : -1.0;
		std::nth_element(av.begin(), av.begin()+av.size()/2, av.end());
		double thr = 1e-6 * av[av.size()/2];

		// interface band = elements incident to at least one pillow node
		std::vector<int> band;
		for(int iel=0; iel<ne; iel++){
			octant_t *e=(octant_t*)sc_array_index(&mesh->elements,iel);
			for(int i=0;i<8;i++) if(pillow_origin.count(e->nodes[i].id)){ band.push_back(iel); break; }
		}

		int nnodes = (int)(coords.size() / 3);

		// Record initial un-collapsed pillow node positions before any shrinking/clamping
		std::vector<double> pillow_init_pos(3 * nnodes, 0.0);
		for (auto &kv : pillow_origin) {
			int p = kv.first;
			pillow_init_pos[3*p+0] = coords[3*p+0];
			pillow_init_pos[3*p+1] = coords[3*p+1];
			pillow_init_pos[3*p+2] = coords[3*p+2];
		}

		// -- Interface clamp --------------------------------------------------
		// The mat0/mat1 interface is the sea floor: no interface node may sit
		// above the sea surface.
		{
			const double SEA_LEVEL = 0.0, SURFACE_EPS = 1.0;
			const double zceil = SEA_LEVEL - SURFACE_EPS;
			std::vector<unsigned char> has0(nnodes,0), has1(nnodes,0);
			for (int iel = 0; iel < ne; iel++) {
				octant_t *e = (octant_t*) sc_array_index(&mesh->elements, iel);
				for (int i = 0; i < 8; i++) { int n = e->nodes[i].id; if (e->n_mat == 0) has0[n]=1; else has1[n]=1; }
			}
			int nclamp = 0;
			for (int n = 0; n < nnodes; n++)
				if (has0[n] && has1[n] && coords[3*n+2] > zceil) { coords[3*n+2] = zceil; nclamp++; }
			printf("     Clamped %d interface nodes below the sea surface\n", nclamp);
		}

		// -- STEP 7: adaptive minimum thickness shrink — GUARANTEE V > 0.
		// Shrink pillow node positions adaptively toward interface origins with a
		// positive minimum thickness floor ALPHA_MIN (5% of initial displacement).
		// This ensures all elements retain positive 3D cell volume (V > 0).
		const double ALPHA_MIN = 0.05; // 5% minimum thickness preserved
		std::vector<double> pillow_alpha(nnodes, 1.0);
		std::vector<char> snapped(nnodes, 0);

		int it2 = 0, nbad2 = 0, nsnap = 0;
		const int MAXIT2 = 4000;
		for (it2 = 0; it2 < MAXIT2; it2++) {
			// Shrink per ELEMENT, not per node: all of an invalid element's own
			// pillow-face nodes are forced to the SAME alpha together. Shrinking
			// them independently (old behaviour) let sibling corners of one hex
			// drift to different depths whenever their OTHER incident element
			// became valid at different iterations, shearing that hex into a
			// wedge instead of just thinning it.
			std::unordered_map<int,double> proposal;         // node id -> smallest proposed alpha this iter
			nbad2 = 0;
			for (int b = 0; b < ne; b++) {                   // scan all elements
				octant_t *e = (octant_t*) sc_array_index(&mesh->elements, b);
				if (hexvol(e) * sgn > 0.0) continue;         // valid positive volume
				nbad2++;
				double minalpha = 1.0;
				std::vector<int> pnodes;
				for (int i = 0; i < 8; i++) {
					int id = e->nodes[i].id;
					if (pillow_origin.count(id) && !snapped[id]) {
						pnodes.push_back(id);
						if (pillow_alpha[id] < minalpha) minalpha = pillow_alpha[id];
					}
				}
				if (pnodes.empty()) continue;                // residual inversion owns no unsnapped pillow node
				double new_alpha = minalpha * 0.5;
				for (int id : pnodes) {
					auto it = proposal.find(id);
					if (it == proposal.end() || new_alpha < it->second) proposal[id] = new_alpha;
				}
			}
			if (nbad2 == 0) break;
			if (proposal.empty()) break;
			for (auto &kv : proposal) {
				int p = kv.first;
				int o = pillow_origin[p];
				pillow_alpha[p] = kv.second;
				if (pillow_alpha[p] <= ALPHA_MIN) {
					pillow_alpha[p] = ALPHA_MIN;
					snapped[p] = 1;
					nsnap++;
				}
				coords[3*p+0] = coords[3*o+0] + pillow_alpha[p] * (pillow_init_pos[3*p+0] - coords[3*o+0]);
				coords[3*p+1] = coords[3*o+1] + pillow_alpha[p] * (pillow_init_pos[3*p+1] - coords[3*o+1]);
				coords[3*p+2] = coords[3*o+2] + pillow_alpha[p] * (pillow_init_pos[3*p+2] - coords[3*o+2]);
			}
		}
		int n_inv_final = 0, n_inv_pillow = 0, n_inv_nonpillow = 0;
		for (int b = 0; b < ne; b++) {
			octant_t *e = (octant_t*) sc_array_index(&mesh->elements, b);
			if (hexvol(e) * sgn > 0.0) continue;
			n_inv_final++;
			bool hasp = false;
			for (int i = 0; i < 8; i++) if (pillow_origin.count(e->nodes[i].id)) { hasp = true; break; }
			if (hasp) n_inv_pillow++; else n_inv_nonpillow++;
		}
		printf("     Pillow validity repair: adaptive shrink (alpha_min=%.2f) adjusted %d pillow "
				"nodes (%d iters); final invalid=%d (pillow %d, non-pillow %d)\n",
				ALPHA_MIN, nsnap, it2, n_inv_final, n_inv_pillow, n_inv_nonpillow);

		// -- Corner-Jacobian cleanup ------------------------------------------
		// ParaView and the solver validate hexes by the 8 corner Jacobians (in
		// the h5 output node order), not by the tet-decomposition volume above —
		// the two disagree on twisted hexes. Two surgical passes:
		//   pass 1: elements listed inside-out (ALL 8 corners negative) are fixed
		//           by relabeling only (swap bottom/top node slots); geometry and
		//           face-set conformity untouched.
		//   pass 2: elements folded at a few corners get a local node nudge —
		//           only INTERIOR nodes move (never interface, sea-surface,
		//           lateral-wall or bottom nodes), and a move is accepted only if
		//           the worst corner Jacobian over ALL elements touching that
		//           node strictly improves (monotone: nothing valid can break).
		{
			// worst (minimum) corner Jacobian in the output node order
			auto hexworst = [&](octant_t *e)->double {
				static const int ord[8] = {4,5,6,7,0,1,2,3};
				static const int nb[8][3] = {{1,3,4},{2,0,5},{3,1,6},{0,2,7},{7,5,0},{4,6,1},{5,7,2},{6,4,3}};
				double X[8],Y[8],Z[8];
				for(int i=0;i<8;i++){int id=e->nodes[ord[i]].id; X[i]=coords[3*id];Y[i]=coords[3*id+1];Z[i]=coords[3*id+2];}
				double w = 1e300;
				for(int k=0;k<8;k++){
					int a=nb[k][0], b=nb[k][1], d=nb[k][2];
					double ax=X[a]-X[k],ay=Y[a]-Y[k],az=Z[a]-Z[k];
					double bx=X[b]-X[k],by=Y[b]-Y[k],bz=Z[b]-Z[k];
					double cx=X[d]-X[k],cy=Y[d]-Y[k],cz=Z[d]-Z[k];
					double J = ax*(by*cz-bz*cy) - ay*(bx*cz-bz*cx) + az*(bx*cy-by*cx);
					if (sgn*J < w) w = sgn*J;
				}
				return w;
			};
			auto hexbest = [&](octant_t *e)->double {
				static const int ord[8] = {4,5,6,7,0,1,2,3};
				static const int nb[8][3] = {{1,3,4},{2,0,5},{3,1,6},{0,2,7},{7,5,0},{4,6,1},{5,7,2},{6,4,3}};
				double X[8],Y[8],Z[8];
				for(int i=0;i<8;i++){int id=e->nodes[ord[i]].id; X[i]=coords[3*id];Y[i]=coords[3*id+1];Z[i]=coords[3*id+2];}
				double w = -1e300;
				for(int k=0;k<8;k++){
					int a=nb[k][0], b=nb[k][1], d=nb[k][2];
					double ax=X[a]-X[k],ay=Y[a]-Y[k],az=Z[a]-Z[k];
					double bx=X[b]-X[k],by=Y[b]-Y[k],bz=Z[b]-Z[k];
					double cx=X[d]-X[k],cy=Y[d]-Y[k],cz=Z[d]-Z[k];
					double J = ax*(by*cz-bz*cy) - ay*(bx*cz-bz*cx) + az*(bx*cy-by*cx);
					if (sgn*J > w) w = sgn*J;
				}
				return w;
			};

			// -- pass 1: mirror-relabel inside-out elements (best corner still
			// negative => ALL corners negative). Swapping the bottom/top node
			// slots flips the listing orientation; node ids/coords untouched, so
			// shared faces (node sets) and neighbors are unaffected.
			int nflip = 0;
			for (int b = 0; b < ne; b++) {
				octant_t *e = (octant_t*) sc_array_index(&mesh->elements, b);
				if (hexbest(e) >= 0.0) continue;
				for (int i = 0; i < 4; i++) {
					octant_node_t tmpn = e->nodes[i];
					e->nodes[i] = e->nodes[i+4];
					e->nodes[i+4] = tmpn;
				}
				nflip++;
				if (hexworst(e) <= 0.0)                    // paranoia: relabel must fix it
					printf("     WARNING: mirror relabel left element %d invalid\n", b);
			}

			// -- pass 2: local nudge of interior nodes only --------------------
			// Disabled: blind 26-direction pattern search that only maximizes a
			// scalar corner-Jacobian metric with no anchor to the input geometry --
			// produced crumpled/self-intersecting patches near the coastline in gmsh.
			// eligibility: never interface (touches both materials), never at or
			// above the sea clamp, never on the lateral walls or the bottom
			if (false) {
			std::vector<unsigned char> m0(nnodes,0), m1(nnodes,0);
			for (int iel = 0; iel < ne; iel++) {
				octant_t *e = (octant_t*) sc_array_index(&mesh->elements, iel);
				for (int i = 0; i < 8; i++) { int n = e->nodes[i].id; if (e->n_mat == 0) m0[n]=1; else m1[n]=1; }
			}
			// Wall planes are the DOMINANT extreme planes (median of the 6000 most
			// extreme values), NOT the bbox: a handful of nodes sit outside the
			// domain (deep pillow nodes extrapolated past the wall), and a
			// bbox-based guard would "protect" only those outliers while leaving
			// the true wall nodes free to move.
			auto wallplane = [&](int ax, bool low)->double {
				std::vector<double> v; v.reserve(nnodes);
				for (int n = 0; n < nnodes; n++) v.push_back(coords[3*n+ax]);
				int k = nnodes < 6000 ? nnodes : 6000;
				if (low) { std::nth_element(v.begin(), v.begin()+k, v.end());
					std::nth_element(v.begin(), v.begin()+k/2, v.begin()+k); return v[k/2]; }
				std::nth_element(v.begin(), v.end()-k, v.end());
				std::nth_element(v.end()-k, v.end()-k/2, v.end()); return *(v.end()-k/2);
			};
			const double wxlo = wallplane(0,true), wxhi = wallplane(0,false);
			const double wylo = wallplane(1,true), wyhi = wallplane(1,false);
			const double wzlo = wallplane(2,true);
			const double wtol = 1e-3, ZCEIL = -1.0;
			auto eligible = [&](int n)->bool {
				if (m0[n] && m1[n]) return false;                     // interface: fixed
				if (coords[3*n+2] >= ZCEIL - 1e-9) return false;      // sea-surface side: fixed
				if (coords[3*n+0] < wxlo+wtol || coords[3*n+0] > wxhi-wtol) return false; // on/beyond walls
				if (coords[3*n+1] < wylo+wtol || coords[3*n+1] > wyhi-wtol) return false;
				if (coords[3*n+2] < wzlo+wtol) return false;          // bottom: fixed
				return true;
			};

			// free nodes = eligible nodes of corner-invalid elements; star map
			std::vector<char> isfree(nnodes, 0);
			int n_badcorner = 0;
			for (int b = 0; b < ne; b++) {
				octant_t *e = (octant_t*) sc_array_index(&mesh->elements, b);
				if (hexworst(e) > 0.0) continue;
				n_badcorner++;
				for (int i = 0; i < 8; i++) if (eligible(e->nodes[i].id)) isfree[e->nodes[i].id] = 1;
			}
			std::unordered_map<int, std::vector<octant_t*>> star;
			for (int iel = 0; iel < ne; iel++) {
				octant_t *e = (octant_t*) sc_array_index(&mesh->elements, iel);
				for (int i = 0; i < 8; i++) if (isfree[e->nodes[i].id]) star[e->nodes[i].id].push_back(e);
			}
			std::vector<int> freev;
			for (int n = 0; n < nnodes; n++) if (isfree[n]) freev.push_back(n);

			auto starmin = [&](int p)->double {
				double m = 1e300;
				for (octant_t* e : star[p]) { double w = hexworst(e); if (w < m) m = w; }
				return m;
			};

			int nmoves = 0, sweep = 0, nstuck = 0;
			const int MAXSWEEP = 100;
			for (sweep = 0; sweep < MAXSWEEP; sweep++) {
				bool improved = false;
				nstuck = 0;
				for (int p : freev) {
					double base = starmin(p);
					if (base > 0.0) continue;
					nstuck++;
					// step scale from the largest element of the star
					double sc = 0.0;
					for (octant_t* e : star[p]) { double v = fabs(hexvol(e)); if (v > sc) sc = v; }
					sc = 0.25 * cbrt(sc + 1.0);
					double px = coords[3*p+0], py = coords[3*p+1], pz = coords[3*p+2];
					bool moved = false;
					for (double s = sc; s >= 1e-3*sc && !moved; s *= 0.5) {
						double best = base, bx = px, by = py, bz = pz;
						for (int dx = -1; dx <= 1; dx++)
						for (int dy = -1; dy <= 1; dy++)
						for (int dz = -1; dz <= 1; dz++) {
							if (!dx && !dy && !dz) continue;
							double dn = sqrt((double)(dx*dx+dy*dy+dz*dz));
							double tx = px + s*dx/dn, ty = py + s*dy/dn, tz = pz + s*dz/dn;
							if (tz >= ZCEIL) tz = ZCEIL - 1e-9;       // never at/above the sea clamp
							// never cross a domain wall or the bottom
							if (tx < wxlo+wtol || tx > wxhi-wtol ||
								ty < wylo+wtol || ty > wyhi-wtol || tz < wzlo+wtol) continue;
							coords[3*p+0] = tx;
							coords[3*p+1] = ty;
							coords[3*p+2] = tz;
							double v = starmin(p);
							if (v > best) { best = v; bx = coords[3*p+0]; by = coords[3*p+1]; bz = tz; }
						}
						coords[3*p+0] = bx; coords[3*p+1] = by; coords[3*p+2] = bz;
						if (best > base) { moved = true; improved = true; nmoves++; }
						else { coords[3*p+0] = px; coords[3*p+1] = py; coords[3*p+2] = pz; }
					}
				}
				if (nstuck == 0 || !improved) break;
			}
			int n_final = 0;
			for (int b = 0; b < ne; b++) {
				octant_t *e = (octant_t*) sc_array_index(&mesh->elements, b);
				if (hexworst(e) <= 0.0) n_final++;
			}
			printf("     Corner-Jacobian cleanup: %d inside-out elements relabeled; nudge pass on %d "
					"corner-invalid elements (%d interior nodes moved, %d sweeps); corner-invalid %d -> %d\n",
					nflip, n_badcorner, nmoves, sweep, n_badcorner + nflip, n_final);
			} // pass 2 disabled
			int n_badcorner_report = 0;
			for (int b = 0; b < ne; b++) {
				octant_t *e = (octant_t*) sc_array_index(&mesh->elements, b);
				if (hexworst(e) <= 0.0) n_badcorner_report++;
			}
			printf("     Corner-Jacobian cleanup: %d inside-out elements relabeled (pass 2 nudge disabled); "
					"corner-invalid elements remaining: %d\n", nflip, n_badcorner_report);
		}
	}

	// -- debug summary --------------------------------------------------------
	if (deb) {
		char ff[80];
		sprintf(ff, "pillow_%04d_%04d.txt", mesh->mpi_size, mesh->mpi_rank);
		FILE *pf = fopen(ff, "w");
		fprintf(pf, "Global face-based pillowing\n");
		fprintf(pf, "Interface faces: %d\n", (int)ifaces.size());
		fprintf(pf, "Pillow nodes:    %d\n", (int)pillow_map.size());
		fprintf(pf, "Pillow elements: %d\n", (int)pillow_elems.size());
		fclose(pf);
	}

	// -- post-pillow conformity check (authoritative) -------------------------
	if (deb) {
		char fn[80];
		sprintf(fn, "pillow_conformity_%04d_%04d.txt", mesh->mpi_size, mesh->mpi_rank);
		RunTopologyCheck(mesh, fn, "post-pillow");
	}

	// -- finalise: replace mesh->nodes with hash_nodes content ----------------
	sc_array_reset(&mesh->nodes);
	sc_hash_array_rip(hash_nodes, &mesh->nodes);
	sc_hash_array_destroy(hash_b_mat);
}

void SurfaceIdentification(hexa_tree_t *mesh, std::vector<double> &coords)
{

	// assign a color for the nodes...

	// free node
	// color= 0 && fixed = 0;

	// exterior surface fixed nodes
	// fixed = 1 && color = 1

	// exterior global nodes
	// color && fixed = -1;
	//  x- = 2 x+ = 3
	//  y- = 4 y+ = 5
	//  z- = 6 z+ = 7

	// exterior local nodes
	// color && fixed = -1;
	//  x- = -2 x+ = -3
	//  y- = -4 y+ = -5
	//  z- = -6 z+ = -7

	// exterior global edges
	// color && fixed = -1;
	//  z-y- = 0+11
	//  z-x+ = 1+11
	//  z-y+ = 2+11
	//  z-x- = 3+11

	// x-y- = 4+11
	// x+y- = 5+11
	// x+y+ = 6+11
	// x-y+ = 7+11

	// z+y- = 8+11
	// z+x+ = 9+11
	// z+y+ = 10+11
	// z+x- = 11+11

	// exterior global vertex
	// color && fixed = -1;
	// x-y-z- = 30
	// x+y-z- = 31
	// x+y+z- = 32
	// x-y+z- = 33

	// x-y-z+ = 34
	// x+y-z+ = 35
	// x+y+z+ = 36
	// x-y+z+ = 37

	// exterior local edges
	// color && fixed = -1;
	//  z-y- = -0-11
	//  z-x+ = -1-11
	//  z-y+ = -2-11
	//  z-x- = -3-11

	// x-y- = -4-11
	// x+y- = -5-11
	// x+y+ = -6-11
	// x-y+ = -7-11

	// z+y- = -8-11
	// z+x+ = -9-11
	// z+y+ = -10-11
	// z+x- = -11-11

	// exterior local vertex
	// color && fixed = -1;
	// x-y-z- = -30
	// x+y-z- = -31
	// x+y+z- = -32
	// x-y+z- = -33

	// x-y-z+ = -34
	// x+y-z+ = -35
	// x+y+z+ = -36
	// x-y+z+ = -37

	bool deb = false;
	bool clamped = true;
	// vertex hash
	sc_hash_array_t *vertex_hash = (sc_hash_array_t *)sc_hash_array_new(sizeof(octant_vertex_t), vertex_hash_id, vertex_equal_id, &clamped);
	// fazendo vertex hash & assign the color to the local nodes
	for (int iel = 0; iel < mesh->elements.elem_count; iel++)
	{
		size_t position;
		octant_t *elem = (octant_t *)sc_array_index(&mesh->elements, iel);
		for (int ino = 0; ino < 8; ino++)
		{
			// build the hash
			octant_vertex_t key;
			key.id = elem->nodes[ino].id;
			octant_vertex_t *vert = (octant_vertex_t *)sc_hash_array_insert_unique(vertex_hash, &key, &position);
			if (vert != NULL)
			{
				vert->id = elem->nodes[ino].id;
				vert->list_elem = 1;
				vert->elem[vert->list_elem - 1] = elem->id;
			}
			else
			{
				vert = (octant_vertex_t *)sc_array_index(&vertex_hash->a, position);
				vert->elem[vert->list_elem] = elem->id;
				vert->list_elem++;
			}

			// setting free the free nodes
			if (elem->nodes[ino].fixed == 0)
				elem->nodes[ino].color = 0;

			// assign the color for the local nodes...

			// surface and edges
			if (elem->nodes[ino].x == mesh->x_start)
			{
				elem->nodes[ino].fixed = -1;
				elem->nodes[ino].color = -2;
			}

			if (elem->nodes[ino].x == mesh->x_end)
			{
				elem->nodes[ino].fixed = -1;
				elem->nodes[ino].color = -3;
			}

			if (elem->nodes[ino].y == mesh->y_start)
			{
				elem->nodes[ino].fixed = -1;
				elem->nodes[ino].color = -4;
			}

			if (elem->nodes[ino].y == mesh->y_end)
			{
				elem->nodes[ino].fixed = -1;
				elem->nodes[ino].color = -5;
			}

			if (elem->nodes[ino].z == 0)
			{
				elem->nodes[ino].fixed = -1;
				elem->nodes[ino].color = -6;
			}

			if (elem->nodes[ino].z == 3 * mesh->max_z)
			{
				elem->nodes[ino].fixed = -1;
				elem->nodes[ino].color = -7;
			}

			// z-y- = -0-11
			if (elem->nodes[ino].z == 0 && elem->nodes[ino].y == mesh->y_start)
			{
				elem->nodes[ino].fixed = -1;
				elem->nodes[ino].color = -0 - 11;
			}
			// z-x+ = -1-11
			if (elem->nodes[ino].z == 0 && elem->nodes[ino].x == mesh->x_end)
			{
				elem->nodes[ino].fixed = -1;
				elem->nodes[ino].color = -1 - 11;
			}
			// z-y+ = -2-11
			if (elem->nodes[ino].z == 0 && elem->nodes[ino].y == mesh->y_end)
			{
				elem->nodes[ino].fixed = -1;
				elem->nodes[ino].color = -2 - 11;
			}
			// z-x- = -3-11
			if (elem->nodes[ino].z == 0 && elem->nodes[ino].x == mesh->x_start)
			{
				elem->nodes[ino].fixed = -1;
				elem->nodes[ino].color = -3 - 11;
			}

			// x-y- = -4-11
			if (elem->nodes[ino].x == mesh->x_start && elem->nodes[ino].y == mesh->y_start)
			{
				elem->nodes[ino].fixed = -1;
				elem->nodes[ino].color = -4 - 11;
			}
			// x+y- = -5-11
			if (elem->nodes[ino].x == mesh->x_end && elem->nodes[ino].y == mesh->y_start)
			{
				elem->nodes[ino].fixed = -1;
				elem->nodes[ino].color = -5 - 11;
			}
			// x+y+ = -6-11
			if (elem->nodes[ino].x == mesh->x_end && elem->nodes[ino].y == mesh->y_end)
			{
				elem->nodes[ino].fixed = -1;
				elem->nodes[ino].color = -6 - 11;
			}
			// x-y+ = -7-11
			if (elem->nodes[ino].x == mesh->x_start && elem->nodes[ino].y == mesh->y_end)
			{
				elem->nodes[ino].fixed = -1;
				elem->nodes[ino].color = -7 - 11;
			}

			// z+y- = -8-11
			if (elem->nodes[ino].z == 3 * mesh->max_z && elem->nodes[ino].y == mesh->y_start)
			{
				elem->nodes[ino].fixed = -1;
				elem->nodes[ino].color = -8 - 11;
			}
			// z+x+ = -9-11
			if (elem->nodes[ino].z == 3 * mesh->max_z && elem->nodes[ino].x == mesh->x_end)
			{
				elem->nodes[ino].fixed = -1;
				elem->nodes[ino].color = -9 - 11;
			}
			// z+y+ = -10-11
			if (elem->nodes[ino].z == 3 * mesh->max_z && elem->nodes[ino].y == mesh->y_end)
			{
				elem->nodes[ino].fixed = -1;
				elem->nodes[ino].color = -10 - 11;
			}
			// z+x- = -11-11
			if (elem->nodes[ino].z == 3 * mesh->max_z && elem->nodes[ino].x == mesh->x_start)
			{
				elem->nodes[ino].fixed = -1;
				elem->nodes[ino].color = -11 - 11;
			}

			// x-y-z- = -30
			if (elem->nodes[ino].x == mesh->x_start && elem->nodes[ino].y == mesh->y_start && elem->nodes[ino].z == 0)
			{
				elem->nodes[ino].fixed = -1;
				elem->nodes[ino].color = -30;
			}
			// x+y-z- = -31
			if (elem->nodes[ino].x == mesh->x_end && elem->nodes[ino].y == mesh->y_start && elem->nodes[ino].z == 0)
			{
				elem->nodes[ino].fixed = -1;
				elem->nodes[ino].color = -31;
			}
			// x+y+z- = -32
			if (elem->nodes[ino].x == mesh->x_end && elem->nodes[ino].y == mesh->y_end && elem->nodes[ino].z == 0)
			{
				elem->nodes[ino].fixed = -1;
				elem->nodes[ino].color = -32;
			}
			// x-y+z- = -33
			if (elem->nodes[ino].x == mesh->x_start && elem->nodes[ino].y == mesh->y_end && elem->nodes[ino].z == 0)
			{
				elem->nodes[ino].fixed = -1;
				elem->nodes[ino].color = -33;
			}
			// x-y-z+ = -34
			if (elem->nodes[ino].x == mesh->x_start && elem->nodes[ino].y == mesh->y_start && elem->nodes[ino].z == 3 * mesh->max_z)
			{
				elem->nodes[ino].fixed = -1;
				elem->nodes[ino].color = -34;
			}
			// x+y-z+ = -35
			if (elem->nodes[ino].x == mesh->x_end && elem->nodes[ino].y == mesh->y_start && elem->nodes[ino].z == 3 * mesh->max_z)
			{
				elem->nodes[ino].fixed = -1;
				elem->nodes[ino].color = -35;
			}
			// x+y+z+ = -36
			if (elem->nodes[ino].x == mesh->x_end && elem->nodes[ino].y == mesh->y_end && elem->nodes[ino].z == 3 * mesh->max_z)
			{
				elem->nodes[ino].fixed = -1;
				elem->nodes[ino].color = -36;
			}
			// x-y+z+ = -37
			if (elem->nodes[ino].x == mesh->x_start && elem->nodes[ino].y == mesh->y_end && elem->nodes[ino].z == 3 * mesh->max_z)
			{
				elem->nodes[ino].fixed = -1;
				elem->nodes[ino].color = -37;
			}
		}
	}

	sc_array_init(&mesh->outsurf, sizeof(octant_t));
	// id global exterior surface
	for (int iel = 0; iel < mesh->elements.elem_count; iel++)
	{
		octant_t *elem = (octant_t *)sc_array_index(&mesh->elements, iel);
		// octant_t * elem = (octant_t*) sc_array_push(&toto);
		// hexa_element_copy(elemOrig,elem);

		int isurf;
		isurf = 0;
		elem->surf[isurf].ext = false;
		if (elem->nodes[FaceNodesMap[isurf][0]].x == 0 && elem->nodes[FaceNodesMap[isurf][1]].x == 0 &&
			elem->nodes[FaceNodesMap[isurf][2]].x == 0 && elem->nodes[FaceNodesMap[isurf][3]].x == 0)
		{
			elem->surf[isurf].ext = true;
			if (deb)
				for (int ino = 0; ino < 4; ino++)
					mesh->part_nodes[elem->nodes[FaceNodesMap[isurf][ino]].id] = isurf + 20;
		}

		isurf = 1;
		elem->surf[isurf].ext = false;
		if (elem->nodes[FaceNodesMap[isurf][0]].x == 3 * mesh->ncellx && elem->nodes[FaceNodesMap[isurf][1]].x == 3 * mesh->ncellx &&
			elem->nodes[FaceNodesMap[isurf][2]].x == 3 * mesh->ncellx && elem->nodes[FaceNodesMap[isurf][3]].x == 3 * mesh->ncellx)
		{
			elem->surf[isurf].ext = true;
			if (deb)
				for (int ino = 0; ino < 4; ino++)
					mesh->part_nodes[elem->nodes[FaceNodesMap[isurf][ino]].id] = isurf + 20;
		}

		isurf = 2;
		elem->surf[isurf].ext = false;
		if (elem->nodes[FaceNodesMap[isurf][0]].y == 0 && elem->nodes[FaceNodesMap[isurf][1]].y == 0 &&
			elem->nodes[FaceNodesMap[isurf][2]].y == 0 && elem->nodes[FaceNodesMap[isurf][3]].y == 0)
		{
			elem->surf[isurf].ext = true;
			if (deb)
				for (int ino = 0; ino < 4; ino++)
					mesh->part_nodes[elem->nodes[FaceNodesMap[isurf][ino]].id] = isurf + 20;
		}

		isurf = 3;
		elem->surf[isurf].ext = false;
		if (elem->nodes[FaceNodesMap[isurf][0]].y == 3 * mesh->ncelly && elem->nodes[FaceNodesMap[isurf][1]].y == 3 * mesh->ncelly &&
			elem->nodes[FaceNodesMap[isurf][2]].y == 3 * mesh->ncelly && elem->nodes[FaceNodesMap[isurf][3]].y == 3 * mesh->ncelly)
		{
			elem->surf[isurf].ext = true;
			if (deb)
				for (int ino = 0; ino < 4; ino++)
					mesh->part_nodes[elem->nodes[FaceNodesMap[isurf][ino]].id] = isurf + 20;
		}

		isurf = 4;
		elem->surf[isurf].ext = false;
		if (elem->nodes[FaceNodesMap[isurf][0]].z == 0 && elem->nodes[FaceNodesMap[isurf][1]].z == 0 &&
			elem->nodes[FaceNodesMap[isurf][2]].z == 0 && elem->nodes[FaceNodesMap[isurf][3]].z == 0)
		{
			elem->surf[isurf].ext = true;
			if (deb)
				for (int ino = 0; ino < 4; ino++)
					mesh->part_nodes[elem->nodes[FaceNodesMap[isurf][ino]].id] = isurf + 20;
		}

		isurf = 5;
		elem->surf[isurf].ext = false;
		if (elem->nodes[FaceNodesMap[isurf][0]].z == 3 * mesh->max_z && elem->nodes[FaceNodesMap[isurf][1]].z == 3 * mesh->max_z &&
			elem->nodes[FaceNodesMap[isurf][2]].z == 3 * mesh->max_z && elem->nodes[FaceNodesMap[isurf][3]].z == 3 * mesh->max_z)
		{
			elem->surf[isurf].ext = true;
			if (deb)
				for (int ino = 0; ino < 4; ino++)
					mesh->part_nodes[elem->nodes[FaceNodesMap[isurf][ino]].id] = isurf + 20;
		}

		bool aux = false;
		for (int isurf = 0; isurf < 6; isurf++)
			if (elem->surf[isurf].ext)
				aux = true;
		if (aux)
		{
			octant_t *elem1 = (octant_t *)sc_array_push(&mesh->outsurf);

			elem1->level = -1;
			elem1->id = elem->id;
			elem1->tem = elem->tem;
			elem1->pad = elem->pad;
			elem1->n_mat = elem->n_mat;
			elem1->pml_id = elem->pml_id;
			elem1->father = elem->id;
			elem1->boundary = elem->boundary;
			elem1->x = elem->x;
			elem1->y = elem->y;
			elem1->z = elem->z;

			for (int ino = 0; ino < 8; ino++)
			{
				elem1->nodes[ino].color = elem->nodes[ino].color;
				elem1->nodes[ino].fixed = elem->nodes[ino].fixed;

				elem1->nodes[ino].id = elem->nodes[ino].id;
				elem1->nodes[ino].x = elem->nodes[ino].x;
				elem1->nodes[ino].y = elem->nodes[ino].y;
				elem1->nodes[ino].z = elem->nodes[ino].z;
			}

			for (int iedge = 0; iedge < 12; iedge++)
			{
				elem1->edge[iedge].coord[0] = elem->edge[iedge].coord[0];
				elem1->edge[iedge].coord[1] = elem->edge[iedge].coord[1];
				elem1->edge[iedge].id = elem->edge[iedge].id;
				elem1->edge[iedge].ref = false;
			}

			for (int isurf = 0; isurf < 6; isurf++)
				elem1->surf[isurf].ext = elem->surf[isurf].ext;
		}
	}

	// id global exterior edges
	for (int iel = 0; iel < mesh->outsurf.elem_count; iel++)
	{
		octant_t *elem = (octant_t *)sc_array_index(&mesh->outsurf, iel);

		// edge 0
		int iedge;
		iedge = 0;
		if (elem->nodes[EdgeVerticesMap[iedge][0]].y == 0 && elem->nodes[EdgeVerticesMap[iedge][1]].y == 0)
		{
			if (elem->nodes[EdgeVerticesMap[iedge][0]].z == 0 && elem->nodes[EdgeVerticesMap[iedge][1]].z == 0)
			{
				elem->edge[iedge].ref = true;
				if (deb)
					mesh->part_nodes[elem->nodes[EdgeVerticesMap[iedge][0]].id] = iedge;
				if (deb)
					mesh->part_nodes[elem->nodes[EdgeVerticesMap[iedge][1]].id] = iedge;
			}
		}

		iedge = 1;
		if (elem->nodes[EdgeVerticesMap[iedge][0]].x == 3 * mesh->ncellx && elem->nodes[EdgeVerticesMap[iedge][1]].x == 3 * mesh->ncellx)
		{
			if (elem->nodes[EdgeVerticesMap[iedge][0]].z == 0 && elem->nodes[EdgeVerticesMap[iedge][1]].z == 0)
			{
				elem->edge[iedge].ref = true;
				if (deb)
					mesh->part_nodes[elem->nodes[EdgeVerticesMap[iedge][0]].id] = iedge;
				if (deb)
					mesh->part_nodes[elem->nodes[EdgeVerticesMap[iedge][1]].id] = iedge;
			}
		}

		iedge = 2;
		if (elem->nodes[EdgeVerticesMap[iedge][0]].y == 3 * mesh->ncelly && elem->nodes[EdgeVerticesMap[iedge][1]].y == 3 * mesh->ncelly)
		{
			if (elem->nodes[EdgeVerticesMap[iedge][0]].z == 0 && elem->nodes[EdgeVerticesMap[iedge][1]].z == 0)
			{
				elem->edge[iedge].ref = true;
				if (deb)
					mesh->part_nodes[elem->nodes[EdgeVerticesMap[iedge][0]].id] = iedge;
				if (deb)
					mesh->part_nodes[elem->nodes[EdgeVerticesMap[iedge][1]].id] = iedge;
			}
		}

		iedge = 3;
		if (elem->nodes[EdgeVerticesMap[iedge][0]].x == 0 && elem->nodes[EdgeVerticesMap[iedge][1]].x == 0)
		{
			if (elem->nodes[EdgeVerticesMap[iedge][0]].z == 0 && elem->nodes[EdgeVerticesMap[iedge][1]].z == 0)
			{
				elem->edge[iedge].ref = true;
				if (deb)
					mesh->part_nodes[elem->nodes[EdgeVerticesMap[iedge][0]].id] = iedge;
				if (deb)
					mesh->part_nodes[elem->nodes[EdgeVerticesMap[iedge][1]].id] = iedge;
			}
		}

		iedge = 4;
		if (elem->nodes[EdgeVerticesMap[iedge][0]].x == 0 && elem->nodes[EdgeVerticesMap[iedge][1]].x == 0)
		{
			if (elem->nodes[EdgeVerticesMap[iedge][0]].y == 0 && elem->nodes[EdgeVerticesMap[iedge][1]].y == 0)
			{
				elem->edge[iedge].ref = true;
				if (deb)
					mesh->part_nodes[elem->nodes[EdgeVerticesMap[iedge][0]].id] = iedge;
				if (deb)
					mesh->part_nodes[elem->nodes[EdgeVerticesMap[iedge][1]].id] = iedge;
			}
		}

		iedge = 5;
		if (elem->nodes[EdgeVerticesMap[iedge][0]].x == 3 * mesh->ncellx && elem->nodes[EdgeVerticesMap[iedge][1]].x == 3 * mesh->ncellx)
		{
			if (elem->nodes[EdgeVerticesMap[iedge][0]].y == 0 && elem->nodes[EdgeVerticesMap[iedge][1]].y == 0)
			{
				elem->edge[iedge].ref = true;
				if (deb)
					mesh->part_nodes[elem->nodes[EdgeVerticesMap[iedge][0]].id] = iedge;
				if (deb)
					mesh->part_nodes[elem->nodes[EdgeVerticesMap[iedge][1]].id] = iedge;
			}
		}

		iedge = 6;
		if (elem->nodes[EdgeVerticesMap[iedge][0]].y == 3 * mesh->ncelly && elem->nodes[EdgeVerticesMap[iedge][1]].y == 3 * mesh->ncelly)
		{
			if (elem->nodes[EdgeVerticesMap[iedge][0]].x == 3 * mesh->ncellx && elem->nodes[EdgeVerticesMap[iedge][1]].x == 3 * mesh->ncellx)
			{
				elem->edge[iedge].ref = true;
				if (deb)
					mesh->part_nodes[elem->nodes[EdgeVerticesMap[iedge][0]].id] = iedge;
				if (deb)
					mesh->part_nodes[elem->nodes[EdgeVerticesMap[iedge][1]].id] = iedge;
			}
		}

		iedge = 7;
		if (elem->nodes[EdgeVerticesMap[iedge][0]].x == 0 && elem->nodes[EdgeVerticesMap[iedge][1]].x == 0)
		{
			if (elem->nodes[EdgeVerticesMap[iedge][0]].y == 3 * mesh->ncelly && elem->nodes[EdgeVerticesMap[iedge][1]].y == 3 * mesh->ncelly)
			{
				elem->edge[iedge].ref = true;
				if (deb)
					mesh->part_nodes[elem->nodes[EdgeVerticesMap[iedge][0]].id] = iedge;
				if (deb)
					mesh->part_nodes[elem->nodes[EdgeVerticesMap[iedge][1]].id] = iedge;
			}
		}

		iedge = 8;
		if (elem->nodes[EdgeVerticesMap[iedge][0]].y == 0 && elem->nodes[EdgeVerticesMap[iedge][1]].y == 0)
		{
			if (elem->nodes[EdgeVerticesMap[iedge][0]].z == 3 * mesh->max_z && elem->nodes[EdgeVerticesMap[iedge][1]].z == 3 * mesh->max_z)
			{
				elem->edge[iedge].ref = true;
				if (deb)
					mesh->part_nodes[elem->nodes[EdgeVerticesMap[iedge][0]].id] = iedge;
				if (deb)
					mesh->part_nodes[elem->nodes[EdgeVerticesMap[iedge][1]].id] = iedge;
			}
		}

		iedge = 9;
		if (elem->nodes[EdgeVerticesMap[iedge][0]].x == 3 * mesh->ncellx && elem->nodes[EdgeVerticesMap[iedge][1]].x == 3 * mesh->ncellx)
		{
			if (elem->nodes[EdgeVerticesMap[iedge][0]].z == 3 * mesh->max_z && elem->nodes[EdgeVerticesMap[iedge][1]].z == 3 * mesh->max_z)
			{
				elem->edge[iedge].ref = true;
				if (deb)
					mesh->part_nodes[elem->nodes[EdgeVerticesMap[iedge][0]].id] = iedge;
				if (deb)
					mesh->part_nodes[elem->nodes[EdgeVerticesMap[iedge][1]].id] = iedge;
			}
		}

		iedge = 10;
		if (elem->nodes[EdgeVerticesMap[iedge][0]].y == 3 * mesh->ncelly && elem->nodes[EdgeVerticesMap[iedge][1]].y == 3 * mesh->ncelly)
		{
			if (elem->nodes[EdgeVerticesMap[iedge][0]].z == 3 * mesh->max_z && elem->nodes[EdgeVerticesMap[iedge][1]].z == 3 * mesh->max_z)
			{
				elem->edge[iedge].ref = true;
				if (deb)
					mesh->part_nodes[elem->nodes[EdgeVerticesMap[iedge][0]].id] = iedge;
				if (deb)
					mesh->part_nodes[elem->nodes[EdgeVerticesMap[iedge][1]].id] = iedge;
			}
		}

		iedge = 11;
		if (elem->nodes[EdgeVerticesMap[iedge][0]].x == 0 && elem->nodes[EdgeVerticesMap[iedge][1]].x == 0)
		{
			if (elem->nodes[EdgeVerticesMap[iedge][0]].z == 3 * mesh->max_z && elem->nodes[EdgeVerticesMap[iedge][1]].z == 3 * mesh->max_z)
			{
				elem->edge[iedge].ref = true;
				if (deb)
					mesh->part_nodes[elem->nodes[EdgeVerticesMap[iedge][0]].id] = iedge;
				if (deb)
					mesh->part_nodes[elem->nodes[EdgeVerticesMap[iedge][1]].id] = iedge;
			}
		}
	}

	if (deb)
	{
		// debug
		for (int iel = 0; iel < mesh->outsurf.elem_count; iel++)
		{
			octant_t *elem = (octant_t *)sc_array_index(&mesh->outsurf, iel);
			for (int ino = 0; ino < 8; ino++)
			{
				if (elem->nodes[ino].fixed == 1)
					mesh->part_nodes[elem->nodes[ino].id] = 1;
			}
		}

		for (int iel = 0; iel < mesh->outsurf.elem_count; iel++)
		{
			octant_t *elem = (octant_t *)sc_array_index(&mesh->outsurf, iel);

			printf("Sou o elemento: %d\n", elem->id);
			printf("Surface: \n", elem->id);
			for (int isurf = 0; isurf < 6; isurf++)
				printf("%s ", elem->surf[isurf].ext ? "true" : "false");

			printf("\nEdge: \n", elem->id);
			for (int isurf = 0; isurf < 12; isurf++)
				printf("%d ", elem->edge[isurf].ref);

			printf("\nNode: \n", elem->id);
			for (int isurf = 0; isurf < 8; isurf++)
				printf("%d ", elem->nodes[isurf].fixed);
			printf("\n", elem->id);
		}
	}
}

void PillowingInterface(hexa_tree_t *mesh, std::vector<double> &coords, std::vector<int> &nodes_b_mat)
{
	int elem_old = mesh->elements.elem_count;
	int nodes_old = mesh->nodes.elem_count;
	bool clamped = true;
	fprintf(mesh->profile, "Time inside PillowingInterface\n");

	// redo the mapping in the nodes
	auto start = std::chrono::steady_clock::now();
	printf("     Redo Node Mapping\n");
	RedoNodeMapping(mesh);
	auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now() - start);
	fprintf(mesh->profile, "    Redo the mapping in the nodes %lld millisecond(s).\n", elapsed.count());
	// std::cout << "Redo the mapping in the nodes "<< elapsed.count() <<" millisecond(s)."<< std::endl;

	// Make the pillow
	start = std::chrono::steady_clock::now();
	printf("     Pillow Layer\n");
	Pillowing(mesh, coords,nodes_b_mat);
	fprintf(mesh->profile, "    Time in PillowLayer %lld millisecond(s).\n", elapsed.count());
	// std::cout << "Time SurfaceIdentification "<< elapsed.count() <<" millisecond(s)."<< std::endl;

	// update the vectors
	mesh->local_n_elements = mesh->elements.elem_count;
	mesh->local_n_nodes = mesh->nodes.elem_count;
	MPI_Allreduce(&mesh->local_n_elements, &mesh->total_n_elements, 1, MPI_LONG, MPI_SUM, MPI_COMM_WORLD);
	MPI_Allreduce(&mesh->local_n_nodes, &mesh->total_n_nodes, 1, MPI_LONG, MPI_SUM, MPI_COMM_WORLD);
	printf("     %d nodes were created in the pillowing process\n", mesh->nodes.elem_count - nodes_old);
	printf("     %d elements were created in the pillowing process\n", mesh->elements.elem_count - elem_old);

	free(mesh->part_nodes);
	mesh->part_nodes = (int *)malloc(mesh->local_n_nodes * sizeof(int));
	for (int ino = 0; ino < mesh->local_n_nodes; ino++)
	{
		mesh->part_nodes[ino] = mesh->mpi_rank;
	}
	for (int ino = 0; ino < nodes_b_mat.size(); ino++)
	{
		mesh->part_nodes[nodes_b_mat[ino]] = 1;
	}

	// Identify the global and local boundaries
	start = std::chrono::steady_clock::now();
	printf("     Surface Identification\n");
	SurfaceIdentification(mesh, coords);
	fprintf(mesh->profile, "    Time in SurfaceIdentification %lld millisecond(s).\n", elapsed.count());
	// std::cout << "Time SurfaceIdentification "<< elapsed.count() <<" millisecond(s)."<< std::endl;
	// NOTE: the interface clamp (z <= sea surface) now runs at the end of
	// Pillowing(), after the pillow layer is built and before its validity
	// dissolve, so the dissolve sees the final clamped geometry. SurfaceIdentification
	// does not move nodes, so the geometry written out is the post-dissolve state.
}
