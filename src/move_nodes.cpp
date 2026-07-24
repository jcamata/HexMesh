#include <gts.h>
#include <glib.h>
#include <vector>
#include <iostream>
using namespace std;
#include <set>
#include <algorithm>
#include <unordered_map>
#include <unordered_set>
#include <sc.h>
#include <sc_io.h>
#include <sc_containers.h>
#include <mpi.h>

#include "hexa.h"
#include "hilbert.h"

#include <ctime>

unsigned vertex_hash_id(const void *v, const void *u) {
	const octant_vertex_t *q = (const octant_vertex_t*) v;
	uint32_t a, b, c;

	a = (uint32_t) q->id;
	b = (uint32_t) 0;
	c = (uint32_t) 1;
	sc_hash_mix(a, b, c);
	sc_hash_final(a, b, c);
	return (unsigned) c;
}

int vertex_equal_id(const void *v, const void *u, const void *w) {
	const octant_vertex_t *e1 = (const octant_vertex_t*) v;
	const octant_vertex_t *e2 = (const octant_vertex_t*) u;

	return (unsigned) (e1->id==e2->id);
}

unsigned el_hash_id(const void *v, const void *u) {
	const octant_t *q = (const octant_t*) v;
	uint64_t a, b, c;

	a = (uint32_t) q->id;
	b = (uint32_t) 0;
	c = (uint32_t) 0;
	sc_hash_mix(a, b, c);
	sc_hash_final(a, b, c);
	return (unsigned) c;
}

int el_equal_id(const void *v, const void *u, const void *w) {
	const octant_t *e1 = (const octant_t*) v;
	const octant_t *e2 = (const octant_t*) u;

	return (unsigned) (e1->id == e2->id);

}

typedef struct {
	bitmask_t coord[3];
	int id;
} node_in_edge_t;

unsigned no_hash_fn1(const void *v, const void *u) {
	const node_in_edge_t *q = (const node_in_edge_t*) v;
	uint32_t a, b, c;

	a = (uint32_t) q->id;
	b = (uint32_t) 0;
	c = (uint32_t) 0;
	sc_hash_mix(a, b, c);
	sc_hash_final(a, b, c);
	return (unsigned) c;
}

int no_equal_fn1(const void *v, const void *u, const void *w) {
	const node_in_edge_t *e1 = (const node_in_edge_t*) v;
	const node_in_edge_t *e2 = (const node_in_edge_t*) u;

	return (unsigned) ((e1->id == e2->id));

}

static unsigned node_hash_id(const void *v, const void *u)
{
	const node_t *q = (const node_t *) v;
	uint32_t a, b, c;

	a = (uint32_t) q->node_id;
	b = 0;
	c = 0;
	sc_hash_mix(a, b, c);
	sc_hash_final(a, b, c);
	return (unsigned) c;
}

static int node_equal_id(const void *v, const void *u, const void *w)
{
	const node_t *e1 = (const node_t *) v;
	const node_t *e2 = (const node_t *) u;
	return (unsigned) (e1->node_id == e2->node_id);
}

static void InitializeOctreeEdgeInfo(octree_t *oct)
{
	oct->edge_info.n_neighbors = 0;
	oct->edge_info.n_intercepted_edges = 0;

	for (int i = 0; i < 26; i++) {
		oct->edge_info.neighbors[i] = -1;
		oct->edge_info.n_intercepted_edges_by_neighbor[i] = 0;
		for (int j = 0; j < 12; j++) {
			oct->edge_info.intercepted_edges_by_neighbor[i][j] = -1;
		}
	}

	for (int i = 0; i < 12; i++) {
		oct->edge_info.intercepted_edges[i] = -1;
	}
}

static bool HasNeighborEdge(const octree_t *oct, int ineighbor, int iedge)
{
	int n = oct->edge_info.n_intercepted_edges_by_neighbor[ineighbor];
	for (int i = 0; i < n; i++) {
		if (oct->edge_info.intercepted_edges_by_neighbor[ineighbor][i] == iedge) {
			return true;
		}
	}
	return false;
}

static void AddNeighborEdgeInfo(octree_t *oct, int neighbor_octree, int iedge)
{
	int ineighbor = -1;
	for (int i = 0; i < oct->edge_info.n_neighbors; i++) {
		if (oct->edge_info.neighbors[i] == neighbor_octree) {
			ineighbor = i;
			break;
		}
	}

	if (ineighbor == -1) {
		if (oct->edge_info.n_neighbors >= 26) {
			return;
		}
		ineighbor = oct->edge_info.n_neighbors++;
		oct->edge_info.neighbors[ineighbor] = neighbor_octree;
		oct->edge_info.n_intercepted_edges_by_neighbor[ineighbor] = 0;
	}

	if (!HasNeighborEdge(oct, ineighbor, iedge) &&
			oct->edge_info.n_intercepted_edges_by_neighbor[ineighbor] < 12) {
		int n = oct->edge_info.n_intercepted_edges_by_neighbor[ineighbor];
		oct->edge_info.intercepted_edges_by_neighbor[ineighbor][n] = iedge;
		oct->edge_info.n_intercepted_edges_by_neighbor[ineighbor]++;
	}
}

static uint64_t GetOctreeEdgeId(hexa_tree_t *mesh, octree_t *oct, int iedge)
{
	int iel0 = EdgeElemOctMap[iedge][0];
	int iel1 = EdgeElemOctMap[iedge][1];

	if (oct->id[iel0] != -1) {
		octant_t *elem = (octant_t *) sc_array_index(&mesh->elements, oct->id[iel0]);
		return elem->edge[iedge].id;
	}
	if (oct->id[iel1] != -1) {
		octant_t *elem = (octant_t *) sc_array_index(&mesh->elements, oct->id[iel1]);
		return elem->edge[iedge].id;
	}
	return 0;
}

static void BuildOctreeNeighborEdgeInfo(hexa_tree_t *mesh)
{
	std::unordered_map<uint64_t, std::vector<std::pair<int, int> > > edge_to_octree;

	for (int ioc = 0; ioc < mesh->oct.elem_count; ioc++) {
		octree_t *oct = (octree_t *) sc_array_index(&mesh->oct, ioc);
		InitializeOctreeEdgeInfo(oct);

		for (int iedge = 0; iedge < 12; iedge++) {
			if (!oct->edge[iedge]) {
				continue;
			}

			if (oct->edge_info.n_intercepted_edges < 12) {
				int n = oct->edge_info.n_intercepted_edges++;
				oct->edge_info.intercepted_edges[n] = iedge;
			}

			uint64_t edge_id = GetOctreeEdgeId(mesh, oct, iedge);
			if (edge_id != 0) {
				edge_to_octree[edge_id].push_back(std::make_pair(ioc, iedge));
			}
		}
	}

	for (std::unordered_map<uint64_t, std::vector<std::pair<int, int> > >::iterator it = edge_to_octree.begin();
			it != edge_to_octree.end(); ++it) {
		std::vector<std::pair<int, int> > &owners = it->second;
		for (size_t i = 0; i < owners.size(); i++) {
			for (size_t j = i + 1; j < owners.size(); j++) {
				int oct_i = owners[i].first;
				int oct_j = owners[j].first;
				int edge_i = owners[i].second;
				int edge_j = owners[j].second;

				octree_t *oi = (octree_t *) sc_array_index(&mesh->oct, oct_i);
				octree_t *oj = (octree_t *) sc_array_index(&mesh->oct, oct_j);

				AddNeighborEdgeInfo(oi, oct_j, edge_i);
				AddNeighborEdgeInfo(oj, oct_i, edge_j);
			}
		}
	}
}

static int CountOctreeInterceptedEdges(const octree_t *oct)
{
	int n_edges = 0;
	for (int iedge = 0; iedge < 12; iedge++) {
		if (oct->edge[iedge]) {
			n_edges++;
		}
	}
	return n_edges;
}

static bool IsOctreeCutPatternRegular(const octree_t *oct)
{
	// A single surface separating the 8 corners of a cube into two nonempty
	// groups must cross at least 3 edges (isolating one corner needs cutting
	// all 3 edges meeting at it -- the cube graph's minimum edge cut). 1-2
	// active edges is therefore not a legitimate simple case, it is evidence
	// of an inconsistent corner classification (e.g. ClassifyOctreeCorners'
	// ray-triangle test landing wrong on one corner).
	int n_edges = CountOctreeInterceptedEdges(oct);
	return n_edges >= 3 && n_edges <= 6;
}

struct moved_node_position_t
{
	int node;
	double coord[3];
};

static bool IsCompleteOctree(const octree_t *oct)
{
	for (int i = 0; i < 8; i++) {
		if (oct->id[i] == -1) {
			return false;
		}
	}
	return true;
}

static void AddUniqueNode(std::vector<int> &nodes, int node)
{
	if (std::find(nodes.begin(), nodes.end(), node) == nodes.end()) {
		nodes.push_back(node);
	}
}

static void AddMovedNode(std::vector<moved_node_position_t> &moves, int node, double x, double y, double z)
{
	for (size_t i = 0; i < moves.size(); i++) {
		if (moves[i].node == node) {
			return;
		}
	}

	moved_node_position_t move;
	move.node = node;
	move.coord[0] = x;
	move.coord[1] = y;
	move.coord[2] = z;
	moves.push_back(move);
}

static bool IsNodeMoved(const std::vector<char> &node_moved, int node)
{
	return node >= 0 && node < (int) node_moved.size() && node_moved[node];
}

static void ApplyMovedNodes(std::vector<double> &coords,
		std::vector<int> &nodes_b_mat,
		std::vector<char> &node_moved,
		const std::vector<moved_node_position_t> &moves)
{
	for (size_t i = 0; i < moves.size(); i++) {
		int node = moves[i].node;
		if (node < 0 || node >= (int) node_moved.size()) {
			continue;
		}
		coords[3 * node + 0] = moves[i].coord[0];
		coords[3 * node + 1] = moves[i].coord[1];
		coords[3 * node + 2] = moves[i].coord[2];
		if (!node_moved[node]) {
			nodes_b_mat.push_back(node);
			node_moved[node] = 1;
		}
	}
}

static void CollectIrregularOctreeMovableNodes(hexa_tree_t *mesh, octree_t *oct, std::vector<int> &nodes)
{
	nodes.clear();
	for (int iel = 0; iel < 8; iel++) {
		octant_t *elem = (octant_t *) sc_array_index(&mesh->elements, oct->id[iel]);
		for (int ino = 0; ino < 8; ino++) {
			if (elem->nodes[ino].fixed == 0) {
				AddUniqueNode(nodes, elem->nodes[ino].id);
			}
		}
	}
}

static int AverageAdjacentMovedNodes(hexa_tree_t *mesh,
		octree_t *oct,
		const std::vector<double> &coords,
		const std::vector<char> &node_moved,
		int target_node,
		double avg[3])
{
	std::vector<int> adjacent_nodes;

	for (int iel = 0; iel < 8; iel++) {
		octant_t *elem = (octant_t *) sc_array_index(&mesh->elements, oct->id[iel]);
		for (int iedge = 0; iedge < 12; iedge++) {
			int node0 = elem->nodes[EdgeVerticesMap[iedge][0]].id;
			int node1 = elem->nodes[EdgeVerticesMap[iedge][1]].id;

			if (node0 == target_node && IsNodeMoved(node_moved, node1)) {
				AddUniqueNode(adjacent_nodes, node1);
			}
			if (node1 == target_node && IsNodeMoved(node_moved, node0)) {
				AddUniqueNode(adjacent_nodes, node0);
			}
		}
	}

	avg[0] = 0;
	avg[1] = 0;
	avg[2] = 0;
	for (size_t i = 0; i < adjacent_nodes.size(); i++) {
		int node = adjacent_nodes[i];
		avg[0] += coords[3 * node + 0];
		avg[1] += coords[3 * node + 1];
		avg[2] += coords[3 * node + 2];
	}

	int count = adjacent_nodes.size();
	if (count > 0) {
		avg[0] /= count;
		avg[1] /= count;
		avg[2] /= count;
	}
	return count;
}

static int AverageMovedNodesInOctree(hexa_tree_t *mesh,
		octree_t *oct,
		const std::vector<double> &coords,
		const std::vector<char> &node_moved,
		double avg[3])
{
	std::vector<int> moved_nodes;

	for (int iel = 0; iel < 8; iel++) {
		octant_t *elem = (octant_t *) sc_array_index(&mesh->elements, oct->id[iel]);
		for (int ino = 0; ino < 8; ino++) {
			int node = elem->nodes[ino].id;
			if (IsNodeMoved(node_moved, node)) {
				AddUniqueNode(moved_nodes, node);
			}
		}
	}

	avg[0] = 0;
	avg[1] = 0;
	avg[2] = 0;
	for (size_t i = 0; i < moved_nodes.size(); i++) {
		int node = moved_nodes[i];
		avg[0] += coords[3 * node + 0];
		avg[1] += coords[3 * node + 1];
		avg[2] += coords[3 * node + 2];
	}

	int count = moved_nodes.size();
	if (count > 0) {
		avg[0] /= count;
		avg[1] /= count;
		avg[2] /= count;
	}
	return count;
}

static void RegularizeSkippedOctreeNodes(hexa_tree_t *mesh,
		std::vector<double> &coords,
		std::vector<int> &nodes_b_mat)
{
	int n_nodes = coords.size() / 3;
	std::vector<char> node_moved(n_nodes, 0);
	for (size_t i = 0; i < nodes_b_mat.size(); i++) {
		int node = nodes_b_mat[i];
		if (node >= 0 && node < n_nodes) {
			node_moved[node] = 1;
		}
	}

	int skipped_octrees = 0;
	for (int ioc = 0; ioc < mesh->oct.elem_count; ioc++) {
		octree_t *oct = (octree_t *) sc_array_index(&mesh->oct, ioc);
		if (IsCompleteOctree(oct) && !IsOctreeCutPatternRegular(oct)) {
			skipped_octrees++;
		}
	}

	int regularized_nodes = 0;
	const int max_iters = 8;
	for (int iter = 0; iter < max_iters; iter++) {
		std::vector<moved_node_position_t> moves;

		for (int ioc = 0; ioc < mesh->oct.elem_count; ioc++) {
			octree_t *oct = (octree_t *) sc_array_index(&mesh->oct, ioc);
			if (!IsCompleteOctree(oct) || IsOctreeCutPatternRegular(oct)) {
				continue;
			}

			std::vector<int> candidates;
			CollectIrregularOctreeMovableNodes(mesh, oct, candidates);
			for (size_t inode = 0; inode < candidates.size(); inode++) {
				int node = candidates[inode];
				if (IsNodeMoved(node_moved, node)) {
					continue;
				}

				double avg[3];
				int count = AverageAdjacentMovedNodes(mesh, oct, coords, node_moved, node, avg);
				if (count >= 2) {
					AddMovedNode(moves, node, avg[0], avg[1], avg[2]);
				}
			}
		}

		if (moves.empty()) {
			break;
		}

		regularized_nodes += moves.size();
		ApplyMovedNodes(coords, nodes_b_mat, node_moved, moves);
	}

	std::vector<moved_node_position_t> fallback_moves;
	for (int ioc = 0; ioc < mesh->oct.elem_count; ioc++) {
		octree_t *oct = (octree_t *) sc_array_index(&mesh->oct, ioc);
		if (!IsCompleteOctree(oct) || IsOctreeCutPatternRegular(oct)) {
			continue;
		}

		double oct_avg[3];
		int oct_moved_count = AverageMovedNodesInOctree(mesh, oct, coords, node_moved, oct_avg);
		if (oct_moved_count < 2) {
			continue;
		}

		std::vector<int> candidates;
		CollectIrregularOctreeMovableNodes(mesh, oct, candidates);
		for (size_t inode = 0; inode < candidates.size(); inode++) {
			int node = candidates[inode];
			if (!IsNodeMoved(node_moved, node)) {
				AddMovedNode(fallback_moves, node, oct_avg[0], oct_avg[1], oct_avg[2]);
			}
		}
	}

	regularized_nodes += fallback_moves.size();
	ApplyMovedNodes(coords, nodes_b_mat, node_moved, fallback_moves);

	if (skipped_octrees > 0) {
		printf(" Regularized %d nodes from %d skipped octrees\n",
				regularized_nodes, skipped_octrees);
	}
}

static bool IsOctreeCutTopologicallyValid(const octree_t *oct)
{
	// 2-coloring of the 8 hex corners: cut edges must connect corners of different colors,
	// non-cut edges must connect corners of the same color.
	// Inconsistency means the cut cannot be explained by a single planar cut.
	static const int adj[8][3][2] = {
		{{0,1},{3,3},{4,4}},
		{{0,0},{1,2},{5,5}},
		{{1,1},{2,3},{6,6}},
		{{2,2},{3,0},{7,7}},
		{{4,0},{8,5},{11,7}},
		{{5,1},{8,4},{9,6}},
		{{6,2},{9,5},{10,7}},
		{{7,3},{10,6},{11,4}}
	};
	int color[8];
	for (int i = 0; i < 8; i++) color[i] = -1;
	color[0] = 0;
	int queue[8], head = 0, tail = 0;
	queue[tail++] = 0;
	while (head < tail) {
		int c = queue[head++];
		for (int i = 0; i < 3; i++) {
			int edge = adj[c][i][0];
			int nb   = adj[c][i][1];
			int expected = oct->edge[edge] ? (1 - color[c]) : color[c];
			if (color[nb] == -1) {
				color[nb] = expected;
				queue[tail++] = nb;
			} else if (color[nb] != expected) {
				return false;
			}
		}
	}
	return true;
}

static void FixOctreesCutTopology(hexa_tree_t *mesh)
{
	// Hexahedron edges: each entry is {corner_a, corner_b}
	static const int hex_edges[12][2] = {
		{0,1},{1,2},{2,3},{3,0},
		{0,4},{1,5},{2,6},{3,7},
		{4,5},{5,6},{6,7},{7,4}
	};

	int total_removed = 0;
	for (int ioc = 0; ioc < mesh->oct.elem_count; ioc++) {
		octree_t *oct = (octree_t *) sc_array_index(&mesh->oct, ioc);
		if (IsOctreeCutTopologicallyValid(oct)) continue;

		// Brute-force search over all 128 vertex bipartitions.
		// Find the one that removes the fewest existing cut edges.
		int best_removed = 13;
		int best_mask = 0;

		for (int mask = 1; mask < 128; mask++) {
			int removed = 0;
			bool feasible = true;
			for (int e = 0; e < 12; e++) {
				int va = hex_edges[e][0];
				int vb = hex_edges[e][1];
				bool should_cut = (((mask >> va) & 1) != ((mask >> vb) & 1));
				if (oct->edge[e] && !should_cut) {
					removed++;             // cut edge inconsistent with bipartition → must remove
				} else if (!oct->edge[e] && should_cut) {
					feasible = false;      // non-cut edge inconsistent → can't add cuts, skip
					break;
				}
			}
			if (feasible && removed < best_removed) {
				best_removed = removed;
				best_mask = mask;
			}
		}

		// Remove cut edges that are inconsistent with the best bipartition
		for (int e = 0; e < 12; e++) {
			if (!oct->edge[e]) continue;
			int va = hex_edges[e][0];
			int vb = hex_edges[e][1];
			bool should_cut = (((best_mask >> va) & 1) != ((best_mask >> vb) & 1));
			if (!should_cut) {
				oct->edge[e] = false;
				total_removed++;
			}
		}
	}

	// Recompute face flags from the corrected edge flags.
	// oct->face[] were set by IdentifyMovableNodes using the original (pre-fix) edges.
	// Stale face flags cause count==0 divisions in ProjectFreeNodes face center computation.
	for (int ioc = 0; ioc < mesh->oct.elem_count; ioc++) {
		octree_t *oct = (octree_t *) sc_array_index(&mesh->oct, ioc);
		oct->face[0] = oct->edge[4] || oct->edge[11] || oct->edge[7]  || oct->edge[3];
		oct->face[1] = oct->edge[5] || oct->edge[1]  || oct->edge[6]  || oct->edge[9];
		oct->face[2] = oct->edge[0] || oct->edge[5]  || oct->edge[8]  || oct->edge[4];
		oct->face[3] = oct->edge[2] || oct->edge[6]  || oct->edge[10] || oct->edge[7];
		oct->face[4] = oct->edge[8] || oct->edge[9]  || oct->edge[10] || oct->edge[11];
		oct->face[5] = oct->edge[0] || oct->edge[1]  || oct->edge[2]  || oct->edge[3];
	}

	if (total_removed > 0)
		printf("    Fixed %d inconsistent cut edges across topologically invalid octrees\n", total_removed);
}

static void RestrictInvalidOctreeEdgeNodes(hexa_tree_t *mesh)
{
	// For each octree edge, which two elements (iel) own it and which node index is movable
	static const int edge_to_iel[12][2] = {
		{0,1},{1,2},{2,3},{3,0},
		{0,4},{1,5},{2,6},{3,7},
		{4,5},{5,6},{6,7},{7,4}
	};
	// edge_iel_node[e][0] = node index in element edge_to_iel[e][0]
	// edge_iel_node[e][1] = node index in element edge_to_iel[e][1]
	static const int edge_iel_node[12][2] = {
		{1,0},{2,1},{3,2},{0,3},
		{4,0},{5,1},{6,2},{7,3},
		{5,4},{6,5},{7,6},{4,7}
	};

	// Pass 1: collect all edge node IDs that belong to at least one topologically valid octree
	std::unordered_set<int> valid_edge_nodes;
	for (int ioc = 0; ioc < mesh->oct.elem_count; ioc++) {
		octree_t *oct = (octree_t *) sc_array_index(&mesh->oct, ioc);
		if (!IsOctreeCutPatternRegular(oct) || !IsOctreeCutTopologicallyValid(oct)) continue;
		for (int iedge = 0; iedge < 12; iedge++) {
			if (!oct->edge[iedge]) continue;
			for (int k = 0; k < 2; k++) {
				int iel = edge_to_iel[iedge][k];
				int ino = edge_iel_node[iedge][k];
				if (oct->id[iel] == -1) continue;
				octant_t *elem = (octant_t *) sc_array_index(&mesh->elements, oct->id[iel]);
				valid_edge_nodes.insert(elem->nodes[ino].id);
			}
		}
	}

	// Pass 2: for regular but topologically invalid octrees, restrict edge nodes
	// that are NOT shared with any valid octree
	int count = 0;
	for (int ioc = 0; ioc < mesh->oct.elem_count; ioc++) {
		octree_t *oct = (octree_t *) sc_array_index(&mesh->oct, ioc);
		if (!IsOctreeCutPatternRegular(oct)) continue;    // already skipped by ProjectFreeNodes
		if (IsOctreeCutTopologicallyValid(oct)) continue; // valid, no restriction needed
		for (int iedge = 0; iedge < 12; iedge++) {
			if (!oct->edge[iedge]) continue;
			for (int k = 0; k < 2; k++) {
				int iel = edge_to_iel[iedge][k];
				int ino = edge_iel_node[iedge][k];
				if (oct->id[iel] == -1) continue;
				octant_t *elem = (octant_t *) sc_array_index(&mesh->elements, oct->id[iel]);
				int nid = elem->nodes[ino].id;
				if (valid_edge_nodes.count(nid)) continue; // shared with valid octree, keep free
				octant_node_t *gnode = (octant_node_t *) sc_array_index(&mesh->nodes, nid);
				if (gnode->fixed == 0) { gnode->fixed = 2; count++; }
				elem->nodes[ino].fixed = 2;
			}
		}
	}
	printf("    Restricted %d edge nodes from topologically invalid octrees\n", count);
}

static void PropagateElementCutEdges(hexa_tree_t *mesh)
{
	// Each geometric edge has a unique cantor-pair ID from its two node IDs.
	// The surface detection in GetInterceptedElements works per-element-per-edge
	// independently. When the surface crosses an edge shared between two elements
	// from different octrees, the intersection may be detected in one but missed
	// in the other (floating point, tangent cases). This leaves the neighboring
	// octree without the cut → FixOctreesCutTopology removes its remaining cuts →
	// gap in the mesh.
	//
	// Fix: for every geometric edge that is cut in ANY element, mark it as cut
	// in ALL elements that share it. This must run before IdentifyMovableNodes
	// so the octree-level edge flags are built from complete information.

	std::unordered_map<uint64_t, std::vector<std::pair<int, int>>> edge_map;
	edge_map.reserve(mesh->elements.elem_count * 6);

	for (int iel = 0; iel < mesh->elements.elem_count; iel++) {
		octant_t *elem = (octant_t *) sc_array_index(&mesh->elements, iel);
		for (int e = 0; e < 12; e++) {
			uint64_t eid = elem->edge[e].id;
			if (eid != 0) {
				edge_map[eid].emplace_back(iel, e);
			}
		}
	}

	int propagated = 0;
	for (auto &kv : edge_map) {
		auto &refs = kv.second;
		if (refs.size() < 2) continue;

		bool any_ref = false;
		for (auto &p : refs) {
			octant_t *elem = (octant_t *) sc_array_index(&mesh->elements, p.first);
			if (elem->edge[p.second].ref) { any_ref = true; break; }
		}
		if (!any_ref) continue;

		for (auto &p : refs) {
			octant_t *elem = (octant_t *) sc_array_index(&mesh->elements, p.first);
			if (!elem->edge[p.second].ref) {
				elem->edge[p.second].ref = true;
				propagated++;
			}
		}
	}

	if (propagated > 0)
		printf("    Propagated %d cut flags to shared element edges\n", propagated);
}

/*
GtsPoint* FoundInterception(hexa_tree_t* mesh,std::vector<double>& coords,int node1, int node2){
	GtsPoint *point = NULL;
	GtsVertex *v1 = gts_vertex_new(gts_vertex_class(), coords[node1 * 3], coords[node1 * 3 + 1], coords[node1 * 3 + 2]);
	GtsVertex *v2 = gts_vertex_new(gts_vertex_class(), coords[node2 * 3], coords[node2 * 3 + 1], coords[node2 * 3 + 2]);

	GtsSegment *segments = gts_segment_new(gts_segment_class(), v1, v2);
	GtsBBox *sb = gts_bbox_segment(gts_bbox_class(), segments);
	GSList* list = gts_bb_tree_overlap(mesh->gdata.bbt, sb);
	//if (list == NULL) continue;
	while (list) {
		GtsBBox *b = GTS_BBOX(list->data);
		point = SegmentTriangleIntersectionCgal(segments, GTS_TRIANGLE(b->bounded));
		if (point) {
			break;
		}
		list = list->next;
	}
	return point;
}
*/
// Adjacency table for the 8 hex corners: adj[c][k] = {edge_index, neighbor_corner}
static const int oct_corner_adj[8][3][2] = {
	{{0,1},{3,3},{4,4}},
	{{0,0},{1,2},{5,5}},
	{{1,1},{2,3},{6,6}},
	{{2,2},{3,0},{7,7}},
	{{4,0},{8,5},{11,7}},
	{{5,1},{8,4},{9,6}},
	{{6,2},{9,5},{10,7}},
	{{7,3},{10,6},{11,4}}
};

// 2-color the 8 hex corners from oct->edge[]: returns false if inconsistent.
// B[i] = 0 or 1 (relative inside/outside label, arbitrary orientation).
static bool GetOctreeBipartition(const octree_t *oct, int B[8])
{
	for (int i = 0; i < 8; i++) B[i] = -1;
	B[0] = 0;
	int queue[8], head = 0, tail = 0;
	queue[tail++] = 0;
	while (head < tail) {
		int c = queue[head++];
		for (int k = 0; k < 3; k++) {
			int edge = oct_corner_adj[c][k][0];
			int nb   = oct_corner_adj[c][k][1];
			int expected = oct->edge[edge] ? (1 - B[c]) : B[c];
			if (B[nb] == -1) {
				B[nb] = expected;
				queue[tail++] = nb;
			} else if (B[nb] != expected) {
				return false;
			}
		}
	}
	return true;
}

void ProjectFreeNodes(hexa_tree_t* mesh, std::vector<double>& coords, std::vector<int>& nodes_b_mat) {

	// For each face, for each of the 4 edges (FaceEdgesMap order):
	// which octree element position and local node index holds the already-moved edge node.
	static const int FaceEdgeElem[6][4] = {
		{0, 0, 7, 7}, // face 0: edges {3,4,7,11}
		{2, 5, 2, 5}, // face 1: edges {1,5,6,9}
		{0, 0, 5, 5}, // face 2: edges {0,4,5,8}
		{2, 2, 7, 7}, // face 3: edges {2,6,7,10}
		{5, 5, 7, 7}, // face 4: edges {8,9,10,11}
		{0, 2, 2, 0}  // face 5: edges {0,1,2,3}
	};
	static const int FaceEdgeNode[6][4] = {
		{3, 4, 3, 4},
		{1, 1, 6, 6},
		{1, 4, 1, 4},
		{3, 6, 3, 6},
		{4, 6, 6, 4},
		{1, 1, 3, 3}
	};
	// Face center node location (written in pass 2)
	static const int FaceCenterElem[6] = {0, 5, 0, 2, 5, 0};
	static const int FaceCenterNode[6] = {7, 2, 5, 7, 7, 2};
	// Face center node location (read in pass 3)
	static const int CenterFaceElem[6] = {0, 2, 0, 2, 5, 0};
	static const int CenterFaceNode[6] = {7, 5, 5, 7, 7, 2};

	// pending: first-write-wins per node; flushed to coords after each pass
	std::unordered_map<int, std::array<double, 3>> pending;
	std::unordered_set<int> in_b;

	// The mat0/mat1 interface is the sea floor: no projected interface node may
	// sit above the sea surface (z = SEA_LEVEL). The bathy surface (gdata) fills
	// on-land vertices with a high sentinel (~+4996 m), so near the coast a cut
	// edge can intersect a steep sentinel triangle and snap the node far above
	// the sea surface, shearing the flat-water (mat-1) cells into spikes. Clamp
	// every projected node to z <= SEA_LEVEL - SURFACE_EPS.
	//
	// SURFACE_EPS must be a SMALL FIXED length, not a fraction of the element
	// height. At the shallow coast the sea-surface node IS also an interface
	// node, so any downward clamp lowers the water lid too; an eps of ~10% of
	// the element height (~46 m here) caved the coastal surface into visible
	// pits below z=0. A fixed 1 m only flattens the +z spikes flush to the lid
	// and keeps a 1 m positive thickness where the floor meets the lid (no
	// zero-volume collapse), while the coastal dip stays visually negligible.
	const double SEA_LEVEL = 0.0;
	const double SURFACE_EPS = 1.0;   // metres below the sea surface
	double sea_clamp = SEA_LEVEL - SURFACE_EPS;

	// Snapshot of the original (valid lattice) coords, captured before ANY
	// projection so `record` below can tell a domain-lid node from an
	// interface node. Also used at the end of this function by the validity
	// limiter to pull invalid/sheared projected elements back.
	std::vector<double> coords0 = coords;

	auto record = [&](int node, double x, double y, double z) {
		// The domain lid (z = SEA_LEVEL) is a fixed OUTER boundary of the model
		// box, not part of the bathymetry -- every column, land or water, has
		// its own top face there. A lid node can still be the "inner" corner of
		// a cut octree edge/face/body-diagonal purely because a DEEPER sibling
		// in the same octree group is on the interface; the diagonal segment
		// used for the face/body-centre intersection test then spans lid-to-
		// interface and can pick up a spurious hit far from where this node
		// actually belongs (seen concretely: a lid node dragged to -2777 m,
		// i.e. two whole lattice levels down, while its own element's other 3
		// top-face corners stayed correctly at 0). Any node that started
		// exactly at the lid never needs bathymetry conformance, so skip it.
		if (coords0[3*node+2] >= SEA_LEVEL - 1e-6) return;
		if (z > sea_clamp) z = sea_clamp;
		pending.emplace(node, std::array<double, 3>{x, y, z});
	};

	auto flush = [&]() {
		for (auto &kv : pending) {
			int n = kv.first;
			coords[3*n+0] = kv.second[0];
			coords[3*n+1] = kv.second[1];
			coords[3*n+2] = kv.second[2];
			if (in_b.insert(n).second)
				nodes_b_mat.push_back(n);
		}
		pending.clear();
	};

	// --- Projection validity limiter: setup ----------------------------------
	// (coords0 already snapshotted above, before `record` was defined.)
	// Projecting interface nodes onto a steep/curved bathy collapses or inverts
	// the thin interface-layer elements (zext median 12 m vs ~500 m wide), which
	// the pillowing then inherits as zero-volume elements that render as holes.
	// After projection we pull the offending moved nodes back toward their
	// original positions (see end of function).
	auto hexvol = [&](octant_t *e) -> double {
		// Use the SAME node order the h5 writer emits (assign_elem_nodes =
		// {4,5,6,7,0,1,2,3} in hexa_h5.cpp). e->nodes is not the standard hex
		// order assumed by the tet decomposition T and the reorder is not a clean
		// z-flip, so computing T directly over e->nodes yields a geometrically
		// wrong volume — the limiter then misses real interface inversions.
		static const int ord[8] = {4,5,6,7,0,1,2,3};
		double X[8], Y[8], Z[8];
		for (int i = 0; i < 8; i++) {
			int id = e->nodes[ord[i]].id;
			X[i] = coords[3*id]; Y[i] = coords[3*id+1]; Z[i] = coords[3*id+2];
		}
		static const int T[6][4] = {{0,1,2,6},{0,2,3,6},{0,3,7,6},{0,7,4,6},{0,4,5,6},{0,5,1,6}};
		double v = 0.0;
		for (auto &t : T) {
			double ax=X[t[1]]-X[t[0]], ay=Y[t[1]]-Y[t[0]], az=Z[t[1]]-Z[t[0]];
			double bx=X[t[2]]-X[t[0]], by=Y[t[2]]-Y[t[0]], bz=Z[t[2]]-Z[t[0]];
			double cx=X[t[3]]-X[t[0]], cy=Y[t[3]]-Y[t[0]], cz=Z[t[3]]-Z[t[0]];
			v += (ax*(by*cz-bz*cy) - ay*(bx*cz-bz*cx) + az*(bx*cy-by*cx))/6.0;
		}
		return v;
	};
	std::vector<double> ref_vol(mesh->elements.elem_count);
	for (int iel = 0; iel < mesh->elements.elem_count; iel++)
		ref_vol[iel] = hexvol((octant_t*) sc_array_index(&mesh->elements, iel));

	// Pass 1: find surface intersection on each cut octree edge and snap the two
	// adjacent inner nodes to that intersection point.
	for (int ioc = 0; ioc < mesh->oct.elem_count; ioc++) {
		octree_t *oct = (octree_t*) sc_array_index(&mesh->oct, ioc);
		if (!IsCompleteOctree(oct) || !IsOctreeCutPatternRegular(oct)) continue;

		for (int iedge = 0; iedge < 12; iedge++) {
			if (!oct->edge[iedge]) continue;

			int iel0 = EdgeElemOctMap[iedge][0];
			int iel1 = EdgeElemOctMap[iedge][1];
			octant_t *elem0 = (octant_t*) sc_array_index(&mesh->elements, oct->id[iel0]);
			octant_t *elem1 = (octant_t*) sc_array_index(&mesh->elements, oct->id[iel1]);

			int nA = elem0->nodes[iel0].id;   // outer corner of octant a
			int nB = elem1->nodes[iel1].id;   // outer corner of octant b

			double x1 = coords[3*nA], y1 = coords[3*nA+1], z1 = coords[3*nA+2];
			double x2 = coords[3*nB], y2 = coords[3*nB+1], z2 = coords[3*nB+2];
			const double ext = 0.02;
			double dx = x2-x1, dy = y2-y1, dz = z2-z1;

			GtsVertex *v1 = gts_vertex_new(gts_vertex_class(), x1-ext*dx, y1-ext*dy, z1-ext*dz);
			GtsVertex *v2 = gts_vertex_new(gts_vertex_class(), x2+ext*dx, y2+ext*dy, z2+ext*dz);
			GtsSegment *seg = gts_segment_new(gts_segment_class(), v1, v2);
			GtsBBox *bb = gts_bbox_segment(gts_bbox_class(), seg);

			GSList *list = gts_bb_tree_overlap(mesh->gdata.bbt, bb);
			GtsPoint *pt = NULL;
			while (list) {
				GtsBBox *b = GTS_BBOX(list->data);
				pt = mesh->input.CgalUse
					? SegmentTriangleIntersectionCgal(seg, GTS_TRIANGLE(b->bounded))
					: SegmentTriangleIntersection(seg, GTS_TRIANGLE(b->bounded));
				if (pt) break;
				list = list->next;
			}
			if (!pt) continue;

			// Each edge is shared between elem0 (inner node = vertex[1]) and
			// elem1 (inner node = vertex[0]); both get snapped to the same point.
			record(elem0->nodes[EdgeVerticesMap[iedge][1]].id, pt->x, pt->y, pt->z);
			record(elem1->nodes[EdgeVerticesMap[iedge][0]].id, pt->x, pt->y, pt->z);
		}
	}
	flush();

	// Pass 2: place the face-center node at the active face diagonal ∩ surface.
	// Falls back to centroid of edge nodes when no active diagonal or no intersection.
	for (int ioc = 0; ioc < mesh->oct.elem_count; ioc++) {
		octree_t *oct = (octree_t*) sc_array_index(&mesh->oct, ioc);
		if (!IsCompleteOctree(oct) || !IsOctreeCutPatternRegular(oct)) continue;

		octant_t *elems[8];
		for (int i = 0; i < 8; i++)
			elems[i] = (octant_t*) sc_array_index(&mesh->elements, oct->id[i]);

		int B[8];
		if (!GetOctreeBipartition(oct, B)) continue; // inconsistent edges → skip

		for (int iface = 0; iface < 6; iface++) {
			if (!oct->face[iface]) continue;

			int c0 = FaceNodesMap[iface][0], c1 = FaceNodesMap[iface][1];
			int c2 = FaceNodesMap[iface][2], c3 = FaceNodesMap[iface][3];
			int da0 = -1, da1 = -1;
			if (B[c0] != B[c2])      { da0 = c0; da1 = c2; }
			else if (B[c1] != B[c3]) { da0 = c1; da1 = c3; }

			int center = elems[FaceCenterElem[iface]]->nodes[FaceCenterNode[iface]].id;

			auto face_centroid_fallback = [&]() {
				double xx = 0, yy = 0, zz = 0; int count = 0;
				for (int k = 0; k < 4; k++) {
					if (!oct->edge[FaceEdgesMap[iface][k]]) continue;
					int nid = elems[FaceEdgeElem[iface][k]]->nodes[FaceEdgeNode[iface][k]].id;
					xx += coords[3*nid+0]; yy += coords[3*nid+1]; zz += coords[3*nid+2];
					count++;
				}
				if (count > 0) record(center, xx/count, yy/count, zz/count);
			};

			if (da0 == -1) { face_centroid_fallback(); continue; }

			int nA = elems[da0]->nodes[da0].id;
			int nB = elems[da1]->nodes[da1].id;
			double x1 = coords[3*nA], y1 = coords[3*nA+1], z1 = coords[3*nA+2];
			double x2 = coords[3*nB], y2 = coords[3*nB+1], z2 = coords[3*nB+2];
			const double ext = 0.02;
			double dx = x2-x1, dy = y2-y1, dz = z2-z1;

			GtsVertex *v1 = gts_vertex_new(gts_vertex_class(), x1-ext*dx, y1-ext*dy, z1-ext*dz);
			GtsVertex *v2 = gts_vertex_new(gts_vertex_class(), x2+ext*dx, y2+ext*dy, z2+ext*dz);
			GtsSegment *seg = gts_segment_new(gts_segment_class(), v1, v2);
			GtsBBox *bb = gts_bbox_segment(gts_bbox_class(), seg);

			GSList *list = gts_bb_tree_overlap(mesh->gdata.bbt, bb);
			GtsPoint *pt = NULL;
			while (list) {
				GtsBBox *b = GTS_BBOX(list->data);
				pt = mesh->input.CgalUse
					? SegmentTriangleIntersectionCgal(seg, GTS_TRIANGLE(b->bounded))
					: SegmentTriangleIntersection(seg, GTS_TRIANGLE(b->bounded));
				if (pt) break;
				list = list->next;
			}

			if (!pt) { face_centroid_fallback(); continue; }
			record(center, pt->x, pt->y, pt->z);
		}
	}
	flush();

	// Pass 3: place the octree-center node at the active body diagonal ∩ surface.
	// Falls back to centroid of face-center nodes when no intersection is found.
	static const int body_diags[4][2] = {{0,6},{1,7},{2,4},{3,5}};
	for (int ioc = 0; ioc < mesh->oct.elem_count; ioc++) {
		octree_t *oct = (octree_t*) sc_array_index(&mesh->oct, ioc);
		if (!IsCompleteOctree(oct) || !IsOctreeCutPatternRegular(oct)) continue;

		octant_t *elems[8];
		for (int i = 0; i < 8; i++)
			elems[i] = (octant_t*) sc_array_index(&mesh->elements, oct->id[i]);

		int B[8];
		if (!GetOctreeBipartition(oct, B)) continue; // inconsistent edges → skip

		int da0 = -1, da1 = -1;
		for (int d = 0; d < 4 && da0 == -1; d++) {
			int a = body_diags[d][0], b = body_diags[d][1];
			if (B[a] != B[b]) { da0 = a; da1 = b; }
		}
		if (da0 == -1) continue;

		int center_node = elems[0]->nodes[6].id;

		int nA = elems[da0]->nodes[da0].id;
		int nB = elems[da1]->nodes[da1].id;
		double x1 = coords[3*nA], y1 = coords[3*nA+1], z1 = coords[3*nA+2];
		double x2 = coords[3*nB], y2 = coords[3*nB+1], z2 = coords[3*nB+2];
		const double ext = 0.02;
		double dx = x2-x1, dy = y2-y1, dz = z2-z1;

		GtsVertex *v1 = gts_vertex_new(gts_vertex_class(), x1-ext*dx, y1-ext*dy, z1-ext*dz);
		GtsVertex *v2 = gts_vertex_new(gts_vertex_class(), x2+ext*dx, y2+ext*dy, z2+ext*dz);
		GtsSegment *seg = gts_segment_new(gts_segment_class(), v1, v2);
		GtsBBox *bb = gts_bbox_segment(gts_bbox_class(), seg);

		GSList *list = gts_bb_tree_overlap(mesh->gdata.bbt, bb);
		GtsPoint *pt = NULL;
		while (list) {
			GtsBBox *b = GTS_BBOX(list->data);
			pt = mesh->input.CgalUse
				? SegmentTriangleIntersectionCgal(seg, GTS_TRIANGLE(b->bounded))
				: SegmentTriangleIntersection(seg, GTS_TRIANGLE(b->bounded));
			if (pt) break;
			list = list->next;
		}

		if (!pt) {
			double xx = 0, yy = 0, zz = 0; int count = 0;
			for (int iface = 0; iface < 6; iface++) {
				if (!oct->face[iface]) continue;
				int nid = elems[CenterFaceElem[iface]]->nodes[CenterFaceNode[iface]].id;
				xx += coords[3*nid+0]; yy += coords[3*nid+1]; zz += coords[3*nid+2];
				count++;
			}
			if (count > 0) record(center_node, xx/count, yy/count, zz/count);
			continue;
		}
		record(center_node, pt->x, pt->y, pt->z);
	}
	flush();

	RegularizeSkippedOctreeNodes(mesh, coords, nodes_b_mat);

	// --- Projection validity limiter: repair ---------------------------------
	// Pull moved nodes back toward their original lattice positions wherever an
	// incident element collapsed or inverted (volume sign flip vs reference, or
	// magnitude < 5% of its original). Halving toward the known-valid lattice
	// config each iteration is monotone toward validity, so it never inverts a
	// currently-valid element. This keeps the interface-layer elements valid so
	// the pillow layer built on them does not collapse into holes. Cost: only
	// nodes of invalid elements move, so it converges in a few iterations.
	{
		int n_nodes_loc = (int)(coords.size() / 3);
		const int MAXIT = 50;
		int it = 0, nbad = 0, npull_total = 0;
		std::vector<char> pull(n_nodes_loc, 0);
		for (it = 0; it < MAXIT; it++) {
			std::fill(pull.begin(), pull.end(), 0);
			nbad = 0;
			for (int iel = 0; iel < mesh->elements.elem_count; iel++) {
				octant_t *e = (octant_t*) sc_array_index(&mesh->elements, iel);
				double v = hexvol(e), rv = ref_vol[iel];
				double av = v < 0 ? -v : v, arv = rv < 0 ? -rv : rv;
				bool vol_ok = (v * rv > 0.0 && av >= 0.05 * arv);
				// Shear check removed: it could not tell legitimate steep
				// interface conformance (a real cliff genuinely needs a large
				// horizontal offset) from the original cross-octree-group
				// mismatch bug, so it was undoing correct MovingNodes output
				// wherever the terrain was steep -- confirmed by comparing
				// against move_nodes.cpp from commits 8c8fb74/bfc2409, which
				// conform to the interface correctly with only the volume
				// check below.
				if (vol_ok) continue;   // still valid
				nbad++;
				for (int i = 0; i < 8; i++) pull[e->nodes[i].id] = 1;
			}
			if (nbad == 0) break;
			for (int n = 0; n < n_nodes_loc; n++) if (pull[n]) {
				coords[3*n+0] = 0.5 * (coords[3*n+0] + coords0[3*n+0]);
				coords[3*n+1] = 0.5 * (coords[3*n+1] + coords0[3*n+1]);
				coords[3*n+2] = 0.5 * (coords[3*n+2] + coords0[3*n+2]);
				npull_total++;
			}
		}
		printf("    Projection limiter: %d invalid elements after %d iters (%d node pull-backs)\n",
				nbad, it, npull_total);
	}

	// Deduplicate nodes_b_mat via hash and update per-element node fixity.
	bool clamped = true;
	sc_hash_array_t *hash_fixed = sc_hash_array_new(sizeof(node_t), node_hash_id, node_equal_id, &clamped);
	size_t pos;
	node_t key;

	for (int ii = 0; ii < (int)nodes_b_mat.size(); ii++) {
		int n = nodes_b_mat[ii];
		key.node_id = n;
		key.coord[0] = coords[3*n+0];
		key.coord[1] = coords[3*n+1];
		key.coord[2] = coords[3*n+2];
		node_t *r = (node_t*) sc_hash_array_insert_unique(hash_fixed, &key, &pos);
		if (r) *r = key;
	}

	for (int ioc = 0; ioc < mesh->oct.elem_count; ioc++) {
		octree_t *oct = (octree_t*) sc_array_index(&mesh->oct, ioc);
		for (int iel = 0; iel < 8; iel++) {
			octant_t *elem = (octant_t*) sc_array_index(&mesh->elements, oct->id[iel]);
			for (int ino = 0; ino < 8; ino++) {
				int n = elem->nodes[ino].id;
				elem->nodes[ino].fixed = 2;
				key.node_id = n;
				key.coord[0] = coords[3*n+0];
				key.coord[1] = coords[3*n+1];
				key.coord[2] = coords[3*n+2];
				if (sc_hash_array_lookup(hash_fixed, &key, &pos))
					elem->nodes[ino].fixed = 1;
			}
		}
	}

	nodes_b_mat.clear();
	for (int ii = 0; ii < (int)hash_fixed->a.elem_count; ii++) {
		node_t *node = (node_t*) sc_array_index(&hash_fixed->a, ii);
		nodes_b_mat.push_back(node->node_id);
	}
	printf(" Total of %d fixed nodes\n", (int)nodes_b_mat.size());

	// Assign node colors: propagate material-side label across octree cut edges.
	for (int ioc = 0; ioc < mesh->oct.elem_count; ioc++) {
		octree_t *oct = (octree_t*) sc_array_index(&mesh->oct, ioc);
		octant_t *elems[8];
		for (int iel = 0; iel < 8; iel++)
			elems[iel] = (octant_t*) sc_array_index(&mesh->elements, oct->id[iel]);

		for (int ino = 0; ino < 8; ino++) {
			for (int k = 0; k < 3; k++) {
				int nb = OctNeighbourMap[ino][k];
				if (!oct->edge[VertexEdgeMap[ino][k]]) {
					if (elems[ino]->nodes[ino].color == -1)
						elems[ino]->nodes[ino].color = 1;
					elems[nb]->nodes[nb].color = elems[ino]->nodes[ino].color;
				} else {
					if (elems[ino]->nodes[ino].color != 2 &&
							elems[nb]->nodes[nb].color == -1)
						elems[nb]->nodes[nb].color = 2;
				}
			}
		}
	}
}

static void VerifyOctreeCutTemplates(hexa_tree_t *mesh)
{
	int n_standard   = 0;
	int n_topo_fail  = 0;
	int n_complex    = 0;
	int n_degenerate = 0;
	int n_uncut      = 0;

	for (int ioc = 0; ioc < mesh->oct.elem_count; ioc++) {
		octree_t *oct = (octree_t*) sc_array_index(&mesh->oct, ioc);
		if (!IsCompleteOctree(oct)) continue;

		int n = CountOctreeInterceptedEdges(oct);

		if (n == 0)      { n_uncut++;      continue; }
		if (n <= 2)      { n_degenerate++; continue; }
		if (n >= 7)      { n_complex++;    continue; }

		if (!IsOctreeCutTopologicallyValid(oct))
			n_topo_fail++;
		else
			n_standard++;
	}

	printf("    Cut template verification:\n");
	printf("      Standard (3-6 edges, topologically valid): %d\n", n_standard);
	printf("      Topological failure (edges still inconsistent): %d\n", n_topo_fail);
	printf("      Complex (7+ edges, impossible plane cut):   %d\n", n_complex);
	printf("      Degenerate (1-2 edges):                     %d\n", n_degenerate);
	printf("      Uncut octrees (0 edges):                    %d\n", n_uncut);

	if (n_topo_fail > 0)
		printf("    WARNING: %d octrees still topologically inconsistent after fix\n", n_topo_fail);
	if (n_complex > 0)
		printf("    WARNING: %d octrees have 7+ cut edges — will be skipped\n", n_complex);
}

void IdentifyMovableNodes(hexa_tree_t* mesh) {

	// Derive canonical cut edges from the domain bipartition.
	// n_mat for each octant's outer corner was set by ClassifyOctreeCorners
	// (called from MovingNodes before this function).
	for (int ioc = 0; ioc < mesh->oct.elem_count; ioc++) {
		octree_t *oct = (octree_t*) sc_array_index(&mesh->oct, ioc);

		if (!IsCompleteOctree(oct)) {
			memset(oct->edge, 0, sizeof(oct->edge));
			memset(oct->face, 0, sizeof(oct->face));
			continue;
		}

		octant_t *elems[8];
		int B[8];
		bool all_same = true;
		for (int i = 0; i < 8; i++) {
			elems[i] = (octant_t*) sc_array_index(&mesh->elements, oct->id[i]);
			B[i] = (elems[i]->n_mat == 0) ? 0 : 1;
		}
		for (int i = 1; i < 8; i++)
			if (B[i] != B[0]) { all_same = false; break; }

		if (all_same) {
			memset(oct->edge, 0, sizeof(oct->edge));
			memset(oct->face, 0, sizeof(oct->face));
			continue;
		}

		for (int e = 0; e < 12; e++) {
			int a = EdgeElemOctMap[e][0];
			int b = EdgeElemOctMap[e][1];
			oct->edge[e] = (B[a] != B[b]);
		}

		for (int f = 0; f < 6; f++) {
			oct->face[f] = false;
			for (int k = 0; k < 4; k++)
				if (oct->edge[FaceEdgesMap[f][k]]) { oct->face[f] = true; break; }
		}
	}

	for (int ino = 0; ino < mesh->nodes.elem_count; ino++) {
		octant_node_t *node = (octant_node_t*) sc_array_index(&mesh->nodes, ino);
		node->fixed = 0;
	}
	for (int iel = 0; iel < mesh->elements.elem_count; iel++) {
		octant_t *elem = (octant_t*) sc_array_index(&mesh->elements, iel);
		for (int ino = 0; ino < 8; ino++) {
			elem->nodes[ino].fixed = 0;
			elem->nodes[ino].color = -1;
		}
	}

	VerifyOctreeCutTemplates(mesh);
	BuildOctreeNeighborEdgeInfo(mesh);
}

void DoOctree(hexa_tree_t* mesh){

	bool               clamped = true;
	size_t              position;
	sc_hash_array_t    *indep_vertex;
	indep_vertex    = (sc_hash_array_t *)sc_hash_array_new(sizeof (octant_vertex_t), vertex_hash_id, vertex_equal_id, &clamped);

	/////////////////
	// create the vertex structure
	for (int iel = 0; iel < mesh->elements.elem_count; ++iel) {

		octant_t *elem = (octant_t*) sc_array_index(&mesh->elements, iel);

		for (int ino = 0; ino < 8; ino++){

			octant_vertex_t key;
			key.id = elem->nodes[ino].id;

			octant_vertex_t* vert = (octant_vertex_t*) sc_hash_array_insert_unique (indep_vertex, &key, &position);
			if(vert != NULL){
				vert->id = elem->nodes[ino].id;
				vert->list_elem = 1;
				vert->elem[vert->list_elem-1] = elem->id;
			}else{
				vert = (octant_vertex_t*) sc_array_index(&indep_vertex->a, position);
				vert->elem[vert->list_elem] = elem->id;
				vert->list_elem++;
			}
		}
	}

	//creating the octree structure
	sc_array_init(&mesh->oct, sizeof(octree_t));
	sc_array_reset(&mesh->oct);

	//hash for the elements
	sc_hash_array_t * hashOctree = (sc_hash_array_t *)sc_hash_array_new(sizeof (octant_t), el_hash_id, el_equal_id, &clamped);
	octant_t* r;
	bool elem_lookup;
	bool vert_lookup;

	//find the elements...
	for (int iel = 0; iel < mesh->elements.elem_count; ++iel) {
		octant_t *cut = (octant_t*) sc_array_index(&mesh->elements, iel);

		size_t position;
		octant_t key;
		key.id = cut->id;

		elem_lookup = sc_hash_array_lookup(hashOctree, &key, &position);

		int temp_id[8];
		bool temp_cut = false;;

		for(int ii = 0; ii<8 ; ii++){
			temp_id[ii] = -1;
		}

		if(elem_lookup){

		}else{

			octant_vertex_t key1;
			size_t position1;
			key1.id = cut->nodes[6].id;

			vert_lookup = sc_hash_array_lookup(indep_vertex, &key1, &position1);

			if(vert_lookup){
				octant_vertex_t* vert = (octant_vertex_t*) sc_array_index(&indep_vertex->a, position1);

				for(int iele = 0; iele < vert->list_elem; iele++){

					int test = vert->elem[iele];

					octant_t *elem = (octant_t*) sc_array_index(&mesh->elements, test);
					key.id = elem->id;

					if(cut->id==elem->id){
						r = (octant_t*) sc_hash_array_insert_unique(hashOctree,&key,&position);
						temp_id[0] = elem->id;
						if(r!=NULL){
							r->id= elem->id;
						}else{
							r =  (octant_t*) sc_array_index(&hashOctree->a,position);
						}
						if(elem->pad==-1 || cut->pad == -1){
							temp_cut = true;
						}
					}//x+

					if(cut->nodes[1].id == elem->nodes[0].id && cut->nodes[2].id == elem->nodes[3].id &&
							cut->nodes[5].id==elem->nodes[4].id && cut->nodes[6].id==elem->nodes[7].id){
						temp_id[1] = elem->id;
						r = (octant_t*) sc_hash_array_insert_unique(hashOctree,&key,&position);
						if(r!=NULL){
							r->id= elem->id;
						}else{
							r =  (octant_t*) sc_array_index(&hashOctree->a,position);
						}
						if(elem->pad==-1 || cut->pad == -1){
							temp_cut = true;
						}
					}//x+y+

					if(cut->nodes[2].id==elem->nodes[0].id && cut->nodes[6].id==elem->nodes[4].id){
						temp_id[2]= elem->id;
						r = (octant_t*) sc_hash_array_insert_unique(hashOctree,&key,&position);
						if(r!=NULL){
							r->id= elem->id;
						}else{
							r =  (octant_t*) sc_array_index(&hashOctree->a,position);
						}
						if(elem->pad==-1 || cut->pad == -1){
							temp_cut = true;
						}
					}//y+

					if(cut->nodes[3].id==elem->nodes[0].id && cut->nodes[2].id==elem->nodes[1].id &&
							cut->nodes[6].id==elem->nodes[5].id && cut->nodes[7].id==elem->nodes[4].id){
						temp_id[3]= elem->id;
						r = (octant_t*) sc_hash_array_insert_unique(hashOctree,&key,&position);
						if(r!=NULL){
							r->id= elem->id;
						}else{
							r =  (octant_t*) sc_array_index(&hashOctree->a,position);
						}
						if(elem->pad==-1 || cut->pad == -1){
							temp_cut = true;
						}
					}//z-

					if(cut->nodes[4].id==elem->nodes[0].id && cut->nodes[5].id==elem->nodes[1].id &&
							cut->nodes[6].id==elem->nodes[2].id && cut->nodes[7].id==elem->nodes[3].id){
						temp_id[4]= elem->id;
						r = (octant_t*) sc_hash_array_insert_unique(hashOctree,&key,&position);
						if(r!=NULL){
							r->id= elem->id;
						}else{
							r =  (octant_t*) sc_array_index(&hashOctree->a,position);
						}
						if(elem->pad==-1 || cut->pad == -1){
							temp_cut = true;
						}
					}//z- x+

					if(cut->nodes[5].id==elem->nodes[0].id && cut->nodes[6].id==elem->nodes[3].id){
						temp_id[5]= elem->id;
						r = (octant_t*) sc_hash_array_insert_unique(hashOctree,&key,&position);
						if(r!=NULL){
							r->id= elem->id;
						}else{
							r =  (octant_t*) sc_array_index(&hashOctree->a,position);
						}
						if(elem->pad==-1 || cut->pad == -1){
							temp_cut = true;
						}
					}//z- x+ y+

					if(cut->nodes[6].id==elem->nodes[0].id){
						temp_id[6]= elem->id;
						r = (octant_t*) sc_hash_array_insert_unique(hashOctree,&key,&position);
						if(r!=NULL){
							r->id= elem->id;
						}else{
							r =  (octant_t*) sc_array_index(&hashOctree->a,position);
						}
						if(elem->pad==-1 || cut->pad == -1){
							temp_cut = true;
						}
					} //z- y+

					if(cut->nodes[6].id==elem->nodes[1].id && cut->nodes[7].id==elem->nodes[0].id){
						temp_id[7]= elem->id;
						r = (octant_t*) sc_hash_array_insert_unique(hashOctree,&key,&position);
						if(r!=NULL){
							r->id= elem->id;
						}else{
							r =  (octant_t*) sc_array_index(&hashOctree->a,position);
						}
						if(elem->pad==-1 || cut->pad == -1){
							temp_cut = true;
						}
					}
				}
			}
		}

		if(temp_cut){
			//create the octree
			octree_t * oc = (octree_t*) sc_array_push(&mesh->oct);
			// fill with -1
			for(int i = 0; i<8; i++){
				oc->id[i] = temp_id[i];
			}
			oc->cut = temp_cut;

			//initialization of edges and faces as non intercepted
			for(int iedge = 0; iedge<12;iedge++){
				oc->edge[iedge] = false;
			}
			for(int isurf = 0; isurf<6;isurf++){
				oc->face[isurf] = false;
			}
			InitializeOctreeEdgeInfo(oc);
		}
	}

	//for debug
	if(false){
		for(int ioc = 0; ioc< mesh->oct.elem_count; ioc++){
			octree_t *oc = (octree_t*) sc_array_index(&mesh->oct, ioc);
			printf("Octree number:%d, cut?:%d\n",ioc,oc->cut);
			printf("Element ids: ");
			for(int i = 0; i<8; i++){
				printf("%d ",oc->id[i]);
			}
			printf("\n");
			printf("Edges: ");

			for(int i = 0; i<12; i++){
				printf("%d ",oc->edge[i]);
			}
			printf("\n");
			printf("Faces: ");

			for(int i = 0; i<6; i++){
				printf("%d ",oc->face[i]);
			}
			printf("\n");
		}
	}
}

void MovingNodes(hexa_tree_t* mesh, std::vector<double>& coords, std::vector<int>& nodes_b_mat, const char* surface) {

	bool deb = false;

	time_t tstart, tend;

	tstart = time(0);
	printf("    Building the octree structure...\n");
	DoOctree(mesh);
	tend = time(0);
	//	cout << "Time in DoOctree "<< difftime(tend, tstart) <<" second(s)."<< endl;

	tstart = time(0);
	printf("    Classifying octree corners...\n");
	ClassifyOctreeCorners(mesh, coords);
	tend = time(0);

	tstart = time(0);
	printf("    Identifying the movable nodes...\n");
	IdentifyMovableNodes(mesh);
	tend = time(0);
	cout << "blabla bla bla" << endl;

	//cout << "Time in IdentifyMovableNodes "<< difftime(tend, tstart) <<" second(s)."<< endl;

	if(deb){
		char fdname[80];
		sprintf(fdname,"free_node_%04d.txt", mesh->mpi_rank);
		FILE* fnode = fopen(fdname,"w");
		for(int ino = 0; ino < mesh->nodes.elem_count; ino ++){
			octant_node_t* node = (octant_node_t*) sc_array_index (&mesh->nodes, ino);
			if(node->fixed==0){
				int nnode = 3*node->id;
				fprintf(fnode,"%f %f %f\n",coords[nnode+0],coords[nnode+1],coords[nnode+2]);
			}
		}
		fclose(fnode);
	}
	cout << "blabla bla bla" << endl;

	if(deb){
		char fdname[80];
		sprintf(fdname,"debug_edges_%04d.txt", mesh->mpi_rank);
		FILE* dedges = fopen(fdname,"w");
		for(int ioc = 0; ioc < mesh->oct.elem_count; ioc ++){
			octree_t* oct = (octree_t*) sc_array_index (&mesh->oct, ioc);
			for(int iel = 0; iel < 8; iel++){
				if(oct->id[iel] !=-1){
					octant_t* elem = (octant_t*) sc_array_index (&mesh->nodes, oct->id[iel]);
					//printf("%s\n", oct->edge[0] ? "true" : "false")
					fprintf(dedges,"El %d\n", oct->id[iel]);
					for(int iedge = 0; iedge < 12;iedge++){
						fprintf(dedges,"%s ", elem->edge[iedge].ref ? "T" : "F");
					}
					fprintf(dedges,"\n");
				}
			}

		}
		fclose(dedges);
	}
	cout << "blabla bla bla" << endl;

	tstart = time(0);
	printf("    Make the projection of the nodes into the surface...\n");
	nodes_b_mat.clear();
	ProjectFreeNodes(mesh,coords,nodes_b_mat);
	tend = time(0);
	//cout << "Time in ProjectFreeNodes "<< difftime(tend, tstart) <<" second(s)."<< endl;

	if(deb){
		int8_t* flag_nodes = (int8_t*) malloc(sizeof (int8_t) * mesh->local_n_nodes);
		memset(flag_nodes, 0, sizeof (int8_t) * mesh->local_n_nodes);

		for (int i = 0; i < mesh->local_n_nodes; i++) {
			flag_nodes[i]=0;
		}

		//verifica a criacao os nos livres e fixos na criacao dos octrees

		//std::cout <<  mesh->elements.elem_count << " coisinhas para mover" << std::endl;
		for (int iel = 0; iel < mesh->elements.elem_count; ++iel) {
			octant_t *elem = (octant_t*) sc_array_index(&mesh->elements, iel);
			for(int ino = 0; ino<8;ino++){
				if(flag_nodes[elem->nodes[ino].id]==2){

				}else{
					flag_nodes[elem->nodes[ino].id] = elem->nodes[ino].fixed;
				}
			}
		}
		for (int i = 0; i < mesh->local_n_nodes; i++) {
			mesh->part_nodes[i] = flag_nodes[i];
		}
		free(flag_nodes);
		for(int i = 0; i<nodes_b_mat.size();i++){
			mesh->part_nodes[nodes_b_mat[i]] = 1;
		}
	}
}
