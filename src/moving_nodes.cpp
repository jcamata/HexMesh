#include <limits>
#include <cmath>
#include <gts.h>
#include <glib.h>
#include <vector>
#include <iostream>
using namespace std;
#include <set>
#include <map>
#include <algorithm>
#include <cstring>
#include <unordered_map>
#include <unordered_set>
#include <sc.h>
#include <sc_io.h>
#include <sc_containers.h>
#include <mpi.h>

#include "hexa.h"
#include "hilbert.h"
#include "mesh_geom.h"
#include "pillow_geometry.h"
#include "verify_mesh.h"
#include "optimize_mesh.h"

#include <ctime>
#include <cstdlib>
#include <climits>

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

unsigned octree_hash_fn(const void *v, const void *u) {
	const octree_t *oct = (const octree_t *) v;
	uint32_t a = (uint32_t) oct->id[0];
	uint32_t b = (uint32_t) oct->id[1];
	uint32_t c = (uint32_t) oct->id[6];
	sc_hash_mix(a, b, c);
	sc_hash_final(a, b, c);
	return (unsigned) c;
}

int octree_equal_fn(const void *v, const void *u, const void *w) {
	const octree_t *o1 = (const octree_t *) v;
	const octree_t *o2 = (const octree_t *) u;
	for (int i = 0; i < 8; i++) {
		if (o1->id[i] != o2->id[i]) return 0;
	}
	return 1;
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
	int hist[13] = {0};
	for (int ioc = 0; ioc < mesh->oct.elem_count; ioc++) {
		octree_t *oct = (octree_t *) sc_array_index(&mesh->oct, ioc);
		if (!IsCompleteOctree(oct)) continue;
		int ne = CountOctreeInterceptedEdges(oct);
		if (ne > 0) hist[ne]++;
		if (!IsOctreeCutPatternRegular(oct)) skipped_octrees++;
	}
	printf("    Cut octrees by intercepted-edge count:");
	for (int i = 1; i <= 12; i++) if (hist[i]) printf(" %d:%d", i, hist[i]);
	printf("\n");

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

// ---------------------------------------------------------------------------
// Warp the x-y lattice onto the coastline.
//
// The topography is followed well because GetMeshFromSurface spreads the vertical deformation
// over the WHOLE column: the ncellz nodes are stretched between zmax and zmin, nobody takes an
// isolated jump. Horizontally the opposite used to happen -- one interface node was snapped
// sideways onto the wall while its neighbours stayed on the rigid lattice, shearing the element
// between them until the projection limiter dragged it back.
//
// The coastline is a VERTICAL wall, so its intersection sits at the same (x, y) for the whole
// depth of the wall. Moving an entire column horizontally is therefore geometrically
// consistent: columns stay vertical and straight, hexahedra keep their shape, and connectivity
// is untouched -- this is purely a change of coordinates on the integer lattice that
// octant_node_t already carries in (x, y).
//
//   anchors : where a horizontal octree edge cuts the surface, the column is placed exactly
//             there (the user's "coloco o no exatamente onde a surface cortou o octree");
//   field   : that displacement is diffused to the surrounding columns by Laplace smoothing,
//             with the domain border pinned at zero, so it decays with distance and the faces
//             stay fixed;
//   z       : re-sampled per column afterwards, because GetMeshFromSurface sampled it at the
//             old (x, y).
// ---------------------------------------------------------------------------
static std::vector<char>   g_warp_col_state;   // 0 untouched, 1 diffused, 2 anchored, 3 clamped
static std::vector<double> g_warp_col_disp;
static std::vector<double> g_warp_col_ux, g_warp_col_uy;
static double g_warp_hx = 0.0, g_warp_hy = 0.0;

static void WarpLatticeToCoastline(hexa_tree_t *mesh, std::vector<double> &coords)
{
	if (!mesh->gdata.bbt || !mesh->tdata.bbt) return;

	const int nx = mesh->ncellx + 1, ny = mesh->ncelly + 1;
	const size_t ncol = (size_t) nx * ny;
	auto COL = [nx](int i, int j) { return (size_t) j * nx + i; };

	std::vector<double> ax(ncol, 0.0), ay(ncol, 0.0);
	std::vector<int>    an(ncol, 0);

	// --- Phase 1: anchors from horizontal cut octree edges ---------------------
	for (int ioc = 0; ioc < mesh->oct.elem_count; ioc++) {
		octree_t *oct = (octree_t*) sc_array_index(&mesh->oct, ioc);
		if (!IsCompleteOctree(oct)) continue;

		for (int iedge = 0; iedge < 12; iedge++) {
			if (!oct->edge[iedge]) continue;

			int iel0 = EdgeElemOctMap[iedge][0];
			int iel1 = EdgeElemOctMap[iedge][1];
			octant_t *elem0 = (octant_t*) sc_array_index(&mesh->elements, oct->id[iel0]);
			octant_t *elem1 = (octant_t*) sc_array_index(&mesh->elements, oct->id[iel1]);

			octant_node_t *ndA = &elem0->nodes[iel0];
			octant_node_t *ndB = &elem1->nodes[iel1];
			// only horizontal edges see the vertical wall
			if (ndA->z != ndB->z) continue;

			int nA = ndA->id, nB = ndB->id;
			double x1 = coords[3*nA], y1 = coords[3*nA+1], z1 = coords[3*nA+2];
			double x2 = coords[3*nB], y2 = coords[3*nB+1], z2 = coords[3*nB+2];
			const double ext = 0.02;
			double dx = x2-x1, dy = y2-y1, dz = z2-z1;

			GtsVertex *v1 = gts_vertex_new(gts_vertex_class(), x1-ext*dx, y1-ext*dy, z1-ext*dz);
			GtsVertex *v2 = gts_vertex_new(gts_vertex_class(), x2+ext*dx, y2+ext*dy, z2+ext*dz);
			GtsSegment *seg = gts_segment_new(gts_segment_class(), v1, v2);
			GtsBBox *bb = gts_bbox_segment(gts_bbox_class(), seg);

			GSList *list = gts_bb_tree_overlap(mesh->gdata.bbt, bb);
			const double mx = 0.5*(x1+x2), my = 0.5*(y1+y2), mz = 0.5*(z1+z2);
			GtsPoint *pt = NULL; double best = 0.0;
			for (GSList *l = list; l; l = l->next) {
				GtsBBox *b = GTS_BBOX(l->data);
				GtsPoint *q = mesh->input.CgalUse
					? SegmentTriangleIntersectionCgal(seg, GTS_TRIANGLE(b->bounded))
					: SegmentTriangleIntersection(seg, GTS_TRIANGLE(b->bounded));
				if (!q) continue;
				double d = (q->x-mx)*(q->x-mx) + (q->y-my)*(q->y-my) + (q->z-mz)*(q->z-mz);
				if (!pt || d < best) { if (pt) gts_object_destroy(GTS_OBJECT(pt)); pt = q; best = d; }
				else gts_object_destroy(GTS_OBJECT(q));
			}
			if (list) g_slist_free(list);

			gts_object_destroy(GTS_OBJECT(bb));
			gts_object_destroy(GTS_OBJECT(seg));

			if (!pt) continue;

			// the two inner nodes of this octree edge share the midpoint column
			for (int side = 0; side < 2; side++) {
				octant_node_t *nd = side == 0
					? &elem0->nodes[EdgeVerticesMap[iedge][1]]
					: &elem1->nodes[EdgeVerticesMap[iedge][0]];
				if (nd->x < 0 || nd->x >= nx || nd->y < 0 || nd->y >= ny) continue;
				size_t c = COL(nd->x, nd->y);
				ax[c] += pt->x - coords[3*nd->id+0];
				ay[c] += pt->y - coords[3*nd->id+1];
				an[c]++;
			}
			gts_object_destroy(GTS_OBJECT(pt));
		}
	}

	int n_anchor = 0;
	std::vector<char> fixed(ncol, 0);
	std::vector<double> ux(ncol, 0.0), uy(ncol, 0.0);
	for (size_t c = 0; c < ncol; c++) if (an[c]) {
		ux[c] = ax[c] / an[c];
		uy[c] = ay[c] / an[c];
		fixed[c] = 1;
		n_anchor++;
	}
	// the domain faces stay put
	for (int i = 0; i < nx; i++) { fixed[COL(i,0)] = 1; fixed[COL(i,ny-1)] = 1;
	                               ux[COL(i,0)] = uy[COL(i,0)] = 0.0;
	                               ux[COL(i,ny-1)] = uy[COL(i,ny-1)] = 0.0; }
	for (int j = 0; j < ny; j++) { fixed[COL(0,j)] = 1; fixed[COL(nx-1,j)] = 1;
	                               ux[COL(0,j)] = uy[COL(0,j)] = 0.0;
	                               ux[COL(nx-1,j)] = uy[COL(nx-1,j)] = 0.0; }

	if (n_anchor == 0) { printf("    Lattice warp: no coastline anchors, skipped\n"); return; }

	// --- Phase 2: Laplace-diffuse the displacement over the free columns -------
	const int SWEEPS = 300;
	std::vector<double> vx(ux), vy(uy);
	for (int s = 0; s < SWEEPS; s++) {
		for (int j = 1; j < ny-1; j++)
			for (int i = 1; i < nx-1; i++) {
				size_t c = COL(i,j);
				if (fixed[c]) continue;
				vx[c] = 0.25*(ux[COL(i-1,j)] + ux[COL(i+1,j)] + ux[COL(i,j-1)] + ux[COL(i,j+1)]);
				vy[c] = 0.25*(uy[COL(i-1,j)] + uy[COL(i+1,j)] + uy[COL(i,j-1)] + uy[COL(i,j+1)]);
			}
		ux.swap(vx); uy.swap(vy);
	}

	// --- Phase 4: do not fold the lattice -------------------------------------
	double hx = (mesh->tdata.bbox->x2 - mesh->tdata.bbox->x1) / (double) mesh->ncellx;
	double hy = (mesh->tdata.bbox->y2 - mesh->tdata.bbox->y1) / (double) mesh->ncelly;
	// A column may not travel more than ~half a cell or it crosses its neighbour. Clamp each
	// column on its own: a single greedy column must not scale down the whole field.
	const double lim = 0.90 * std::min(hx, hy);
	double umax = 0.0;
	int n_clamped = 0;
	for (size_t c = 0; c < ncol; c++) {
		double m = std::sqrt(ux[c]*ux[c] + uy[c]*uy[c]);
		if (m > umax) umax = m;
		if (m > lim) { ux[c] *= lim/m; uy[c] *= lim/m; n_clamped++; }
	}
	// record per-column state for the inversion map (diagnostic only)
	g_warp_col_state.assign(ncol, 0);
	g_warp_col_disp.assign(ncol, 0.0);
	g_warp_col_ux = ux; g_warp_col_uy = uy;
	g_warp_hx = hx;     g_warp_hy = hy;
	for (size_t c = 0; c < ncol; c++) {
		double m = std::sqrt(ux[c]*ux[c] + uy[c]*uy[c]);
		g_warp_col_disp[c] = m;
		if (m > 1e-9) g_warp_col_state[c] = 1;          // diffused
		if (an[c])    g_warp_col_state[c] = 2;          // anchored
		if (m >= lim - 1e-9 && m > 1e-9) g_warp_col_state[c] = 3;  // clamped at the limit
	}

	// --- Phase 3: apply, then re-sample z -------------------------------------
	double zmin = -mesh->input.z;
	std::vector<double> zmax_col(ncol, 0.0);
	std::vector<char> have_z(ncol, 0);

	for (int i = 0; i < mesh->nodes.elem_count; i++) {
		octant_node_t *n = (octant_node_t*) sc_array_index(&mesh->nodes, i);
		if (n->x < 0 || n->x >= nx || n->y < 0 || n->y >= ny) continue;
		size_t c = COL(n->x, n->y);
		coords[3*n->id+0] += ux[c];
		coords[3*n->id+1] += uy[c];

		if (!have_z[c]) {
			double zmax;
			// Vertical ray-cast: the true surface height under (x, y), via the
			// same helper GetMeshFromSurface uses -- not the nearest-point
			// Euclidean distance, which is wrong on steep slopes and
			// near-vertical walls (e.g. a coastline).
			if (!eval_gts_height(mesh, 1000, coords[3*n->id+0], coords[3*n->id+1], zmax)) {
				GtsPoint *p = gts_point_new(gts_point_class(), coords[3*n->id+0], coords[3*n->id+1], mesh->tdata.bbox->z2);
				double d = gts_bb_tree_point_distance(mesh->tdata.bbt, p, distance, NULL);
				zmax = mesh->tdata.bbox->z2 - d;
				gts_object_destroy(GTS_OBJECT(p));
				printf("    Lattice warp: vertical ray missed topo surface at (%f, %f); "
				       "using nearest-point fallback\n", coords[3*n->id+0], coords[3*n->id+1]);
			}
			zmax_col[c] = zmax;
			have_z[c] = 1;
		}
		double dz = (zmax_col[c] - zmin) / (double) mesh->ncellz;
		coords[3*n->id+2] = zmax_col[c] - n->z * dz;
	}

	printf("    Lattice warp: %d anchored columns of %zu, max displacement %.1f m of a "
	       "%.1f x %.1f m cell, %d columns clamped at %.1f m\n",
	       n_anchor, ncol, umax, hx, hy, n_clamped, lim);
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

	double dom_min_x = 1e300, dom_max_x = -1e300;
	double dom_min_y = 1e300, dom_max_y = -1e300;
	double dom_min_z = 1e300, dom_max_z = -1e300;
	for (size_t i = 0; i < coords0.size() / 3; i++) {
		if (coords0[3*i+0] < dom_min_x) dom_min_x = coords0[3*i+0];
		if (coords0[3*i+0] > dom_max_x) dom_max_x = coords0[3*i+0];
		if (coords0[3*i+1] < dom_min_y) dom_min_y = coords0[3*i+1];
		if (coords0[3*i+1] > dom_max_y) dom_max_y = coords0[3*i+1];
		if (coords0[3*i+2] < dom_min_z) dom_min_z = coords0[3*i+2];
		if (coords0[3*i+2] > dom_max_z) dom_max_z = coords0[3*i+2];
	}

	// Boundary-plane snap tolerance. It only has to absorb float noise on a node that
	// was built exactly on the wall plane, so it must scale with the domain -- a fixed
	// 1.0 m is larger than the WHOLE domain on a unit-scale model (HexMesh_plane.input
	// spans 0.96), where it made every node match both the min and the max test and the
	// max won: the entire interface layer collapsed onto the (x_max, y_max) corner and
	// the pillow layer built on it became 5832 zero-Jacobian pyramids.
	// min(1.0, 1e-3*L) keeps metre-scale models byte-identical to the old behaviour
	// (1e-3*L there is hundreds of metres) and is ~10% of an element edge otherwise.
	const double snap_x = std::min(1.0, 1e-3 * (dom_max_x - dom_min_x));
	const double snap_y = std::min(1.0, 1e-3 * (dom_max_y - dom_min_y));
	const double snap_z = std::min(1.0, 1e-3 * (dom_max_z - dom_min_z));

	// Two different nodes projected onto the exact same point is a collapse: the hexes
	// between them lose a face and turn into prisms/pyramids, which is what renders as a
	// hole downstream. The second claimant keeps its lattice position instead.
	struct pt3_hash {
		size_t operator()(const std::array<double,3> &p) const {
			size_t h = 1469598103934665603ULL;
			for (int i = 0; i < 3; i++) {
				size_t b; double v = p[i]; std::memcpy(&b, &v, sizeof(b));
				h ^= b; h *= 1099511628211ULL;
			}
			return h;
		}
	};
	std::unordered_map<std::array<double,3>, int, pt3_hash> claimed_by;
	mesh->proj_target.assign(coords.size(), std::numeric_limits<double>::quiet_NaN());
	mesh->proj_len.assign(coords.size() / 3, 0.0);
	int n_collapse_rejected = 0;
	int n_collapse_backoff = 0;

	// Boundary snaps and the sea-surface clamp, applied to any candidate position.
	auto constrain = [&](int node, double &x, double &y, double &z) {
		// If the node was on an exterior boundary plane (X+, X-, Y+, Y-, Z-), keep its boundary coordinate locked!
		if (std::fabs(coords0[3*node+0] - dom_min_x) < snap_x) x = dom_min_x;
		else if (std::fabs(coords0[3*node+0] - dom_max_x) < snap_x) x = dom_max_x;
		if (std::fabs(coords0[3*node+1] - dom_min_y) < snap_y) y = dom_min_y;
		else if (std::fabs(coords0[3*node+1] - dom_max_y) < snap_y) y = dom_max_y;
		if (std::fabs(coords0[3*node+2] - dom_min_z) < snap_z) z = dom_min_z;

		// If the node is on the top surface (z >= SEA_LEVEL), snap its (x, y) to the coastline
		// while preserving its top surface elevation (coords0 z) so it doesn't get dragged down.
		if (coords0[3*node+2] >= SEA_LEVEL - 1e-6) z = coords0[3*node+2];
		else if (z > sea_clamp) z = sea_clamp;
	};

	auto record = [&](int node, double x, double y, double z) {
		// Two different nodes projected onto the same point collapse the hexes between
		// them, so the point is claimed first-come. The loser used to be ABANDONED on the
		// lattice -- not moved, and (because this returned before pending.emplace) not even
		// entered into nodes_b_mat, so nothing downstream knew the interface had a hole
		// there. Measured: 29 nodes on hyeres, 18 on kashiwazaki, 2 on mauna_loa_small.
		// Instead, back the loser off along its OWN placement segment (lattice -> target)
		// until it lands on a free point: still near the surface, still a genuine interface
		// node, and distinct by a real distance rather than by an epsilon.
		static const double BACKOFF[] = {1.0, 0.9, 0.75, 0.5, 0.25};
		const double ox = coords0[3*node+0], oy = coords0[3*node+1], oz = coords0[3*node+2];
		if (std::isnan(mesh->proj_target[3*node])) {   // first write wins, as for pending
			double tx = x, ty = y, tz = z;
			constrain(node, tx, ty, tz);
			mesh->proj_target[3*node] = tx; mesh->proj_target[3*node+1] = ty; mesh->proj_target[3*node+2] = tz;
			mesh->proj_len[node] = std::sqrt((tx-ox)*(tx-ox) + (ty-oy)*(ty-oy) + (tz-oz)*(tz-oz));
		}
		for (size_t k = 0; k < sizeof(BACKOFF)/sizeof(BACKOFF[0]); k++) {
			const double a = BACKOFF[k];
			double cx = ox + a*(x-ox), cy = oy + a*(y-oy), cz = oz + a*(z-oz);
			constrain(node, cx, cy, cz);
			std::array<double, 3> target{cx, cy, cz};
			auto claim = claimed_by.emplace(target, node);
			if (claim.second || claim.first->second == node) {
				pending.emplace(node, target);
				if (k) n_collapse_backoff++;
				return;
			}
		}
		// Every backoff was taken too: leave it on the lattice, as before.
		n_collapse_rejected++;
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
	auto report_collapses = [&]() {
		if (n_collapse_backoff)
			printf("    Projection: %d nodes backed off along their own segment to avoid "
			       "collapsing onto an already-projected node\n", n_collapse_backoff);
		if (n_collapse_rejected)
			printf("    Projection: %d node moves rejected (no free point on the segment, "
			       "left on the lattice)\n", n_collapse_rejected);
	};

	// --- Projection validity limiter: setup ----------------------------------
	// (coords0 already snapshotted above, before `record` was defined.)
	// Projecting interface nodes onto a steep/curved bathy collapses or inverts
	// the thin interface-layer elements (zext median 12 m vs ~500 m wide), which
	// the pillowing then inherits as zero-volume elements that render as holes.
	// After projection we pull the offending moved nodes back toward their
	// original positions (see end of function).
	// Use the SAME node order the h5 writer emits (mgeom::H5_ORD =
	// {4,5,6,7,0,1,2,3} in hexa_h5.cpp/assign_elem_nodes). e->nodes is not the
	// standard hex order assumed by the tet decomposition and the reorder is not
	// a clean z-flip, so computing it directly over e->nodes yields a
	// geometrically wrong volume — the limiter then misses real inversions.
	auto elem_xyz = [&](octant_t *e, double X[8], double Y[8], double Z[8]) {
		for (int i = 0; i < 8; i++) {
			int id = e->nodes[mgeom::H5_ORD[i]].id;
			X[i] = coords[3*id]; Y[i] = coords[3*id+1]; Z[i] = coords[3*id+2];
		}
	};
	auto hexvol = [&](octant_t *e) -> double {
		double X[8], Y[8], Z[8];
		elem_xyz(e, X, Y, Z);
		return mgeom::hex_signed_volume(X, Y, Z);
	};
	// Shortest edge of an element, in the same H5_ORD node order elem_xyz emits.
	auto shortest_edge = [](const double X[8], const double Y[8], const double Z[8]) -> double {
		static const int E[12][2] = {
			{0,1},{1,2},{2,3},{3,0},{4,5},{5,6},{6,7},{7,4},{0,4},{1,5},{2,6},{3,7}
		};
		double lo = 1e300;
		for (int k = 0; k < 12; k++) {
			int a = E[k][0], b = E[k][1];
			double dx = X[a]-X[b], dy = Y[a]-Y[b], dz = Z[a]-Z[b];
			double l = std::sqrt(dx*dx + dy*dy + dz*dz);
			if (l < lo) lo = l;
		}
		return lo;
	};
	std::vector<double> ref_vol(mesh->elements.elem_count);
	std::vector<double> ref_edge(mesh->elements.elem_count);
	for (int iel = 0; iel < mesh->elements.elem_count; iel++) {
		octant_t *e = (octant_t*) sc_array_index(&mesh->elements, iel);
		double X[8], Y[8], Z[8];
		elem_xyz(e, X, Y, Z);
		ref_vol[iel]  = mgeom::hex_signed_volume(X, Y, Z);
		ref_edge[iel] = shortest_edge(X, Y, Z);
	}

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

			// Near the coast an octree edge crosses the surface TWICE -- once on the sea
			// floor and once on the vertical wall. Taking whichever triangle the bb-tree
			// happens to return first made the result depend on heap pointer order, i.e. it
			// changed between runs of the same binary on the same input. Evaluate every
			// candidate and keep the one nearest the edge midpoint, which is where both
			// inner nodes sit.
			GtsPoint *pt = NULL;
			{
				const double mx = 0.5*(x1+x2), my = 0.5*(y1+y2), mz = 0.5*(z1+z2);
				double best = 0.0;
				for (size_t k = 0; k < mesh->gdata_vec.size(); k++) {
					if (!mesh->gdata_vec[k].bbt) continue;
					GSList *list = gts_bb_tree_overlap(mesh->gdata_vec[k].bbt, bb);
					for (GSList *l = list; l; l = l->next) {
						GtsBBox *b = GTS_BBOX(l->data);
						GtsPoint *q = mesh->input.CgalUse
							? SegmentTriangleIntersectionCgal(seg, GTS_TRIANGLE(b->bounded))
							: SegmentTriangleIntersection(seg, GTS_TRIANGLE(b->bounded));
						if (!q) continue;
						double d = (q->x-mx)*(q->x-mx) + (q->y-my)*(q->y-my) + (q->z-mz)*(q->z-mz);
						if (!pt || d < best) { if (pt) gts_object_destroy(GTS_OBJECT(pt)); pt = q; best = d; }
						else gts_object_destroy(GTS_OBJECT(q));
					}
					if (list) g_slist_free(list);
				}
			}

			gts_object_destroy(GTS_OBJECT(bb));
			gts_object_destroy(GTS_OBJECT(seg));

			if (!pt) continue;

			// Each edge is shared between elem0 (inner node = vertex[1]) and
			// elem1 (inner node = vertex[0]); both get snapped to the same point.
			record(elem0->nodes[EdgeVerticesMap[iedge][1]].id, pt->x, pt->y, pt->z);
			record(elem1->nodes[EdgeVerticesMap[iedge][0]].id, pt->x, pt->y, pt->z);
			gts_object_destroy(GTS_OBJECT(pt));
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

			// Same rationale as Pass 1: keep the candidate nearest the diagonal's
			// midpoint, not whichever the bb-tree happens to return first.
			GtsPoint *pt = NULL;
			{
				const double mx = 0.5*(x1+x2), my = 0.5*(y1+y2), mz = 0.5*(z1+z2);
				double best = 0.0;
				for (size_t k = 0; k < mesh->gdata_vec.size(); k++) {
					if (!mesh->gdata_vec[k].bbt) continue;
					GSList *list = gts_bb_tree_overlap(mesh->gdata_vec[k].bbt, bb);
					for (GSList *l = list; l; l = l->next) {
						GtsBBox *b = GTS_BBOX(l->data);
						GtsPoint *q = mesh->input.CgalUse
							? SegmentTriangleIntersectionCgal(seg, GTS_TRIANGLE(b->bounded))
							: SegmentTriangleIntersection(seg, GTS_TRIANGLE(b->bounded));
						if (!q) continue;
						double d = (q->x-mx)*(q->x-mx) + (q->y-my)*(q->y-my) + (q->z-mz)*(q->z-mz);
						if (!pt || d < best) { if (pt) gts_object_destroy(GTS_OBJECT(pt)); pt = q; best = d; }
						else gts_object_destroy(GTS_OBJECT(q));
					}
					if (list) g_slist_free(list);
				}
			}

			gts_object_destroy(GTS_OBJECT(bb));
			gts_object_destroy(GTS_OBJECT(seg));

			if (!pt) { face_centroid_fallback(); continue; }
			record(center, pt->x, pt->y, pt->z);
			gts_object_destroy(GTS_OBJECT(pt));
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

		// Same rationale as Pass 1: keep the candidate nearest the diagonal's
		// midpoint, not whichever the bb-tree happens to return first.
		GtsPoint *pt = NULL;
		{
			const double mx = 0.5*(x1+x2), my = 0.5*(y1+y2), mz = 0.5*(z1+z2);
			double best = 0.0;
			for (size_t k = 0; k < mesh->gdata_vec.size(); k++) {
				if (!mesh->gdata_vec[k].bbt) continue;
				GSList *list = gts_bb_tree_overlap(mesh->gdata_vec[k].bbt, bb);
				for (GSList *l = list; l; l = l->next) {
					GtsBBox *b = GTS_BBOX(l->data);
					GtsPoint *q = mesh->input.CgalUse
						? SegmentTriangleIntersectionCgal(seg, GTS_TRIANGLE(b->bounded))
						: SegmentTriangleIntersection(seg, GTS_TRIANGLE(b->bounded));
					if (!q) continue;
					double d = (q->x-mx)*(q->x-mx) + (q->y-my)*(q->y-my) + (q->z-mz)*(q->z-mz);
					if (!pt || d < best) { if (pt) gts_object_destroy(GTS_OBJECT(pt)); pt = q; best = d; }
					else gts_object_destroy(GTS_OBJECT(q));
				}
				if (list) g_slist_free(list);
			}
		}

		gts_object_destroy(GTS_OBJECT(bb));
		gts_object_destroy(GTS_OBJECT(seg));

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
		gts_object_destroy(GTS_OBJECT(pt));
	}
	flush();

	report_collapses();

	RegularizeSkippedOctreeNodes(mesh, coords, nodes_b_mat);
	if(mesh->input.refinement) {
		mesh->refinement_projection_requests.clear();
		for(int ie=0;ie<int(mesh->elements.elem_count);++ie) {
			auto *e=static_cast<octant_t *>(sc_array_index(&mesh->elements,ie));
			pillowgeom::Hex p;
			for(int k=0;k<8;++k)for(int d=0;d<3;++d)p[k][d]=coords[3*e->nodes[mgeom::H5_ORD[k]].id+d];
			if(!pillowgeom::valid(p,1))mesh->refinement_projection_requests.push_back(ie);
		}
		printf("    Adaptive refinement: captured %zu raw projection failures before repair\n",mesh->refinement_projection_requests.size());
	}

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
		// Retreat step. The node is moved a fraction (1 - PULL_STEP) of the way back to its
		// lattice position each round, so after k rounds it keeps PULL_STEP^k of the
		// displacement the surface asked for. At the old 0.5 the very first retreat threw away
		// half of it, and there was no way to sit at 80%: the reachable positions were
		// 100/50/25/12.5%. A small step walks back gently and stops as soon as the element is
		// valid, at the cost of needing proportionally more iterations to reach the same depth
		// -- MAXIT is sized so that PULL_STEP^MAXIT is still below the old 0.5^8.
		const double PULL_STEP = 0.92;
		// Re-enabled. At 0 the projection was exact everywhere, including where the surface
		// crushed an octree cell into a sliver -- a 55 m edge inside a 3 km element at the
		// coastline. Those slivers tangle with their neighbour and no amount of node
		// smoothing downstream can untangle them; they were the entire residual inverted
		// count the untangler could not clear. 0.92^40 = 0.036, so a node can retreat to
		// 3.6% of the displacement the surface asked for before the loop gives up.
		const int MAXIT = 40;
		// How tightly the mesh is allowed to conform to the surface. Every node of an element
		// failing these tests is pulled back toward its lattice position, so the stricter
		// they are, the further the interface ends up from the real coastline/sea floor.
		//
		// Margin, not just sign: a barely-valid element (minSJ ~1e-3) leaves the pillow layer
		// no room and inverts as soon as a buffer node is inserted -- that is why these are
		// not simply 0. Loosened from 0.05 to let the wall follow the coastline more closely;
		// raise them back if pillowing starts producing inverted elements.
		// Min corner scaled Jacobian. The volume and edge tests below are blind to a CONCAVE
		// CORNER: a hex keeps a healthy positive volume, and every edge its full length, with
		// 1-3 of its 8 corner Jacobians negative. That is not a corner case -- measured on
		// Argostoli_ref4, all 236 elements the final VerifyMeshInversion flagged were of
		// exactly that shape (vol > 0, worst minSJ -1.6e-4), so with this test disabled the
		// limiter reported "0 invalid" while handing 236 invalid elements downstream. The
		// limiter and VerifyMeshInversion now use the same definition of valid.
		// 0.0 (validity, no margin) removes all 236 for 0.6 pp of surviving displacement and
		// an unchanged h_min; 0.05 costs 6 pp, 0.20 costs 15 pp. Raise it only if pillowing
		// starts inverting elements again -- it wants a margin, this only guarantees validity.
		// Overridable with PROJ_SJ_MIN for A/B runs (-1.0 restores the old, blind behaviour).
		double SJ_MIN = 0.0;
		if (const char *ev = getenv("PROJ_SJ_MIN")) SJ_MIN = atof(ev);
		const double VOL_MIN = 1e-6;   // min |volume| as a fraction of the reference volume
		// The volume test alone is too weak: a sliver keeps a healthy positive volume while
		// one of its edges collapses (the tangled kefalonia pair had vol +9.6e7 with a 55 m
		// edge against a 3.1 km one). Cap how much of its shortest edge the projection may
		// take from an element. This is the fidelity knob: LOWER lets the mesh hug the
		// coastline more closely and risks slivers, HIGHER keeps elements healthy and leaves
		// the interface further from the GTS. Only elements that violate it are pulled back,
		// and only until they recover.
		// With the corner-Jacobian test above back on, this no longer has to stand in for
		// shape validity -- it only guards against slivers (an element keeping a valid shape
		// while one edge is crushed). Overridable with PROJ_EDGE_MIN_FRAC for A/B runs.
		double EDGE_MIN_FRAC = 0.25;
		if (const char *ev = getenv("PROJ_EDGE_MIN_FRAC")) EDGE_MIN_FRAC = atof(ev);
		// Shear: none of the tests above sees an element the projection turned on its side. At a
		// coastline the two ends of one lattice column are placed by different cuts (top-edge node
		// on the coast wall, face-centre node on the sea floor) and can land on opposite sides of
		// the cell: hawaii ref3 had a 0.5 km tall water cell sheared 11 km sideways, valid by every
		// test above, whose pillow children set the whole mesh's dt_min (2.5e-06 s). dt follows
		// sigma_min(J) (tools/dt_report.py), so cap how much of it the projection may take, relative
		// to the same element before the projection: scale-free, any coastline, any interface.
		// Overridable with PROJ_SIGMA_MIN_FRAC (0 disables).
		double SIGMA_MIN_FRAC = 0.05;   // 0.02/0.05/0.1 on hawaii, kefalonia_small1 ref3: see HANDOFF §5
		if (const char *ev = getenv("PROJ_SIGMA_MIN_FRAC")) SIGMA_MIN_FRAC = atof(ev);
		auto sigma_of = [&](octant_t *e, const std::vector<double> &c) {
			pillowgeom::Hex h;
			for (int k = 0; k < 8; k++) for (int d = 0; d < 3; d++) h[k][d] = c[3*e->nodes[mgeom::H5_ORD[k]].id+d];
			return pillowgeom::size_metric(h, false);
		};
		std::vector<double> ref_sigma;
		if (SIGMA_MIN_FRAC > 0) {
			ref_sigma.resize(mesh->elements.elem_count);
			for (int iel = 0; iel < mesh->elements.elem_count; iel++)
				ref_sigma[iel] = sigma_of((octant_t*) sc_array_index(&mesh->elements, iel), coords0);
		}
		// Stall detector. The retreat is only monotone if ALL nodes of an element retreat
		// together; a partial retreat can invert a neighbour that was fine, so the failing
		// SET can rotate instead of shrinking and the loop burns every iteration without
		// converging. Measured on belle_ile: 40 sweeps, 8711 pull-backs, yet no node
		// retreated more than 7 times -- and the extra sweeps only spent fidelity (74%
		// surviving displacement against 78-80% on every case that converged). Stop once
		// the best count has not improved for STALL_LIMIT sweeps and let the octree
		// untangler finish the job with the displacement still intact.
		const int STALL_LIMIT = 5;
		int best_bad = INT_MAX, stall = 0;
		int it = 0, nbad = 0, npull_total = 0, n_prebad = 0;
		std::vector<char> pull(n_nodes_loc, 0);
		std::vector<int> pull_count(n_nodes_loc, 0);
		std::vector<double> keep(n_nodes_loc, 1.0);

		// Only a node this projection actually displaced can be retreated -- a node still at
		// its lattice position has nowhere to go, and marking all 8 corners of a bad element
		// dragged correctly-placed neighbours off the coastline with it. coords0 is the
		// pre-projection snapshot, so "moved" needs no extra bookkeeping.
		std::vector<char> moved(n_nodes_loc, 0);
		int n_moved = 0;
		for (int n = 0; n < n_nodes_loc; n++) {
			if (coords[3*n+0] != coords0[3*n+0] ||
			    coords[3*n+1] != coords0[3*n+1] ||
			    coords[3*n+2] != coords0[3*n+2]) { moved[n] = 1; n_moved++; }
		}
		for (it = 0; it < MAXIT; it++) {
			std::fill(pull.begin(), pull.end(), 0);
			nbad = 0; n_prebad = 0;
			for (int iel = 0; iel < mesh->elements.elem_count; iel++) {
				octant_t *e = (octant_t*) sc_array_index(&mesh->elements, iel);
				double X[8], Y[8], Z[8];
				elem_xyz(e, X, Y, Z);
				double v = mgeom::hex_signed_volume(X, Y, Z), rv = ref_vol[iel];
				double av = v < 0 ? -v : v, arv = rv < 0 ? -rv : rv;
				// Check for true geometric volume inversion:
				bool vol_ok = (v * rv > 0.0 && av >= VOL_MIN * arv);
				// ... and for an edge the projection has crushed.
				bool edge_ok = (shortest_edge(X, Y, Z) >= EDGE_MIN_FRAC * ref_edge[iel]);
				bool sj_ok = true;
				if (SJ_MIN > -1.0) {
					double sgn = (rv >= 0.0) ? 1.0 : -1.0;
					sj_ok = (mgeom::hex_min_corner_sj(X, Y, Z) * sgn >= SJ_MIN);
				}
				bool sigma_ok = SIGMA_MIN_FRAC <= 0 || sigma_of(e, coords) >= SIGMA_MIN_FRAC * ref_sigma[iel];
				if (vol_ok && edge_ok && sj_ok && sigma_ok) continue;   // still healthy
				// An element none of whose nodes this projection moved was already invalid
				// before it ran; retreating cannot repair it, and counting it would keep the
				// loop spinning to MAXIT. Leave it to the untangler.
				int n_movable = 0;
				for (int i = 0; i < 8; i++) if (moved[e->nodes[i].id]) n_movable++;
				if (n_movable == 0) { n_prebad++; continue; }
				nbad++;
				for (int i = 0; i < 8; i++)
					if (moved[e->nodes[i].id]) pull[e->nodes[i].id] = 1;
			}
			if (nbad == 0) break;
			if (nbad < best_bad) { best_bad = nbad; stall = 0; }
			else if (++stall >= STALL_LIMIT) {
				printf("    Projection limiter: stalled at %d invalid after %d iters, "
				       "handing the rest to the untangler\n", nbad, it + 1);
				break;
			}
			for (int n = 0; n < n_nodes_loc; n++) if (pull[n]) {
				for (int k = 0; k < 3; k++)
					coords[3*n+k] = PULL_STEP * coords[3*n+k] + (1.0 - PULL_STEP) * coords0[3*n+k];
				keep[n] *= PULL_STEP;
				npull_total++;
				pull_count[n]++;
			}
		}
		printf("    Projection limiter: %d invalid elements after %d iters (%d node pull-backs, "
				"%d moved nodes, %d elements already invalid before projection)\n",
				nbad, it, npull_total, n_moved, n_prebad);
		// How much of the projection survived: the fraction of the displacement the surface
		// asked for that the node still carries.
		{
			int npulled = 0, worst = 0;
			double resid = 0.0;
			for (int n = 0; n < n_nodes_loc; n++) if (pull_count[n]) {
				npulled++;
				resid += keep[n];
				if (pull_count[n] > worst) worst = pull_count[n];
			}
			if (npulled)
				printf("    Pulled-back nodes: %d, mean surviving displacement %.1f%% "
				       "(worst node %d retreats)\n", npulled, 100.0 * resid / npulled, worst);
		}

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
		// An incomplete octree carries oct->id[iel] == -1; sc_array_index would then read
		// before the elements array. Every other octree loop in this file already guards.
		if (!IsCompleteOctree(oct)) continue;
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

	// Assign node colors: verified 2-coloring of the octree's cut-edge pattern via
	// GetOctreeBipartition (moving_nodes.cpp:771), the same routine ProjectFreeNodes's own
	// Pass 2/3 already trust. The previous ad-hoc propagation here iterated corners in a
	// fixed 0..7 order and, on every non-cut edge, unconditionally overwrote the neighbour's
	// color with its own -- with no contradiction check, unlike GetOctreeBipartition's BFS
	// (which explicitly rejects an inconsistent cut pattern). It could silently produce a
	// self-inconsistent coloring for a cube-graph cut pattern GetOctreeBipartition would have
	// refused. Apply_material's octree-corner consistency check trusts these colors, so an
	// inconsistent one undermines that check regardless of how it compares mat vs. color.
	for (int ioc = 0; ioc < mesh->oct.elem_count; ioc++) {
		octree_t *oct = (octree_t*) sc_array_index(&mesh->oct, ioc);
		if (!IsCompleteOctree(oct) || !IsOctreeCutPatternRegular(oct)) continue; // matches the 3 sibling GetOctreeBipartition call sites in this file
		int B[8];
		if (!GetOctreeBipartition(oct, B)) continue;   // inconsistent cut pattern -- leave colors unset
		for (int iel = 0; iel < 8; iel++) {
			octant_t *elem = (octant_t*) sc_array_index(&mesh->elements, oct->id[iel]);
			elem->nodes[iel].color = B[iel] ? 2 : 1;
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
	static const int FaceCenterElem_local[6] = {0, 5, 0, 2, 5, 0};
	static const int FaceCenterNode_local[6] = {7, 2, 5, 7, 7, 2};

	// Default all nodes to fixed (1). The 8 outer corners of every octree block
	// will remain fixed (1). Only the 19 internal nodes of cut octrees will be
	// marked as movable (fixed = 0).
	for (int ino = 0; ino < mesh->nodes.elem_count; ino++) {
		octant_node_t *node = (octant_node_t*) sc_array_index(&mesh->nodes, ino);
		node->fixed = 1;
	}
	for (int iel = 0; iel < mesh->elements.elem_count; iel++) {
		octant_t *elem = (octant_t*) sc_array_index(&mesh->elements, iel);
		for (int ino = 0; ino < 8; ino++) {
			elem->nodes[ino].fixed = 1;
			elem->nodes[ino].color = -1;
		}
	}

	// Derive canonical cut edges from the domain bipartition B[0..7].
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
		for (int i = 1; i < 8; i++) {
			if (B[i] != B[0]) { all_same = false; break; }
		}

		if (all_same) {
			memset(oct->edge, 0, sizeof(oct->edge));
			memset(oct->face, 0, sizeof(oct->face));
			continue;
		}

		int n_cut_edges = 0;
		for (int e = 0; e < 12; e++) {
			int a = EdgeElemOctMap[e][0];
			int b = EdgeElemOctMap[e][1];
			oct->edge[e] = (B[a] != B[b]);
			if (oct->edge[e]) n_cut_edges++;
		}

		// TOPOLOGICAL SAFETY SYSTEM:
		// A continuous domain-separating surface cutting a 3D hex cell cannot
		// intersect only 1 or 2 edges (which represents numerical noise or vertex tangency).
		// Zero out cut flags for degenerate (1-2 edge) cuts.
		if (n_cut_edges > 0 && n_cut_edges <= 2) {
			memset(oct->edge, 0, sizeof(oct->edge));
			memset(oct->face, 0, sizeof(oct->face));
			continue;
		}

		// Classify active faces
		for (int f = 0; f < 6; f++) {
			oct->face[f] = false;
			for (int k = 0; k < 4; k++) {
				if (oct->edge[FaceEdgesMap[f][k]]) { oct->face[f] = true; break; }
			}
		}

		// Mark internal movable nodes (fixed = 0) for valid cut octrees:
		// a) 12 edge mid-nodes on cut edges
		for (int e = 0; e < 12; e++) {
			if (oct->edge[e]) {
				int a = EdgeElemOctMap[e][0];
				int b = EdgeElemOctMap[e][1];
				int edge_node_id = elems[a]->nodes[b].id;
				if (edge_node_id >= 0 && edge_node_id < mesh->nodes.elem_count) {
					octant_node_t *node = (octant_node_t*) sc_array_index(&mesh->nodes, edge_node_id);
					node->fixed = 0;
				}
			}
		}

		// b) 6 face center nodes on active faces
		for (int f = 0; f < 6; f++) {
			if (oct->face[f]) {
				int el = FaceCenterElem_local[f];
				int nd = FaceCenterNode_local[f];
				int face_node_id = elems[el]->nodes[nd].id;
				if (face_node_id >= 0 && face_node_id < mesh->nodes.elem_count) {
					octant_node_t *node = (octant_node_t*) sc_array_index(&mesh->nodes, face_node_id);
					node->fixed = 0;
				}
			}
		}

		// c) 1 octree center node
		if (n_cut_edges >= 3) {
			int center_node_id = elems[0]->nodes[6].id;
			if (center_node_id >= 0 && center_node_id < mesh->nodes.elem_count) {
				octant_node_t *node = (octant_node_t*) sc_array_index(&mesh->nodes, center_node_id);
				node->fixed = 0;
			}
		}
	}

	// Sync fixed status to element node local copies
	for (int iel = 0; iel < mesh->elements.elem_count; iel++) {
		octant_t *elem = (octant_t*) sc_array_index(&mesh->elements, iel);
		for (int ino = 0; ino < 8; ino++) {
			int nid = elem->nodes[ino].id;
			if (nid >= 0 && nid < mesh->nodes.elem_count) {
				octant_node_t *gnode = (octant_node_t*) sc_array_index(&mesh->nodes, nid);
				elem->nodes[ino].fixed = gnode->fixed;
			}
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

// ---------------------------------------------------------------------------
// Inversion map (diagnostic). For every inverted element, record WHAT KIND of
// element it is -- parent-octree cut pattern, position in the octree, node
// colors, how far THIS stage moved its nodes, warp column state -- then print
// cross-tabs and dump one CSV row per element so two runs can be diffed.
// Written for the lattice-warp on/off A/B; pure measurement, changes nothing.
// ---------------------------------------------------------------------------

static int count_neg_corners(const double X[8], const double Y[8], const double Z[8], int ref)
{
	int n = 0;
	for (int k = 0; k < 8; k++) {
		int i0=k, i1=mgeom::CORNER_NB[k][0], i2=mgeom::CORNER_NB[k][1], i3=mgeom::CORNER_NB[k][2];
		double ax=X[i1]-X[i0], ay=Y[i1]-Y[i0], az=Z[i1]-Z[i0];
		double bx=X[i2]-X[i0], by=Y[i2]-Y[i0], bz=Z[i2]-Z[i0];
		double cx=X[i3]-X[i0], cy=Y[i3]-Y[i0], cz=Z[i3]-Z[i0];
		double det = ax*(by*cz-bz*cy) - ay*(bx*cz-bz*cx) + az*(bx*cy-by*cx);
		if (det * ref <= 0.0) n++;
	}
	return n;
}

void DumpInversionMap(hexa_tree_t *mesh, const std::vector<double> &coords,
                             const std::vector<double> &prev, const char *tag)
{
	const int ne = (int) mesh->elements.elem_count;
	if (ne == 0) return;
	const bool have_prev = (prev.size() == coords.size());

	// --- element -> parent octree ------------------------------------------
	std::vector<signed char> oct_pos(ne, -1), oct_cut(ne, 0), oct_nedge(ne, -1), oct_nface(ne, -1);
	for (int ioc = 0; ioc < (int) mesh->oct.elem_count; ioc++) {
		octree_t *oct = (octree_t*) sc_array_index(&mesh->oct, ioc);
		int nce = 0, ncf = 0;
		for (int k = 0; k < 12; k++) if (oct->edge[k]) nce++;
		for (int k = 0; k < 6;  k++) if (oct->face[k]) ncf++;
		for (int p = 0; p < 8; p++) {
			int id = (int) oct->id[p];
			if (id < 0 || id >= ne) continue;
			oct_pos[id]   = (signed char) p;
			oct_cut[id]   = oct->cut ? 1 : 0;
			oct_nedge[id] = (signed char) nce;
			oct_nface[id] = (signed char) ncf;
		}
	}

	// --- reference sign (same rule as analyze_mesh) -------------------------
	std::vector<double> vols(ne), sjs(ne);
	double vol_sum = 0.0;
	for (int iel = 0; iel < ne; iel++) {
		double X[8],Y[8],Z[8];
		load_elem_xyz(mesh, coords, iel, X, Y, Z);
		vols[iel] = mgeom::hex_signed_volume(X,Y,Z);
		sjs[iel]  = mgeom::hex_min_corner_sj(X,Y,Z);
		vol_sum  += vols[iel];
	}
	const int ref = mgeom::reference_sign(vol_sum);

	char fname[256];
	snprintf(fname, sizeof(fname), "invmap_%s.csv", tag);
	FILE *f = fopen(fname, "w");
	if (f) fprintf(f, "elem,nneg,minsj,vol,cut,octpos,noctedge,noctface,intercepted,"
	                  "nfree,ncolor1,ncolor2,nmoved,maxdisp,maxdxy,maxdz,"
	                  "zlayer,dzlat,level,dxlat,dylat,warpstate,warpdisp,relu,relu_over_h,cx,cy,cz\n");

	// --- histograms ---------------------------------------------------------
	int h_nneg[9]      = {0};
	int h_octpos[9]    = {0};        // 8 = not in any complete octree
	int h_noctedge[14] = {0};        // 13 = no octree
	int h_cut[3]       = {0};        // not-in-octree / uncut octree / cut octree
	int h_inter[2]     = {0};
	int h_nfree[9]     = {0};
	int h_nmoved[9]    = {0};
	int h_warp[5]      = {0};        // 4 = unknown column
	int h_zlayer[64]   = {0};
	int h_colormix[3]  = {0};        // 0 = uniform color, 1 = mixed, 2 = unset(-1) present
	int n_inv = 0;
	int zlayer_max = 0;
	for (int i = 0; i < (int) mesh->nodes.elem_count; i++) {
		octant_node_t *n = (octant_node_t*) sc_array_index(&mesh->nodes, i);
		if (n->z > zlayer_max) zlayer_max = n->z;
	}
	const int nx = mesh->ncellx + 1;

	// (dxlat, dylat, dzlat, cut) population, so the map reports RATES not just counts
	std::map<std::string, std::pair<int,int> > pop;   // class -> (total, inverted)

	for (int iel = 0; iel < ne; iel++) {
		{
			octant_t *e0 = (octant_t*) sc_array_index(&mesh->elements, iel);
			int ax=1<<30, bx=-(1<<30), ay=1<<30, by=-(1<<30), az=1<<30, bz=-(1<<30);
			for (int k = 0; k < 8; k++) {
				octant_node_t *nd = &e0->nodes[k];
				if (nd->x < ax) ax = nd->x; if (nd->x > bx) bx = nd->x;
				if (nd->y < ay) ay = nd->y; if (nd->y > by) by = nd->y;
				if (nd->z < az) az = nd->z; if (nd->z > bz) bz = nd->z;
			}
			char key[64];
			snprintf(key, sizeof(key), "dx%d dy%d dz%d %s", bx-ax, by-ay, bz-az,
			         oct_pos[iel] < 0 ? "noOct" : (oct_cut[iel] ? "CUT  " : "uncut"));
			std::pair<int,int> &pv = pop[key];
			pv.first++;
			if (mgeom::is_inverted(vols[iel], sjs[iel], ref)) pv.second++;
		}
		if (!mgeom::is_inverted(vols[iel], sjs[iel], ref)) continue;
		n_inv++;
		octant_t *e = (octant_t*) sc_array_index(&mesh->elements, iel);

		double X[8],Y[8],Z[8];
		load_elem_xyz(mesh, coords, iel, X, Y, Z);
		int nneg = count_neg_corners(X, Y, Z, ref);

		int nfree = 0, nc1 = 0, nc2 = 0, ncunset = 0, nmoved = 0, zl = 0;
		double maxdisp = 0.0, maxdxy = 0.0, maxdz = 0.0;
		double cx = 0.0, cy = 0.0, cz = 0.0;
		int wstate = 4; double wdisp = 0.0;
		int lxmin=1<<30, lxmax=-(1<<30), lymin=1<<30, lymax=-(1<<30), lzmin=1<<30, lzmax=-(1<<30);
		double cux[8], cuy[8]; int ncu = 0;
		for (int k = 0; k < 8; k++) {
			octant_node_t *nd = &e->nodes[k];
			int id = nd->id;
			if (nd->fixed == 0) nfree++;
			if (nd->color == 1) nc1++; else if (nd->color == 2) nc2++; else ncunset++;
			if (nd->z > zl) zl = nd->z;
			cx += coords[3*id+0]/8.0; cy += coords[3*id+1]/8.0; cz += coords[3*id+2]/8.0;
			if (have_prev) {
				double dx = coords[3*id+0]-prev[3*id+0];
				double dy = coords[3*id+1]-prev[3*id+1];
				double dz = coords[3*id+2]-prev[3*id+2];
				double d = std::sqrt(dx*dx+dy*dy+dz*dz);
				if (d > 1e-9) nmoved++;
				if (d > maxdisp) maxdisp = d;
				double dxy = std::sqrt(dx*dx+dy*dy);
				if (dxy > maxdxy) maxdxy = dxy;
				if (std::fabs(dz) > maxdz) maxdz = std::fabs(dz);
			}
			if (nd->x < lxmin) lxmin = nd->x; if (nd->x > lxmax) lxmax = nd->x;
			if (nd->y < lymin) lymin = nd->y; if (nd->y > lymax) lymax = nd->y;
			if (nd->z < lzmin) lzmin = nd->z; if (nd->z > lzmax) lzmax = nd->z;
			if (!g_warp_col_state.empty() && nd->x >= 0 && nd->y >= 0) {
				size_t c = (size_t) nd->y * nx + nd->x;
				if (c < g_warp_col_state.size()) {
					if (wstate == 4 || g_warp_col_state[c] > wstate) wstate = g_warp_col_state[c];
					if (g_warp_col_disp[c] > wdisp) wdisp = g_warp_col_disp[c];
					cux[ncu] = g_warp_col_ux[c]; cuy[ncu] = g_warp_col_uy[c]; ncu++;
				}
			}
		}
		// widest disagreement between the 8 corner columns: what actually shears the hex
		double relu = 0.0;
		for (int a = 0; a < ncu; a++) for (int b = a+1; b < ncu; b++) {
			double dx = cux[a]-cux[b], dy = cuy[a]-cuy[b];
			double d = std::sqrt(dx*dx+dy*dy);
			if (d > relu) relu = d;
		}
		double hcell = std::sqrt((double)(lxmax-lxmin)*(lxmax-lxmin)*g_warp_hx*g_warp_hx
		                       + (double)(lymax-lymin)*(lymax-lymin)*g_warp_hy*g_warp_hy);
		double relu_over_h = (hcell > 1e-9) ? relu/hcell : 0.0;
		int inter = (e->pad == -1) ? 1 : 0;

		h_nneg[nneg]++;
		h_octpos[oct_pos[iel] < 0 ? 8 : oct_pos[iel]]++;
		h_noctedge[oct_nedge[iel] < 0 ? 13 : oct_nedge[iel]]++;
		h_cut[oct_pos[iel] < 0 ? 0 : (oct_cut[iel] ? 2 : 1)]++;
		h_inter[inter]++;
		h_nfree[nfree]++;
		h_nmoved[nmoved]++;
		h_warp[wstate > 4 ? 4 : wstate]++;
		if (zl < 64) h_zlayer[zl]++;
		h_colormix[ncunset ? 2 : ((nc1 && nc2) ? 1 : 0)]++;

		if (f) fprintf(f, "%d,%d,%.6e,%.6e,%d,%d,%d,%d,%d,%d,%d,%d,%d,%.4f,%.4f,%.4f,%d,%d,%d,%d,%d,%d,%.4f,%.4f,%.4f,%.2f,%.2f,%.2f\n",
		               iel, nneg, sjs[iel]*ref, vols[iel]*ref,
		               (int) (oct_pos[iel] < 0 ? -1 : oct_cut[iel]), (int) oct_pos[iel],
		               (int) oct_nedge[iel], (int) oct_nface[iel], inter,
		               nfree, nc1, nc2, nmoved, maxdisp, maxdxy, maxdz,
		               zl, lzmax-lzmin, (int) e->level, lxmax-lxmin, lymax-lymin,
		               wstate, wdisp, relu, relu_over_h, cx, cy, cz);
	}
	if (f) fclose(f);

	printf(" =========================================================\n");
	printf("   INVERSION MAP [%s]  (ref sign %+d, %d inverted of %d)\n", tag, ref, n_inv, ne);
	printf(" =========================================================\n");
	{
		static const int E[12][2] = {
			{0,1},{1,2},{2,3},{3,0},{4,5},{5,6},{6,7},{7,4},{0,4},{1,5},{2,6},{3,7}
		};
		int b20 = 0, b50 = 0, b100 = 0;
		for (int iel = 0; iel < ne; iel++) {
			double X[8],Y[8],Z[8];
			load_elem_xyz(mesh, coords, iel, X, Y, Z);
			double lo = 1e300;
			for (int k = 0; k < 12; k++) {
				int a = E[k][0], c = E[k][1];
				double dx=X[a]-X[c], dy=Y[a]-Y[c], dz=Z[a]-Z[c];
				double l = std::sqrt(dx*dx+dy*dy+dz*dz);
				if (l < lo) lo = l;
			}
			if (lo < 20.0)  b20++;
			if (lo < 50.0)  b50++;
			if (lo < 100.0) b100++;
		}
		printf("    thin elements (shortest edge): <20 m: %d  <50 m: %d  <100 m: %d\n", b20, b50, b100);
	}
	if (n_inv == 0) { printf("    nothing inverted\n =========================================================\n\n"); return; }
	printf("    negative corners (of 8):   ");
	for (int k = 1; k <= 8; k++) if (h_nneg[k]) printf("%d/8:%d  ", k, h_nneg[k]);
	printf("\n    parent octree:             none:%d  uncut:%d  CUT:%d\n", h_cut[0], h_cut[1], h_cut[2]);
	printf("    position in octree:        ");
	for (int k = 0; k < 8; k++) if (h_octpos[k]) printf("p%d:%d  ", k, h_octpos[k]);
	if (h_octpos[8]) printf("none:%d", h_octpos[8]);
	printf("\n    cut edges of the octree:   ");
	for (int k = 0; k <= 12; k++) if (h_noctedge[k]) printf("%d:%d  ", k, h_noctedge[k]);
	if (h_noctedge[13]) printf("noOct:%d", h_noctedge[13]);
	printf("\n    intercepted (pad==-1):     yes:%d  no:%d\n", h_inter[1], h_inter[0]);
	printf("    corner colors:             uniform:%d  mixed(1&2):%d  unset:%d\n",
	       h_colormix[0], h_colormix[1], h_colormix[2]);
	printf("    free nodes (fixed==0):     ");
	for (int k = 0; k <= 8; k++) if (h_nfree[k]) printf("%d:%d  ", k, h_nfree[k]);
	printf("\n    nodes moved by this stage: ");
	for (int k = 0; k <= 8; k++) if (h_nmoved[k]) printf("%d:%d  ", k, h_nmoved[k]);
	printf("\n    warp column state:         untouched:%d  diffused:%d  anchored:%d  clamped:%d  n/a:%d\n",
	       h_warp[0], h_warp[1], h_warp[2], h_warp[3], h_warp[4]);
	printf("    z layer (0 = top/surface): ");
	for (int k = 0; k < 64; k++) if (h_zlayer[k]) printf("%d:%d  ", k, h_zlayer[k]);
	printf("\n    (zlayer_max = %d)   CSV: %s\n", zlayer_max, fname);
	printf("    population by lattice extent x octree class (inverted / total, rate):\n");
	for (std::map<std::string, std::pair<int,int> >::const_iterator it = pop.begin(); it != pop.end(); ++it)
		if (it->second.second)
			printf("      %-22s %6d / %8d   %7.3f %%\n", it->first.c_str(),
			       it->second.second, it->second.first,
			       100.0 * it->second.second / it->second.first);
	printf("    (classes with 0 inverted omitted; full population:)\n");
	for (std::map<std::string, std::pair<int,int> >::const_iterator it = pop.begin(); it != pop.end(); ++it)
		printf("      %-22s %8d\n", it->first.c_str(), it->second.first);
	printf(" =========================================================\n\n");
}

// ---------------------------------------------------------------------------
// Octree-interior untangling.
//
// ProjectFreeNodes' pull-back can only slide a node back along the single line
// from its lattice position to the projected point. A PARTIAL retreat on one
// node can invert a neighbour that was fine, so the loop oscillates instead of
// converging -- measured on belle_ile: 3 elements still bad after all 40
// sweeps, 8711 pull-backs, yet no node retreated more than 7 times, i.e. the
// failing set kept changing rather than being ground down.
//
// The 3x3x3 node grid of an octree is 8 corners plus 19 interior nodes (12 edge
// midpoints, 6 face centres, 1 body centre), and those 19 are exactly what the
// projection moves -- elem->nodes[iel] at octree position iel is the octant's
// OUTER corner, the other 7 are interior. Pinning the corners and freeing the
// interior gives every hex 7 movable nodes of 8, with the lattice configuration
// guaranteed to lie inside that set, so a monotone relaxer cannot get stuck the
// way the 1-D pull-back does. Interface nodes keep their gts_surface_id, so the
// untangler slides them ALONG the surface: validity at no fidelity cost.
//
// The corners stay pinned on purpose: a corner is shared by up to 8 blocks,
// most of them not cut at all, so moving one spreads deformation into regions
// that were already fine.
// ---------------------------------------------------------------------------
static void UntangleOctreeInterior(hexa_tree_t *mesh, std::vector<double> &coords,
                                   std::vector<int> &nodes_b_mat)
{
	const int nn = (int)(coords.size() / 3);
	if (nn == 0 || mesh->oct.elem_count == 0) return;

	MeshAnalysis a0 = analyze_mesh(mesh, coords);
	if (a0.n_inverted == 0) {
		printf("    Octree untangler: nothing inverted after projection, skipped\n");
		return;
	}

	// Corner vs interior, straight from the octree layout.
	std::vector<char> is_corner(nn, 0), is_interior(nn, 0);
	for (int ioc = 0; ioc < (int) mesh->oct.elem_count; ioc++) {
		octree_t *oct = (octree_t*) sc_array_index(&mesh->oct, ioc);
		if (!IsCompleteOctree(oct)) continue;
		for (int iel = 0; iel < 8; iel++) {
			octant_t *e = (octant_t*) sc_array_index(&mesh->elements, oct->id[iel]);
			for (int ino = 0; ino < 8; ino++) {
				int id = e->nodes[ino].id;
				if (id < 0 || id >= nn) continue;
				if (ino == iel) is_corner[id] = 1; else is_interior[id] = 1;
			}
		}
	}

	std::vector<uint8_t> wall_lock;
	std::vector<NodeConstraint> cons =
		classify_node_constraints(mesh, coords, nodes_b_mat, &wall_lock);

	// Everything that is not an octree interior node is frozen. Interior nodes keep
	// whatever wall lock classify_node_constraints gave them -- a node on a domain
	// wall must not be pushed off it.
	const uint8_t FULL = mgeom::LOCK_X | mgeom::LOCK_Y | mgeom::LOCK_Z;
	int n_free = 0, n_free_surf = 0;
	for (int i = 0; i < nn; i++) {
		if (is_interior[i] && !is_corner[i]) {
			n_free++;
			if (cons[i].gts_surface_id >= 0) n_free_surf++;
		} else {
			cons[i].lock_mask = FULL;
		}
	}

	printf("    Octree untangler: %d inverted in, %d free interior nodes of %d "
	       "(%d of them on a GTS surface)\n", a0.n_inverted, n_free, nn, n_free_surf);

	int remaining = untangle_inversions(mesh, coords, cons, wall_lock, a0.reference_sign);

	MeshAnalysis a1 = analyze_mesh(mesh, coords);
	printf("    Octree untangler: %d inverted out (untangler reported %d), "
	       "h_min %.6e -> %.6e\n", a1.n_inverted, remaining, a0.h_min, a1.h_min);
}

// Projection/untangling used to move shared vertices independently on each
// rank. Agree on the lowest-rank copy before certifying the mesh, and agree on
// the retreat target too. The subsequent collective retreat then preserves
// identical positions on every side of a partition boundary.
static void SynchronizeProjectedNodes(hexa_tree_t *mesh, std::vector<double> &coords,
                                      std::vector<double> &reference) {
	struct SharedPosition {
		int xyz[3], rank, id;
		double position[3], target[3];
	};
	std::array<int,2> low={INT_MAX,INT_MAX}, high={INT_MIN,INT_MIN};
	for (size_t id=0; id<mesh->nodes.elem_count; ++id) {
		auto *n=static_cast<octant_node_t *>(sc_array_index(&mesh->nodes,id));
		low[0]=std::min(low[0],n->x); high[0]=std::max(high[0],n->x);
		low[1]=std::min(low[1],n->y); high[1]=std::max(high[1],n->y);
	}
	std::vector<SharedPosition> send;
	for (size_t id=0; id<mesh->nodes.elem_count; ++id) {
		auto *n=static_cast<octant_node_t *>(sc_array_index(&mesh->nodes,id));
		if (n->x!=low[0] && n->x!=high[0] && n->y!=low[1] && n->y!=high[1]) continue;
		SharedPosition p{};
		p.xyz[0]=n->x; p.xyz[1]=n->y; p.xyz[2]=n->z;
		p.rank=mesh->mpi_rank; p.id=static_cast<int>(id);
		for (int d=0; d<3; ++d) { p.position[d]=coords[3*id+d]; p.target[d]=reference[3*id+d]; }
		send.push_back(p);
	}
	int bytes=static_cast<int>(send.size()*sizeof(SharedPosition)), total=0;
	std::vector<int> sizes(mesh->mpi_size), offsets(mesh->mpi_size);
	MPI_Allgather(&bytes,1,MPI_INT,sizes.data(),1,MPI_INT,MPI_COMM_WORLD);
	for (int r=0; r<mesh->mpi_size; ++r) { offsets[r]=total; total+=sizes[r]; }
	std::vector<SharedPosition> positions(total/sizeof(SharedPosition));
	MPI_Allgatherv(send.data(),bytes,MPI_BYTE,positions.data(),sizes.data(),offsets.data(),MPI_BYTE,MPI_COMM_WORLD);
	std::sort(positions.begin(),positions.end(),[](const SharedPosition &a,const SharedPosition &b) {
		for (int d=0; d<3; ++d) if (a.xyz[d]!=b.xyz[d]) return a.xyz[d]<b.xyz[d];
		return a.rank<b.rank;
	});
	int changed=0;
	double maximum=0;
	for (size_t i=0; i<positions.size();) {
		size_t j=i+1;
		while (j<positions.size() && std::equal(positions[i].xyz,positions[i].xyz+3,positions[j].xyz)) ++j;
		for (size_t k=i; k<j; ++k) if (positions[k].rank==mesh->mpi_rank) {
			int id=positions[k].id;
			double gap=0;
			for (int d=0; d<3; ++d) {
				double delta=coords[3*id+d]-positions[i].position[d]; gap+=delta*delta;
				coords[3*id+d]=positions[i].position[d];
				reference[3*id+d]=positions[i].target[d];
			}
			if (gap>0) ++changed;
			maximum=std::max(maximum,std::sqrt(gap));
		}
		i=j;
	}
	printf("    Projection MPI agreement: rank %d, %d shared nodes synchronized (maximum gap %.6g m)\n",
	       mesh->mpi_rank,changed,maximum);
}

// The corner-only untangler can leave a Jacobian sign change inside a hex.
// Retreat only nodes incident to those cells toward the saved pre-projection
// mesh, then recheck every affected neighbour with the full-volume certificate.
// This uses the same projection retreat as ProjectFreeNodes; material labels
// have not been assigned yet. Pillow insertion must receive a valid host mesh.
static void ValidateProjectedJacobians(hexa_tree_t *mesh, std::vector<double> &coords,
                                      const std::vector<double> &reference_input) {
	std::vector<double> synchronized_reference;
	if (mesh->mpi_size>1) {
		synchronized_reference=reference_input;
		SynchronizeProjectedNodes(mesh,coords,synchronized_reference);
	}
	const auto &reference=mesh->mpi_size>1 ? synchronized_reference : reference_input;
	const int ne=mesh->elements.elem_count, nn=mesh->nodes.elem_count;
	std::vector<unsigned char> bad(ne,0);
	std::vector<std::vector<int>> incidence(nn);
	auto valid=[&](int ie) {
		auto *e=static_cast<octant_t *>(sc_array_index(&mesh->elements,ie));
		pillowgeom::Hex p;
		for (int k=0; k<8; ++k) for (int d=0; d<3; ++d)
			p[k][d]=coords[3*e->nodes[mgeom::H5_ORD[k]].id+d];
		return pillowgeom::valid(p,1);
	};
	int count=0;
	for (int ie=0; ie<ne; ++ie) {
		bad[ie]=!valid(ie); count+=bad[ie];
		auto *e=static_cast<octant_t *>(sc_array_index(&mesh->elements,ie));
		for (auto &n : e->nodes) incidence[n.id].push_back(ie);
	}
	const int initial=count;
	std::map<std::array<int,3>,int> lattice;
	if (mesh->mpi_size>1) for (int id=0; id<nn; ++id) {
		auto *n=static_cast<octant_node_t *>(sc_array_index(&mesh->nodes,id));
		lattice[{n->x,n->y,n->z}]=id;
	}
	int retreats=0, iteration=0;
	for (; iteration<80; ++iteration) {
		int global_count=count;
		if (mesh->mpi_size>1) MPI_Allreduce(&count,&global_count,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
		if (!global_count) break;
		std::set<int> move, dirty;
		for (int ie=0; ie<ne; ++ie) if (bad[ie]) {
			auto *e=static_cast<octant_t *>(sc_array_index(&mesh->elements,ie));
			for (auto &n : e->nodes)
				if (!std::equal(coords.begin()+3*n.id,coords.begin()+3*n.id+3,reference.begin()+3*n.id))
					move.insert(n.id);
		}
		// A shared projected node retreats on both ranks in the same iteration.
		if (mesh->mpi_size>1) {
			std::vector<int> send, sizes(mesh->mpi_size), offsets(mesh->mpi_size);
			for (int id : move) {
				auto *n=static_cast<octant_node_t *>(sc_array_index(&mesh->nodes,id));
				send.insert(send.end(),{n->x,n->y,n->z});
			}
			int sz=send.size(), total=0;
			MPI_Allgather(&sz,1,MPI_INT,sizes.data(),1,MPI_INT,MPI_COMM_WORLD);
			for (int r=0; r<mesh->mpi_size; ++r) { offsets[r]=total; total+=sizes[r]; }
			std::vector<int> received(total);
			MPI_Allgatherv(send.data(),sz,MPI_INT,received.data(),sizes.data(),offsets.data(),MPI_INT,MPI_COMM_WORLD);
			for (int i=0; i<total; i+=3) {
				auto it=lattice.find({received[i],received[i+1],received[i+2]});
				if (it!=lattice.end()) move.insert(it->second);
			}
		}
		for (int id : move) {
			for (int d=0; d<3; ++d) coords[3*id+d]=0.8*coords[3*id+d]+0.2*reference[3*id+d];
			dirty.insert(incidence[id].begin(),incidence[id].end());
			++retreats;
		}
		for (int ie : dirty) {
			count-=bad[ie]; bad[ie]=!valid(ie); count+=bad[ie];
		}
	}
	printf("    Full-Jacobian projection check: %d invalid before, %d after (%d local node retreats, %d rounds)\n",
	       initial,count,retreats,iteration);
}

void MovingNodes(hexa_tree_t* mesh, std::vector<double>& coords, std::vector<int>& nodes_b_mat) {

	bool deb = false;

	printf("    Building the octree structure...\n");
	DoOctree(mesh);

	printf("    Classifying octree corners...\n");
	ClassifyOctreeCorners(mesh, coords);

	printf("    Identifying the movable nodes...\n");
	IdentifyMovableNodes(mesh);

	// if(deb){
	// 	char fdname[80];
	// 	sprintf(fdname,"free_node_%04d.txt", mesh->mpi_rank);
	// 	FILE* fnode = fopen(fdname,"w");
	// 	for(int ino = 0; ino < mesh->nodes.elem_count; ino ++){
	// 		octant_node_t* node = (octant_node_t*) sc_array_index (&mesh->nodes, ino);
	// 		if(node->fixed==0){
	// 			int nnode = 3*node->id;
	// 			fprintf(fnode,"%f %f %f\n",coords[nnode+0],coords[nnode+1],coords[nnode+2]);
	// 		}
	// 	}
	// 	fclose(fnode);
	// }

	// if(deb){
	// 	char fdname[80];
	// 	sprintf(fdname,"debug_edges_%04d.txt", mesh->mpi_rank);
	// 	FILE* dedges = fopen(fdname,"w");
	// 	for(int ioc = 0; ioc < mesh->oct.elem_count; ioc ++){
	// 		octree_t* oct = (octree_t*) sc_array_index (&mesh->oct, ioc);
	// 		for(int iel = 0; iel < 8; iel++){
	// 			if(oct->id[iel] !=-1){
	// 				octant_t* elem = (octant_t*) sc_array_index (&mesh->nodes, oct->id[iel]);
	// 				fprintf(dedges,"El %d\n", oct->id[iel]);
	// 				for(int iedge = 0; iedge < 12;iedge++){
	// 					fprintf(dedges,"%s ", elem->edge[iedge].ref ? "T" : "F");
	// 				}
	// 				fprintf(dedges,"\n");
	// 			}
	// 		}
	// 	}
	// 	fclose(dedges);
	// }

	printf("    Make the projection of the nodes into the surface...\n");
	nodes_b_mat.clear();
	// LATTICE_WARP=1 enables WarpLatticeToCoastline (off by default) so the same
	// binary can be run as the warp-on / warp-off A/B.
	const bool use_warp = (getenv("LATTICE_WARP") != NULL && atoi(getenv("LATTICE_WARP")) != 0);
	printf("    MovingNodes: lattice warp %s\n", use_warp ? "ON" : "OFF");

	std::vector<double> coords_in = coords;
	printf("    MovingNodes inversion check: before lattice warp\n");
	VerifyMeshInversion(mesh, &coords);
	DumpInversionMap(mesh, coords, coords_in, use_warp ? "warpON_0pre" : "warpOFF_0pre");

	if (use_warp) WarpLatticeToCoastline(mesh, coords);
	printf("    MovingNodes inversion check: after lattice warp\n");
	VerifyMeshInversion(mesh, &coords);
	DumpInversionMap(mesh, coords, coords_in, use_warp ? "warpON_1warp" : "warpOFF_1warp");

	std::vector<double> coords_prewarp = coords;
	ProjectFreeNodes(mesh,coords,nodes_b_mat);
	printf("    MovingNodes inversion check: after surface projection\n");
	VerifyMeshInversion(mesh, &coords);
	DumpInversionMap(mesh, coords, coords_prewarp, use_warp ? "warpON_2proj" : "warpOFF_2proj");

	// OCT_UNTANGLE=0 skips it, for A/B against the projection-only mesh.
	if (!(getenv("OCT_UNTANGLE") && atoi(getenv("OCT_UNTANGLE")) == 0)) {
		UntangleOctreeInterior(mesh, coords, nodes_b_mat);
		printf("    MovingNodes inversion check: after octree untangling\n");
		VerifyMeshInversion(mesh, &coords);
	}
	ValidateProjectedJacobians(mesh,coords,coords_prewarp);
	DumpInversionMap(mesh, coords, coords_in,      use_warp ? "warpON_3total" : "warpOFF_3total");

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

// ---------------------------------------------------------------------------------------------
// pillowFirst (log/plans/2026-10-12_pillow-before-projection.md). The mesh was pillowed while it
// was still the unmoved staircase (every host a near-box, half-cell shells on the interface
// faces); now the staircase interface is glued to its GTS, the shells taking up the deformation:
//  a. each interface node goes to the closest point of its GTS (a vertical ray alone would send
//     both rows of a stair riser to the same z and collapse the riser face);
//  b. a few sweeps of tangential smoothing on the GTS (neighbour average, re-projected) spread the
//     riser rows along the slope, so the staircase quads become regular quads on the surface;
//  c. the core of every touched host is rebuilt from the host's new geometry and insets, so the
//     shells keep their proportions;
//  d. limiter: an element incident to a moved node that is not valid over its whole volume, or
//     lost more than PROJ_SIGMA_MIN_FRAC of its sigma_min, halves the move of its interface nodes
//     (repeated; at the end they go back). The intended point is kept in proj_target (meshOpt).
// Fixed: nodes on the domain boundary (walls, bottom, topography / sea lid -> the coastline) and
// nodes whose interface faces involve more than two materials (triple junctions). Only nodes within
// 1.5 local edges of a GTS move: material boundaries that are not on a GTS (e.g. the zcuts) stay.
// Env for A/B: PF_NOPROJ=1 skips all of it, PF_LIMITER=0 skips d, PF_SMOOTH=<sweeps> (default 5).
// closest point of triangle abc to p (Ericson, Real-Time Collision Detection 5.1.5)
static void closest_on_triangle(const double p[3], const double a[3], const double b[3], const double c[3], double q[3]) {
	auto dot = [](const double *u, const double *v) { return u[0]*v[0] + u[1]*v[1] + u[2]*v[2]; };
	double ab[3], ac[3], ap[3], bp[3], cp[3];
	for (int d = 0; d < 3; d++) { ab[d] = b[d]-a[d]; ac[d] = c[d]-a[d]; ap[d] = p[d]-a[d]; bp[d] = p[d]-b[d]; cp[d] = p[d]-c[d]; }
	auto set = [&](const double *o, double s, const double *u, double t, const double *v) { for (int d = 0; d < 3; d++) q[d] = o[d] + s*u[d] + t*v[d]; };
	const double d1 = dot(ab, ap), d2 = dot(ac, ap);
	if (d1 <= 0 && d2 <= 0) { set(a, 0, ab, 0, ac); return; }
	const double d3 = dot(ab, bp), d4 = dot(ac, bp);
	if (d3 >= 0 && d4 <= d3) { set(b, 0, ab, 0, ac); return; }
	const double vc = d1*d4 - d3*d2;
	if (vc <= 0 && d1 >= 0 && d3 <= 0) { set(a, d1/(d1-d3), ab, 0, ac); return; }
	const double d5 = dot(ab, cp), d6 = dot(ac, cp);
	if (d6 >= 0 && d5 <= d6) { set(c, 0, ab, 0, ac); return; }
	const double vb = d5*d2 - d1*d6;
	if (vb <= 0 && d2 >= 0 && d6 <= 0) { set(a, 0, ab, d2/(d2-d6), ac); return; }
	const double va = d3*d6 - d5*d4;
	if (va <= 0 && (d4-d3) >= 0 && (d5-d6) >= 0) {
		const double w = (d4-d3)/((d4-d3)+(d5-d6)), bc[3] = {c[0]-b[0], c[1]-b[1], c[2]-b[2]};
		set(b, w, bc, 0, ac); return;
	}
	const double den = 1.0/(va+vb+vc);
	set(a, vb*den, ab, vc*den, ac);
}

// closest point of GTS surface k to x. GTS's own gts_point_triangle_closest aborts here (stack
// check inside the library), so the bb-tree only gives the candidate triangles.
bool closest_on_gts(hexa_tree_t *mesh, int k, const double x[3], double out[3], double *dist) {
	if (k < 0 || k >= (int)mesh->gdata_vec.size() || !mesh->gdata_vec[k].bbt) return false;
	GtsPoint *p = gts_point_new(gts_point_class(), x[0], x[1], x[2]);
	GSList *list = gts_bb_tree_point_closest_bboxes(mesh->gdata_vec[k].bbt, p);
	gts_object_destroy(GTS_OBJECT(p));
	double best = std::numeric_limits<double>::infinity();
	for (GSList *l = list; l; l = l->next) {
		GtsTriangle *t = GTS_TRIANGLE(GTS_BBOX(l->data)->bounded);
		GtsVertex *v1, *v2, *v3;
		gts_triangle_vertices(t, &v1, &v2, &v3);
		const double a[3] = {GTS_POINT(v1)->x, GTS_POINT(v1)->y, GTS_POINT(v1)->z},
		             b[3] = {GTS_POINT(v2)->x, GTS_POINT(v2)->y, GTS_POINT(v2)->z},
		             c[3] = {GTS_POINT(v3)->x, GTS_POINT(v3)->y, GTS_POINT(v3)->z};
		double q[3];
		closest_on_triangle(x, a, b, c, q);
		const double d2 = (q[0]-x[0])*(q[0]-x[0]) + (q[1]-x[1])*(q[1]-x[1]) + (q[2]-x[2])*(q[2]-x[2]);
		if (d2 < best) { best = d2; for (int d = 0; d < 3; d++) out[d] = q[d]; }
	}
	g_slist_free(list);
	if (!std::isfinite(best)) return false;
	if (dist) *dist = std::sqrt(best);
	return true;
}

void ProjectPillowedInterface(hexa_tree_t *mesh, std::vector<double> &coords) {
	if (getenv("PF_NOPROJ") && atoi(getenv("PF_NOPROJ"))) { printf("    ProjectPillowedInterface: skipped (PF_NOPROJ)\n"); return; }
	if (mesh->mpi_size > 1) { printf("    ProjectPillowedInterface: serial only for now\n"); return; }
	const int nn = mesh->nodes.elem_count, ne = mesh->elements.elem_count;
	auto ids_of = [&](int e, int id[8]) {
		octant_t *el = (octant_t *) sc_array_index(&mesh->elements, e);
		for (int k = 0; k < 8; k++) id[k] = el->nodes[mgeom::H5_ORD[k]].id;
	};
	auto mat_of = [&](int e) { return ((octant_t *) sc_array_index(&mesh->elements, e))->n_mat; };
	auto hex_of = [&](int e) {
		int id[8]; ids_of(e, id);
		pillowgeom::Hex h;
		for (int k = 0; k < 8; k++) for (int d = 0; d < 3; d++) h[k][d] = coords[3*id[k]+d];
		return h;
	};
	double vol = 0;
	for (int e = 0; e < ne; e++) vol += pillowgeom::jacobian(hex_of(e), 0.5, 0.5, 0.5);
	const int orient = vol >= 0 ? 1 : -1;

	// faces: interface (two elements, two materials) and boundary (one element)
	struct F { std::array<int,4> key, ring; int e; };
	std::vector<F> faces;
	faces.reserve(6 * ne);
	for (int e = 0; e < ne; e++) {
		int id[8]; ids_of(e, id);
		for (const auto &f : pillowgeom::faces) {
			F x; x.e = e;
			for (int k = 0; k < 4; k++) x.ring[k] = x.key[k] = id[f[k]];
			std::sort(x.key.begin(), x.key.end());
			faces.push_back(x);
		}
	}
	std::sort(faces.begin(), faces.end(), [](const F &a, const F &b) { return a.key < b.key; });
	std::vector<char> on_boundary(nn, 0);
	std::vector<std::vector<int>> node_mats(nn), nbr(nn);
	std::vector<std::vector<int>> inc(nn);
	for (int e = 0; e < ne; e++) { int id[8]; ids_of(e, id); for (int k = 0; k < 8; k++) inc[id[k]].push_back(e); }
	size_t n_iface_faces = 0;
	for (size_t i = 0; i < faces.size();) {
		size_t j = i + 1;
		while (j < faces.size() && faces[j].key == faces[i].key) j++;
		if (j - i == 1) for (int n : faces[i].key) on_boundary[n] = 1;
		else if (j - i == 2 && mat_of(faces[i].e) != mat_of(faces[i+1].e)) {
			n_iface_faces++;
			const int ma = mat_of(faces[i].e), mb = mat_of(faces[i+1].e);
			for (int k = 0; k < 4; k++) {
				const int n = faces[i].ring[k];
				for (int m : {ma, mb}) if (std::find(node_mats[n].begin(), node_mats[n].end(), m) == node_mats[n].end()) node_mats[n].push_back(m);
				for (int nb : {faces[i].ring[(k+1)%4], faces[i].ring[(k+3)%4]})
					if (std::find(nbr[n].begin(), nbr[n].end(), nb) == nbr[n].end()) nbr[n].push_back(nb);
			}
		}
		i = j;
	}
	// local cell size: longest interface edge at the node (lattice edges of the host faces; the
	// pillow's own edges are shorter, and the cells are often much taller than wide)
	auto h_of = [&](int n) {
		double h = 0;
		for (int nb : nbr[n]) {
			double d2 = 0;
			for (int d = 0; d < 3; d++) { const double t = coords[3*n+d] - coords[3*nb+d]; d2 += t*t; }
			h = std::max(h, std::sqrt(d2));
		}
		return h;
	};
	// movers: interface nodes off the boundary, two materials, a GTS within 1.5 local edges
	const double FAR = getenv("PF_FAR") ? atof(getenv("PF_FAR")) : 1.5;   // A/B
	const std::vector<double> x0 = coords;
	std::vector<int> movers, surf_of(nn, -1);
	int n_iface = 0, n_bnd = 0, n_triple = 0, n_far = 0;
	for (int n = 0; n < nn; n++) {
		if (node_mats[n].empty()) continue;
		n_iface++;
		if (on_boundary[n]) { n_bnd++; continue; }
		if (node_mats[n].size() > 2) { n_triple++; continue; }
		double best = 1e300, c[3], d;
		for (int k = 0; k < (int)mesh->gdata_vec.size(); k++)
			if (closest_on_gts(mesh, k, &coords[3*n], c, &d) && d < best) { best = d; surf_of[n] = k; }
		if (surf_of[n] < 0 || best > FAR * h_of(n)) {
			if (getenv("PF_DEBUG")) printf("      far node %d mats %d/%d z %.1f dist %.1f h %.1f\n", n, node_mats[n][0], node_mats[n][1], coords[3*n+2], best, h_of(n));
			surf_of[n] = -1; n_far++; continue;
		}
		movers.push_back(n);
	}
	printf("    ProjectPillowedInterface: %zu interface faces, %d interface nodes: %zu move, %d on the boundary (coast/walls), %d triple junctions, %d off any GTS\n",
	       n_iface_faces, n_iface, movers.size(), n_bnd, n_triple, n_far);
	std::vector<char> is_mover(nn, 0);
	for (int n : movers) is_mover[n] = 1;

	// a. closest point
	for (int n : movers) closest_on_gts(mesh, surf_of[n], &x0[3*n], &coords[3*n], nullptr);
	// b. tangential smoothing on the GTS (Jacobi, relaxed); fixed interface nodes are anchors
	int SWEEPS = getenv("PF_SMOOTH") ? atoi(getenv("PF_SMOOTH")) : 5;
	for (int sw = 0; sw < SWEEPS; sw++) {
		std::vector<double> nxt(3 * movers.size());
		for (size_t i = 0; i < movers.size(); i++) {
			const int n = movers[i];
			double a[3] = {0, 0, 0};
			for (int nb : nbr[n]) for (int d = 0; d < 3; d++) a[d] += coords[3*nb+d] / nbr[n].size();
			double y[3];
			for (int d = 0; d < 3; d++) y[d] = 0.5 * coords[3*n+d] + 0.5 * a[d];
			if (!closest_on_gts(mesh, surf_of[n], y, &nxt[3*i], nullptr)) for (int d = 0; d < 3; d++) nxt[3*i+d] = coords[3*n+d];
		}
		for (size_t i = 0; i < movers.size(); i++) for (int d = 0; d < 3; d++) coords[3*movers[i]+d] = nxt[3*i+d];
	}
	// intended points, for meshOpt closeness
	mesh->proj_target.assign(coords.size(), std::numeric_limits<double>::quiet_NaN());
	mesh->proj_len.assign(nn, 0.0);
	for (int n : movers) {
		double l2 = 0;
		for (int d = 0; d < 3; d++) { mesh->proj_target[3*n+d] = coords[3*n+d]; l2 += (coords[3*n+d]-x0[3*n+d])*(coords[3*n+d]-x0[3*n+d]); }
		mesh->proj_len[n] = std::sqrt(l2);
	}
	const std::vector<double> target = coords;

	// columns: the original nodes sharing a lattice (x, y), bottom to top. The displacement of the
	// interface nodes is spread along their column (option 1 of the plan's results), linearly
	// between anchors: the column's ends, every interface or boundary node that does not move
	// (displacement 0) and the movers (their own) -- as GetMeshFromSurface stretches the columns to
	// the topography. The shells alone absorb about half a cell; most movers need one or more.
	const bool STRETCH = !(getenv("PF_STRETCH") && atoi(getenv("PF_STRETCH")) == 0);
	int first_new = nn;
	for (const auto &r : mesh->pillow_hosts) for (int k = 0; k < 8; k++) first_new = std::min(first_new, r.core[k]);
	std::vector<std::vector<int>> columns;
	std::vector<int> col_of(nn, -1);
	{
		std::map<std::pair<int,int>, int> idx;
		for (int n = 0; n < first_new; n++) {
			const auto *nd = (octant_node_t *) sc_array_index(&mesh->nodes, n);
			auto it = idx.emplace(std::make_pair(nd->x, nd->y), (int)columns.size()).first;
			if (it->second == (int)columns.size()) columns.emplace_back();
			columns[it->second].push_back(n);
			col_of[n] = it->second;
		}
		for (auto &c : columns) std::sort(c.begin(), c.end(), [&](int a, int b) { return x0[3*a+2] < x0[3*b+2]; });
	}
	std::vector<std::vector<int>> col_movers(columns.size());
	for (int n : movers) col_movers[col_of[n]].push_back(n);
	std::vector<char> anchor(nn, 0);
	for (int n = 0; n < first_new; n++) anchor[n] = is_mover[n] || on_boundary[n] || !node_mats[n].empty();

	// c. cores follow their hosts
	std::vector<std::vector<int>> hosts_of(nn);
	for (size_t i = 0; i < mesh->pillow_hosts.size(); i++) for (int k = 0; k < 8; k++) hosts_of[mesh->pillow_hosts[i].outer[k]].push_back(i);
	auto rebuild_all = [&]() {
		for (const auto &r : mesh->pillow_hosts) {
			pillowgeom::Hex outer, core;
			for (int k = 0; k < 8; k++) for (int d = 0; d < 3; d++) outer[k][d] = coords[3*r.outer[k]+d];
			pillowgeom::core_from_insets(outer, r.t, core);
			for (int k = 0; k < 8; k++) for (int d = 0; d < 3; d++) coords[3*r.core[k]+d] = core[k][d];
		}
	};
	// place everything for the share frac[n] of each mover's move
	std::vector<double> frac(nn, 1.0);
	auto apply = [&]() {
		for (int n : movers) for (int d = 0; d < 3; d++) coords[3*n+d] = x0[3*n+d] + frac[n] * (target[3*n+d] - x0[3*n+d]);
		if (STRETCH) for (size_t c = 0; c < columns.size(); c++) {
			if (col_movers[c].empty()) continue;
			const auto &col = columns[c];
			int a = 0;   // last anchor (the bottom end is one)
			for (size_t i = 1; i < col.size(); i++) {
				if (!anchor[col[i]] && i + 1 < col.size()) continue;
				const int na = col[a], nb = col[i];
				const double za = x0[3*na+2], zb = x0[3*nb+2];
				for (size_t j = a + 1; j < i; j++) {
					const int n = col[j];
					const double w = zb > za ? (x0[3*n+2] - za) / (zb - za) : 0.0;
					for (int d = 0; d < 3; d++)
						coords[3*n+d] = x0[3*n+d] + (1-w) * (coords[3*na+d] - x0[3*na+d]) + w * (coords[3*nb+d] - x0[3*nb+d]);
				}
				a = i;
			}
		}
		rebuild_all();
	};
	apply();

	// elements to watch: every element with a node that moved
	std::vector<int> watch;
	for (int e = 0; e < ne; e++) {
		int id[8]; ids_of(e, id);
		for (int k = 0; k < 8; k++) if (coords[3*id[k]] != x0[3*id[k]] || coords[3*id[k]+1] != x0[3*id[k]+1] || coords[3*id[k]+2] != x0[3*id[k]+2]) { watch.push_back(e); break; }
	}
	double FRAC = 0.05;
	if (const char *ev = getenv("PROJ_SIGMA_MIN_FRAC")) FRAC = atof(ev);
	std::unordered_map<int,double> sigma0;
	{   // reference sigma on the pre-move geometry
		std::vector<double> cur = coords; coords = x0;
		for (int e : watch) sigma0[e] = pillowgeom::size_metric(hex_of(e), false);
		coords = cur;
	}
	auto bad = [&](int e) {
		const auto h = hex_of(e);
		return !pillowgeom::valid(h, orient) || (FRAC > 0 && pillowgeom::size_metric(h, false) < FRAC * sigma0[e]);
	};
	int n_bad0 = 0;
	for (int e : watch) n_bad0 += bad(e);
	printf("    ProjectPillowedInterface: after projection + %d smoothing sweeps%s: %d of %zu watched elements fail (invalid or sigma < %.2f x start)\n",
	       SWEEPS, STRETCH ? " + column stretch" : "", n_bad0, watch.size(), FRAC);
	if (getenv("PF_DEBUG")) {   // who fails, and how far the movers go relative to their cell
		std::vector<char> is_core(ne, 0), is_shell(ne, 0);
		const int first_shell = ne - 6 * (int)mesh->pillow_hosts.size();
		for (int e = first_shell; e < ne; e++) is_shell[e] = 1;
		std::unordered_set<int> core_nodes;
		for (const auto &r : mesh->pillow_hosts) for (int k = 0; k < 8; k++) core_nodes.insert(r.core[k]);
		for (int e = 0; e < first_shell; e++) { int id[8]; ids_of(e, id); is_core[e] = core_nodes.count(id[0]) > 0; }
		int c_core = 0, c_shell = 0, c_plain = 0, c_inv = 0;
		for (int e : watch) if (bad(e)) {
			(is_shell[e] ? c_shell : is_core[e] ? c_core : c_plain)++;
			c_inv += !pillowgeom::valid(hex_of(e), orient);
		}
		printf("      failing: %d shells, %d cores, %d unpillowed; %d of them invalid (the rest too thin)\n", c_shell, c_core, c_plain, c_inv);
		int hist[6] = {0};
		for (int n : movers) {
			double m2 = 0; for (int d = 0; d < 3; d++) m2 += (coords[3*n+d]-x0[3*n+d])*(coords[3*n+d]-x0[3*n+d]);
			const double r = std::sqrt(m2) / std::max(1e-9, h_of(n));
			hist[r < 0.1 ? 0 : r < 0.25 ? 1 : r < 0.5 ? 2 : r < 0.75 ? 3 : r < 1 ? 4 : 5]++;
		}
		printf("      move / cell size: <0.1 %d, <0.25 %d, <0.5 %d, <0.75 %d, <1 %d, >=1 %d\n", hist[0], hist[1], hist[2], hist[3], hist[4], hist[5]);
	}

	// d. limiter: a failing element halves the moves of the movers of the columns through its nodes
	// (a core node has no column: its host's outer nodes stand for it)
	const bool LIMIT = !(getenv("PF_LIMITER") && atoi(getenv("PF_LIMITER")) == 0);
	int rounds = 0, n_bad = n_bad0;
	if (LIMIT) {
		std::vector<std::vector<int>> host_of_core(nn);
		for (size_t i = 0; i < mesh->pillow_hosts.size(); i++) for (int k = 0; k < 8; k++) host_of_core[mesh->pillow_hosts[i].core[k]].push_back(i);
		auto add_col = [&](int n, std::unordered_set<int> &pull) {
			if (col_of[n] < 0) return;
			if (STRETCH) for (int m : col_movers[col_of[n]]) pull.insert(m);
			else if (is_mover[n]) pull.insert(n);
		};
		for (; rounds < 40 && n_bad > 0; rounds++) {
			std::unordered_set<int> pull;
			for (int e : watch) if (bad(e)) {
				int id[8]; ids_of(e, id);
				for (int k = 0; k < 8; k++) {
					add_col(id[k], pull);
					for (int h : host_of_core[id[k]]) for (int j = 0; j < 8; j++) add_col(mesh->pillow_hosts[h].outer[j], pull);
				}
			}
			if (pull.empty()) break;
			for (int n : pull) frac[n] = rounds < 30 ? 0.5 * frac[n] : 0.0;
			apply();
			n_bad = 0;
			for (int e : watch) n_bad += bad(e);
		}
		double kept = 0;
		for (int n : movers) kept += frac[n];
		printf("    ProjectPillowedInterface: limiter %d rounds, %d elements still failing, movers keep %.1f%% of their move on average\n",
		       rounds, n_bad, movers.empty() ? 100.0 : 100.0 * kept / movers.size());
	}
}
