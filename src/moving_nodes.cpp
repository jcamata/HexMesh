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
				if (vol_ok && edge_ok && sj_ok) continue;   // still healthy
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
