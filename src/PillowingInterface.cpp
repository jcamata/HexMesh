#include <gts.h>
#include <glib.h>
#include <vector>
#include <iostream>
using namespace std;
#include <set>
#include <algorithm>
#include <sc.h>
#include <sc_io.h>
#include <sc_containers.h>
#include <mpi.h>
#include <cassert>

#include "hexa.h"
#include "hilbert.h"

#include <ctime>
#include <chrono>

unsigned edget_id_hash(const void *v, const void *u) {
	const edge_t *q = (const edge_t*) v;
	uint64_t a, b, c;

	a = (uint64_t) q->id;
	b = (uint64_t) 1;
	c = (uint64_t) 1;
	sc_hash_mix(a, b, c);
	sc_hash_final(a, b, c);
	return (unsigned) c;
}

int edget_id_equal(const void *v, const void *u, const void *w) {
	const edge_t *e1 = (const edge_t*) v;
	const edge_t *e2 = (const edge_t*) u;

	return (unsigned) ((e1->id == e2->id));

}

unsigned id_hash(const void *v, const void *u) {
	const octant_edge_t *q = (const octant_edge_t*) v;
	uint64_t a, b, c;

	a = (uint64_t) q->id;
	b = (uint64_t) 1;
	c = (uint64_t) 1;
	sc_hash_mix(a, b, c);
	sc_hash_final(a, b, c);
	return (unsigned) c;
}

int id_equal(const void *v, const void *u, const void *w) {
	const octant_edge_t *e1 = (const octant_edge_t*) v;
	const octant_edge_t *e2 = (const octant_edge_t*) u;

	return (unsigned) ((e1->id == e2->id));

}

void edge_add(uint64_t id, sc_hash_array_t* hash_edge_ref) {
	size_t position;
	octant_edge_t *r;
	octant_edge_t key;
	key.id = id;

	r = (octant_edge_t*) sc_hash_array_insert_unique(hash_edge_ref, &key, &position);
	if(r != NULL){
		r->id = key.id;
		r->ref = true;
	}
}

void BuildHash(hexa_tree_t* mesh, std::vector<double>& coords, std::vector<int>& nodes_b_mat,
		sc_hash_array_t* hash_nodes, sc_hash_array_t* hash_b_mat,sc_hash_array_t*vertex_hash,
		sc_hash_array_t*hash_edge_ref){

	//hash nodes
	for(int n = 0; n < mesh->nodes.elem_count; n++){
		size_t position;
		node_t *r;
		node_t key;
		octant_node_t* node = (octant_node_t*) sc_array_index (&mesh->nodes, n);
		key.coord[0] = coords[3*node->id+0];
		key.coord[1] = coords[3*node->id+1];
		key.coord[2] = coords[3*node->id+2];
		key.node_id = node->id;

		r = (node_t*) sc_hash_array_insert_unique(hash_nodes, &key, &position);
		if(r!=NULL){
			r->coord[0] = coords[3*node->id+0];
			r->coord[1] = coords[3*node->id+1];
			r->coord[2] = coords[3*node->id+2];
			r->node_id = node->id;
		}else{
			printf("Verificar o no numero %d\n",node->id);
		}
	}

	//fazendo hash nos nodes_b_mat
	for(int n = 0; n < nodes_b_mat.size(); n++){
		size_t position;
		node_t *r;
		node_t key;
		int node = nodes_b_mat[n];
		key.coord[0] = coords[3*node+0];
		key.coord[1] = coords[3*node+1];
		key.coord[2] = coords[3*node+2];
		key.node_id = node;

		r = (node_t*) sc_hash_array_insert_unique(hash_b_mat, &key, &position);
		if(r!=NULL){
			r->coord[0] = coords[3*node+0];
			r->coord[1] = coords[3*node+1];
			r->coord[2] = coords[3*node+2];
			r->node_id = node;
			r->flag = true;
		}else{
			printf("Verificar o no numero %d\n",node);
		}
	}

	//fazendo vertex hash
	for(int ioc = 0 ; ioc < mesh->oct.elem_count; ioc++){
		octree_t * oc = (octree_t*) sc_array_index (&mesh->oct, ioc);
		for(int iel = 0; iel < 8 ; iel++){
			size_t  position;
			octant_t *elem = (octant_t*) sc_array_index(&mesh->elements, oc->id[iel]);

			for (int ino = 0; ino < 8; ino++){
				octant_vertex_t key;
				key.id = elem->nodes[ino].id;
				octant_vertex_t* vert = (octant_vertex_t*) sc_hash_array_insert_unique (vertex_hash, &key, &position);
				if(vert != NULL){
					vert->id = elem->nodes[ino].id;
					vert->list_elem = 1;
					vert->elem[vert->list_elem-1] = elem->id;
				}else{
					vert = (octant_vertex_t*) sc_array_index(&vertex_hash->a, position);
					vert->elem[vert->list_elem] = elem->id;
					vert->list_elem++;
				}
			}
		}
	}

	//fazendo edge hash
	//mesh->oct.elem_count
	for(int ioc = 0 ; ioc < mesh->oct.elem_count; ioc++){
		octree_t * oc = (octree_t*) sc_array_index (&mesh->oct, ioc);
		for(int iel = 0; iel < 8 ; iel++){
			octant_t *elem = (octant_t*) sc_array_index(&mesh->elements, oc->id[iel]);
			size_t  position;
			int ncount = 0;
			int aux[8];
			for(int ino = 0; ino < 8; ino++){
				node_t key;
				int node = elem->nodes[ino].id;
				key.coord[0] = coords[3*node+0];
				key.coord[1] = coords[3*node+1];
				key.coord[2] = coords[3*node+2];
				key.node_id = node;

				bool nodel = sc_hash_array_lookup(hash_b_mat, &key, &position);
				if(nodel){
					aux[ncount] = ino;
					//VertexEdgeMap[ino][0];
					//VertexEdgeMap[ino][1];
					//VertexEdgeMap[ino][2];
					ncount++;
				}
			}

			int id;
			if(ncount == 1 && false){

				elem->edge[VertexEdgeMap[aux[0]][0]].ref = true;
				elem->edge[VertexEdgeMap[aux[0]][1]].ref = true;
				elem->edge[VertexEdgeMap[aux[0]][2]].ref = true;

				id = elem->edge[VertexEdgeMap[aux[0]][0]].id;
				edge_add(id, hash_edge_ref);

				id = elem->edge[VertexEdgeMap[aux[0]][1]].id;
				edge_add(id, hash_edge_ref);

				id = elem->edge[VertexEdgeMap[aux[0]][2]].id;
				edge_add(id, hash_edge_ref);
			}

			if(ncount != 1 || true){
				std::vector<int> listaux;
				std::vector<int> list;
				//printf("vetor auxiliar: ");
				for(int ive = 0; ive < ncount; ive++){
					for(int iedge = 0; iedge < 3; iedge++){
						list.push_back(VertexEdgeMap[aux[ive]][iedge]);
						//printf("%d ", VertexEdgeMap[aux[ive]][iedge]);
					}
				}
				//printf("\n");

				for(int i = 0; i < list.size(); i++){
					for(int j = i + 1; j < list.size(); j++){
						if(list[i] == list[j]){
							listaux.push_back(list[i]);
						}
					}
				}

				//printf("vetor: ");
				//for(int i = 0; i < list.size(); i++){
				//	printf("%d ", list[i]);
				//}
				//printf("\n");

				for(int i = 0; i < list.size(); i++){
					bool aux = true;
					for(int j = 0; j < listaux.size(); j++){
						if(list[i]==listaux[j]){
							aux = false;
						}
					}

					if(aux){
						elem->edge[list[i]].ref = true;
						id = elem->edge[list[i]].id;
						edge_add(id, hash_edge_ref);
					}
				}
			}

			//printf("numero de nos do elemento: %d eh %d\n",elem->id,ncount);

		}

	}

}


void IdentifyTemplate(hexa_tree_t* mesh, std::vector<int>& elements_ids, sc_hash_array_t* hash_edge_ref){

	int el_1 = 0;
	int el_2 = 0;
	int el_3 = 0;
	int el_4 = 0;
	int el_5 = 0;
	int el_6 = 0;
	int el_7 = 0;
	int el_8 = 0;
	int el_9 = 0;
	int el_10 = 0;
	int el_11 = 0;

	//update the element information
	for (int iel = 0; iel < elements_ids.size(); ++iel) {
		size_t position;
		octant_t *elem = (octant_t*) sc_array_index(&mesh->elements, elements_ids[iel]);
		for(int edge=0; edge<12; edge++){
			bool out = sc_hash_array_lookup(hash_edge_ref, &elem->edge[edge].id, &position);
			elem->edge[edge].ref = out;
		}
	}

	for (int iel = 0; iel < elements_ids.size(); ++iel) {

		octant_t *elem = (octant_t*) sc_array_index(&mesh->elements, elements_ids[iel]);
		int ed_cont = 0;

		for (int edge = 0; edge < 12; ++edge) {
			if(elem->edge[edge].ref){
				ed_cont++;
			}
		}

		//template 11
		elem->pad = 144;
		elem->tem = 11;

		// element verification
		// template 1
		if(ed_cont == 1){
			if ( (elem->edge[0].ref) && (!elem->edge[1].ref) && (!elem->edge[2].ref) && (!elem->edge[3].ref) &&
					(!elem->edge[4].ref) && (!elem->edge[5].ref) && (!elem->edge[6].ref) && (!elem->edge[7].ref) &&
					(!elem->edge[8].ref) && (!elem->edge[9].ref) && (!elem->edge[10].ref) && (!elem->edge[11].ref)) {

				elem->pad = 10;
				elem->tem = 1;

			}
			else if ( (!elem->edge[0].ref) && (elem->edge[1].ref) && (!elem->edge[2].ref) && (!elem->edge[3].ref) &&
					(!elem->edge[4].ref) && (!elem->edge[5].ref) && (!elem->edge[6].ref) && (!elem->edge[7].ref) &&
					(!elem->edge[8].ref) && (!elem->edge[9].ref) && (!elem->edge[10].ref) && (!elem->edge[11].ref)) {

				elem->pad = 11;
				elem->tem = 1;

			}
			else if ( (!elem->edge[0].ref) && (!elem->edge[1].ref) && (elem->edge[2].ref) && (!elem->edge[3].ref) &&
					(!elem->edge[4].ref) && (!elem->edge[5].ref) && (!elem->edge[6].ref) && (!elem->edge[7].ref) &&
					(!elem->edge[8].ref) && (!elem->edge[9].ref) && (!elem->edge[10].ref) && (!elem->edge[11].ref)) {

				elem->pad = 12;
				elem->tem = 1;

			}
			else if ( (!elem->edge[0].ref) && (!elem->edge[1].ref) && (!elem->edge[2].ref) && (elem->edge[3].ref) &&
					(!elem->edge[4].ref) && (!elem->edge[5].ref) && (!elem->edge[6].ref) && (!elem->edge[7].ref) &&
					(!elem->edge[8].ref) && (!elem->edge[9].ref) && (!elem->edge[10].ref) && (!elem->edge[11].ref)) {

				elem->pad = 13;
				elem->tem = 1;

			}
			else if ( (!elem->edge[0].ref) && (!elem->edge[1].ref) && (!elem->edge[2].ref) && (!elem->edge[3].ref) &&
					(elem->edge[4].ref) && (!elem->edge[5].ref) && (!elem->edge[6].ref) && (!elem->edge[7].ref) &&
					(!elem->edge[8].ref) && (!elem->edge[9].ref) && (!elem->edge[10].ref) && (!elem->edge[11].ref)) {

				elem->pad = 14;
				elem->tem = 1;

			}
			else if ( (!elem->edge[0].ref) && (!elem->edge[1].ref) && (!elem->edge[2].ref) && (!elem->edge[3].ref) &&
					(!elem->edge[4].ref) && (elem->edge[5].ref) && (!elem->edge[6].ref) && (!elem->edge[7].ref) &&
					(!elem->edge[8].ref) && (!elem->edge[9].ref) && (!elem->edge[10].ref) && (!elem->edge[11].ref)) {

				elem->pad = 15;
				elem->tem = 1;

			}
			else if ( (!elem->edge[0].ref) && (!elem->edge[1].ref) && (!elem->edge[2].ref) && (!elem->edge[3].ref) &&
					(!elem->edge[4].ref) && (!elem->edge[5].ref) && (elem->edge[6].ref) && (!elem->edge[7].ref) &&
					(!elem->edge[8].ref) && (!elem->edge[9].ref) && (!elem->edge[10].ref) && (!elem->edge[11].ref)) {

				elem->pad = 16;
				elem->tem = 1;

			}
			else if ( (!elem->edge[0].ref) && (!elem->edge[1].ref) && (!elem->edge[2].ref) && (!elem->edge[3].ref) &&
					(!elem->edge[4].ref) && (!elem->edge[5].ref) && (!elem->edge[6].ref) && (elem->edge[7].ref) &&
					(!elem->edge[8].ref) && (!elem->edge[9].ref) && (!elem->edge[10].ref) && (!elem->edge[11].ref)) {

				elem->pad = 17;
				elem->tem = 1;

			}
			else if ( (!elem->edge[0].ref) && (!elem->edge[1].ref) && (!elem->edge[2].ref) && (!elem->edge[3].ref) &&
					(!elem->edge[4].ref) && (!elem->edge[5].ref) && (!elem->edge[6].ref) && (!elem->edge[7].ref) &&
					(elem->edge[8].ref) && (!elem->edge[9].ref) && (!elem->edge[10].ref) && (!elem->edge[11].ref)) {

				elem->pad = 18;
				elem->tem = 1;

			}
			else if ( (!elem->edge[0].ref) && (!elem->edge[1].ref) && (!elem->edge[2].ref) && (!elem->edge[3].ref) &&
					(!elem->edge[4].ref) && (!elem->edge[5].ref) && (!elem->edge[6].ref) && (!elem->edge[7].ref) &&
					(!elem->edge[8].ref) && (elem->edge[9].ref) && (!elem->edge[10].ref) && (!elem->edge[11].ref)) {

				elem->pad = 19;
				elem->tem = 1;

			}
			else if ( (!elem->edge[0].ref) && (!elem->edge[1].ref) && (!elem->edge[2].ref) && (!elem->edge[3].ref) &&
					(!elem->edge[4].ref) && (!elem->edge[5].ref) && (!elem->edge[6].ref) && (!elem->edge[7].ref) &&
					(!elem->edge[8].ref) && (!elem->edge[9].ref) && (elem->edge[10].ref) && (!elem->edge[11].ref)) {

				elem->pad = 20;
				elem->tem = 1;

			}
			else if ( (!elem->edge[0].ref) && (!elem->edge[1].ref) && (!elem->edge[2].ref) && (!elem->edge[3].ref) &&
					(!elem->edge[4].ref) && (!elem->edge[5].ref) && (!elem->edge[6].ref) && (!elem->edge[7].ref) &&
					(!elem->edge[8].ref) && (!elem->edge[9].ref) && (!elem->edge[10].ref) && (elem->edge[11].ref)) {

				elem->pad = 21;
				elem->tem = 1;

			}
		}
		//template 2
		if(ed_cont == 4){
			if ( (elem->edge[0].ref) && (!elem->edge[1].ref) && (elem->edge[2].ref) && (!elem->edge[3].ref) &&
					(!elem->edge[4].ref) && (!elem->edge[5].ref) && (!elem->edge[6].ref) && (!elem->edge[7].ref) &&
					(elem->edge[8].ref) && (!elem->edge[9].ref) && (elem->edge[10].ref) && (!elem->edge[11].ref)) {

				elem->pad = 22;
				elem->tem = 2;

			}
			else if ( (!elem->edge[0].ref) && (elem->edge[1].ref) && (!elem->edge[2].ref) && (elem->edge[3].ref) &&
					(!elem->edge[4].ref) && (!elem->edge[5].ref) && (!elem->edge[6].ref) && (!elem->edge[7].ref) &&
					(!elem->edge[8].ref) && (elem->edge[9].ref) && (!elem->edge[10].ref) && (elem->edge[11].ref)) {

				elem->pad = 23;
				elem->tem = 2;

			}
			else if ( (!elem->edge[0].ref) && (!elem->edge[1].ref) && (!elem->edge[2].ref) && (!elem->edge[3].ref) &&
					(elem->edge[4].ref) && (elem->edge[5].ref) && (elem->edge[6].ref) && (elem->edge[7].ref) &&
					(!elem->edge[8].ref) && (!elem->edge[9].ref) && (!elem->edge[10].ref) && (!elem->edge[11].ref)) {

				elem->pad = 24;
				elem->tem = 2;

			}
		}
		//template 3
		if(ed_cont == 2){
			if ( (elem->edge[0].ref) && (!elem->edge[1].ref) && (elem->edge[2].ref) && (!elem->edge[3].ref) &&
					(!elem->edge[4].ref) && (!elem->edge[5].ref) && (!elem->edge[6].ref) && (!elem->edge[7].ref) &&
					(!elem->edge[8].ref) && (!elem->edge[9].ref) && (!elem->edge[10].ref) && (!elem->edge[11].ref)) {

				elem->pad = 25;
				elem->tem = 3;

			}
			else if ( (!elem->edge[0].ref) && (elem->edge[1].ref) && (!elem->edge[2].ref) && (elem->edge[3].ref) &&
					(!elem->edge[4].ref) && (!elem->edge[5].ref) && (!elem->edge[6].ref) && (!elem->edge[7].ref) &&
					(!elem->edge[8].ref) && (!elem->edge[9].ref) && (!elem->edge[10].ref) && (!elem->edge[11].ref)) {

				elem->pad = 26;
				elem->tem = 3;

			}
			else if ( (!elem->edge[0].ref) && (!elem->edge[1].ref) && (!elem->edge[2].ref) && (!elem->edge[3].ref) &&
					(!elem->edge[4].ref) && (!elem->edge[5].ref) && (!elem->edge[6].ref) && (!elem->edge[7].ref) &&
					(elem->edge[8].ref) && (!elem->edge[9].ref) && (elem->edge[10].ref) && (!elem->edge[11].ref)) {

				elem->pad = 27;
				elem->tem = 3;

			}
			else if ( (!elem->edge[0].ref) && (!elem->edge[1].ref) && (!elem->edge[2].ref) && (!elem->edge[3].ref) &&
					(!elem->edge[4].ref) && (!elem->edge[5].ref) && (!elem->edge[6].ref) && (!elem->edge[7].ref) &&
					(!elem->edge[8].ref) && (elem->edge[9].ref) && (!elem->edge[10].ref) && (elem->edge[11].ref)) {

				elem->pad = 28;
				elem->tem = 3;

			}
			else if ( (!elem->edge[0].ref) && (!elem->edge[1].ref) && (!elem->edge[2].ref) && (!elem->edge[3].ref) &&
					(elem->edge[4].ref) && (elem->edge[5].ref) && (!elem->edge[6].ref) && (!elem->edge[7].ref) &&
					(!elem->edge[8].ref) && (!elem->edge[9].ref) && (!elem->edge[10].ref) && (!elem->edge[11].ref)) {

				elem->pad = 29;
				elem->tem = 3;

			}
			else if ( (!elem->edge[0].ref) && (!elem->edge[1].ref) && (!elem->edge[2].ref) && (!elem->edge[3].ref) &&
					(!elem->edge[4].ref) && (elem->edge[5].ref) && (elem->edge[6].ref) && (!elem->edge[7].ref) &&
					(!elem->edge[8].ref) && (!elem->edge[9].ref) && (!elem->edge[10].ref) && (!elem->edge[11].ref)) {

				elem->pad = 30;
				elem->tem = 3;

			}
			else if ( (!elem->edge[0].ref) && (!elem->edge[1].ref) && (!elem->edge[2].ref) && (!elem->edge[3].ref) &&
					(!elem->edge[4].ref) && (!elem->edge[5].ref) && (elem->edge[6].ref) && (elem->edge[7].ref) &&
					(!elem->edge[8].ref) && (!elem->edge[9].ref) && (!elem->edge[10].ref) && (!elem->edge[11].ref)) {

				elem->pad = 31;
				elem->tem = 3;

			}
			else if ( (!elem->edge[0].ref) && (!elem->edge[1].ref) && (!elem->edge[2].ref) && (!elem->edge[3].ref) &&
					(elem->edge[4].ref) && (!elem->edge[5].ref) && (!elem->edge[6].ref) && (elem->edge[7].ref) &&
					(!elem->edge[8].ref) && (!elem->edge[9].ref) && (!elem->edge[10].ref) && (!elem->edge[11].ref)) {

				elem->pad = 32;
				elem->tem = 3;

			}
			else if ( (elem->edge[0].ref) && (!elem->edge[1].ref) && (!elem->edge[2].ref) && (!elem->edge[3].ref) &&
					(!elem->edge[4].ref) && (!elem->edge[5].ref) && (!elem->edge[6].ref) && (!elem->edge[7].ref) &&
					(elem->edge[8].ref) && (!elem->edge[9].ref) && (!elem->edge[10].ref) && (!elem->edge[11].ref)) {

				elem->pad = 33;
				elem->tem = 3;

			}
			else if ( (!elem->edge[0].ref) && (elem->edge[1].ref) && (!elem->edge[2].ref) && (!elem->edge[3].ref) &&
					(!elem->edge[4].ref) && (!elem->edge[5].ref) && (!elem->edge[6].ref) && (!elem->edge[7].ref) &&
					(!elem->edge[8].ref) && (elem->edge[9].ref) && (!elem->edge[10].ref) && (!elem->edge[11].ref)) {

				elem->pad = 34;
				elem->tem = 3;

			}
			else if ( (!elem->edge[0].ref) && (!elem->edge[1].ref) && (elem->edge[2].ref) && (!elem->edge[3].ref) &&
					(!elem->edge[4].ref) && (!elem->edge[5].ref) && (!elem->edge[6].ref) && (!elem->edge[7].ref) &&
					(!elem->edge[8].ref) && (!elem->edge[9].ref) && (elem->edge[10].ref) && (!elem->edge[11].ref)) {

				elem->pad = 35;
				elem->tem = 3;

			}
			else if ( (!elem->edge[0].ref) && (!elem->edge[1].ref) && (!elem->edge[2].ref) && (elem->edge[3].ref) &&
					(!elem->edge[4].ref) && (!elem->edge[5].ref) && (!elem->edge[6].ref) && (!elem->edge[7].ref) &&
					(!elem->edge[8].ref) && (!elem->edge[9].ref) && (!elem->edge[10].ref) && (elem->edge[11].ref)) {

				elem->pad = 36;
				elem->tem = 3;

			}
		}
		//template 4
		if(ed_cont == 4){
			if ( (elem->edge[0].ref) && (elem->edge[1].ref) && (!elem->edge[2].ref) && (!elem->edge[3].ref) &&
					(!elem->edge[4].ref) && (!elem->edge[5].ref) && (!elem->edge[6].ref) && (!elem->edge[7].ref) &&
					(elem->edge[8].ref) && (elem->edge[9].ref) && (!elem->edge[10].ref) && (!elem->edge[11].ref)) {

				elem->pad = 37;
				elem->tem = 4;

			}
			else if ( (!elem->edge[0].ref) && (elem->edge[1].ref) && (elem->edge[2].ref) && (!elem->edge[3].ref) &&
					(!elem->edge[4].ref) && (!elem->edge[5].ref) && (!elem->edge[6].ref) && (!elem->edge[7].ref) &&
					(!elem->edge[8].ref) && (elem->edge[9].ref) && (elem->edge[10].ref) && (!elem->edge[11].ref)) {

				elem->pad = 38;
				elem->tem = 4;

			}
			else if ( (!elem->edge[0].ref) && (!elem->edge[1].ref) && (elem->edge[2].ref) && (elem->edge[3].ref) &&
					(!elem->edge[4].ref) && (!elem->edge[5].ref) && (!elem->edge[6].ref) && (!elem->edge[7].ref) &&
					(!elem->edge[8].ref) && (!elem->edge[9].ref) && (elem->edge[10].ref) && (elem->edge[11].ref)) {

				elem->pad = 39;
				elem->tem = 4;

			}
			else if ( (elem->edge[0].ref) && (!elem->edge[1].ref) && (!elem->edge[2].ref) && (elem->edge[3].ref) &&
					(!elem->edge[4].ref) && (!elem->edge[5].ref) && (!elem->edge[6].ref) && (!elem->edge[7].ref) &&
					(elem->edge[8].ref) && (!elem->edge[9].ref) && (!elem->edge[10].ref) && (elem->edge[11].ref)) {

				elem->pad = 40;
				elem->tem = 4;

			}
			else if ( (elem->edge[0].ref) && (!elem->edge[1].ref) && (elem->edge[2].ref) && (!elem->edge[3].ref) &&
					(elem->edge[4].ref) && (!elem->edge[5].ref) && (!elem->edge[6].ref) && (elem->edge[7].ref) &&
					(!elem->edge[8].ref) && (!elem->edge[9].ref) && (!elem->edge[10].ref) && (!elem->edge[11].ref)) {

				elem->pad = 41;
				elem->tem = 4;

			}
			else if ( (elem->edge[0].ref) && (!elem->edge[1].ref) && (elem->edge[2].ref) && (!elem->edge[3].ref) &&
					(!elem->edge[4].ref) && (elem->edge[5].ref) && (elem->edge[6].ref) && (!elem->edge[7].ref) &&
					(!elem->edge[8].ref) && (!elem->edge[9].ref) && (!elem->edge[10].ref) && (!elem->edge[11].ref)) {

				elem->pad = 42;
				elem->tem = 4;

			}
			else if ( (!elem->edge[0].ref) && (!elem->edge[1].ref) && (!elem->edge[2].ref) && (!elem->edge[3].ref) &&
					(!elem->edge[4].ref) && (elem->edge[5].ref) && (elem->edge[6].ref) && (!elem->edge[7].ref) &&
					(elem->edge[8].ref) && (!elem->edge[9].ref) && (elem->edge[10].ref) && (!elem->edge[11].ref)) {

				elem->pad = 43;
				elem->tem = 4;

			}
			else if ( (!elem->edge[0].ref) && (!elem->edge[1].ref) && (!elem->edge[2].ref) && (!elem->edge[3].ref) &&
					(elem->edge[4].ref) && (!elem->edge[5].ref) && (!elem->edge[6].ref) && (elem->edge[7].ref) &&
					(elem->edge[8].ref) && (!elem->edge[9].ref) && (elem->edge[10].ref) && (!elem->edge[11].ref)) {

				elem->pad = 44;
				elem->tem = 4;

			}
			else if ( (!elem->edge[0].ref) && (elem->edge[1].ref) && (!elem->edge[2].ref) && (elem->edge[3].ref) &&
					(elem->edge[4].ref) && (elem->edge[5].ref) && (!elem->edge[6].ref) && (!elem->edge[7].ref) &&
					(!elem->edge[8].ref) && (!elem->edge[9].ref) && (!elem->edge[10].ref) && (!elem->edge[11].ref)) {

				elem->pad = 45;
				elem->tem = 4;

			}
			else if ( (!elem->edge[0].ref) && (elem->edge[1].ref) && (!elem->edge[2].ref) && (elem->edge[3].ref) &&
					(!elem->edge[4].ref) && (!elem->edge[5].ref) && (elem->edge[6].ref) && (elem->edge[7].ref) &&
					(!elem->edge[8].ref) && (!elem->edge[9].ref) && (!elem->edge[10].ref) && (!elem->edge[11].ref)) {

				elem->pad = 46;
				elem->tem = 4;

			}
			else if ( (!elem->edge[0].ref) && (!elem->edge[1].ref) && (!elem->edge[2].ref) && (!elem->edge[3].ref) &&
					(!elem->edge[4].ref) && (!elem->edge[5].ref) && (elem->edge[6].ref) && (elem->edge[7].ref) &&
					(!elem->edge[8].ref) && (elem->edge[9].ref) && (!elem->edge[10].ref) && (elem->edge[11].ref)) {

				elem->pad = 47;
				elem->tem = 4;

			}
			else if ( (!elem->edge[0].ref) && (!elem->edge[1].ref) && (!elem->edge[2].ref) && (!elem->edge[3].ref) &&
					(elem->edge[4].ref) && (elem->edge[5].ref) && (!elem->edge[6].ref) && (!elem->edge[7].ref) &&
					(!elem->edge[8].ref) && (elem->edge[9].ref) && (!elem->edge[10].ref) && (elem->edge[11].ref)) {

				elem->pad = 48;
				elem->tem = 4;

			}
		}
		//template 5
		if(ed_cont == 4 || ed_cont == 3 || ed_cont == 2){
			//template 5 com 4
			if ( (elem->edge[0].ref) && (elem->edge[1].ref) && (elem->edge[2].ref) && (elem->edge[3].ref) &&
					(!elem->edge[4].ref) && (!elem->edge[5].ref) && (!elem->edge[6].ref) && (!elem->edge[7].ref) &&
					(!elem->edge[8].ref) && (!elem->edge[9].ref) && (!elem->edge[10].ref) && (!elem->edge[11].ref)) {

				elem->pad = 49;
				elem->tem = 5;

			}
			else if ( (!elem->edge[0].ref) && (!elem->edge[1].ref) && (!elem->edge[2].ref) && (!elem->edge[3].ref) &&
					(!elem->edge[4].ref) && (!elem->edge[5].ref) && (!elem->edge[6].ref) && (!elem->edge[7].ref) &&
					(elem->edge[8].ref) && (elem->edge[9].ref) && (elem->edge[10].ref) && (elem->edge[11].ref)) {

				elem->pad = 50;
				elem->tem = 5;

			}
			else if ( (elem->edge[0].ref) && (!elem->edge[1].ref) && (!elem->edge[2].ref) && (!elem->edge[3].ref) &&
					(elem->edge[4].ref) && (elem->edge[5].ref) && (!elem->edge[6].ref) && (!elem->edge[7].ref) &&
					(elem->edge[8].ref) && (!elem->edge[9].ref) && (!elem->edge[10].ref) && (!elem->edge[11].ref)) {

				elem->pad = 51;
				elem->tem = 5;

			}
			else if ( (!elem->edge[0].ref) && (elem->edge[1].ref) && (!elem->edge[2].ref) && (!elem->edge[3].ref) &&
					(!elem->edge[4].ref) && (elem->edge[5].ref) && (elem->edge[6].ref) && (!elem->edge[7].ref) &&
					(!elem->edge[8].ref) && (elem->edge[9].ref) && (!elem->edge[10].ref) && (!elem->edge[11].ref)) {

				elem->pad = 52;
				elem->tem = 5;

			}
			else if ( (!elem->edge[0].ref) && (!elem->edge[1].ref) && (elem->edge[2].ref) && (!elem->edge[3].ref) &&
					(!elem->edge[4].ref) && (!elem->edge[5].ref) && (elem->edge[6].ref) && (elem->edge[7].ref) &&
					(!elem->edge[8].ref) && (!elem->edge[9].ref) && (elem->edge[10].ref) && (!elem->edge[11].ref)) {

				elem->pad = 53;
				elem->tem = 5;

			}
			else if ( (!elem->edge[0].ref) && (!elem->edge[1].ref) && (!elem->edge[2].ref) && (elem->edge[3].ref) &&
					(elem->edge[4].ref) && (!elem->edge[5].ref) && (!elem->edge[6].ref) && (elem->edge[7].ref) &&
					(!elem->edge[8].ref) && (!elem->edge[9].ref) && (!elem->edge[10].ref) && (elem->edge[11].ref)) {

				elem->pad = 54;
				elem->tem = 5;

			}


			//template 5 com 3
			else if ( (elem->edge[0].ref) && (elem->edge[1].ref) && (!elem->edge[2].ref) && (elem->edge[3].ref) &&
					(!elem->edge[4].ref) && (!elem->edge[5].ref) && (!elem->edge[6].ref) && (!elem->edge[7].ref) &&
					(!elem->edge[8].ref) && (!elem->edge[9].ref) && (!elem->edge[10].ref) && (!elem->edge[11].ref)) {

				elem->pad = 55;
				elem->tem = 5;

			}
			else if ( (elem->edge[0].ref) && (elem->edge[1].ref) && (elem->edge[2].ref) && (!elem->edge[3].ref) &&
					(!elem->edge[4].ref) && (!elem->edge[5].ref) && (!elem->edge[6].ref) && (!elem->edge[7].ref) &&
					(!elem->edge[8].ref) && (!elem->edge[9].ref) && (!elem->edge[10].ref) && (!elem->edge[11].ref)) {

				elem->pad = 56;
				elem->tem = 5;

			}
			else if ( (!elem->edge[0].ref) && (elem->edge[1].ref) && (elem->edge[2].ref) && (elem->edge[3].ref) &&
					(!elem->edge[4].ref) && (!elem->edge[5].ref) && (!elem->edge[6].ref) && (!elem->edge[7].ref) &&
					(!elem->edge[8].ref) && (!elem->edge[9].ref) && (!elem->edge[10].ref) && (!elem->edge[11].ref)) {

				elem->pad = 57;
				elem->tem = 5;

			}
			else if ( (elem->edge[0].ref) && (!elem->edge[1].ref) && (elem->edge[2].ref) && (elem->edge[3].ref) &&
					(!elem->edge[4].ref) && (!elem->edge[5].ref) && (!elem->edge[6].ref) && (!elem->edge[7].ref) &&
					(!elem->edge[8].ref) && (!elem->edge[9].ref) && (!elem->edge[10].ref) && (!elem->edge[11].ref)) {

				elem->pad = 58;
				elem->tem = 5;

			}
			else if ( (!elem->edge[0].ref) && (!elem->edge[1].ref) && (!elem->edge[2].ref) && (!elem->edge[3].ref) &&
					(!elem->edge[4].ref) && (!elem->edge[5].ref) && (!elem->edge[6].ref) && (!elem->edge[7].ref) &&
					(elem->edge[8].ref) && (elem->edge[9].ref) && (elem->edge[10].ref) && (!elem->edge[11].ref)) {

				elem->pad = 59;
				elem->tem = 5;

			}
			else if ( (!elem->edge[0].ref) && (!elem->edge[1].ref) && (!elem->edge[2].ref) && (!elem->edge[3].ref) &&
					(!elem->edge[4].ref) && (!elem->edge[5].ref) && (!elem->edge[6].ref) && (!elem->edge[7].ref) &&
					(elem->edge[8].ref) && (elem->edge[9].ref) && (!elem->edge[10].ref) && (elem->edge[11].ref)) {

				elem->pad = 60;
				elem->tem = 5;

			}
			else if ( (!elem->edge[0].ref) && (!elem->edge[1].ref) && (!elem->edge[2].ref) && (!elem->edge[3].ref) &&
					(!elem->edge[4].ref) && (!elem->edge[5].ref) && (!elem->edge[6].ref) && (!elem->edge[7].ref) &&
					(!elem->edge[8].ref) && (elem->edge[9].ref) && (elem->edge[10].ref) && (elem->edge[11].ref)) {

				elem->pad = 61;
				elem->tem = 5;

			}
			else if ( (!elem->edge[0].ref) && (!elem->edge[1].ref) && (!elem->edge[2].ref) && (!elem->edge[3].ref) &&
					(!elem->edge[4].ref) && (!elem->edge[5].ref) && (!elem->edge[6].ref) && (!elem->edge[7].ref) &&
					(elem->edge[8].ref) && (!elem->edge[9].ref) && (elem->edge[10].ref) && (elem->edge[11].ref)) {

				elem->pad = 62;
				elem->tem = 5;

			}
			else if ( (elem->edge[0].ref) && (!elem->edge[1].ref) && (!elem->edge[2].ref) && (!elem->edge[3].ref) &&
					(elem->edge[4].ref) && (elem->edge[5].ref) && (!elem->edge[6].ref) && (!elem->edge[7].ref) &&
					(!elem->edge[8].ref) && (!elem->edge[9].ref) && (!elem->edge[10].ref) && (!elem->edge[11].ref)) {

				elem->pad = 63;
				elem->tem = 5;

			}
			else if ( (elem->edge[0].ref) && (!elem->edge[1].ref) && (!elem->edge[2].ref) && (!elem->edge[3].ref) &&
					(!elem->edge[4].ref) && (elem->edge[5].ref) && (!elem->edge[6].ref) && (!elem->edge[7].ref) &&
					(elem->edge[8].ref) && (!elem->edge[9].ref) && (!elem->edge[10].ref) && (!elem->edge[11].ref)) {

				elem->pad = 64;
				elem->tem = 5;

			}
			else if ( (!elem->edge[0].ref) && (!elem->edge[1].ref) && (!elem->edge[2].ref) && (!elem->edge[3].ref) &&
					(elem->edge[4].ref) && (elem->edge[5].ref) && (!elem->edge[6].ref) && (!elem->edge[7].ref) &&
					(elem->edge[8].ref) && (!elem->edge[9].ref) && (!elem->edge[10].ref) && (!elem->edge[11].ref)) {

				elem->pad = 65;
				elem->tem = 5;

			}
			else if ( (elem->edge[0].ref) && (!elem->edge[1].ref) && (!elem->edge[2].ref) && (!elem->edge[3].ref) &&
					(elem->edge[4].ref) && (!elem->edge[5].ref) && (!elem->edge[6].ref) && (!elem->edge[7].ref) &&
					(elem->edge[8].ref) && (!elem->edge[9].ref) && (!elem->edge[10].ref) && (!elem->edge[11].ref)) {

				elem->pad = 66;
				elem->tem = 5;

			}
			else if ( (!elem->edge[0].ref) && (elem->edge[1].ref) && (!elem->edge[2].ref) && (!elem->edge[3].ref) &&
					(!elem->edge[4].ref) && (elem->edge[5].ref) && (elem->edge[6].ref) && (!elem->edge[7].ref) &&
					(!elem->edge[8].ref) && (!elem->edge[9].ref) && (!elem->edge[10].ref) && (!elem->edge[11].ref)) {

				elem->pad = 67;
				elem->tem = 5;

			}
			else if ( (!elem->edge[0].ref) && (elem->edge[1].ref) && (!elem->edge[2].ref) && (!elem->edge[3].ref) &&
					(!elem->edge[4].ref) && (elem->edge[5].ref) && (!elem->edge[6].ref) && (!elem->edge[7].ref) &&
					(!elem->edge[8].ref) && (elem->edge[9].ref) && (!elem->edge[10].ref) && (!elem->edge[11].ref)) {

				elem->pad = 68;
				elem->tem = 5;

			}
			else if ( (!elem->edge[0].ref) && (!elem->edge[1].ref) && (!elem->edge[2].ref) && (!elem->edge[3].ref) &&
					(!elem->edge[4].ref) && (elem->edge[5].ref) && (elem->edge[6].ref) && (!elem->edge[7].ref) &&
					(!elem->edge[8].ref) && (elem->edge[9].ref) && (!elem->edge[10].ref) && (!elem->edge[11].ref)) {

				elem->pad = 69;
				elem->tem = 5;

			}
			else if ( (!elem->edge[0].ref) && (elem->edge[1].ref) && (!elem->edge[2].ref) && (!elem->edge[3].ref) &&
					(!elem->edge[4].ref) && (!elem->edge[5].ref) && (elem->edge[6].ref) && (!elem->edge[7].ref) &&
					(!elem->edge[8].ref) && (elem->edge[9].ref) && (!elem->edge[10].ref) && (!elem->edge[11].ref)) {

				elem->pad = 70;
				elem->tem = 5;

			}
			else if ( (!elem->edge[0].ref) && (!elem->edge[1].ref) && (elem->edge[2].ref) && (!elem->edge[3].ref) &&
					(!elem->edge[4].ref) && (!elem->edge[5].ref) && (elem->edge[6].ref) && (elem->edge[7].ref) &&
					(!elem->edge[8].ref) && (!elem->edge[9].ref) && (!elem->edge[10].ref) && (!elem->edge[11].ref)) {

				elem->pad = 71;
				elem->tem = 5;

			}
			else if ( (!elem->edge[0].ref) && (!elem->edge[1].ref) && (elem->edge[2].ref) && (!elem->edge[3].ref) &&
					(!elem->edge[4].ref) && (!elem->edge[5].ref) && (elem->edge[6].ref) && (!elem->edge[7].ref) &&
					(!elem->edge[8].ref) && (!elem->edge[9].ref) && (elem->edge[10].ref) && (!elem->edge[11].ref)) {

				elem->pad = 72;
				elem->tem = 5;

			}
			else if ( (!elem->edge[0].ref) && (!elem->edge[1].ref) && (!elem->edge[2].ref) && (!elem->edge[3].ref) &&
					(!elem->edge[4].ref) && (!elem->edge[5].ref) && (elem->edge[6].ref) && (elem->edge[7].ref) &&
					(!elem->edge[8].ref) && (!elem->edge[9].ref) && (elem->edge[10].ref) && (!elem->edge[11].ref)) {

				elem->pad = 73;
				elem->tem = 5;

			}
			else if ( (!elem->edge[0].ref) && (!elem->edge[1].ref) && (elem->edge[2].ref) && (!elem->edge[3].ref) &&
					(!elem->edge[4].ref) && (!elem->edge[5].ref) && (!elem->edge[6].ref) && (elem->edge[7].ref) &&
					(!elem->edge[8].ref) && (!elem->edge[9].ref) && (elem->edge[10].ref) && (!elem->edge[11].ref)) {

				elem->pad = 74;
				elem->tem = 5;

			}
			else if ( (!elem->edge[0].ref) && (!elem->edge[1].ref) && (!elem->edge[2].ref) && (elem->edge[3].ref) &&
					(elem->edge[4].ref) && (!elem->edge[5].ref) && (!elem->edge[6].ref) && (elem->edge[7].ref) &&
					(!elem->edge[8].ref) && (!elem->edge[9].ref) && (!elem->edge[10].ref) && (!elem->edge[11].ref)) {

				elem->pad = 75;
				elem->tem = 5;

			}
			else if ( (!elem->edge[0].ref) && (!elem->edge[1].ref) && (!elem->edge[2].ref) && (elem->edge[3].ref) &&
					(elem->edge[4].ref) && (!elem->edge[5].ref) && (!elem->edge[6].ref) && (!elem->edge[7].ref) &&
					(!elem->edge[8].ref) && (!elem->edge[9].ref) && (!elem->edge[10].ref) && (elem->edge[11].ref)) {

				elem->pad = 76;
				elem->tem = 5;

			}
			else if ( (!elem->edge[0].ref) && (!elem->edge[1].ref) && (!elem->edge[2].ref) && (!elem->edge[3].ref) &&
					(elem->edge[4].ref) && (!elem->edge[5].ref) && (!elem->edge[6].ref) && (elem->edge[7].ref) &&
					(!elem->edge[8].ref) && (!elem->edge[9].ref) && (!elem->edge[10].ref) && (elem->edge[11].ref)) {

				elem->pad = 77;
				elem->tem = 5;

			}
			else if ( (!elem->edge[0].ref) && (!elem->edge[1].ref) && (!elem->edge[2].ref) && (elem->edge[3].ref) &&
					(!elem->edge[4].ref) && (!elem->edge[5].ref) && (!elem->edge[6].ref) && (elem->edge[7].ref) &&
					(!elem->edge[8].ref) && (!elem->edge[9].ref) && (!elem->edge[10].ref) && (elem->edge[11].ref)) {

				elem->pad = 78;
				elem->tem = 5;

			}

			// templete 5 com 2
			else if ( (elem->edge[0].ref) && (elem->edge[1].ref) && (!elem->edge[2].ref) && (!elem->edge[3].ref) &&
					(!elem->edge[4].ref) && (!elem->edge[5].ref) && (!elem->edge[6].ref) && (!elem->edge[7].ref) &&
					(!elem->edge[8].ref) && (!elem->edge[9].ref) && (!elem->edge[10].ref) && (!elem->edge[11].ref)) {

				elem->pad = 79;
				elem->tem = 5;

			}
			else if ( (!elem->edge[0].ref) && (elem->edge[1].ref) && (elem->edge[2].ref) && (!elem->edge[3].ref) &&
					(!elem->edge[4].ref) && (!elem->edge[5].ref) && (!elem->edge[6].ref) && (!elem->edge[7].ref) &&
					(!elem->edge[8].ref) && (!elem->edge[9].ref) && (!elem->edge[10].ref) && (!elem->edge[11].ref)) {

				elem->pad = 80;
				elem->tem = 5;

			}
			else if ( (!elem->edge[0].ref) && (!elem->edge[1].ref) && (elem->edge[2].ref) && (elem->edge[3].ref) &&
					(!elem->edge[4].ref) && (!elem->edge[5].ref) && (!elem->edge[6].ref) && (!elem->edge[7].ref) &&
					(!elem->edge[8].ref) && (!elem->edge[9].ref) && (!elem->edge[10].ref) && (!elem->edge[11].ref)) {

				elem->pad = 81;
				elem->tem = 5;

			}
			else if ( (elem->edge[0].ref) && (!elem->edge[1].ref) && (!elem->edge[2].ref) && (elem->edge[3].ref) &&
					(!elem->edge[4].ref) && (!elem->edge[5].ref) && (!elem->edge[6].ref) && (!elem->edge[7].ref) &&
					(!elem->edge[8].ref) && (!elem->edge[9].ref) && (!elem->edge[10].ref) && (!elem->edge[11].ref)) {

				elem->pad = 82;
				elem->tem = 5;

			}
			else if ( (!elem->edge[0].ref) && (!elem->edge[1].ref) && (!elem->edge[2].ref) && (!elem->edge[3].ref) &&
					(!elem->edge[4].ref) && (!elem->edge[5].ref) && (!elem->edge[6].ref) && (!elem->edge[7].ref) &&
					(elem->edge[8].ref) && (elem->edge[9].ref) && (!elem->edge[10].ref) && (!elem->edge[11].ref)) {

				elem->pad = 83;
				elem->tem = 5;

			}
			else if ( (!elem->edge[0].ref) && (!elem->edge[1].ref) && (!elem->edge[2].ref) && (!elem->edge[3].ref) &&
					(!elem->edge[4].ref) && (!elem->edge[5].ref) && (!elem->edge[6].ref) && (!elem->edge[7].ref) &&
					(!elem->edge[8].ref) && (elem->edge[9].ref) && (elem->edge[10].ref) && (!elem->edge[11].ref)) {

				elem->pad = 84;
				elem->tem = 5;

			}
			else if ( (!elem->edge[0].ref) && (!elem->edge[1].ref) && (!elem->edge[2].ref) && (!elem->edge[3].ref) &&
					(!elem->edge[4].ref) && (!elem->edge[5].ref) && (!elem->edge[6].ref) && (!elem->edge[7].ref) &&
					(!elem->edge[8].ref) && (!elem->edge[9].ref) && (elem->edge[10].ref) && (elem->edge[11].ref)) {

				elem->pad = 85;
				elem->tem = 5;

			}
			else if ( (!elem->edge[0].ref) && (!elem->edge[1].ref) && (!elem->edge[2].ref) && (!elem->edge[3].ref) &&
					(!elem->edge[4].ref) && (!elem->edge[5].ref) && (!elem->edge[6].ref) && (!elem->edge[7].ref) &&
					(elem->edge[8].ref) && (!elem->edge[9].ref) && (!elem->edge[10].ref) && (elem->edge[11].ref)) {

				elem->pad = 86;
				elem->tem = 5;

			}
			else if ( (elem->edge[0].ref) && (!elem->edge[1].ref) && (!elem->edge[2].ref) && (!elem->edge[3].ref) &&
					(elem->edge[4].ref) && (!elem->edge[5].ref) && (!elem->edge[6].ref) && (!elem->edge[7].ref) &&
					(!elem->edge[8].ref) && (!elem->edge[9].ref) && (!elem->edge[10].ref) && (!elem->edge[11].ref)) {

				elem->pad = 87;
				elem->tem = 5;

			}
			else if ( (elem->edge[0].ref) && (!elem->edge[1].ref) && (!elem->edge[2].ref) && (!elem->edge[3].ref) &&
					(!elem->edge[4].ref) && (elem->edge[5].ref) && (!elem->edge[6].ref) && (!elem->edge[7].ref) &&
					(!elem->edge[8].ref) && (!elem->edge[9].ref) && (!elem->edge[10].ref) && (!elem->edge[11].ref)) {

				elem->pad = 88;
				elem->tem = 5;

			}
			else if ( (!elem->edge[0].ref) && (!elem->edge[1].ref) && (!elem->edge[2].ref) && (!elem->edge[3].ref) &&
					(!elem->edge[4].ref) && (elem->edge[5].ref) && (!elem->edge[6].ref) && (!elem->edge[7].ref) &&
					(elem->edge[8].ref) && (!elem->edge[9].ref) && (!elem->edge[10].ref) && (!elem->edge[11].ref)) {

				elem->pad = 89;
				elem->tem = 5;

			}
			else if ( (!elem->edge[0].ref) && (!elem->edge[1].ref) && (!elem->edge[2].ref) && (!elem->edge[3].ref) &&
					(elem->edge[4].ref) && (!elem->edge[5].ref) && (!elem->edge[6].ref) && (!elem->edge[7].ref) &&
					(elem->edge[8].ref) && (!elem->edge[9].ref) && (!elem->edge[10].ref) && (!elem->edge[11].ref)) {

				elem->pad = 90;
				elem->tem = 5;

			}
			else if ( (!elem->edge[0].ref) && (!elem->edge[1].ref) && (!elem->edge[2].ref) && (!elem->edge[3].ref) &&
					(!elem->edge[4].ref) && (elem->edge[5].ref) && (!elem->edge[6].ref) && (!elem->edge[7].ref) &&
					(!elem->edge[8].ref) && (elem->edge[9].ref) && (!elem->edge[10].ref) && (!elem->edge[11].ref)) {

				elem->pad = 91;
				elem->tem = 5;

			}
			else if ( (!elem->edge[0].ref) && (!elem->edge[1].ref) && (!elem->edge[2].ref) && (!elem->edge[3].ref) &&
					(!elem->edge[4].ref) && (!elem->edge[5].ref) && (elem->edge[6].ref) && (!elem->edge[7].ref) &&
					(!elem->edge[8].ref) && (elem->edge[9].ref) && (!elem->edge[10].ref) && (!elem->edge[11].ref)) {

				elem->pad = 92;
				elem->tem = 5;

			}
			else if ( (!elem->edge[0].ref) && (elem->edge[1].ref) && (!elem->edge[2].ref) && (!elem->edge[3].ref) &&
					(!elem->edge[4].ref) && (!elem->edge[5].ref) && (elem->edge[6].ref) && (!elem->edge[7].ref) &&
					(!elem->edge[8].ref) && (!elem->edge[9].ref) && (!elem->edge[10].ref) && (!elem->edge[11].ref)) {

				elem->pad = 93;
				elem->tem = 5;

			}
			else if ( (!elem->edge[0].ref) && (elem->edge[1].ref) && (!elem->edge[2].ref) && (!elem->edge[3].ref) &&
					(!elem->edge[4].ref) && (elem->edge[5].ref) && (!elem->edge[6].ref) && (!elem->edge[7].ref) &&
					(!elem->edge[8].ref) && (!elem->edge[9].ref) && (!elem->edge[10].ref) && (!elem->edge[11].ref)) {

				elem->pad = 94;
				elem->tem = 5;

			}
			else if ( (!elem->edge[0].ref) && (!elem->edge[1].ref) && (elem->edge[2].ref) && (!elem->edge[3].ref) &&
					(!elem->edge[4].ref) && (!elem->edge[5].ref) && (elem->edge[6].ref) && (!elem->edge[7].ref) &&
					(!elem->edge[8].ref) && (!elem->edge[9].ref) && (!elem->edge[10].ref) && (!elem->edge[11].ref)) {

				elem->pad = 95;
				elem->tem = 5;

			}
			else if ( (!elem->edge[0].ref) && (!elem->edge[1].ref) && (!elem->edge[2].ref) && (!elem->edge[3].ref) &&
					(!elem->edge[4].ref) && (!elem->edge[5].ref) && (elem->edge[6].ref) && (!elem->edge[7].ref) &&
					(!elem->edge[8].ref) && (!elem->edge[9].ref) && (elem->edge[10].ref) && (!elem->edge[11].ref)) {

				elem->pad = 96;
				elem->tem = 5;

			}
			else if ( (!elem->edge[0].ref) && (!elem->edge[1].ref) && (!elem->edge[2].ref) && (!elem->edge[3].ref) &&
					(!elem->edge[4].ref) && (!elem->edge[5].ref) && (!elem->edge[6].ref) && (elem->edge[7].ref) &&
					(!elem->edge[8].ref) && (!elem->edge[9].ref) && (elem->edge[10].ref) && (!elem->edge[11].ref)) {

				elem->pad = 97;
				elem->tem = 5;

			}
			else if ( (!elem->edge[0].ref) && (!elem->edge[1].ref) && (elem->edge[2].ref) && (!elem->edge[3].ref) &&
					(!elem->edge[4].ref) && (!elem->edge[5].ref) && (!elem->edge[6].ref) && (elem->edge[7].ref) &&
					(!elem->edge[8].ref) && (!elem->edge[9].ref) && (!elem->edge[10].ref) && (!elem->edge[11].ref)) {

				elem->pad = 98;
				elem->tem = 5;

			}
			else if ( (!elem->edge[0].ref) && (!elem->edge[1].ref) && (!elem->edge[2].ref) && (elem->edge[3].ref) &&
					(elem->edge[4].ref) && (!elem->edge[5].ref) && (!elem->edge[6].ref) && (!elem->edge[7].ref) &&
					(!elem->edge[8].ref) && (!elem->edge[9].ref) && (!elem->edge[10].ref) && (!elem->edge[11].ref)) {

				elem->pad = 99;
				elem->tem = 5;

			}
			else if ( (!elem->edge[0].ref) && (!elem->edge[1].ref) && (!elem->edge[2].ref) && (!elem->edge[3].ref) &&
					(elem->edge[4].ref) && (!elem->edge[5].ref) && (elem->edge[6].ref) && (!elem->edge[7].ref) &&
					(!elem->edge[8].ref) && (!elem->edge[9].ref) && (!elem->edge[10].ref) && (!elem->edge[11].ref)) {

				elem->pad = 100;
				elem->tem = 5;

			}
			else if ( (!elem->edge[0].ref) && (!elem->edge[1].ref) && (!elem->edge[2].ref) && (!elem->edge[3].ref) &&
					(!elem->edge[4].ref) && (!elem->edge[5].ref) && (elem->edge[6].ref) && (elem->edge[7].ref) &&
					(!elem->edge[8].ref) && (!elem->edge[9].ref) && (!elem->edge[10].ref) && (!elem->edge[11].ref)) {

				elem->pad = 101;
				elem->tem = 5;

			}
			else if ( (!elem->edge[0].ref) && (!elem->edge[1].ref) && (!elem->edge[2].ref) && (elem->edge[3].ref) &&
					(!elem->edge[4].ref) && (!elem->edge[5].ref) && (!elem->edge[6].ref) && (elem->edge[7].ref) &&
					(!elem->edge[8].ref) && (!elem->edge[9].ref) && (!elem->edge[10].ref) && (!elem->edge[11].ref)) {

				elem->pad = 102;
				elem->tem = 5;

			}
		}
		//template 6
		if(ed_cont == 6){
			if ( (elem->edge[0].ref) && (elem->edge[1].ref) && (!elem->edge[2].ref) && (elem->edge[3].ref) &&
					(!elem->edge[4].ref) && (!elem->edge[5].ref) && (!elem->edge[6].ref) && (!elem->edge[7].ref) &&
					(elem->edge[8].ref) && (elem->edge[9].ref) && (!elem->edge[10].ref) && (elem->edge[11].ref)) {

				elem->pad = 103;
				elem->tem = 6;

			}
			else if ( (elem->edge[0].ref) && (elem->edge[1].ref) && (elem->edge[2].ref) && (!elem->edge[3].ref) &&
					(!elem->edge[4].ref) && (!elem->edge[5].ref) && (!elem->edge[6].ref) && (!elem->edge[7].ref) &&
					(elem->edge[8].ref) && (elem->edge[9].ref) && (elem->edge[10].ref) && (!elem->edge[11].ref)) {

				elem->pad = 104;
				elem->tem = 6;

			}
			else if ( (!elem->edge[0].ref) && (elem->edge[1].ref) && (elem->edge[2].ref) && (elem->edge[3].ref) &&
					(!elem->edge[4].ref) && (!elem->edge[5].ref) && (!elem->edge[6].ref) && (!elem->edge[7].ref) &&
					(!elem->edge[8].ref) && (elem->edge[9].ref) && (elem->edge[10].ref) && (elem->edge[11].ref)) {

				elem->pad = 105;
				elem->tem = 6;

			}
			else if ( (elem->edge[0].ref) && (!elem->edge[1].ref) && (elem->edge[2].ref) && (elem->edge[3].ref) &&
					(!elem->edge[4].ref) && (!elem->edge[5].ref) && (!elem->edge[6].ref) && (!elem->edge[7].ref) &&
					(elem->edge[8].ref) && (!elem->edge[9].ref) && (elem->edge[10].ref) && (elem->edge[11].ref)) {

				elem->pad = 106;
				elem->tem = 6;

			}
			else if ( (elem->edge[0].ref) && (!elem->edge[1].ref) && (elem->edge[2].ref) && (!elem->edge[3].ref) &&
					(!elem->edge[4].ref) && (elem->edge[5].ref) && (elem->edge[6].ref) && (!elem->edge[7].ref) &&
					(elem->edge[8].ref) && (!elem->edge[9].ref) && (elem->edge[10].ref) && (!elem->edge[11].ref)) {

				elem->pad = 107;
				elem->tem = 6;

			}
			else if ( (!elem->edge[0].ref) && (!elem->edge[1].ref) && (!elem->edge[2].ref) && (!elem->edge[3].ref) &&
					(elem->edge[4].ref) && (elem->edge[5].ref) && (elem->edge[6].ref) && (elem->edge[7].ref) &&
					(elem->edge[8].ref) && (!elem->edge[9].ref) && (elem->edge[10].ref) && (!elem->edge[11].ref)) {

				elem->pad = 108;
				elem->tem = 6;

			}
			else if ( (elem->edge[0].ref) && (!elem->edge[1].ref) && (elem->edge[2].ref) && (!elem->edge[3].ref) &&
					(elem->edge[4].ref) && (!elem->edge[5].ref) && (!elem->edge[6].ref) && (elem->edge[7].ref) &&
					(elem->edge[8].ref) && (!elem->edge[9].ref) && (elem->edge[10].ref) && (!elem->edge[11].ref)) {

				elem->pad = 109;
				elem->tem = 6;

			}
			else if ( (elem->edge[0].ref) && (!elem->edge[1].ref) && (elem->edge[2].ref) && (!elem->edge[3].ref) &&
					(elem->edge[4].ref) && (elem->edge[5].ref) && (elem->edge[6].ref) && (elem->edge[7].ref) &&
					(!elem->edge[8].ref) && (!elem->edge[9].ref) && (!elem->edge[10].ref) && (!elem->edge[11].ref)) {

				elem->pad = 110;
				elem->tem = 6;

			}
			else if ( (!elem->edge[0].ref) && (elem->edge[1].ref) && (!elem->edge[2].ref) && (elem->edge[3].ref) &&
					(elem->edge[4].ref) && (elem->edge[5].ref) && (elem->edge[6].ref) && (elem->edge[7].ref) &&
					(!elem->edge[8].ref) && (!elem->edge[9].ref) && (!elem->edge[10].ref) && (!elem->edge[11].ref)) {

				elem->pad = 111;
				elem->tem = 6;

			}
			else if ( (!elem->edge[0].ref) && (elem->edge[1].ref) && (!elem->edge[2].ref) && (elem->edge[3].ref) &&
					(elem->edge[4].ref) && (elem->edge[5].ref) && (!elem->edge[6].ref) && (!elem->edge[7].ref) &&
					(!elem->edge[8].ref) && (elem->edge[9].ref) && (!elem->edge[10].ref) && (elem->edge[11].ref)) {

				elem->pad = 112;
				elem->tem = 6;

			}
			else if ( (!elem->edge[0].ref) && (!elem->edge[1].ref) && (!elem->edge[2].ref) && (!elem->edge[3].ref) &&
					(elem->edge[4].ref) && (elem->edge[5].ref) && (elem->edge[6].ref) && (elem->edge[7].ref) &&
					(!elem->edge[8].ref) && (elem->edge[9].ref) && (!elem->edge[10].ref) && (elem->edge[11].ref)) {

				elem->pad = 113;
				elem->tem = 6;

			}
			else if ( (!elem->edge[0].ref) && (elem->edge[1].ref) && (!elem->edge[2].ref) && (elem->edge[3].ref) &&
					(!elem->edge[4].ref) && (!elem->edge[5].ref) && (elem->edge[6].ref) && (elem->edge[7].ref) &&
					(!elem->edge[8].ref) && (elem->edge[9].ref) && (!elem->edge[10].ref) && (elem->edge[11].ref)) {

				elem->pad = 114;
				elem->tem = 6;

			}
		}
		//template 7
		if(ed_cont == 3){
			if ( (elem->edge[0].ref) && (!elem->edge[1].ref) && (!elem->edge[2].ref) && (elem->edge[3].ref) &&
					(elem->edge[4].ref) && (!elem->edge[5].ref) && (!elem->edge[6].ref) && (!elem->edge[7].ref) &&
					(!elem->edge[8].ref) && (!elem->edge[9].ref) && (!elem->edge[10].ref) && (!elem->edge[11].ref)) {

				elem->pad = 115;
				elem->tem = 7;

			}
			else if ( (elem->edge[0].ref) && (elem->edge[1].ref) && (!elem->edge[2].ref) && (!elem->edge[3].ref) &&
					(!elem->edge[4].ref) && (elem->edge[5].ref) && (!elem->edge[6].ref) && (!elem->edge[7].ref) &&
					(!elem->edge[8].ref) && (!elem->edge[9].ref) && (!elem->edge[10].ref) && (!elem->edge[11].ref)) {

				elem->pad = 116;
				elem->tem = 7;

			}
			else if ( (!elem->edge[0].ref) && (elem->edge[1].ref) && (elem->edge[2].ref) && (!elem->edge[3].ref) &&
					(!elem->edge[4].ref) && (!elem->edge[5].ref) && (elem->edge[6].ref) && (!elem->edge[7].ref) &&
					(!elem->edge[8].ref) && (!elem->edge[9].ref) && (!elem->edge[10].ref) && (!elem->edge[11].ref)) {

				elem->pad = 117;
				elem->tem = 7;

			}
			else if ( (!elem->edge[0].ref) && (!elem->edge[1].ref) && (elem->edge[2].ref) && (elem->edge[3].ref) &&
					(!elem->edge[4].ref) && (!elem->edge[5].ref) && (!elem->edge[6].ref) && (elem->edge[7].ref) &&
					(!elem->edge[8].ref) && (!elem->edge[9].ref) && (!elem->edge[10].ref) && (!elem->edge[11].ref)) {

				elem->pad = 118;
				elem->tem = 7;

			}
			else if ( (!elem->edge[0].ref) && (!elem->edge[1].ref) && (!elem->edge[2].ref) && (!elem->edge[3].ref) &&
					(elem->edge[4].ref) && (!elem->edge[5].ref) && (!elem->edge[6].ref) && (!elem->edge[7].ref) &&
					(elem->edge[8].ref) && (!elem->edge[9].ref) && (!elem->edge[10].ref) && (elem->edge[11].ref)) {

				elem->pad = 119;
				elem->tem = 7;

			}
			else if ( (!elem->edge[0].ref) && (!elem->edge[1].ref) && (!elem->edge[2].ref) && (!elem->edge[3].ref) &&
					(!elem->edge[4].ref) && (elem->edge[5].ref) && (!elem->edge[6].ref) && (!elem->edge[7].ref) &&
					(elem->edge[8].ref) && (elem->edge[9].ref) && (!elem->edge[10].ref) && (!elem->edge[11].ref)) {

				elem->pad = 120;
				elem->tem = 7;

			}
			else if ( (!elem->edge[0].ref) && (!elem->edge[1].ref) && (!elem->edge[2].ref) && (!elem->edge[3].ref) &&
					(!elem->edge[4].ref) && (!elem->edge[5].ref) && (elem->edge[6].ref) && (!elem->edge[7].ref) &&
					(!elem->edge[8].ref) && (elem->edge[9].ref) && (elem->edge[10].ref) && (!elem->edge[11].ref)) {

				elem->pad = 121;
				elem->tem = 7;

			}
			else if ( (!elem->edge[0].ref) && (!elem->edge[1].ref) && (!elem->edge[2].ref) && (!elem->edge[3].ref) &&
					(!elem->edge[4].ref) && (!elem->edge[5].ref) && (!elem->edge[6].ref) && (elem->edge[7].ref) &&
					(!elem->edge[8].ref) && (!elem->edge[9].ref) && (elem->edge[10].ref) && (elem->edge[11].ref)) {

				elem->pad = 122;
				elem->tem = 7;

			}
		}
		//template 8
		if(ed_cont == 8){
			if ( (elem->edge[0].ref) && (elem->edge[1].ref) && (elem->edge[2].ref) && (elem->edge[3].ref) &&
					(!elem->edge[4].ref) && (!elem->edge[5].ref) && (!elem->edge[6].ref) && (!elem->edge[7].ref) &&
					(elem->edge[8].ref) && (elem->edge[9].ref) && (elem->edge[10].ref) && (elem->edge[11].ref)) {

				elem->pad = 123;
				elem->tem = 8;

			}
			else if ( (elem->edge[0].ref) && (!elem->edge[1].ref) && (elem->edge[2].ref) && (!elem->edge[3].ref) &&
					(elem->edge[4].ref) && (elem->edge[5].ref) && (elem->edge[6].ref) && (elem->edge[7].ref) &&
					(elem->edge[8].ref) && (!elem->edge[9].ref) && (elem->edge[10].ref) && (!elem->edge[11].ref)) {

				elem->pad = 124;
				elem->tem = 8;

			}
			else if ( (!elem->edge[0].ref) && (elem->edge[1].ref) && (!elem->edge[2].ref) && (elem->edge[3].ref) &&
					(elem->edge[4].ref) && (elem->edge[5].ref) && (elem->edge[6].ref) && (elem->edge[7].ref) &&
					(!elem->edge[8].ref) && (elem->edge[9].ref) && (!elem->edge[10].ref) && (elem->edge[11].ref)) {

				elem->pad = 125;
				elem->tem = 8;

			}
		}
		//template 9
		if(ed_cont == 8){
			if ( (elem->edge[0].ref) && (elem->edge[1].ref) && (elem->edge[2].ref) && (elem->edge[3].ref) &&
					(elem->edge[4].ref) && (elem->edge[5].ref) && (elem->edge[6].ref) && (elem->edge[7].ref) &&
					(!elem->edge[8].ref) && (!elem->edge[9].ref) && (!elem->edge[10].ref) && (!elem->edge[11].ref)) {

				elem->pad = 126;
				elem->tem = 9;

			}
			else if ( (!elem->edge[0].ref) && (!elem->edge[1].ref) && (!elem->edge[2].ref) && (!elem->edge[3].ref) &&
					(elem->edge[4].ref) && (elem->edge[5].ref) && (elem->edge[6].ref) && (elem->edge[7].ref) &&
					(elem->edge[8].ref) && (elem->edge[9].ref) && (elem->edge[10].ref) && (elem->edge[11].ref)) {

				elem->pad = 127;
				elem->tem = 9;

			}
			else if ( (elem->edge[0].ref) && (elem->edge[1].ref) && (!elem->edge[2].ref) && (elem->edge[3].ref) &&
					(elem->edge[4].ref) && (elem->edge[5].ref) && (!elem->edge[6].ref) && (!elem->edge[7].ref) &&
					(elem->edge[8].ref) && (elem->edge[9].ref) && (!elem->edge[10].ref) && (elem->edge[11].ref)) {

				elem->pad = 128;
				elem->tem = 9;

			}
			else if ( (!elem->edge[0].ref) && (elem->edge[1].ref) && (elem->edge[2].ref) && (elem->edge[3].ref) &&
					(!elem->edge[4].ref) && (!elem->edge[5].ref) && (elem->edge[6].ref) && (elem->edge[7].ref) &&
					(!elem->edge[8].ref) && (elem->edge[9].ref) && (elem->edge[10].ref) && (elem->edge[11].ref)) {

				elem->pad = 129;
				elem->tem = 9;

			}
			else if ( (elem->edge[0].ref) && (elem->edge[1].ref) && (elem->edge[2].ref) && (!elem->edge[3].ref) &&
					(!elem->edge[4].ref) && (elem->edge[5].ref) && (elem->edge[6].ref) && (!elem->edge[7].ref) &&
					(elem->edge[8].ref) && (elem->edge[9].ref) && (elem->edge[10].ref) && (!elem->edge[11].ref)) {

				elem->pad = 130;
				elem->tem = 9;

			}
			else if ( (elem->edge[0].ref) && (!elem->edge[1].ref) && (elem->edge[2].ref) && (elem->edge[3].ref) &&
					(elem->edge[4].ref) && (!elem->edge[5].ref) && (!elem->edge[6].ref) && (elem->edge[7].ref) &&
					(elem->edge[8].ref) && (!elem->edge[9].ref) && (elem->edge[10].ref) && (elem->edge[11].ref)) {

				elem->pad = 131;
				elem->tem = 9;

			}
		}
		//template 10
		if(ed_cont == 5){
			if ( (elem->edge[0].ref) && (elem->edge[1].ref) && (!elem->edge[2].ref) && (elem->edge[3].ref) &&
					(elem->edge[4].ref) && (elem->edge[5].ref) && (!elem->edge[6].ref) && (!elem->edge[7].ref) &&
					(!elem->edge[8].ref) && (!elem->edge[9].ref) && (!elem->edge[10].ref) && (!elem->edge[11].ref)) {

				elem->pad = 132;
				elem->tem = 10;

			}
			else if ( (elem->edge[0].ref) && (elem->edge[1].ref) && (elem->edge[2].ref) && (!elem->edge[3].ref) &&
					(!elem->edge[4].ref) && (elem->edge[5].ref) && (elem->edge[6].ref) && (!elem->edge[7].ref) &&
					(!elem->edge[8].ref) && (!elem->edge[9].ref) && (!elem->edge[10].ref) && (!elem->edge[11].ref)) {

				elem->pad = 133;
				elem->tem = 10;

			}
			else if ( (!elem->edge[0].ref) && (elem->edge[1].ref) && (elem->edge[2].ref) && (elem->edge[3].ref) &&
					(!elem->edge[4].ref) && (!elem->edge[5].ref) && (elem->edge[6].ref) && (elem->edge[7].ref) &&
					(!elem->edge[8].ref) && (!elem->edge[9].ref) && (!elem->edge[10].ref) && (!elem->edge[11].ref)) {

				elem->pad = 134;
				elem->tem = 10;

			}
			else if ( (elem->edge[0].ref) && (!elem->edge[1].ref) && (elem->edge[2].ref) && (elem->edge[3].ref) &&
					(elem->edge[4].ref) && (!elem->edge[5].ref) && (!elem->edge[6].ref) && (elem->edge[7].ref) &&
					(!elem->edge[8].ref) && (!elem->edge[9].ref) && (!elem->edge[10].ref) && (!elem->edge[11].ref)) {

				elem->pad = 135;
				elem->tem = 10;

			}
			else if ( (elem->edge[0].ref) && (!elem->edge[1].ref) && (!elem->edge[2].ref) && (elem->edge[3].ref) &&
					(elem->edge[4].ref) && (!elem->edge[5].ref) && (!elem->edge[6].ref) && (!elem->edge[7].ref) &&
					(elem->edge[8].ref) && (!elem->edge[9].ref) && (!elem->edge[10].ref) && (elem->edge[11].ref)) {

				elem->pad = 136;
				elem->tem = 10;

			}
			else if ( (elem->edge[0].ref) && (elem->edge[1].ref) && (!elem->edge[2].ref) && (!elem->edge[3].ref) &&
					(!elem->edge[4].ref) && (elem->edge[5].ref) && (!elem->edge[6].ref) && (!elem->edge[7].ref) &&
					(elem->edge[8].ref) && (elem->edge[9].ref) && (!elem->edge[10].ref) && (!elem->edge[11].ref)) {

				elem->pad = 137;
				elem->tem = 10;

			}
			else if ( (!elem->edge[0].ref) && (elem->edge[1].ref) && (elem->edge[2].ref) && (!elem->edge[3].ref) &&
					(!elem->edge[4].ref) && (!elem->edge[5].ref) && (elem->edge[6].ref) && (!elem->edge[7].ref) &&
					(!elem->edge[8].ref) && (elem->edge[9].ref) && (elem->edge[10].ref) && (!elem->edge[11].ref)) {

				elem->pad = 138;
				elem->tem = 10;

			}
			else if ( (!elem->edge[0].ref) && (!elem->edge[1].ref) && (elem->edge[2].ref) && (elem->edge[3].ref) &&
					(!elem->edge[4].ref) && (!elem->edge[5].ref) && (!elem->edge[6].ref) && (elem->edge[7].ref) &&
					(!elem->edge[8].ref) && (!elem->edge[9].ref) && (elem->edge[10].ref) && (elem->edge[11].ref)) {

				elem->pad = 139;
				elem->tem = 10;

			}
			else if ( (!elem->edge[0].ref) && (!elem->edge[1].ref) && (!elem->edge[2].ref) && (!elem->edge[3].ref) &&
					(elem->edge[4].ref) && (elem->edge[5].ref) && (!elem->edge[6].ref) && (!elem->edge[7].ref) &&
					(elem->edge[8].ref) && (elem->edge[9].ref) && (!elem->edge[10].ref) && (elem->edge[11].ref)) {

				elem->pad = 140;
				elem->tem = 10;

			}
			else if ( (!elem->edge[0].ref) && (!elem->edge[1].ref) && (!elem->edge[2].ref) && (!elem->edge[3].ref) &&
					(!elem->edge[4].ref) && (elem->edge[5].ref) && (elem->edge[6].ref) && (!elem->edge[7].ref) &&
					(elem->edge[8].ref) && (elem->edge[9].ref) && (elem->edge[10].ref) && (!elem->edge[11].ref)) {

				elem->pad = 141;
				elem->tem = 10;

			}
			else if ( (!elem->edge[0].ref) && (!elem->edge[1].ref) && (!elem->edge[2].ref) && (!elem->edge[3].ref) &&
					(!elem->edge[4].ref) && (!elem->edge[5].ref) && (elem->edge[6].ref) && (elem->edge[7].ref) &&
					(!elem->edge[8].ref) && (elem->edge[9].ref) && (elem->edge[10].ref) && (elem->edge[11].ref)) {

				elem->pad = 142;
				elem->tem = 10;

			}
			else if ( (!elem->edge[0].ref) && (!elem->edge[1].ref) && (!elem->edge[2].ref) && (!elem->edge[3].ref) &&
					(elem->edge[4].ref) && (!elem->edge[5].ref) && (!elem->edge[6].ref) && (elem->edge[7].ref) &&
					(elem->edge[8].ref) && (!elem->edge[9].ref) && (elem->edge[10].ref) && (elem->edge[11].ref)) {

				elem->pad = 143;
				elem->tem = 10;

			}
		}
		if(elem->tem== 1){el_1++;}
		if(elem->tem== 2){el_2++;}
		if(elem->tem== 3){el_3++;}
		if(elem->tem== 4){el_4++;}
		if(elem->tem== 5){el_5++;}
		if(elem->tem== 6){el_6++;}
		if(elem->tem== 7){el_7++;}
		if(elem->tem== 8){el_8++;}
		if(elem->tem== 9){el_9++;}
		if(elem->tem==10){el_10++;}
		if(elem->tem==11){el_11++;}
	}

	int su = 0;
	su = el_1+el_2+el_3+el_4+el_5+el_6+el_7+el_8+el_9+el_10+el_11;

#ifdef HEXA_DEBUG_
	fprintf(mesh->fdbg,"case_1: %d\n",el_1);
	fprintf(mesh->fdbg,"case_2: %d\n",el_2);
	fprintf(mesh->fdbg,"case_3: %d\n",el_3);
	fprintf(mesh->fdbg,"case_4: %d\n",el_4);
	fprintf(mesh->fdbg,"case_5: %d\n",el_5);
	fprintf(mesh->fdbg,"case_6: %d\n",el_6);
	fprintf(mesh->fdbg,"case_7: %d\n",el_7);
	fprintf(mesh->fdbg,"case_8: %d\n",el_8);
	fprintf(mesh->fdbg,"case_9: %d\n",el_9);
	fprintf(mesh->fdbg,"case_10: %d\n",el_10);
	fprintf(mesh->fdbg,"case_11: %d\n",el_11);
	fprintf(mesh->fdbg,"sum: %d\n",su);
#endif
}

void Edge_identification(hexa_tree_t* mesh, std::vector<int>& elements_ids, sc_hash_array_t* hash_edge_ref) {

	for (int iel = 0; iel < elements_ids.size(); ++iel) {

		octant_t *elem = (octant_t*) sc_array_index(&mesh->elements, elements_ids[iel]);

		//template 11
		if(elem->pad==144){
			for (int edge = 0; edge < 12; ++edge) {
				edge_add(elem->edge[edge].id, hash_edge_ref );
			}
		}
		//template 10
		else if(elem->pad==143){
			edge_add(elem->edge[4].id, hash_edge_ref );
			edge_add(elem->edge[7].id, hash_edge_ref );
			edge_add(elem->edge[8].id, hash_edge_ref );
			edge_add(elem->edge[10].id, hash_edge_ref );
			edge_add(elem->edge[11].id, hash_edge_ref );
		}else if(elem->pad==142){
			edge_add(elem->edge[6].id, hash_edge_ref );
			edge_add(elem->edge[7].id, hash_edge_ref );
			edge_add(elem->edge[9].id, hash_edge_ref );
			edge_add(elem->edge[10].id, hash_edge_ref );
			edge_add(elem->edge[11].id, hash_edge_ref );
		}else if(elem->pad==141){
			edge_add(elem->edge[5].id, hash_edge_ref );
			edge_add(elem->edge[6].id, hash_edge_ref );
			edge_add(elem->edge[8].id, hash_edge_ref );
			edge_add(elem->edge[9].id, hash_edge_ref );
			edge_add(elem->edge[10].id, hash_edge_ref );
		}else if(elem->pad==140){
			edge_add(elem->edge[4].id, hash_edge_ref );
			edge_add(elem->edge[5].id, hash_edge_ref );
			edge_add(elem->edge[8].id, hash_edge_ref );
			edge_add(elem->edge[9].id, hash_edge_ref );
			edge_add(elem->edge[11].id, hash_edge_ref );
		}else if(elem->pad==139){
			edge_add(elem->edge[2].id, hash_edge_ref );
			edge_add(elem->edge[3].id, hash_edge_ref );
			edge_add(elem->edge[7].id, hash_edge_ref );
			edge_add(elem->edge[10].id, hash_edge_ref );
			edge_add(elem->edge[11].id, hash_edge_ref );
		}else if(elem->pad==138){
			edge_add(elem->edge[1].id, hash_edge_ref );
			edge_add(elem->edge[2].id, hash_edge_ref );
			edge_add(elem->edge[6].id, hash_edge_ref );
			edge_add(elem->edge[9].id, hash_edge_ref );
			edge_add(elem->edge[10].id, hash_edge_ref );
		}else if(elem->pad==137){
			edge_add(elem->edge[0].id, hash_edge_ref );
			edge_add(elem->edge[1].id, hash_edge_ref );
			edge_add(elem->edge[5].id, hash_edge_ref );
			edge_add(elem->edge[8].id, hash_edge_ref );
			edge_add(elem->edge[9].id, hash_edge_ref );
		}else if(elem->pad==136){
			edge_add(elem->edge[0].id, hash_edge_ref );
			edge_add(elem->edge[3].id, hash_edge_ref );
			edge_add(elem->edge[4].id, hash_edge_ref );
			edge_add(elem->edge[8].id, hash_edge_ref );
			edge_add(elem->edge[11].id, hash_edge_ref );
		}else if(elem->pad==135){
			edge_add(elem->edge[0].id, hash_edge_ref );
			edge_add(elem->edge[2].id, hash_edge_ref );
			edge_add(elem->edge[3].id, hash_edge_ref );
			edge_add(elem->edge[4].id, hash_edge_ref );
			edge_add(elem->edge[7].id, hash_edge_ref );
		}else if(elem->pad==134){
			edge_add(elem->edge[1].id, hash_edge_ref );
			edge_add(elem->edge[2].id, hash_edge_ref );
			edge_add(elem->edge[3].id, hash_edge_ref );
			edge_add(elem->edge[6].id, hash_edge_ref );
			edge_add(elem->edge[7].id, hash_edge_ref );
		}else if(elem->pad==133){
			edge_add(elem->edge[0].id, hash_edge_ref );
			edge_add(elem->edge[1].id, hash_edge_ref );
			edge_add(elem->edge[2].id, hash_edge_ref );
			edge_add(elem->edge[5].id, hash_edge_ref );
			edge_add(elem->edge[6].id, hash_edge_ref );
		}else if(elem->pad==132){
			edge_add(elem->edge[0].id, hash_edge_ref );
			edge_add(elem->edge[1].id, hash_edge_ref );
			edge_add(elem->edge[3].id, hash_edge_ref );
			edge_add(elem->edge[4].id, hash_edge_ref );
			edge_add(elem->edge[5].id, hash_edge_ref );
		}
		//template 9
		else if(elem->pad==131){
			edge_add(elem->edge[0].id, hash_edge_ref );
			edge_add(elem->edge[2].id, hash_edge_ref );
			edge_add(elem->edge[3].id, hash_edge_ref );
			edge_add(elem->edge[4].id, hash_edge_ref );
			edge_add(elem->edge[7].id, hash_edge_ref );
			edge_add(elem->edge[8].id, hash_edge_ref );
			edge_add(elem->edge[10].id, hash_edge_ref );
			edge_add(elem->edge[11].id, hash_edge_ref );
		}else if(elem->pad==130){
			edge_add(elem->edge[0].id, hash_edge_ref );
			edge_add(elem->edge[1].id, hash_edge_ref );
			edge_add(elem->edge[3].id, hash_edge_ref );
			edge_add(elem->edge[5].id, hash_edge_ref );
			edge_add(elem->edge[6].id, hash_edge_ref );
			edge_add(elem->edge[8].id, hash_edge_ref );
			edge_add(elem->edge[9].id, hash_edge_ref );
			edge_add(elem->edge[11].id, hash_edge_ref );
		}else if(elem->pad==129){
			edge_add(elem->edge[1].id, hash_edge_ref );
			edge_add(elem->edge[2].id, hash_edge_ref );
			edge_add(elem->edge[3].id, hash_edge_ref );
			edge_add(elem->edge[6].id, hash_edge_ref );
			edge_add(elem->edge[7].id, hash_edge_ref );
			edge_add(elem->edge[9].id, hash_edge_ref );
			edge_add(elem->edge[10].id, hash_edge_ref );
			edge_add(elem->edge[11].id, hash_edge_ref );
		}else if(elem->pad==128){
			edge_add(elem->edge[0].id, hash_edge_ref );
			edge_add(elem->edge[1].id, hash_edge_ref );
			edge_add(elem->edge[3].id, hash_edge_ref );
			edge_add(elem->edge[5].id, hash_edge_ref );
			edge_add(elem->edge[8].id, hash_edge_ref );
			edge_add(elem->edge[9].id, hash_edge_ref );
			edge_add(elem->edge[9].id, hash_edge_ref );
			edge_add(elem->edge[11].id, hash_edge_ref );
		}else if(elem->pad==127){
			edge_add(elem->edge[4].id, hash_edge_ref );
			edge_add(elem->edge[5].id, hash_edge_ref );
			edge_add(elem->edge[6].id, hash_edge_ref );
			edge_add(elem->edge[7].id, hash_edge_ref );
			edge_add(elem->edge[8].id, hash_edge_ref );
			edge_add(elem->edge[9].id, hash_edge_ref );
			edge_add(elem->edge[10].id, hash_edge_ref );
			edge_add(elem->edge[11].id, hash_edge_ref );
		}else if(elem->pad==126){
			edge_add(elem->edge[0].id, hash_edge_ref );
			edge_add(elem->edge[1].id, hash_edge_ref );
			edge_add(elem->edge[2].id, hash_edge_ref );
			edge_add(elem->edge[3].id, hash_edge_ref );
			edge_add(elem->edge[4].id, hash_edge_ref );
			edge_add(elem->edge[5].id, hash_edge_ref );
			edge_add(elem->edge[6].id, hash_edge_ref );
			edge_add(elem->edge[7].id, hash_edge_ref );
		}
		//template 8
		else if(elem->pad==125){
			edge_add(elem->edge[1].id, hash_edge_ref );
			edge_add(elem->edge[3].id, hash_edge_ref );
			edge_add(elem->edge[4].id, hash_edge_ref );
			edge_add(elem->edge[5].id, hash_edge_ref );
			edge_add(elem->edge[6].id, hash_edge_ref );
			edge_add(elem->edge[7].id, hash_edge_ref );
			edge_add(elem->edge[9].id, hash_edge_ref );
			edge_add(elem->edge[11].id, hash_edge_ref );
		}else if(elem->pad==124){
			edge_add(elem->edge[0].id, hash_edge_ref );
			edge_add(elem->edge[2].id, hash_edge_ref );
			edge_add(elem->edge[4].id, hash_edge_ref );
			edge_add(elem->edge[5].id, hash_edge_ref );
			edge_add(elem->edge[6].id, hash_edge_ref );
			edge_add(elem->edge[7].id, hash_edge_ref );
			edge_add(elem->edge[8].id, hash_edge_ref );
			edge_add(elem->edge[10].id, hash_edge_ref );;
		}else if(elem->pad==123){
			edge_add(elem->edge[0].id, hash_edge_ref );
			edge_add(elem->edge[1].id, hash_edge_ref );
			edge_add(elem->edge[2].id, hash_edge_ref );
			edge_add(elem->edge[3].id, hash_edge_ref );
			edge_add(elem->edge[8].id, hash_edge_ref );
			edge_add(elem->edge[9].id, hash_edge_ref );
			edge_add(elem->edge[10].id, hash_edge_ref );
			edge_add(elem->edge[11].id, hash_edge_ref );
		}
		//template 7
		else if(elem->pad==122){
			edge_add(elem->edge[7].id, hash_edge_ref );
			edge_add(elem->edge[10].id, hash_edge_ref );
			edge_add(elem->edge[11].id, hash_edge_ref );
		}else if(elem->pad==121){
			edge_add(elem->edge[6].id, hash_edge_ref );
			edge_add(elem->edge[9].id, hash_edge_ref );
			edge_add(elem->edge[10].id, hash_edge_ref );
		}else if(elem->pad==120){
			edge_add(elem->edge[5].id, hash_edge_ref );
			edge_add(elem->edge[8].id, hash_edge_ref );
			edge_add(elem->edge[9].id, hash_edge_ref );
		}else if(elem->pad==119){
			edge_add(elem->edge[4].id, hash_edge_ref );
			edge_add(elem->edge[8].id, hash_edge_ref );
			edge_add(elem->edge[11].id, hash_edge_ref );
		}else if(elem->pad==118){
			edge_add(elem->edge[2].id, hash_edge_ref );
			edge_add(elem->edge[3].id, hash_edge_ref );
			edge_add(elem->edge[7].id, hash_edge_ref );
		}else if(elem->pad==117){
			edge_add(elem->edge[1].id, hash_edge_ref );
			edge_add(elem->edge[2].id, hash_edge_ref );
			edge_add(elem->edge[6].id, hash_edge_ref );
		}else if(elem->pad==116){
			edge_add(elem->edge[0].id, hash_edge_ref );
			edge_add(elem->edge[1].id, hash_edge_ref );
			edge_add(elem->edge[5].id, hash_edge_ref );
		}else if(elem->pad==115){
			edge_add(elem->edge[0].id, hash_edge_ref );
			edge_add(elem->edge[3].id, hash_edge_ref );
			edge_add(elem->edge[4].id, hash_edge_ref );
		}
		//template 6
		else if(elem->pad==114){
			edge_add(elem->edge[1].id, hash_edge_ref );
			edge_add(elem->edge[3].id, hash_edge_ref );
			edge_add(elem->edge[6].id, hash_edge_ref );
			edge_add(elem->edge[7].id, hash_edge_ref );
			edge_add(elem->edge[9].id, hash_edge_ref );
			edge_add(elem->edge[11].id, hash_edge_ref );
		}else if(elem->pad==113){
			edge_add(elem->edge[4].id, hash_edge_ref );
			edge_add(elem->edge[5].id, hash_edge_ref );
			edge_add(elem->edge[6].id, hash_edge_ref );
			edge_add(elem->edge[7].id, hash_edge_ref );
			edge_add(elem->edge[9].id, hash_edge_ref );
			edge_add(elem->edge[11].id, hash_edge_ref );
		}else if(elem->pad==112){
			edge_add(elem->edge[1].id, hash_edge_ref );
			edge_add(elem->edge[3].id, hash_edge_ref );
			edge_add(elem->edge[4].id, hash_edge_ref );
			edge_add(elem->edge[5].id, hash_edge_ref );
			edge_add(elem->edge[8].id, hash_edge_ref );
			edge_add(elem->edge[11].id, hash_edge_ref );
		}else if(elem->pad==111){
			edge_add(elem->edge[1].id, hash_edge_ref );
			edge_add(elem->edge[3].id, hash_edge_ref );
			edge_add(elem->edge[4].id, hash_edge_ref );
			edge_add(elem->edge[5].id, hash_edge_ref );
			edge_add(elem->edge[6].id, hash_edge_ref );
			edge_add(elem->edge[7].id, hash_edge_ref );
		}else if(elem->pad==110){
			edge_add(elem->edge[0].id, hash_edge_ref );
			edge_add(elem->edge[2].id, hash_edge_ref );
			edge_add(elem->edge[4].id, hash_edge_ref );
			edge_add(elem->edge[5].id, hash_edge_ref );
			edge_add(elem->edge[6].id, hash_edge_ref );
			edge_add(elem->edge[7].id, hash_edge_ref );
		}else if(elem->pad==109){
			edge_add(elem->edge[0].id, hash_edge_ref );
			edge_add(elem->edge[2].id, hash_edge_ref );
			edge_add(elem->edge[4].id, hash_edge_ref );
			edge_add(elem->edge[7].id, hash_edge_ref );
			edge_add(elem->edge[8].id, hash_edge_ref );
			edge_add(elem->edge[10].id, hash_edge_ref );
		}else if(elem->pad==108){
			edge_add(elem->edge[4].id, hash_edge_ref );
			edge_add(elem->edge[5].id, hash_edge_ref );
			edge_add(elem->edge[6].id, hash_edge_ref );
			edge_add(elem->edge[7].id, hash_edge_ref );
			edge_add(elem->edge[8].id, hash_edge_ref );
			edge_add(elem->edge[10].id, hash_edge_ref );
		}else if(elem->pad==107){
			edge_add(elem->edge[0].id, hash_edge_ref );
			edge_add(elem->edge[2].id, hash_edge_ref );
			edge_add(elem->edge[5].id, hash_edge_ref );
			edge_add(elem->edge[6].id, hash_edge_ref );
			edge_add(elem->edge[8].id, hash_edge_ref );
			edge_add(elem->edge[10].id, hash_edge_ref );
		}else if(elem->pad==106){
			edge_add(elem->edge[0].id, hash_edge_ref );
			edge_add(elem->edge[2].id, hash_edge_ref );
			edge_add(elem->edge[3].id, hash_edge_ref );
			edge_add(elem->edge[8].id, hash_edge_ref );
			edge_add(elem->edge[10].id, hash_edge_ref );
			edge_add(elem->edge[11].id, hash_edge_ref );
		}else if(elem->pad==105){
			edge_add(elem->edge[1].id, hash_edge_ref );
			edge_add(elem->edge[2].id, hash_edge_ref );
			edge_add(elem->edge[3].id, hash_edge_ref );
			edge_add(elem->edge[9].id, hash_edge_ref );
			edge_add(elem->edge[10].id, hash_edge_ref );
			edge_add(elem->edge[11].id, hash_edge_ref );;
		}else if(elem->pad==104){
			edge_add(elem->edge[0].id, hash_edge_ref );
			edge_add(elem->edge[1].id, hash_edge_ref );
			edge_add(elem->edge[2].id, hash_edge_ref );
			edge_add(elem->edge[8].id, hash_edge_ref );
			edge_add(elem->edge[9].id, hash_edge_ref );
			edge_add(elem->edge[10].id, hash_edge_ref );
		}else if(elem->pad==103){
			edge_add(elem->edge[0].id, hash_edge_ref );
			edge_add(elem->edge[1].id, hash_edge_ref );
			edge_add(elem->edge[3].id, hash_edge_ref );
			edge_add(elem->edge[8].id, hash_edge_ref );
			edge_add(elem->edge[9].id, hash_edge_ref );
			edge_add(elem->edge[11].id, hash_edge_ref );
		}
		//template 5
		else if(elem->pad==102){
			edge_add(elem->edge[3].id, hash_edge_ref );
			edge_add(elem->edge[7].id, hash_edge_ref );
		}else if(elem->pad==101){
			edge_add(elem->edge[6].id, hash_edge_ref );
			edge_add(elem->edge[7].id, hash_edge_ref );
		}else if(elem->pad==100){
			edge_add(elem->edge[4].id, hash_edge_ref );
			edge_add(elem->edge[6].id, hash_edge_ref );
		}else if(elem->pad==99){
			edge_add(elem->edge[3].id, hash_edge_ref );
			edge_add(elem->edge[4].id, hash_edge_ref );
		}else if(elem->pad==98){
			edge_add(elem->edge[2].id, hash_edge_ref );
			edge_add(elem->edge[7].id, hash_edge_ref );
		}else if(elem->pad==97){
			edge_add(elem->edge[7].id, hash_edge_ref );
			edge_add(elem->edge[10].id, hash_edge_ref );
		}else if(elem->pad==96){
			edge_add(elem->edge[6].id, hash_edge_ref );
			edge_add(elem->edge[10].id, hash_edge_ref );
		}else if(elem->pad==95){
			edge_add(elem->edge[2].id, hash_edge_ref );
			edge_add(elem->edge[6].id, hash_edge_ref );;
		}else if(elem->pad==94){
			edge_add(elem->edge[1].id, hash_edge_ref );
			edge_add(elem->edge[5].id, hash_edge_ref );
		}else if(elem->pad==93){
			edge_add(elem->edge[1].id, hash_edge_ref );
			edge_add(elem->edge[6].id, hash_edge_ref );
		}else if(elem->pad==92){
			edge_add(elem->edge[6].id, hash_edge_ref );
			edge_add(elem->edge[9].id, hash_edge_ref );
		}else if(elem->pad==91){
			edge_add(elem->edge[5].id, hash_edge_ref );
			edge_add(elem->edge[9].id, hash_edge_ref );
		}else if(elem->pad==90){
			edge_add(elem->edge[0].id, hash_edge_ref );
			edge_add(elem->edge[4].id, hash_edge_ref );
		}else if(elem->pad==89){
			edge_add(elem->edge[0].id, hash_edge_ref );
			edge_add(elem->edge[5].id, hash_edge_ref );
		}else if(elem->pad==88){
			edge_add(elem->edge[5].id, hash_edge_ref );
			edge_add(elem->edge[8].id, hash_edge_ref );;
		}else if(elem->pad==87){
			edge_add(elem->edge[8].id, hash_edge_ref );
			edge_add(elem->edge[4].id, hash_edge_ref );
		}else if(elem->pad==86){
			edge_add(elem->edge[8].id, hash_edge_ref );
			edge_add(elem->edge[11].id, hash_edge_ref );
		}else if(elem->pad==85){
			edge_add(elem->edge[10].id, hash_edge_ref );
			edge_add(elem->edge[11].id, hash_edge_ref );
		}else if(elem->pad==84){
			edge_add(elem->edge[9].id, hash_edge_ref );
			edge_add(elem->edge[10].id, hash_edge_ref );
		}else if(elem->pad==83){
			edge_add(elem->edge[8].id, hash_edge_ref );
			edge_add(elem->edge[9].id, hash_edge_ref );
		}else if(elem->pad==82){
			edge_add(elem->edge[0].id, hash_edge_ref );
			edge_add(elem->edge[3].id, hash_edge_ref );
		}else if(elem->pad==81){
			edge_add(elem->edge[2].id, hash_edge_ref );
			edge_add(elem->edge[3].id, hash_edge_ref );
		}else if(elem->pad==80){
			edge_add(elem->edge[1].id, hash_edge_ref );
			edge_add(elem->edge[2].id, hash_edge_ref );
		}else if(elem->pad==79){
			edge_add(elem->edge[0].id, hash_edge_ref );
			edge_add(elem->edge[1].id, hash_edge_ref );
		}else if(elem->pad==78){
			edge_add(elem->edge[3].id, hash_edge_ref );
			edge_add(elem->edge[7].id, hash_edge_ref );
			edge_add(elem->edge[11].id, hash_edge_ref );
		}else if(elem->pad==77){
			edge_add(elem->edge[4].id, hash_edge_ref );
			edge_add(elem->edge[7].id, hash_edge_ref );
			edge_add(elem->edge[11].id, hash_edge_ref );
		}else if(elem->pad==76){
			edge_add(elem->edge[3].id, hash_edge_ref );
			edge_add(elem->edge[4].id, hash_edge_ref );
			edge_add(elem->edge[11].id, hash_edge_ref );
		}else if(elem->pad==75){
			edge_add(elem->edge[3].id, hash_edge_ref );
			edge_add(elem->edge[4].id, hash_edge_ref );
			edge_add(elem->edge[7].id, hash_edge_ref );
		}else if(elem->pad==74){
			edge_add(elem->edge[2].id, hash_edge_ref );
			edge_add(elem->edge[7].id, hash_edge_ref );
			edge_add(elem->edge[10].id, hash_edge_ref );
		}else if(elem->pad==73){
			edge_add(elem->edge[6].id, hash_edge_ref );
			edge_add(elem->edge[7].id, hash_edge_ref );
			edge_add(elem->edge[10].id, hash_edge_ref );
		}else if(elem->pad==72){
			edge_add(elem->edge[2].id, hash_edge_ref );
			edge_add(elem->edge[6].id, hash_edge_ref );
			edge_add(elem->edge[10].id, hash_edge_ref );
		}else if(elem->pad==71){
			edge_add(elem->edge[2].id, hash_edge_ref );
			edge_add(elem->edge[6].id, hash_edge_ref );
			edge_add(elem->edge[7].id, hash_edge_ref );
		}else if(elem->pad==70){
			edge_add(elem->edge[1].id, hash_edge_ref );
			edge_add(elem->edge[6].id, hash_edge_ref );
			edge_add(elem->edge[9].id, hash_edge_ref );
		}else if(elem->pad==69){
			edge_add(elem->edge[5].id, hash_edge_ref );
			edge_add(elem->edge[6].id, hash_edge_ref );
			edge_add(elem->edge[9].id, hash_edge_ref );
		}else if(elem->pad==68){
			edge_add(elem->edge[1].id, hash_edge_ref );
			edge_add(elem->edge[5].id, hash_edge_ref );
			edge_add(elem->edge[9].id, hash_edge_ref );
		}else if(elem->pad==67){
			edge_add(elem->edge[1].id, hash_edge_ref );
			edge_add(elem->edge[5].id, hash_edge_ref );
			edge_add(elem->edge[6].id, hash_edge_ref );
		}else if(elem->pad==66){
			edge_add(elem->edge[0].id, hash_edge_ref );
			edge_add(elem->edge[4].id, hash_edge_ref );
			edge_add(elem->edge[8].id, hash_edge_ref );
		}else if(elem->pad==65){
			edge_add(elem->edge[4].id, hash_edge_ref );
			edge_add(elem->edge[5].id, hash_edge_ref );
			edge_add(elem->edge[8].id, hash_edge_ref );
		}else if(elem->pad==64){
			edge_add(elem->edge[0].id, hash_edge_ref );
			edge_add(elem->edge[5].id, hash_edge_ref );
			edge_add(elem->edge[8].id, hash_edge_ref );
		}else if(elem->pad==63){
			edge_add(elem->edge[0].id, hash_edge_ref );
			edge_add(elem->edge[4].id, hash_edge_ref );
			edge_add(elem->edge[5].id, hash_edge_ref );
		}else if(elem->pad==62){
			edge_add(elem->edge[8].id, hash_edge_ref );
			edge_add(elem->edge[10].id, hash_edge_ref );
			edge_add(elem->edge[11].id, hash_edge_ref );
		}else if(elem->pad==61){
			edge_add(elem->edge[9].id, hash_edge_ref );
			edge_add(elem->edge[10].id, hash_edge_ref );
			edge_add(elem->edge[11].id, hash_edge_ref );
		}else if(elem->pad==60){
			edge_add(elem->edge[8].id, hash_edge_ref );
			edge_add(elem->edge[9].id, hash_edge_ref );
			edge_add(elem->edge[11].id, hash_edge_ref );
		}else if(elem->pad==59){
			edge_add(elem->edge[8].id, hash_edge_ref );
			edge_add(elem->edge[9].id, hash_edge_ref );
			edge_add(elem->edge[10].id, hash_edge_ref );
		}else if(elem->pad==58){
			edge_add(elem->edge[0].id, hash_edge_ref );
			edge_add(elem->edge[2].id, hash_edge_ref );
			edge_add(elem->edge[11].id, hash_edge_ref );
		}else if(elem->pad==57){
			edge_add(elem->edge[1].id, hash_edge_ref );
			edge_add(elem->edge[2].id, hash_edge_ref );
			edge_add(elem->edge[3].id, hash_edge_ref );
		}else if(elem->pad==56){
			edge_add(elem->edge[0].id, hash_edge_ref );
			edge_add(elem->edge[1].id, hash_edge_ref );
			edge_add(elem->edge[2].id, hash_edge_ref );
		}else if(elem->pad==55){
			edge_add(elem->edge[0].id, hash_edge_ref );
			edge_add(elem->edge[1].id, hash_edge_ref );
			edge_add(elem->edge[3].id, hash_edge_ref );
		}else if(elem->pad==54){
			edge_add(elem->edge[3].id, hash_edge_ref );
			edge_add(elem->edge[4].id, hash_edge_ref );
			edge_add(elem->edge[7].id, hash_edge_ref );
			edge_add(elem->edge[11].id, hash_edge_ref );
		}else if(elem->pad==53){
			edge_add(elem->edge[2].id, hash_edge_ref );
			edge_add(elem->edge[6].id, hash_edge_ref );
			edge_add(elem->edge[7].id, hash_edge_ref );
			edge_add(elem->edge[10].id, hash_edge_ref );
		}else if(elem->pad==52){
			edge_add(elem->edge[1].id, hash_edge_ref );
			edge_add(elem->edge[5].id, hash_edge_ref );
			edge_add(elem->edge[6].id, hash_edge_ref );
			edge_add(elem->edge[9].id, hash_edge_ref );
		}else if(elem->pad==51){
			edge_add(elem->edge[0].id, hash_edge_ref );
			edge_add(elem->edge[4].id, hash_edge_ref );
			edge_add(elem->edge[5].id, hash_edge_ref );
			edge_add(elem->edge[8].id, hash_edge_ref );
		}else if(elem->pad==50){
			edge_add(elem->edge[8].id, hash_edge_ref );
			edge_add(elem->edge[9].id, hash_edge_ref );
			edge_add(elem->edge[10].id, hash_edge_ref );
			edge_add(elem->edge[11].id, hash_edge_ref );
		}else if(elem->pad==49){
			edge_add(elem->edge[0].id, hash_edge_ref );
			edge_add(elem->edge[1].id, hash_edge_ref );
			edge_add(elem->edge[2].id, hash_edge_ref );
			edge_add(elem->edge[3].id, hash_edge_ref );
		}
		//template 4
		else if(elem->pad==48){
			edge_add(elem->edge[4].id, hash_edge_ref );
			edge_add(elem->edge[5].id, hash_edge_ref );
			edge_add(elem->edge[9].id, hash_edge_ref );
			edge_add(elem->edge[11].id, hash_edge_ref );
		}else if(elem->pad==47){
			edge_add(elem->edge[6].id, hash_edge_ref );
			edge_add(elem->edge[7].id, hash_edge_ref );
			edge_add(elem->edge[9].id, hash_edge_ref );
			edge_add(elem->edge[11].id, hash_edge_ref );
		}else if(elem->pad==46){
			edge_add(elem->edge[1].id, hash_edge_ref );
			edge_add(elem->edge[3].id, hash_edge_ref );
			edge_add(elem->edge[6].id, hash_edge_ref );
			edge_add(elem->edge[7].id, hash_edge_ref );
		}else if(elem->pad==45){
			edge_add(elem->edge[1].id, hash_edge_ref );
			edge_add(elem->edge[3].id, hash_edge_ref );
			edge_add(elem->edge[4].id, hash_edge_ref );
			edge_add(elem->edge[5].id, hash_edge_ref );
		}else if(elem->pad==44){
			edge_add(elem->edge[4].id, hash_edge_ref );
			edge_add(elem->edge[7].id, hash_edge_ref );
			edge_add(elem->edge[8].id, hash_edge_ref );
			edge_add(elem->edge[10].id, hash_edge_ref );
		}else if(elem->pad==43){
			edge_add(elem->edge[5].id, hash_edge_ref );
			edge_add(elem->edge[6].id, hash_edge_ref );
			edge_add(elem->edge[8].id, hash_edge_ref );
			edge_add(elem->edge[10].id, hash_edge_ref );
		}else if(elem->pad==42){
			edge_add(elem->edge[0].id, hash_edge_ref );
			edge_add(elem->edge[2].id, hash_edge_ref );
			edge_add(elem->edge[5].id, hash_edge_ref );
			edge_add(elem->edge[6].id, hash_edge_ref );
		}else if(elem->pad==41){
			edge_add(elem->edge[0].id, hash_edge_ref );
			edge_add(elem->edge[2].id, hash_edge_ref );
			edge_add(elem->edge[4].id, hash_edge_ref );
			edge_add(elem->edge[7].id, hash_edge_ref );
		}else if(elem->pad==40){
			edge_add(elem->edge[0].id, hash_edge_ref );
			edge_add(elem->edge[3].id, hash_edge_ref );
			edge_add(elem->edge[8].id, hash_edge_ref );
			edge_add(elem->edge[11].id, hash_edge_ref );
		}else if(elem->pad==39){
			edge_add(elem->edge[2].id, hash_edge_ref );
			edge_add(elem->edge[3].id, hash_edge_ref );
			edge_add(elem->edge[10].id, hash_edge_ref );
			edge_add(elem->edge[11].id, hash_edge_ref );
		}else if(elem->pad==38){
			edge_add(elem->edge[1].id, hash_edge_ref );
			edge_add(elem->edge[2].id, hash_edge_ref );
			edge_add(elem->edge[9].id, hash_edge_ref );
			edge_add(elem->edge[10].id, hash_edge_ref );
		}else if(elem->pad==37){
			edge_add(elem->edge[0].id, hash_edge_ref );
			edge_add(elem->edge[1].id, hash_edge_ref );
			edge_add(elem->edge[8].id, hash_edge_ref );
			edge_add(elem->edge[9].id, hash_edge_ref );
		}
		//template 3
		else if(elem->pad==36){
			edge_add(elem->edge[3].id, hash_edge_ref );
			edge_add(elem->edge[11].id, hash_edge_ref );
		}else if(elem->pad==35){
			edge_add(elem->edge[1].id, hash_edge_ref );
			edge_add(elem->edge[10].id, hash_edge_ref );
		}else if(elem->pad==34){
			edge_add(elem->edge[1].id, hash_edge_ref );
			edge_add(elem->edge[9].id, hash_edge_ref );
		}else if(elem->pad==33){
			edge_add(elem->edge[0].id, hash_edge_ref );
			edge_add(elem->edge[8].id, hash_edge_ref );
		}else if(elem->pad==32){
			edge_add(elem->edge[4].id, hash_edge_ref );
			edge_add(elem->edge[7].id, hash_edge_ref );
		}else if(elem->pad==31){
			edge_add(elem->edge[6].id, hash_edge_ref );
			edge_add(elem->edge[7].id, hash_edge_ref );
		}else if(elem->pad==30){
			edge_add(elem->edge[5].id, hash_edge_ref );
			edge_add(elem->edge[6].id, hash_edge_ref );
		}else if(elem->pad==29){
			edge_add(elem->edge[4].id, hash_edge_ref );
			edge_add(elem->edge[5].id, hash_edge_ref );
		}else if(elem->pad==28){
			edge_add(elem->edge[9].id, hash_edge_ref );
			edge_add(elem->edge[11].id, hash_edge_ref );
		}else if(elem->pad==27){
			edge_add(elem->edge[8].id, hash_edge_ref );
			edge_add(elem->edge[10].id, hash_edge_ref );
		}else if(elem->pad==26){
			edge_add(elem->edge[1].id, hash_edge_ref );
			edge_add(elem->edge[3].id, hash_edge_ref );
		}else if(elem->pad==25){
			edge_add(elem->edge[0].id, hash_edge_ref );
			edge_add(elem->edge[2].id, hash_edge_ref );
		}
		//template 2
		else if(elem->pad==24){
			edge_add(elem->edge[4].id, hash_edge_ref );
			edge_add(elem->edge[5].id, hash_edge_ref );
			edge_add(elem->edge[6].id, hash_edge_ref );
			edge_add(elem->edge[7].id, hash_edge_ref );
		}else if(elem->pad==23){
			edge_add(elem->edge[1].id, hash_edge_ref );
			edge_add(elem->edge[3].id, hash_edge_ref );
			edge_add(elem->edge[9].id, hash_edge_ref );
			edge_add(elem->edge[11].id, hash_edge_ref );
		}else if(elem->pad==22){
			edge_add(elem->edge[0].id, hash_edge_ref );
			edge_add(elem->edge[2].id, hash_edge_ref );
			edge_add(elem->edge[8].id, hash_edge_ref );
			edge_add(elem->edge[10].id, hash_edge_ref );
		}
		//template 1
		else if(elem->pad==21){
			edge_add(elem->edge[11].id, hash_edge_ref );
		}else if(elem->pad==20){
			edge_add(elem->edge[10].id, hash_edge_ref );
		}else if(elem->pad==19){
			edge_add(elem->edge[9].id, hash_edge_ref );
		}else if(elem->pad==18){
			edge_add(elem->edge[8].id, hash_edge_ref );
		}else if(elem->pad==17){
			edge_add(elem->edge[7].id, hash_edge_ref );
		}else if(elem->pad==16){
			edge_add(elem->edge[6].id, hash_edge_ref );
		}else if(elem->pad==15){
			edge_add(elem->edge[5].id, hash_edge_ref );
		}else if(elem->pad==14){
			edge_add(elem->edge[4].id, hash_edge_ref );
		}else if(elem->pad==13){
			edge_add(elem->edge[3].id, hash_edge_ref );
		}else if(elem->pad==12){
			edge_add(elem->edge[2].id, hash_edge_ref );
		}else if(elem->pad==11){
			edge_add(elem->edge[1].id, hash_edge_ref );
		}else if(elem->pad==10){
			edge_add(elem->edge[0].id, hash_edge_ref );
		}
	}

	size_t position;
	for (int iedge = 0; iedge< mesh->shared_edges.elem_count; iedge++){

		octant_edge_t* shared_ed = (octant_edge_t*) sc_array_index(&mesh->shared_edges, iedge);

		bool out =  sc_hash_array_lookup(hash_edge_ref, &shared_ed->id, &position);
		shared_ed->ref = out;
	}

#ifdef HEXA_DEBUG_
	if(0){
		fprintf(mesh->fdbg ,"shared edges status\n");
		fprintf(mesh->fdbg ,"shared edges number:%d\n",mesh->shared_edges.elem_count);
		for (int iedge = 0; iedge< mesh->shared_edges.elem_count; iedge++){
			octant_edge_t* shared_ed = (octant_edge_t*) sc_array_index(&mesh->shared_edges, iedge);
			fprintf(mesh->fdbg ,"id: %d  ref:%d  rank:%d\n",shared_ed->id,shared_ed->ref,mesh->mpi_rank);
		}
	}
#endif

}

void Edge_comunication(hexa_tree_t* mesh, std::vector<int>& elements_ids, sc_hash_array_t* hash_edge_ref){

	size_t position;
	bool clamped=true;
	MPI_Request *requests;
	MPI_Status  *statuses;
	long long   *recvbuf;
	long long   *sendbuf;
	int n_requests;

	n_requests = mesh->comm_map_edge.nrequests;
	recvbuf    = (long long*)malloc(2*mesh->comm_map_edge.max_recvbuf_size*sizeof(long long));
	sendbuf    = (long long*)malloc(2*mesh->comm_map_edge.max_sendbuf_size*sizeof(long long));

	requests = (MPI_Request*) malloc (n_requests*sizeof(MPI_Request));
	statuses = (MPI_Status*)  malloc (n_requests*sizeof(MPI_Status));

	// rank[n+1] send to rank[n]
	int c = 0;

	int max_recvbuf_size = 2*mesh->comm_map_edge.max_recvbuf_size;
	int max_sendbuf_size = 2*mesh->comm_map_edge.max_sendbuf_size;

	int offset = 0;

	// post all non-blocking receives
	for(int i = 0; i < mesh->comm_map_edge.RecvFrom.elem_count; ++i) {
		message_t *m = (message_t*) sc_array_index(&mesh->comm_map_edge.RecvFrom, i);
		MPI_Irecv(&recvbuf[offset], 2*m->idxs.elem_count, MPI_LONG_LONG, m->rank,0,MPI_COMM_WORLD, &requests[c]);
		offset += 2*m->idxs.elem_count;
		c++;
	}
	assert(offset == max_recvbuf_size);


	offset = 0;
	for(int i=0; i < mesh->comm_map_edge.SendTo.elem_count; i++){
		message_t* m = (message_t*) sc_array_index(&mesh->comm_map_edge.SendTo, i);
		for(int j = 0; j< m->idxs.elem_count;j++){
			long long* mm = (long long*) sc_array_index(&m->idxs, j);

			octant_edge_t key;
			key.id = *mm;
			long long out = sc_hash_array_lookup(hash_edge_ref, &key, &position);
			sendbuf[2*j+offset] = (long long) *mm;
			sendbuf[2*j+offset+1] = (long long) out;
		}
		MPI_Isend(&sendbuf[offset], 2*m->idxs.elem_count, MPI_LONG_LONG, m->rank,0,MPI_COMM_WORLD, &requests[c]);
		offset += 2*m->idxs.elem_count;
		c++;
	}

	assert(offset == max_sendbuf_size);

	assert(c == n_requests);

	MPI_Waitall(n_requests,requests,statuses);


	offset =0;
	//spread the information
	for(int i = 0; i < mesh->comm_map_edge.RecvFrom.elem_count; ++i) {
		message_t *m = (message_t*) sc_array_index(&mesh->comm_map_edge.RecvFrom, i);
		for(int j = 0; j < m->idxs.elem_count; ++j){
			if(recvbuf[2*j+offset+1]){
				edge_add(recvbuf[2*j+offset], hash_edge_ref );
			}
		}
		offset += 2*m->idxs.elem_count;
	}

	free(recvbuf);
	free(sendbuf);
	free(requests);
	free(statuses);

	// rank[n] send to rank[n+1]
	recvbuf    = (long long*)malloc(2*mesh->comm_map_edge.max_sendbuf_size*sizeof(long long));
	sendbuf    = (long long*)malloc(2*mesh->comm_map_edge.max_recvbuf_size*sizeof(long long));

	requests = (MPI_Request*) malloc (n_requests*sizeof(MPI_Request));
	statuses = (MPI_Status*)  malloc (n_requests*sizeof(MPI_Status));

	max_recvbuf_size = 2*mesh->comm_map_edge.max_sendbuf_size;
	max_sendbuf_size = 2*mesh->comm_map_edge.max_recvbuf_size;

	c = 0;
	offset = 0;

	// post all non-blocking receives
	for(int i = 0; i < mesh->comm_map_edge.SendTo.elem_count; ++i) {
		message_t *m = (message_t*) sc_array_index(&mesh->comm_map_edge.SendTo, i);
		MPI_Irecv(&recvbuf[offset], 2*m->idxs.elem_count, MPI_LONG_LONG, m->rank,0,MPI_COMM_WORLD, &requests[c]);
		offset += 2*m->idxs.elem_count;
		c++;
	}
	assert(offset == max_recvbuf_size);

	offset = 0;
	for(int i=0; i < mesh->comm_map_edge.RecvFrom.elem_count; i++){
		message_t* m = (message_t*) sc_array_index(&mesh->comm_map_edge.RecvFrom, i);
		for(int j = 0; j< m->idxs.elem_count;j++){
			long long* mm = (long long*) sc_array_index(&m->idxs, j);
			octant_edge_t key;
			key.id = *mm;
			int out =  sc_hash_array_lookup(hash_edge_ref, &key, &position);
			sendbuf[2*j+offset] = (long long) *mm;
			sendbuf[2*j+offset+1] = (long long)   out;
		}
		MPI_Isend(&sendbuf[offset], 2*m->idxs.elem_count, MPI_LONG_LONG, m->rank,0,MPI_COMM_WORLD, &requests[c]);
		offset += 2*m->idxs.elem_count;
		c++;
	}

	assert(offset == max_sendbuf_size);

	assert(c == n_requests);

	MPI_Waitall(n_requests,requests,statuses);

	offset =0;
	//add the information in the hash_edge_ref
	for(int i = 0; i < mesh->comm_map_edge.SendTo.elem_count; ++i) {
		message_t *m = (message_t*) sc_array_index(&mesh->comm_map_edge.SendTo, i);
		for(int j = 0; j < m->idxs.elem_count; ++j){
			if(recvbuf[2*j+offset+1]){
				edge_add(recvbuf[2*j+offset], hash_edge_ref );
			}
		}
		offset += 2*m->idxs.elem_count;
	}

	free(&recvbuf[0]);
	free(&sendbuf[0]);
	free(requests);
	free(statuses);

}

void Edge_propagationOld(hexa_tree_t* mesh, std::vector<int>& elements_ids, sc_hash_array_t* hash_edge_ref) {

	size_t position;

	for(int iel= 0; iel < mesh->elements.elem_count; iel++ ){
		octant_t *elem = (octant_t*) sc_array_index(&mesh->elements, iel);

		for (int edge = 0; edge < 12; ++edge){
			bool out = false;
			out =  sc_hash_array_lookup(hash_edge_ref, &elem->edge[edge].id, &position);

			if(out){
				elem->edge[edge].ref=true;
				elem->pad = -1;
				elements_ids.push_back(iel);
			}
		}
	}
	//cleaning the element vector
	std::sort( elements_ids.begin(), elements_ids.end() );
	elements_ids.erase( std::unique( elements_ids.begin(), elements_ids.end() ), elements_ids.end() );
}

void Edge_propagation(hexa_tree_t* mesh, std::vector<int>& elements_ids, sc_hash_array_t* hash_edge_ref, sc_hash_array_t* hash_edge, sc_hash_array_t* hash_elem) {

	/*
	for(int iedgeref = 0; iedgeref < hash_edge_ref->a.elem_count; iedgeref++){
		octant_edge_t *edgeref = (octant_edge_t*) sc_array_index(&hash_edge_ref->a, iedgeref);

		size_t position;
		edge_t key;

		key.id = edgeref->id;
		edgeref->ref = true;
		bool out =  sc_hash_array_lookup(hash_edge, &key, &position);
		if(out){
			edge_t *edge = (edge_t*) sc_array_index(&hash_edge->a, position);
			edge->ref = true;
			for(int iel = 0; iel< edge->list_elem; iel++){
				octant_t *elem = (octant_t*) sc_array_index(&mesh->elements, edge->elem[iel]);

				for(int iedge = 0; iedge < 12; iedge++){
					if(elem->edge[iedge].id == edgeref->id){
						elem->edge[iedge].ref = true;
					}
				}

				elem->pad = -1;
				octant_t *r;
				key.id = elem->id;
				r = (octant_t*) sc_hash_array_insert_unique(hash_elem, &key, &position);
				if(r!=NULL){
					r->id = elem->id;
				}
			}
		}
	}


	//elements_ids.clear();
	for(int iel= elements_ids.size(); iel < hash_elem->a.elem_count; iel++ ){
		octant_t *elem = (octant_t*) sc_array_index(&hash_elem->a, iel);
		elements_ids.push_back(elem->id);
	}
	 */
	size_t position;
	octant_t key;

	for(int iel= 0; iel < mesh->elements.elem_count; iel++ ){
		octant_t *elem = (octant_t*) sc_array_index(&mesh->elements, iel);

		for (int edge = 0; edge < 12; ++edge){
			bool out = false;
			out =  sc_hash_array_lookup(hash_edge_ref, &elem->edge[edge].id, &position);

			if(out){
				octant_t *r;
				key.id = elem->id;
				r = (octant_t*) sc_hash_array_insert_unique(hash_elem, &key, &position);
				if(r!=NULL){
					r->id = elem->id;
				}
				elem->edge[edge].ref=true;
				elem->pad = -1;
				//elements_ids.push_back(iel);
			}
		}
	}

	//elements_ids.clear();
	for(int iel= elements_ids.size(); iel < hash_elem->a.elem_count; iel++ ){
		octant_t *elem = (octant_t*) sc_array_index(&hash_elem->a, iel);
		elements_ids.push_back(elem->id);
	}

}

void CheckOctreeTemplate(hexa_tree_t* mesh,sc_hash_array_t*hash_edge_ref) {

	bool clamped = true;
	//criando a hash de elementos afetados para evitar nos duplicados
	sc_hash_array_t*   hash_elem  = (sc_hash_array_t *)sc_hash_array_new(sizeof(octant_t), el_hash_id, el_equal_id, &clamped);
	size_t position;
	octant_t key;
	std::vector<int> elements_ids;

	for(int ioc = 0; ioc <mesh->oct.elem_count; ioc++){
		octree_t * oc = (octree_t*) sc_array_index (&mesh->oct, ioc);
		for(int iel = 0; iel<8; iel++){
			elements_ids.push_back(oc->id[iel]);
			octant_t *r;
			key.id = elements_ids[iel];
			r = (octant_t*) sc_hash_array_insert_unique(hash_elem, &key, &position);
			if(r!=NULL){
				r->id = elements_ids[iel];
			}else{
				printf("Verificar o no numero %d\n",elements_ids[iel]);
			}
		}
	}

	//auto start = std::chrono::steady_clock::now( );



	for(int iel= 0; iel < elements_ids.size(); iel++ ){

	}
	//auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>( std::chrono::steady_clock::now( ) - start );
	//std::cout << "Time in hash elements_ids "<< elapsed.count() <<" millisecond(s)."<< std::endl;

	//TODO trocar o hashid e fn...
	sc_hash_array_t* hash_edge  = (sc_hash_array_t *)sc_hash_array_new(sizeof(edge_t), edget_id_hash, edget_id_equal, &clamped);
	/////////////////
	// create the edge structure
	//mesh->elements.elem_count
	for (int iel = 0; iel < 0; ++iel) {

		octant_t *elem = (octant_t*) sc_array_index(&mesh->elements, iel);

		for (int iedge = 0; iedge < 12; iedge++){
			edge_t key;
			key.id = elem->edge[iedge].id;
			edge_t* edge = (edge_t*) sc_hash_array_insert_unique (hash_edge, &key, &position);
			if(edge != NULL){
				edge->id = elem->edge[iedge].id;
				edge->list_elem = 1;
				edge->elem[edge->list_elem-1] = elem->id;
				edge->ref = elem->edge[iedge].ref;
			}else{
				edge = (edge_t*) sc_array_index(&hash_edge->a, position);
				edge->elem[edge->list_elem] = elem->id;
				edge->list_elem++;
			}
		}
	}


	int iter_count = 0;
	int diff = 50;
	while(iter_count < 1000 && diff != 0){
		printf("Numero de elementos:%d numero de arestas:%d\n",elements_ids.size(),hash_edge_ref->a.elem_count);
		auto start = std::chrono::steady_clock::now( );
		int edgecount = hash_edge_ref->a.elem_count;
		IdentifyTemplate(mesh, elements_ids, hash_edge_ref);
		auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>( std::chrono::steady_clock::now( ) - start );
		std::cout << "Time in IdentifyTemplate "<< elapsed.count() <<" millisecond(s)."<< std::endl;

		start = std::chrono::steady_clock::now( );
		Edge_identification( mesh, elements_ids, hash_edge_ref);
		elapsed = std::chrono::duration_cast<std::chrono::milliseconds>( std::chrono::steady_clock::now( ) - start );
		std::cout << "Time in Edge_identification "<< elapsed.count() <<" millisecond(s)."<< std::endl;
		start = std::chrono::steady_clock::now( );
		Edge_propagation (mesh, elements_ids, hash_edge_ref,hash_edge,hash_elem);
		elapsed = std::chrono::duration_cast<std::chrono::milliseconds>( std::chrono::steady_clock::now( ) - start );
		std::cout << "Time in Edge_propagation "<< elapsed.count() <<" millisecond(s)."<< std::endl;
		start = std::chrono::steady_clock::now( );
		Edge_comunication(mesh, elements_ids, hash_edge_ref);
		elapsed = std::chrono::duration_cast<std::chrono::milliseconds>( std::chrono::steady_clock::now( ) - start );
		std::cout << "Time in Edge_comunication "<< elapsed.count() <<" millisecond(s)."<< std::endl;
		diff = hash_edge_ref->a.elem_count - edgecount;

		if((iter_count % 20) == 0){
			printf("         Iteration number:%d; diff equals to:%d; number of edges in the hash:%d\n",iter_count,diff,hash_edge_ref->a.elem_count);
		}

		if(diff == 0){
			printf("         %d iteractions to propagate the edge contamination\n",iter_count);
		}
		iter_count++;
	}

	IdentifyTemplate(mesh, elements_ids, hash_edge_ref);

	//sc_hash_array_rip(hash_edge_ref,&mesh->edges_ref);

}

vector<int> RotateHex(int* rot, int* sym){

	vector<int> order;
	int aux[8];

	for(int k = 0;k<8;k++){
		order.push_back(k);
	}

	if(rot[0]==-1){

		for(int i=0;i<8;i++){
			aux[i]=order[i];
		}

		//order = {3,2,6,7,0,1,5,4};
		order.clear();

		order.push_back(aux[3]);
		order.push_back(aux[2]);
		order.push_back(aux[6]);
		order.push_back(aux[7]);

		order.push_back(aux[0]);
		order.push_back(aux[1]);
		order.push_back(aux[5]);
		order.push_back(aux[4]);

	}

	if(rot[0]==1){

		for(int i=0;i<8;i++){
			aux[i]=order[i];
		}

		//order = {4,0,7,3,5,1,2,6};
		order.clear();

		order.push_back(aux[4]);
		order.push_back(aux[5]);
		order.push_back(aux[1]);
		order.push_back(aux[0]);

		order.push_back(aux[7]);
		order.push_back(aux[6]);
		order.push_back(aux[2]);
		order.push_back(aux[3]);

	}

	if(rot[1]==1){

		for(int i=0;i<8;i++){
			aux[i]=order[i];
		}

		//order = {1,5,6,2,0,4,7,3};
		order.clear();

		order.push_back(aux[1]);
		order.push_back(aux[5]);
		order.push_back(aux[6]);
		order.push_back(aux[2]);

		order.push_back(aux[0]);
		order.push_back(aux[4]);
		order.push_back(aux[7]);
		order.push_back(aux[3]);


	}

	if(rot[1]==-1){

		for(int i=0;i<8;i++){
			aux[i]=order[i];
		}

		//order = {4,5,1,0,7,6,2,3};
		order.clear();

		order.push_back(aux[4]);
		order.push_back(aux[0]);
		order.push_back(aux[3]);
		order.push_back(aux[7]);

		order.push_back(aux[5]);
		order.push_back(aux[1]);
		order.push_back(aux[2]);
		order.push_back(aux[6]);

	}

	if(rot[2]==1){

		for(int i=0;i<8;i++){
			aux[i]=order[i];
		}

		//order = {3,0,1,2,7,4,5,6};
		order.clear();

		order.push_back(aux[3]);
		order.push_back(aux[0]);
		order.push_back(aux[1]);
		order.push_back(aux[2]);

		order.push_back(aux[7]);
		order.push_back(aux[4]);
		order.push_back(aux[5]);
		order.push_back(aux[6]);

	}

	if(rot[2]==-1){

		for(int i=0;i<8;i++){
			aux[i]=order[i];
		}

		//order = {1,2,3,0,5,6,7,4};
		order.clear();

		order.push_back(aux[1]);
		order.push_back(aux[2]);
		order.push_back(aux[3]);
		order.push_back(aux[0]);

		order.push_back(aux[5]);
		order.push_back(aux[6]);
		order.push_back(aux[7]);
		order.push_back(aux[4]);

	}

	if(sym[0]==1){

		for(int i=0;i<8;i++){
			aux[i]=order[i];
		}

		//order = {aux[1],aux[0],aux[3],aux[2],aux[5],aux[4],aux[7],aux[6]};
		order.clear();

		order.push_back(aux[1]);
		order.push_back(aux[0]);
		order.push_back(aux[3]);
		order.push_back(aux[2]);

		order.push_back(aux[5]);
		order.push_back(aux[4]);
		order.push_back(aux[7]);
		order.push_back(aux[6]);
	}

	if(sym[1]==1){

		for(int i=0;i<8;i++){
			aux[i]=order[i];
		}

		//order = {aux[3],aux[2],aux[6],aux[7],aux[0],aux[1],aux[5],aux[6]};
		order.clear();

		order.push_back(aux[3]);
		order.push_back(aux[2]);
		order.push_back(aux[1]);
		order.push_back(aux[0]);

		order.push_back(aux[7]);
		order.push_back(aux[6]);
		order.push_back(aux[5]);
		order.push_back(aux[4]);
	}

	if(sym[2]==1){

		for(int i=0;i<8;i++){
			aux[i]=order[i];
		}

		//order = {aux[4],aux[5],aux[6],aux[7],aux[0],aux[1],aux[2],aux[3]};
		order.clear();

		order.push_back(aux[4]);
		order.push_back(aux[5]);
		order.push_back(aux[6]);
		order.push_back(aux[7]);

		order.push_back(aux[0]);
		order.push_back(aux[1]);
		order.push_back(aux[2]);
		order.push_back(aux[3]);
	}

	return order;
}

unsigned node_shared_hash_fn(const void *v, const void *u) {
	const shared_node_t *q = (const shared_node_t*) v;
	uint32_t a, b, c;

	a = (uint32_t) q->id;
	b = (uint32_t) 0;
	c = (uint32_t) 0;
	sc_hash_mix(a, b, c);
	sc_hash_final(a, b, c);
	return (unsigned) c;
}

int node_shared_equal_fn(const void *v, const void *u, const void *w) {
	const shared_node_t *e1 = (const shared_node_t*) v;
	const shared_node_t *e2 = (const shared_node_t*) u;

	return (unsigned) (e1->id == e2->id);
}

void CopyPropEl(hexa_tree_t* mesh, int id, octant_t *elem1){

	octant_t *elem = (octant_t*) sc_array_index(&mesh->elements, id);

	elem1->level = -1;
	elem1->tem = elem->tem;
	elem1->pad = elem->pad;
	elem1->n_mat = elem->n_mat;
	elem1->pml_id = elem->pml_id;
	elem1->father = id;
	elem1->boundary = elem->boundary;

	//TODO try to correct the index of the nodes and elements because of the 27-tree
	for(int i=0; i<8; i++){
		elem1->nodes[i].color = elem->nodes[i].color;
		elem1->nodes[i].x = elem->nodes[i].x;
		elem1->nodes[i].y = elem->nodes[i].y;
		elem1->nodes[i].z = elem->nodes[i].z;
	}
	elem1->x=elem->x;
	elem1->y=elem->y;
	elem1->z=elem->z;
}

void ApplyTemplate1(hexa_tree_t* mesh, std::vector<double>& coords, std::vector<int>& elements_ids, int iel, sc_hash_array_t* hash_nodes, double step){

	sc_array_t toto;
	sc_array_init(&toto, sizeof(octant_t));

	octant_t* elemOrig = (octant_t*) sc_array_index(&mesh->elements, elements_ids[iel]);
	octant_t * elem = (octant_t*) sc_array_push(&toto);
	hexa_element_copy(elemOrig,elem);

	int id = elements_ids[iel];

	double cord_in_ref[3];
	cord_in_ref[0] = 0;
	cord_in_ref[1] = 0;
	cord_in_ref[2] = 0;

	//reference element edge 0
	double local_ref[5][8][3];
	if(true){
		//element 0
		local_ref[0][0][0] = -1;
		local_ref[0][0][1] = -1;
		local_ref[0][0][2] = -1;

		local_ref[0][1][0] = -1+step;
		local_ref[0][1][1] = -1;
		local_ref[0][1][2] = -1;

		local_ref[0][2][0] = -1+step;
		local_ref[0][2][1] = -1+2*step;
		local_ref[0][2][2] = -1;

		local_ref[0][3][0] = -1;
		local_ref[0][3][1] = 1;
		local_ref[0][3][2] = -1;

		local_ref[0][4][0] = -1;
		local_ref[0][4][1] = -1;
		local_ref[0][4][2] = 1;

		local_ref[0][5][0] = -1+step;
		local_ref[0][5][1] = -1;
		local_ref[0][5][2] = -1+2*step;

		local_ref[0][6][0] = -1+step;
		local_ref[0][6][1] = -1+2*step;
		local_ref[0][6][2] = -1+2*step;

		local_ref[0][7][0] = -1;
		local_ref[0][7][1] = 1;
		local_ref[0][7][2] = 1;

		//element 1
		local_ref[1][0][0] = -1+step;
		local_ref[1][0][1] = -1;
		local_ref[1][0][2] = -1;

		local_ref[1][1][0] = -1+2*step;
		local_ref[1][1][1] = -1;
		local_ref[1][1][2] = -1;

		local_ref[1][2][0] = -1+2*step;
		local_ref[1][2][1] = -1+2*step;
		local_ref[1][2][2] = -1;

		local_ref[1][3][0] = -1+step;
		local_ref[1][3][1] = -1+2*step;
		local_ref[1][3][2] = -1;

		local_ref[1][4][0] = -1+step;
		local_ref[1][4][1] = -1;
		local_ref[1][4][2] = -1+2*step;

		local_ref[1][5][0] = -1+2*step;
		local_ref[1][5][1] = -1;
		local_ref[1][5][2] = -1+2*step;

		local_ref[1][6][0] = -1+2*step;
		local_ref[1][6][1] = -1+2*step;
		local_ref[1][6][2] = -1+2*step;

		local_ref[1][7][0] = -1+step;
		local_ref[1][7][1] = -1+2*step;
		local_ref[1][7][2] = -1+2*step;


		//element 2
		local_ref[2][0][0] = -1+2*step;
		local_ref[2][0][1] = -1;
		local_ref[2][0][2] = -1;

		local_ref[2][1][0] =  1;
		local_ref[2][1][1] = -1;
		local_ref[2][1][2] = -1;

		local_ref[2][2][0] =  1;
		local_ref[2][2][1] =  1;
		local_ref[2][2][2] = -1;

		local_ref[2][3][0] = -1+2*step;
		local_ref[2][3][1] = -1+2*step;
		local_ref[2][3][2] = -1;

		local_ref[2][4][0] = -1+2*step;
		local_ref[2][4][1] = -1;
		local_ref[2][4][2] = -1+2*step;

		local_ref[2][5][0] =  1;
		local_ref[2][5][1] = -1;
		local_ref[2][5][2] =  1;

		local_ref[2][6][0] = 1;
		local_ref[2][6][1] = 1;
		local_ref[2][6][2] = 1;

		local_ref[2][7][0] = -1+2*step;
		local_ref[2][7][1] = -1+2*step;
		local_ref[2][7][2] = -1+2*step;


		//element 3
		local_ref[3][0][0] = -1+1*step;
		local_ref[3][0][1] = -1+2*step;
		local_ref[3][0][2] = -1;

		local_ref[3][1][0] = -1+2*step;
		local_ref[3][1][1] = -1+2*step;
		local_ref[3][1][2] = -1;

		local_ref[3][2][0] =  1;
		local_ref[3][2][1] =  1;
		local_ref[3][2][2] = -1;

		local_ref[3][3][0] = -1;
		local_ref[3][3][1] = 1;
		local_ref[3][3][2] = -1;

		local_ref[3][4][0] = -1+1*step;
		local_ref[3][4][1] = -1+2*step;
		local_ref[3][4][2] = -1+2*step;

		local_ref[3][5][0] = -1+2*step;
		local_ref[3][5][1] = -1+2*step;
		local_ref[3][5][2] = -1+2*step;

		local_ref[3][6][0] = 1;
		local_ref[3][6][1] = 1;
		local_ref[3][6][2] = 1;

		local_ref[3][7][0] = -1;
		local_ref[3][7][1] = 1;
		local_ref[3][7][2] = 1;

		//element 4
		local_ref[4][0][0] = -1+step;
		local_ref[4][0][1] = -1;
		local_ref[4][0][2] = -1+2*step;

		local_ref[4][1][0] = -1+2*step;
		local_ref[4][1][1] = -1;
		local_ref[4][1][2] = -1+2*step;

		local_ref[4][2][0] = -1+2*step;
		local_ref[4][2][1] = -1+2*step;
		local_ref[4][2][2] = -1+2*step;

		local_ref[4][3][0] = -1+step;
		local_ref[4][3][1] = -1+2*step;
		local_ref[4][3][2] = -1+2*step;

		local_ref[4][4][0] = -1;
		local_ref[4][4][1] = -1;
		local_ref[4][4][2] = 1;

		local_ref[4][5][0] = 1;
		local_ref[4][5][1] = -1;
		local_ref[4][5][2] = 1;

		local_ref[4][6][0] = 1;
		local_ref[4][6][1] = 1;
		local_ref[4][6][2] = 1;

		local_ref[4][7][0] = -1;
		local_ref[4][7][1] = 1;
		local_ref[4][7][2] = 1;
	}

	//define the rotation of the reference element
	int rot[3];
	int sym[3];
	int id_node[8];
	vector<int> ord;

	if(elem->pad==10){
		//edge 0
		rot[0] = 0;
		rot[1] = 0;
		rot[2] = 0;

		sym[0] = 0;
		sym[1] = 0;
		sym[2] = 0;

		ord = RotateHex(rot,sym);

		for(int node_id=0;node_id<8;node_id++){
			id_node[node_id] = elem->nodes[ord[node_id]].id;
		}

	}else if(elem->pad==11){
		//edge 1
		rot[0] = 0;
		rot[1] = 0;
		rot[2] = -1;

		sym[0] = 0;
		sym[1] = 0;
		sym[2] = 0;

		ord = RotateHex(rot,sym);

		for(int node_id=0;node_id<8;node_id++){
			id_node[node_id] = elem->nodes[ord[node_id]].id;
		}

	}else if(elem->pad==12){
		//edge 2
		rot[0] = 0;
		rot[1] = 0;
		rot[2] = 0;

		sym[0] = 0;
		sym[1] = 1;
		sym[2] = 0;

		ord = RotateHex(rot,sym);

		for(int node_id=0;node_id<8;node_id++){
			id_node[node_id] = elem->nodes[ord[node_id]].id;
		}

	}else if(elem->pad==13){
		//edge 3
		rot[0] = 0;
		rot[1] = 0;
		rot[2] = 1;

		sym[0] = 0;
		sym[1] = 0;
		sym[2] = 0;

		ord = RotateHex(rot,sym);

		for(int node_id=0;node_id<8;node_id++){
			id_node[node_id] = elem->nodes[ord[node_id]].id;
		}

	}else if(elem->pad==14){
		//edge 4
		rot[0] = 0;
		rot[1] = -1;
		rot[2] = 0;

		sym[0] = 0;
		sym[1] = 0;
		sym[2] = 0;

		ord = RotateHex(rot,sym);

		for(int node_id=0;node_id<8;node_id++){
			id_node[node_id] = elem->nodes[ord[node_id]].id;
		}

	}else if(elem->pad==15){
		//edge 5
		rot[0] = 0;
		rot[1] = 1;
		rot[2] = 0;

		sym[0] = 0;
		sym[1] = 0;
		sym[2] = 0;

		ord = RotateHex(rot,sym);

		for(int node_id=0;node_id<8;node_id++){
			id_node[node_id] = elem->nodes[ord[node_id]].id;
		}

	}else if(elem->pad==16){
		//edge 6
		rot[0] = 0;
		rot[1] = 1;
		rot[2] = 0;

		sym[0] = 0;
		sym[1] = 1;
		sym[2] = 0;

		ord = RotateHex(rot,sym);

		for(int node_id=0;node_id<8;node_id++){
			id_node[node_id] = elem->nodes[ord[node_id]].id;
		}

	}else if(elem->pad==17){
		//edge 7
		rot[0] = 0;
		rot[1] = -1;
		rot[2] = 0;

		sym[0] = 0;
		sym[1] = 1;
		sym[2] = 0;

		ord = RotateHex(rot,sym);

		for(int node_id=0;node_id<8;node_id++){
			id_node[node_id] = elem->nodes[ord[node_id]].id;
		}

	}else if(elem->pad==18){
		//edge 8
		rot[0] = 0;
		rot[1] = 0;
		rot[2] = 0;

		sym[0] = 0;
		sym[1] = 0;
		sym[2] = 1;

		ord = RotateHex(rot,sym);

		for(int node_id=0;node_id<8;node_id++){
			id_node[node_id] = elem->nodes[ord[node_id]].id;
		}

	}else if(elem->pad==19){
		//edge 9
		rot[0] = 0;
		rot[1] = 0;
		rot[2] = -1;

		sym[0] = 0;
		sym[1] = 0;
		sym[2] = 1;

		ord = RotateHex(rot,sym);

		for(int node_id=0;node_id<8;node_id++){
			id_node[node_id] = elem->nodes[ord[node_id]].id;
		}

	}else if(elem->pad==20){
		//edge 10
		rot[0] = -1;
		rot[1] = 0;
		rot[2] = 0;

		sym[0] = 0;
		sym[1] = 1;
		sym[2] = 0;

		ord = RotateHex(rot,sym);

		for(int node_id=0;node_id<8;node_id++){
			id_node[node_id] = elem->nodes[ord[node_id]].id;
		}

	}else if(elem->pad==21){
		//edge 11
		rot[0] = 0;
		rot[1] = 0;
		rot[2] = 1;

		sym[0] = 0;
		sym[1] = 0;
		sym[2] = 1;

		ord = RotateHex(rot,sym);

		for(int node_id=0;node_id<8;node_id++){
			id_node[node_id] = elem->nodes[ord[node_id]].id;
		}

	}else{
		printf ("Error in template 1\n");
		exit (EXIT_FAILURE);
	}


	for(int i=0;i<5;i++){
		int conn_p[8];
		GtsPoint* point[8]={NULL};

		double cord_in_x[8],cord_in_y[8],cord_in_z[8];
		//add the nodes in the coord vector
		for (int ii = 0; ii < 8; ii++){
			cord_in_x[ii]=coords[3*id_node[ii]] ;
			cord_in_y[ii]=coords[3*id_node[ii]+1] ;
			cord_in_z[ii]=coords[3*id_node[ii]+2] ;
			//fprintf(mesh->fdbg,"coord in: %f, %f, %f, in the node: %d\n",cord_in_x[ii],cord_in_y[ii],cord_in_z[ii],elem->nodes[ii].id);
		}

		//add the new nodes in the
		for(int ii=0;ii<8;ii++){
			conn_p[ii] = 0;

			if((local_ref[i][ii][0]==1 || local_ref[i][ii][0]==-1) &&
					(local_ref[i][ii][1]==1 || local_ref[i][ii][1]==-1) &&
					(local_ref[i][ii][2]==1 || local_ref[i][ii][2]==-1)){

				if(local_ref[i][ii][0]==-1 && local_ref[i][ii][1]==-1 && local_ref[i][ii][2]==-1){
					conn_p[ii] = id_node[0];
				}else if(local_ref[i][ii][0]==1 && local_ref[i][ii][1]==-1 && local_ref[i][ii][2]==-1){
					conn_p[ii] = id_node[1];
				}else if(local_ref[i][ii][0]==1 && local_ref[i][ii][1]==1 && local_ref[i][ii][2]==-1){
					conn_p[ii] = id_node[2];
				}else if(local_ref[i][ii][0]==-1 && local_ref[i][ii][1]==1 && local_ref[i][ii][2]==-1){
					conn_p[ii] = id_node[3];
				}else if(local_ref[i][ii][0]==-1 && local_ref[i][ii][1]==-1 && local_ref[i][ii][2]==1){
					conn_p[ii] = id_node[4];
				}else if(local_ref[i][ii][0]==1 && local_ref[i][ii][1]==-1 && local_ref[i][ii][2]==1){
					conn_p[ii] = id_node[5];
				}else if(local_ref[i][ii][0]==1 && local_ref[i][ii][1]==1 && local_ref[i][ii][2]==1){
					conn_p[ii] = id_node[6];
				}else if(local_ref[i][ii][0]==-1 && local_ref[i][ii][1]==1 && local_ref[i][ii][2]==1){
					conn_p[ii] = id_node[7];
				}
				//fprintf(mesh->fdbg,"coord out: %f, %f, %f, in the node: %d\n",coords[3*conn_p[ii]],coords[3*conn_p[ii]+1],coords[3*conn_p[ii]+2],conn_p[ii]);
			}else{
				cord_in_ref[0] = local_ref[i][ii][0];
				cord_in_ref[1] = local_ref[i][ii][1];
				cord_in_ref[2] = local_ref[i][ii][2];

				point[ii] = LinearMapHex(cord_in_ref, cord_in_x,cord_in_y,cord_in_z);
				double var[3];
				var[0] = point[ii]->x;
				var[1] = point[ii]->y;
				var[2] = point[ii]->z;
				conn_p[ii] = AddPoint( mesh, hash_nodes, point[ii] , coords);
				//fprintf(mesh->fdbg,"coord out: %f, %f, %f, in the node: %d\n",var[0],var[1],var[2],conn_p[ii]);
			}
		}

		int aux[8];
		for(int ino=0;ino<8;ino++){
			aux[ino] = conn_p[ino];
		}
		for(int ino=0;ino<8;ino++){
			conn_p[ord[ino]] = aux[ino];
		}
		//add the new elements and make the connectivity
		if(i==0){
			octant_t *elem1 = (octant_t*) sc_array_index(&mesh->elements, elements_ids[iel]);

			elem1->nodes[0].id = conn_p[0];
			elem1->nodes[1].id = conn_p[1];
			elem1->nodes[2].id = conn_p[2];
			elem1->nodes[3].id = conn_p[3];

			elem1->nodes[4].id = conn_p[4];
			elem1->nodes[5].id = conn_p[5];
			elem1->nodes[6].id = conn_p[6];
			elem1->nodes[7].id = conn_p[7];

			CopyPropEl(mesh,id,elem1);

		}else{
			octant_t* elem2 = (octant_t*) sc_array_push(&mesh->elements);

			elem2->nodes[0].id = conn_p[0];
			elem2->nodes[1].id = conn_p[1];
			elem2->nodes[2].id = conn_p[2];
			elem2->nodes[3].id = conn_p[3];

			elem2->nodes[4].id = conn_p[4];
			elem2->nodes[5].id = conn_p[5];
			elem2->nodes[6].id = conn_p[6];
			elem2->nodes[7].id = conn_p[7];

			CopyPropEl(mesh,id,elem2);

		}
	}

	sc_array_reset(&toto);

}

void ApplyTemplate2(hexa_tree_t* mesh, std::vector<double>& coords, std::vector<int>& elements_ids, int iel, sc_hash_array_t* hash_nodes, double step){

	sc_array_t toto;
	sc_array_init(&toto, sizeof(octant_t));

	octant_t* elemOrig = (octant_t*) sc_array_index(&mesh->elements, elements_ids[iel]);
	octant_t * elem = (octant_t*) sc_array_push(&toto);
	hexa_element_copy(elemOrig,elem);

	int id = elements_ids[iel];

	double cord_in_ref[3];
	cord_in_ref[0] = 0;
	cord_in_ref[1] = 0;
	cord_in_ref[2] = 0;

	//reference element edge 0
	double local_ref[3][8][3];
	if(true){
		//element 0
		local_ref[0][0][0] = -1;
		local_ref[0][0][1] = -1;
		local_ref[0][0][2] = -1;

		local_ref[0][1][0] = 1;
		local_ref[0][1][1] = -1;
		local_ref[0][1][2] = -1;

		local_ref[0][2][0] = 1;
		local_ref[0][2][1] = 1;
		local_ref[0][2][2] = -1;

		local_ref[0][3][0] = -1;
		local_ref[0][3][1] = 1;
		local_ref[0][3][2] = -1;

		local_ref[0][4][0] = -1;
		local_ref[0][4][1] = -1;
		local_ref[0][4][2] = -1+step;

		local_ref[0][5][0] = 1;
		local_ref[0][5][1] = -1;
		local_ref[0][5][2] = -1+step;

		local_ref[0][6][0] = 1;
		local_ref[0][6][1] = 1;
		local_ref[0][6][2] = -1+step;

		local_ref[0][7][0] = -1;
		local_ref[0][7][1] = 1;
		local_ref[0][7][2] = -1+step;

		//element 1
		local_ref[1][0][0] = -1;
		local_ref[1][0][1] = -1;
		local_ref[1][0][2] = -1+step;

		local_ref[1][1][0] = 1;
		local_ref[1][1][1] = -1;
		local_ref[1][1][2] = -1+step;

		local_ref[1][2][0] = 1;
		local_ref[1][2][1] = 1;
		local_ref[1][2][2] = -1+step;

		local_ref[1][3][0] = -1;
		local_ref[1][3][1] = 1;
		local_ref[1][3][2] = -1+step;

		local_ref[1][4][0] = -1;
		local_ref[1][4][1] = -1;
		local_ref[1][4][2] = -1+2*step;

		local_ref[1][5][0] = 1;
		local_ref[1][5][1] = -1;
		local_ref[1][5][2] = -1+2*step;

		local_ref[1][6][0] = 1;
		local_ref[1][6][1] = 1;
		local_ref[1][6][2] = -1+2*step;

		local_ref[1][7][0] = -1;
		local_ref[1][7][1] = 1;
		local_ref[1][7][2] = -1+2*step;


		//element 2
		local_ref[2][0][0] = -1;
		local_ref[2][0][1] = -1;
		local_ref[2][0][2] = -1+2*step;

		local_ref[2][1][0] =  1;
		local_ref[2][1][1] = -1;
		local_ref[2][1][2] = -1+2*step;

		local_ref[2][2][0] =  1;
		local_ref[2][2][1] =  1;
		local_ref[2][2][2] = -1+2*step;

		local_ref[2][3][0] = -1;
		local_ref[2][3][1] = 1;
		local_ref[2][3][2] = -1+2*step;

		local_ref[2][4][0] = -1;
		local_ref[2][4][1] = -1;
		local_ref[2][4][2] = 1;

		local_ref[2][5][0] =  1;
		local_ref[2][5][1] = -1;
		local_ref[2][5][2] =  1;

		local_ref[2][6][0] = 1;
		local_ref[2][6][1] = 1;
		local_ref[2][6][2] = 1;

		local_ref[2][7][0] = -1;
		local_ref[2][7][1] = 1;
		local_ref[2][7][2] = 1;
	}

	//define the rotation of the reference element
	int rot[3];
	int sym[3];
	int id_node[8];
	vector<int> ord;

	if(elem->pad==24){
		//edge 4 5 6 7
		rot[0] = 0;
		rot[1] = 0;
		rot[2] = 0;

		sym[0] = 0;
		sym[1] = 0;
		sym[2] = 0;

		ord = RotateHex(rot,sym);

		for(int node_id=0;node_id<8;node_id++){
			id_node[node_id] = elem->nodes[ord[node_id]].id;
		}

	}else if(elem->pad==23){
		//edge 1 3 9 11
		rot[0] = 1;
		rot[1] = 0;
		rot[2] = 0;

		sym[0] = 0;
		sym[1] = 0;
		sym[2] = 0;

		ord = RotateHex(rot,sym);

		for(int node_id=0;node_id<8;node_id++){
			id_node[node_id] = elem->nodes[ord[node_id]].id;
		}

	}else if(elem->pad==22){
		//edge 0 2 8 10
		rot[0] = 0;
		rot[1] = 1;
		rot[2] = 0;

		sym[0] = 0;
		sym[1] = 0;
		sym[2] = 0;

		ord = RotateHex(rot,sym);

		for(int node_id=0;node_id<8;node_id++){
			id_node[node_id] = elem->nodes[ord[node_id]].id;
		}

	}else{
		printf ("Error in template 2\n");
		exit (EXIT_FAILURE);
	}


	for(int i=0;i<3;i++){
		int conn_p[8];
		GtsPoint* point[8]={NULL};

		double cord_in_x[8],cord_in_y[8],cord_in_z[8];
		//add the nodes in the coord vector
		for (int ii = 0; ii < 8; ii++){
			cord_in_x[ii]=coords[3*id_node[ii]] ;
			cord_in_y[ii]=coords[3*id_node[ii]+1] ;
			cord_in_z[ii]=coords[3*id_node[ii]+2] ;
			//fprintf(mesh->fdbg,"coord in: %f, %f, %f, in the node: %d\n",cord_in_x[ii],cord_in_y[ii],cord_in_z[ii],elem->nodes[ii].id);
		}

		//add the new nodes in the
		for(int ii=0;ii<8;ii++){
			conn_p[ii] = 0;

			if((local_ref[i][ii][0]==1 || local_ref[i][ii][0]==-1) &&
					(local_ref[i][ii][1]==1 || local_ref[i][ii][1]==-1) &&
					(local_ref[i][ii][2]==1 || local_ref[i][ii][2]==-1)){

				if(local_ref[i][ii][0]==-1 && local_ref[i][ii][1]==-1 && local_ref[i][ii][2]==-1){
					conn_p[ii] = id_node[0];
				}else if(local_ref[i][ii][0]==1 && local_ref[i][ii][1]==-1 && local_ref[i][ii][2]==-1){
					conn_p[ii] = id_node[1];
				}else if(local_ref[i][ii][0]==1 && local_ref[i][ii][1]==1 && local_ref[i][ii][2]==-1){
					conn_p[ii] = id_node[2];
				}else if(local_ref[i][ii][0]==-1 && local_ref[i][ii][1]==1 && local_ref[i][ii][2]==-1){
					conn_p[ii] = id_node[3];
				}else if(local_ref[i][ii][0]==-1 && local_ref[i][ii][1]==-1 && local_ref[i][ii][2]==1){
					conn_p[ii] = id_node[4];
				}else if(local_ref[i][ii][0]==1 && local_ref[i][ii][1]==-1 && local_ref[i][ii][2]==1){
					conn_p[ii] = id_node[5];
				}else if(local_ref[i][ii][0]==1 && local_ref[i][ii][1]==1 && local_ref[i][ii][2]==1){
					conn_p[ii] = id_node[6];
				}else if(local_ref[i][ii][0]==-1 && local_ref[i][ii][1]==1 && local_ref[i][ii][2]==1){
					conn_p[ii] = id_node[7];
				}
				//fprintf(mesh->fdbg,"coord out: %f, %f, %f, in the node: %d\n",coords[3*conn_p[ii]],coords[3*conn_p[ii]+1],coords[3*conn_p[ii]+2],conn_p[ii]);
			}else{
				cord_in_ref[0] = local_ref[i][ii][0];
				cord_in_ref[1] = local_ref[i][ii][1];
				cord_in_ref[2] = local_ref[i][ii][2];

				point[ii] = LinearMapHex(cord_in_ref, cord_in_x,cord_in_y,cord_in_z);
				double var[3];
				var[0] = point[ii]->x;
				var[1] = point[ii]->y;
				var[2] = point[ii]->z;
				conn_p[ii] = AddPoint( mesh, hash_nodes, point[ii] , coords);
				//fprintf(mesh->fdbg,"coord out: %f, %f, %f, in the node: %d\n",var[0],var[1],var[2],conn_p[ii]);

			}
		}

		int aux[8];
		for(int ino=0;ino<8;ino++){
			aux[ino] = conn_p[ino];
		}
		for(int ino=0;ino<8;ino++){
			conn_p[ord[ino]] = aux[ino];
		}
		//add the new elements and make the connectivity
		if(i==0){
			octant_t *elem1 = (octant_t*) sc_array_index(&mesh->elements, elements_ids[iel]);

			elem1->nodes[0].id = conn_p[0];
			elem1->nodes[1].id = conn_p[1];
			elem1->nodes[2].id = conn_p[2];
			elem1->nodes[3].id = conn_p[3];

			elem1->nodes[4].id = conn_p[4];
			elem1->nodes[5].id = conn_p[5];
			elem1->nodes[6].id = conn_p[6];
			elem1->nodes[7].id = conn_p[7];

			CopyPropEl(mesh,id,elem1);

		}else{
			octant_t* elem2 = (octant_t*) sc_array_push(&mesh->elements);

			elem2->nodes[0].id = conn_p[0];
			elem2->nodes[1].id = conn_p[1];
			elem2->nodes[2].id = conn_p[2];
			elem2->nodes[3].id = conn_p[3];

			elem2->nodes[4].id = conn_p[4];
			elem2->nodes[5].id = conn_p[5];
			elem2->nodes[6].id = conn_p[6];
			elem2->nodes[7].id = conn_p[7];

			CopyPropEl(mesh,id,elem2);

		}
	}
	sc_array_reset(&toto);
}

void ApplyTemplate3(hexa_tree_t* mesh, std::vector<double>& coords, std::vector<int>& elements_ids, int iel, sc_hash_array_t* hash_nodes, double step){

	sc_array_t toto;
	sc_array_init(&toto, sizeof(octant_t));

	octant_t* elemOrig = (octant_t*) sc_array_index(&mesh->elements, elements_ids[iel]);
	octant_t * elem = (octant_t*) sc_array_push(&toto);
	hexa_element_copy(elemOrig,elem);

	int id = elements_ids[iel];

	double cord_in_ref[3];
	cord_in_ref[0] = 0;
	cord_in_ref[1] = 0;
	cord_in_ref[2] = 0;

	//reference element edge 0
	double local_ref[4][8][3];
	if(true){
		//element 0
		local_ref[0][0][0] = -1;
		local_ref[0][0][1] = -1;
		local_ref[0][0][2] = -1;

		local_ref[0][1][0] = 1;
		local_ref[0][1][1] = -1;
		local_ref[0][1][2] = -1;

		local_ref[0][2][0] = 1;
		local_ref[0][2][1] = -1+1*step;
		local_ref[0][2][2] = -1;

		local_ref[0][3][0] = -1;
		local_ref[0][3][1] = -1+1*step;
		local_ref[0][3][2] = -1;

		local_ref[0][4][0] = -1;
		local_ref[0][4][1] = -1;
		local_ref[0][4][2] = 1;

		local_ref[0][5][0] = 1;
		local_ref[0][5][1] = -1;
		local_ref[0][5][2] = 1;

		local_ref[0][6][0] = 1;
		local_ref[0][6][1] = -1+1*step;
		local_ref[0][6][2] = -1+2*step;

		local_ref[0][7][0] = -1;
		local_ref[0][7][1] = -1+1*step;
		local_ref[0][7][2] = -1+2*step;

		//element 1
		local_ref[1][0][0] = -1;
		local_ref[1][0][1] = -1+step;
		local_ref[1][0][2] = -1;

		local_ref[1][1][0] = 1;
		local_ref[1][1][1] = -1+step;
		local_ref[1][1][2] = -1;

		local_ref[1][2][0] = 1;
		local_ref[1][2][1] = -1+2*step;
		local_ref[1][2][2] = -1;

		local_ref[1][3][0] = -1;
		local_ref[1][3][1] = -1+2*step;
		local_ref[1][3][2] = -1;

		local_ref[1][4][0] = -1;
		local_ref[1][4][1] = -1+step;
		local_ref[1][4][2] = -1+2*step;

		local_ref[1][5][0] = 1;
		local_ref[1][5][1] = -1+step;
		local_ref[1][5][2] = -1+2*step;

		local_ref[1][6][0] = 1;
		local_ref[1][6][1] = -1+2*step;
		local_ref[1][6][2] = -1+2*step;

		local_ref[1][7][0] = -1;
		local_ref[1][7][1] = -1+2*step;
		local_ref[1][7][2] = -1+2*step;


		//element 2
		local_ref[2][0][0] = -1;
		local_ref[2][0][1] = -1+2*step;
		local_ref[2][0][2] = -1;

		local_ref[2][1][0] =  1;
		local_ref[2][1][1] = -1+2*step;
		local_ref[2][1][2] = -1;

		local_ref[2][2][0] =  1;
		local_ref[2][2][1] =  1;
		local_ref[2][2][2] = -1;

		local_ref[2][3][0] = -1;
		local_ref[2][3][1] =  1;
		local_ref[2][3][2] = -1;

		local_ref[2][4][0] = -1;
		local_ref[2][4][1] = -1+2*step;
		local_ref[2][4][2] = -1+2*step;

		local_ref[2][5][0] =  1;
		local_ref[2][5][1] = -1+2*step;
		local_ref[2][5][2] =  -1+2*step;

		local_ref[2][6][0] = 1;
		local_ref[2][6][1] = 1;
		local_ref[2][6][2] = 1;

		local_ref[2][7][0] = -1;
		local_ref[2][7][1] = 1;
		local_ref[2][7][2] = 1;

		//element 3
		local_ref[3][0][0] = -1;
		local_ref[3][0][1] = -1+step;
		local_ref[3][0][2] = -1+2*step;

		local_ref[3][1][0] =  1;
		local_ref[3][1][1] = -1+step;
		local_ref[3][1][2] = -1+2*step;

		local_ref[3][2][0] =  1;
		local_ref[3][2][1] = -1+2*step;
		local_ref[3][2][2] = -1+2*step;

		local_ref[3][3][0] = -1;
		local_ref[3][3][1] = -1+2*step;
		local_ref[3][3][2] = -1+2*step;

		local_ref[3][4][0] = -1;
		local_ref[3][4][1] = -1;
		local_ref[3][4][2] = 1;

		local_ref[3][5][0] =  1;
		local_ref[3][5][1] = -1;
		local_ref[3][5][2] =  1;

		local_ref[3][6][0] = 1;
		local_ref[3][6][1] = 1;
		local_ref[3][6][2] = 1;

		local_ref[3][7][0] = -1;
		local_ref[3][7][1] = 1;
		local_ref[3][7][2] = 1;
	}

	//define the rotation of the reference element
	int rot[3];
	int sym[3];
	int id_node[8];
	vector<int> ord;

	if(elem->pad==25){
		//edge 0 2
		rot[0] = 0;
		rot[1] = 0;
		rot[2] = 1;

		sym[0] = 0;
		sym[1] = 0;
		sym[2] = 0;

		ord = RotateHex(rot,sym);

		for(int node_id=0;node_id<8;node_id++){
			id_node[node_id] = elem->nodes[ord[node_id]].id;
		}

	}else if(elem->pad==26){
		//edge 1 3
		rot[0] = 0;
		rot[1] = 0;
		rot[2] = 0;

		sym[0] = 0;
		sym[1] = 0;
		sym[2] = 0;

		ord = RotateHex(rot,sym);

		for(int node_id=0;node_id<8;node_id++){
			id_node[node_id] = elem->nodes[ord[node_id]].id;
		}

	}else if(elem->pad==27){
		//edge 8 10
		rot[0] = 0;
		rot[1] = 0;
		rot[2] = 0;

		sym[0] = 0;
		sym[1] = 0;
		sym[2] = 1;

		ord = RotateHex(rot,sym);

		for(int node_id=0;node_id<8;node_id++){
			id_node[node_id] = elem->nodes[ord[node_id]].id;
		}

	}else if(elem->pad==28){
		//edge 9 11
		rot[0] = 0;
		rot[1] = 0;
		rot[2] = 1;

		sym[0] = 0;
		sym[1] = 0;
		sym[2] = 1;

		ord = RotateHex(rot,sym);

		for(int node_id=0;node_id<8;node_id++){
			id_node[node_id] = elem->nodes[ord[node_id]].id;
		}

	}else if(elem->pad==29){
		//edge 4 5
		rot[0] = 1;
		rot[1] = 0;
		rot[2] = 0;

		sym[0] = 0;
		sym[1] = 0;
		sym[2] = 0;

		ord = RotateHex(rot,sym);

		for(int node_id=0;node_id<8;node_id++){
			id_node[node_id] = elem->nodes[ord[node_id]].id;
		}

	}else if(elem->pad==30){
		//edge 5 6
		rot[0] = -1;
		rot[1] = 1;
		rot[2] = 0;

		sym[0] = 0;
		sym[1] = 0;
		sym[2] = 0;

		ord = RotateHex(rot,sym);

		for(int node_id=0;node_id<8;node_id++){
			id_node[node_id] = elem->nodes[ord[node_id]].id;
		}

	}else if(elem->pad==31){
		//edge 6 7
		rot[0] = -1;
		rot[1] = 0;
		rot[2] = 0;

		sym[0] = 0;
		sym[1] = 0;
		sym[2] = 0;

		ord = RotateHex(rot,sym);

		for(int node_id=0;node_id<8;node_id++){
			id_node[node_id] = elem->nodes[ord[node_id]].id;
		}

	}else if(elem->pad==32){
		//edge 7 4
		rot[0] = -1;
		rot[1] = -1;
		rot[2] = 0;

		sym[0] = 0;
		sym[1] = 0;
		sym[2] = 0;

		ord = RotateHex(rot,sym);

		for(int node_id=0;node_id<8;node_id++){
			id_node[node_id] = elem->nodes[ord[node_id]].id;
		}

	}else if(elem->pad==33){
		//edge 0 8
		rot[0] = 1;
		rot[1] = 0;
		rot[2] = 1;

		sym[0] = 0;
		sym[1] = 0;
		sym[2] = 0;

		ord = RotateHex(rot,sym);

		for(int node_id=0;node_id<8;node_id++){
			id_node[node_id] = elem->nodes[ord[node_id]].id;
		}

	}else if(elem->pad==34){
		//edge 1 9
		rot[0] = 0;
		rot[1] = 1;
		rot[2] = 0;

		sym[0] = 0;
		sym[1] = 0;
		sym[2] = 0;

		ord = RotateHex(rot,sym);

		for(int node_id=0;node_id<8;node_id++){
			id_node[node_id] = elem->nodes[ord[node_id]].id;
		}

	}else if(elem->pad==35){
		//edge 2 10
		rot[0] = -1;
		rot[1] = 0;
		rot[2] = 1;

		sym[0] = 0;
		sym[1] = 0;
		sym[2] = 0;

		ord = RotateHex(rot,sym);

		for(int node_id=0;node_id<8;node_id++){
			id_node[node_id] = elem->nodes[ord[node_id]].id;
		}

	}else if(elem->pad==36){
		//edge 3 11
		rot[0] = 0;
		rot[1] = -1;
		rot[2] = 0;

		sym[0] = 0;
		sym[1] = 0;
		sym[2] = 0;

		ord = RotateHex(rot,sym);

		for(int node_id=0;node_id<8;node_id++){
			id_node[node_id] = elem->nodes[ord[node_id]].id;
		}

	}else{
		printf ("Error in template 3\n");
		exit (EXIT_FAILURE);
	}


	for(int i=0;i<4;i++){
		int conn_p[8];
		GtsPoint* point[8]={NULL};

		double cord_in_x[8],cord_in_y[8],cord_in_z[8];
		//add the nodes in the coord vector
		for (int ii = 0; ii < 8; ii++){
			cord_in_x[ii]=coords[3*id_node[ii]] ;
			cord_in_y[ii]=coords[3*id_node[ii]+1] ;
			cord_in_z[ii]=coords[3*id_node[ii]+2] ;
			//fprintf(mesh->fdbg,"coord in: %f, %f, %f, in the node: %d\n",cord_in_x[ii],cord_in_y[ii],cord_in_z[ii],elem->nodes[ii].id);
		}

		//add the new nodes in the
		for(int ii=0;ii<8;ii++){
			conn_p[ii] = 0;

			if((local_ref[i][ii][0]==1 || local_ref[i][ii][0]==-1) &&
					(local_ref[i][ii][1]==1 || local_ref[i][ii][1]==-1) &&
					(local_ref[i][ii][2]==1 || local_ref[i][ii][2]==-1)){

				if(local_ref[i][ii][0]==-1 && local_ref[i][ii][1]==-1 && local_ref[i][ii][2]==-1){
					conn_p[ii] = id_node[0];
				}else if(local_ref[i][ii][0]==1 && local_ref[i][ii][1]==-1 && local_ref[i][ii][2]==-1){
					conn_p[ii] = id_node[1];
				}else if(local_ref[i][ii][0]==1 && local_ref[i][ii][1]==1 && local_ref[i][ii][2]==-1){
					conn_p[ii] = id_node[2];
				}else if(local_ref[i][ii][0]==-1 && local_ref[i][ii][1]==1 && local_ref[i][ii][2]==-1){
					conn_p[ii] = id_node[3];
				}else if(local_ref[i][ii][0]==-1 && local_ref[i][ii][1]==-1 && local_ref[i][ii][2]==1){
					conn_p[ii] = id_node[4];
				}else if(local_ref[i][ii][0]==1 && local_ref[i][ii][1]==-1 && local_ref[i][ii][2]==1){
					conn_p[ii] = id_node[5];
				}else if(local_ref[i][ii][0]==1 && local_ref[i][ii][1]==1 && local_ref[i][ii][2]==1){
					conn_p[ii] = id_node[6];
				}else if(local_ref[i][ii][0]==-1 && local_ref[i][ii][1]==1 && local_ref[i][ii][2]==1){
					conn_p[ii] = id_node[7];
				}
				//fprintf(mesh->fdbg,"coord out: %f, %f, %f, in the node: %d\n",coords[3*conn_p[ii]],coords[3*conn_p[ii]+1],coords[3*conn_p[ii]+2],conn_p[ii]);
			}else{
				cord_in_ref[0] = local_ref[i][ii][0];
				cord_in_ref[1] = local_ref[i][ii][1];
				cord_in_ref[2] = local_ref[i][ii][2];

				point[ii] = LinearMapHex(cord_in_ref, cord_in_x,cord_in_y,cord_in_z);
				double var[3];
				var[0] = point[ii]->x;
				var[1] = point[ii]->y;
				var[2] = point[ii]->z;
				conn_p[ii] = AddPoint( mesh, hash_nodes, point[ii] , coords);
				//fprintf(mesh->fdbg,"coord out: %f, %f, %f, in the node: %d\n",var[0],var[1],var[2],conn_p[ii]);

			}
		}

		int aux[8];
		for(int ino=0;ino<8;ino++){
			aux[ino] = conn_p[ino];
		}
		for(int ino=0;ino<8;ino++){
			conn_p[ord[ino]] = aux[ino];
		}

		//add the new elements and make the connectivity
		if(i==0){
			octant_t *elem1 = (octant_t*) sc_array_index(&mesh->elements, elements_ids[iel]);

			elem1->nodes[0].id = conn_p[0];
			elem1->nodes[1].id = conn_p[1];
			elem1->nodes[2].id = conn_p[2];
			elem1->nodes[3].id = conn_p[3];

			elem1->nodes[4].id = conn_p[4];
			elem1->nodes[5].id = conn_p[5];
			elem1->nodes[6].id = conn_p[6];
			elem1->nodes[7].id = conn_p[7];

			CopyPropEl(mesh,id,elem1);

		}else{
			octant_t* elem2 = (octant_t*) sc_array_push(&mesh->elements);

			elem2->nodes[0].id = conn_p[0];
			elem2->nodes[1].id = conn_p[1];
			elem2->nodes[2].id = conn_p[2];
			elem2->nodes[3].id = conn_p[3];

			elem2->nodes[4].id = conn_p[4];
			elem2->nodes[5].id = conn_p[5];
			elem2->nodes[6].id = conn_p[6];
			elem2->nodes[7].id = conn_p[7];

			CopyPropEl(mesh,id,elem2);

		}
	}
	sc_array_reset(&toto);
}

void ApplyTemplate4(hexa_tree_t* mesh, std::vector<double>& coords, std::vector<int>& elements_ids, int iel, sc_hash_array_t* hash_nodes, double step){

	sc_array_t toto;
	sc_array_init(&toto, sizeof(octant_t));

	octant_t* elemOrig = (octant_t*) sc_array_index(&mesh->elements, elements_ids[iel]);
	octant_t * elem = (octant_t*) sc_array_push(&toto);
	hexa_element_copy(elemOrig,elem);

	int id = elements_ids[iel];

	double cord_in_ref[3];
	cord_in_ref[0] = 0;
	cord_in_ref[1] = 0;
	cord_in_ref[2] = 0;

	//reference element edge 0
	double local_ref[5][8][3];
	if(true){
		//element 0
		local_ref[0][0][0] = -1;
		local_ref[0][0][1] = -1;
		local_ref[0][0][2] = -1;

		local_ref[0][1][0] = -1+step;
		local_ref[0][1][1] = -1;
		local_ref[0][1][2] = -1;

		local_ref[0][2][0] = -1+step;
		local_ref[0][2][1] = -1+step;
		local_ref[0][2][2] = -1;

		local_ref[0][3][0] = -1;
		local_ref[0][3][1] = -1+step;
		local_ref[0][3][2] = -1;

		local_ref[0][4][0] = -1;
		local_ref[0][4][1] = -1;
		local_ref[0][4][2] =  1;

		local_ref[0][5][0] = -1+step;
		local_ref[0][5][1] = -1;
		local_ref[0][5][2] = 1;

		local_ref[0][6][0] = -1+step;
		local_ref[0][6][1] = -1+step;
		local_ref[0][6][2] = 1;

		local_ref[0][7][0] = -1;
		local_ref[0][7][1] = -1+step;
		local_ref[0][7][2] =  1;

		//element 1
		local_ref[1][0][0] = -1+step;
		local_ref[1][0][1] = -1;
		local_ref[1][0][2] = -1;

		local_ref[1][1][0] = -1+2*step;
		local_ref[1][1][1] = -1;
		local_ref[1][1][2] = -1;

		local_ref[1][2][0] = -1+2*step;
		local_ref[1][2][1] = -1+2*step;
		local_ref[1][2][2] = -1;

		local_ref[1][3][0] = -1+step;
		local_ref[1][3][1] = -1+1*step;
		local_ref[1][3][2] = -1;

		local_ref[1][4][0] = -1+step;
		local_ref[1][4][1] = -1;
		local_ref[1][4][2] = 1;

		local_ref[1][5][0] = -1+2*step;
		local_ref[1][5][1] = -1;
		local_ref[1][5][2] = 1;

		local_ref[1][6][0] = -1+2*step;
		local_ref[1][6][1] = -1+2*step;
		local_ref[1][6][2] = 1;

		local_ref[1][7][0] = -1+step;
		local_ref[1][7][1] = -1+1*step;
		local_ref[1][7][2] = 1;


		//element 2
		local_ref[2][0][0] = -1;
		local_ref[2][0][1] = -1+step;
		local_ref[2][0][2] = -1;

		local_ref[2][1][0] = -1+step;
		local_ref[2][1][1] = -1+step;
		local_ref[2][1][2] = -1;

		local_ref[2][2][0] = -1+2*step;
		local_ref[2][2][1] = -1+2*step;
		local_ref[2][2][2] = -1;

		local_ref[2][3][0] = -1;
		local_ref[2][3][1] = -1+2*step;
		local_ref[2][3][2] = -1;

		local_ref[2][4][0] = -1;
		local_ref[2][4][1] = -1+step;
		local_ref[2][4][2] = 1;

		local_ref[2][5][0] =  -1+step;
		local_ref[2][5][1] = -1+step;
		local_ref[2][5][2] =  1;

		local_ref[2][6][0] = -1+2*step;
		local_ref[2][6][1] = -1+2*step;
		local_ref[2][6][2] = 1;

		local_ref[2][7][0] = -1;
		local_ref[2][7][1] = -1+2*step;
		local_ref[2][7][2] = 1;


		//element 3
		local_ref[3][0][0] = -1+2*step;
		local_ref[3][0][1] = -1;
		local_ref[3][0][2] = -1;

		local_ref[3][1][0] =  1;
		local_ref[3][1][1] = -1;
		local_ref[3][1][2] = -1;

		local_ref[3][2][0] =  1;
		local_ref[3][2][1] =  1;
		local_ref[3][2][2] = -1;

		local_ref[3][3][0] = -1+2*step;
		local_ref[3][3][1] = -1+2*step;
		local_ref[3][3][2] = -1;

		local_ref[3][4][0] = -1+2*step;
		local_ref[3][4][1] = -1;
		local_ref[3][4][2] =  1;

		local_ref[3][5][0] = 1;
		local_ref[3][5][1] = -1;
		local_ref[3][5][2] = 1;

		local_ref[3][6][0] = 1;
		local_ref[3][6][1] = 1;
		local_ref[3][6][2] = 1;

		local_ref[3][7][0] = -1+2*step;
		local_ref[3][7][1] = -1+2*step;
		local_ref[3][7][2] = 1;

		//element 4
		local_ref[4][0][0] = -1;
		local_ref[4][0][1] = -1+2*step;
		local_ref[4][0][2] = -1;

		local_ref[4][1][0] = -1+2*step;
		local_ref[4][1][1] = -1+2*step;
		local_ref[4][1][2] = -1;

		local_ref[4][2][0] = 1;
		local_ref[4][2][1] = 1;
		local_ref[4][2][2] = -1;

		local_ref[4][3][0] = -1;
		local_ref[4][3][1] = 1;
		local_ref[4][3][2] = -1;

		local_ref[4][4][0] = -1;
		local_ref[4][4][1] = -1+2*step;
		local_ref[4][4][2] = 1;

		local_ref[4][5][0] = -1+2*step;
		local_ref[4][5][1] = -1+2*step;
		local_ref[4][5][2] = 1;

		local_ref[4][6][0] = 1;
		local_ref[4][6][1] = 1;
		local_ref[4][6][2] = 1;

		local_ref[4][7][0] = -1;
		local_ref[4][7][1] = 1;
		local_ref[4][7][2] = 1;
	}

	//define the rotation of the reference element
	int rot[3];
	int sym[3];
	int id_node[8];
	vector<int> ord;

	if(elem->pad== 37){
		//edge 0 1 8 9
		rot[0] = 0;
		rot[1] = 0;
		rot[2] = 0;

		sym[0] = 1;
		sym[1] = 0;
		sym[2] = 0;

		ord = RotateHex(rot,sym);

		for(int node_id=0;node_id<8;node_id++){
			id_node[node_id] = elem->nodes[ord[node_id]].id;
		}

	}else if(elem->pad==38){
		//edge 1 2 9 10
		rot[0] = 0;
		rot[1] = 0;
		rot[2] = -1;

		sym[0] = 1;
		sym[1] = 0;
		sym[2] = 0;

		ord = RotateHex(rot,sym);

		for(int node_id=0;node_id<8;node_id++){
			id_node[node_id] = elem->nodes[ord[node_id]].id;
		}

	}else if(elem->pad==39){
		//edge 2 3 10 11
		rot[0] = 0;
		rot[1] = 0;
		rot[2] = 1;

		sym[0] = 0;
		sym[1] = 0;
		sym[2] = 0;

		ord = RotateHex(rot,sym);

		for(int node_id=0;node_id<8;node_id++){
			id_node[node_id] = elem->nodes[ord[node_id]].id;
		}

	}else if(elem->pad==40){
		//edge 0 3 8 11
		rot[0] = 0;
		rot[1] = 0;
		rot[2] = 0;

		sym[0] = 0;
		sym[1] = 0;
		sym[2] = 0;

		ord = RotateHex(rot,sym);

		for(int node_id=0;node_id<8;node_id++){
			id_node[node_id] = elem->nodes[ord[node_id]].id;
		}

	}else if(elem->pad==41){
		//edge 0 4 2 7
		rot[0] = 1;
		rot[1] = 0;
		rot[2] = 0;

		sym[0] = 0;
		sym[1] = 1;
		sym[2] = 0;

		ord = RotateHex(rot,sym);

		for(int node_id=0;node_id<8;node_id++){
			id_node[node_id] = elem->nodes[ord[node_id]].id;
		}

	}else if(elem->pad==42){
		//edge 0 5 2 6
		rot[0] = 1;
		rot[1] = 0;
		rot[2] = 0;

		sym[0] = 1;
		sym[1] = 1;
		sym[2] = 0;

		ord = RotateHex(rot,sym);

		for(int node_id=0;node_id<8;node_id++){
			id_node[node_id] = elem->nodes[ord[node_id]].id;
		}

	}else if(elem->pad==43){
		//edge 5 8 6 10
		rot[0] = 1;
		rot[1] = 0;
		rot[2] = 0;

		sym[0] = 1;
		sym[1] = 0;
		sym[2] = 0;

		ord = RotateHex(rot,sym);

		for(int node_id=0;node_id<8;node_id++){
			id_node[node_id] = elem->nodes[ord[node_id]].id;
		}

	}else if(elem->pad==44){
		//edge 4 8 7 10
		rot[0] = 1;
		rot[1] = 0;
		rot[2] = 0;

		sym[0] = 0;
		sym[1] = 0;
		sym[2] = 0;

		ord = RotateHex(rot,sym);

		for(int node_id=0;node_id<8;node_id++){
			id_node[node_id] = elem->nodes[ord[node_id]].id;
		}

	}else if(elem->pad==45){
		//edge 3 4  1  5
		rot[0] = 0;
		rot[1] = 1;
		rot[2] = 0;

		sym[0] = 0;
		sym[1] = 0;
		sym[2] = 0;

		ord = RotateHex(rot,sym);

		for(int node_id=0;node_id<8;node_id++){
			id_node[node_id] = elem->nodes[ord[node_id]].id;
		}

	}else if(elem->pad==46){
		//edge 1 6 3 7
		rot[0] = 0;
		rot[1] = 1;
		rot[2] = 0;

		sym[0] = 0;
		sym[1] = 1;
		sym[2] = 0;

		ord = RotateHex(rot,sym);

		for(int node_id=0;node_id<8;node_id++){
			id_node[node_id] = elem->nodes[ord[node_id]].id;
		}

	}else if(elem->pad==47){
		//edge 9 6  7 11
		rot[0] = 0;
		rot[1] = 1;
		rot[2] = 0;

		sym[0] = 1;
		sym[1] = 1;
		sym[2] = 0;

		ord = RotateHex(rot,sym);

		for(int node_id=0;node_id<8;node_id++){
			id_node[node_id] = elem->nodes[ord[node_id]].id;
		}

	}else if(elem->pad==48){
		//edge 4 5 9 11
		rot[0] = 0;
		rot[1] = 1;
		rot[2] = 0;

		sym[0] = 1;
		sym[1] = 0;
		sym[2] = 0;

		ord = RotateHex(rot,sym);

		for(int node_id=0;node_id<8;node_id++){
			id_node[node_id] = elem->nodes[ord[node_id]].id;
		}

	}else{
		printf ("Error in template 4\n");
		exit (EXIT_FAILURE);
	}


	for(int i=0;i<5;i++){
		int conn_p[8];
		GtsPoint* point[8]={NULL};

		double cord_in_x[8],cord_in_y[8],cord_in_z[8];
		//add the nodes in the coord vector
		for (int ii = 0; ii < 8; ii++){
			cord_in_x[ii]=coords[3*id_node[ii]] ;
			cord_in_y[ii]=coords[3*id_node[ii]+1] ;
			cord_in_z[ii]=coords[3*id_node[ii]+2] ;
			//fprintf(mesh->fdbg,"coord in: %f, %f, %f, in the node: %d\n",cord_in_x[ii],cord_in_y[ii],cord_in_z[ii],elem->nodes[ii].id);
		}

		//add the new nodes in the
		for(int ii=0;ii<8;ii++){
			conn_p[ii] = 0;

			if((local_ref[i][ii][0]==1 || local_ref[i][ii][0]==-1) &&
					(local_ref[i][ii][1]==1 || local_ref[i][ii][1]==-1) &&
					(local_ref[i][ii][2]==1 || local_ref[i][ii][2]==-1)){

				if(local_ref[i][ii][0]==-1 && local_ref[i][ii][1]==-1 && local_ref[i][ii][2]==-1){
					conn_p[ii] = id_node[0];
				}else if(local_ref[i][ii][0]==1 && local_ref[i][ii][1]==-1 && local_ref[i][ii][2]==-1){
					conn_p[ii] = id_node[1];
				}else if(local_ref[i][ii][0]==1 && local_ref[i][ii][1]==1 && local_ref[i][ii][2]==-1){
					conn_p[ii] = id_node[2];
				}else if(local_ref[i][ii][0]==-1 && local_ref[i][ii][1]==1 && local_ref[i][ii][2]==-1){
					conn_p[ii] = id_node[3];
				}else if(local_ref[i][ii][0]==-1 && local_ref[i][ii][1]==-1 && local_ref[i][ii][2]==1){
					conn_p[ii] = id_node[4];
				}else if(local_ref[i][ii][0]==1 && local_ref[i][ii][1]==-1 && local_ref[i][ii][2]==1){
					conn_p[ii] = id_node[5];
				}else if(local_ref[i][ii][0]==1 && local_ref[i][ii][1]==1 && local_ref[i][ii][2]==1){
					conn_p[ii] = id_node[6];
				}else if(local_ref[i][ii][0]==-1 && local_ref[i][ii][1]==1 && local_ref[i][ii][2]==1){
					conn_p[ii] = id_node[7];
				}
				//fprintf(mesh->fdbg,"coord out: %f, %f, %f, in the node: %d\n",coords[3*conn_p[ii]],coords[3*conn_p[ii]+1],coords[3*conn_p[ii]+2],conn_p[ii]);
			}else{
				cord_in_ref[0] = local_ref[i][ii][0];
				cord_in_ref[1] = local_ref[i][ii][1];
				cord_in_ref[2] = local_ref[i][ii][2];

				point[ii] = LinearMapHex(cord_in_ref, cord_in_x,cord_in_y,cord_in_z);
				conn_p[ii] = AddPoint( mesh, hash_nodes, point[ii] , coords);
				//fprintf(mesh->fdbg,"coord out: %f, %f, %f, in the node: %d\n",var[0],var[1],var[2],conn_p[ii]);

			}
		}

		int aux[8];
		for(int ino=0;ino<8;ino++){
			aux[ino] = conn_p[ino];
		}
		for(int ino=0;ino<8;ino++){
			conn_p[ord[ino]] = aux[ino];
		}

		//add the new elements and make the connectivity
		if(i==0){
			octant_t *elem1 = (octant_t*) sc_array_index(&mesh->elements, elements_ids[iel]);

			elem1->nodes[0].id = conn_p[0];
			elem1->nodes[1].id = conn_p[1];
			elem1->nodes[2].id = conn_p[2];
			elem1->nodes[3].id = conn_p[3];

			elem1->nodes[4].id = conn_p[4];
			elem1->nodes[5].id = conn_p[5];
			elem1->nodes[6].id = conn_p[6];
			elem1->nodes[7].id = conn_p[7];

			CopyPropEl(mesh,id,elem1);

		}else{
			octant_t* elem2 = (octant_t*) sc_array_push(&mesh->elements);

			elem2->nodes[0].id = conn_p[0];
			elem2->nodes[1].id = conn_p[1];
			elem2->nodes[2].id = conn_p[2];
			elem2->nodes[3].id = conn_p[3];

			elem2->nodes[4].id = conn_p[4];
			elem2->nodes[5].id = conn_p[5];
			elem2->nodes[6].id = conn_p[6];
			elem2->nodes[7].id = conn_p[7];

			CopyPropEl(mesh,id,elem2);

		}
	}
	sc_array_reset(&toto);
}

void ApplyTemplate5(hexa_tree_t* mesh, std::vector<double>& coords, std::vector<int>& elements_ids, int iel, sc_hash_array_t* hash_nodes, double step){

	sc_array_t toto;
	sc_array_init(&toto, sizeof(octant_t));

	octant_t* elemOrig = (octant_t*) sc_array_index(&mesh->elements, elements_ids[iel]);
	octant_t * elem = (octant_t*) sc_array_push(&toto);
	hexa_element_copy(elemOrig,elem);

	int id = elements_ids[iel];

	double cord_in_ref[3];
	cord_in_ref[0] = 0;
	cord_in_ref[1] = 0;
	cord_in_ref[2] = 0;

	//reference element edge 0
	double local_ref[13][8][3];
	if(true){
		//element 0
		local_ref[0][0][0] = -1;
		local_ref[0][0][1] = -1;
		local_ref[0][0][2] = -1;

		local_ref[0][1][0] = 1;
		local_ref[0][1][1] = -1;
		local_ref[0][1][2] = -1;

		local_ref[0][2][0] = 1;
		local_ref[0][2][1] = 1;
		local_ref[0][2][2] = -1;

		local_ref[0][3][0] = -1;
		local_ref[0][3][1] = 1;
		local_ref[0][3][2] = -1;

		local_ref[0][4][0] = -1;
		local_ref[0][4][1] = -1+step;
		local_ref[0][4][2] = -1+step;

		local_ref[0][5][0] = 1;
		local_ref[0][5][1] = -1+step;
		local_ref[0][5][2] = -1+step;

		local_ref[0][6][0] = 1;
		local_ref[0][6][1] = -1+2*step;
		local_ref[0][6][2] = -1+step;

		local_ref[0][7][0] = -1;
		local_ref[0][7][1] = -1+2*step;
		local_ref[0][7][2] = -1+step;

		//element 1
		local_ref[1][0][0] = -1;
		local_ref[1][0][1] = -1+step;
		local_ref[1][0][2] = -1+step;

		local_ref[1][1][0] = 1;
		local_ref[1][1][1] = -1+step;
		local_ref[1][1][2] = -1+step;

		local_ref[1][2][0] = 1;
		local_ref[1][2][1] = -1+2*step;
		local_ref[1][2][2] = -1+step;

		local_ref[1][3][0] = -1;
		local_ref[1][3][1] = -1+2*step;
		local_ref[1][3][2] = -1+step;

		local_ref[1][4][0] = -1+step;
		local_ref[1][4][1] = -1+step;
		local_ref[1][4][2] = -1+2*step;

		local_ref[1][5][0] = -1+2*step;
		local_ref[1][5][1] = -1+step;
		local_ref[1][5][2] = -1+2*step;

		local_ref[1][6][0] = -1+2*step;
		local_ref[1][6][1] = -1+2*step;
		local_ref[1][6][2] = -1+2*step;

		local_ref[1][7][0] = -1+step;
		local_ref[1][7][1] = -1+2*step;
		local_ref[1][7][2] = -1+2*step;

		//element 2
		local_ref[2][0][0] = -1+step;
		local_ref[2][0][1] = -1+step;
		local_ref[2][0][2] = -1+2*step;

		local_ref[2][1][0] = -1+2*step;
		local_ref[2][1][1] = -1+step;
		local_ref[2][1][2] = -1+2*step;

		local_ref[2][2][0] = -1+2*step;
		local_ref[2][2][1] = -1+2*step;
		local_ref[2][2][2] = -1+2*step;

		local_ref[2][3][0] = -1+step;
		local_ref[2][3][1] = -1+2*step;
		local_ref[2][3][2] = -1+2*step;

		local_ref[2][4][0] = -1+step;
		local_ref[2][4][1] = -1+step;
		local_ref[2][4][2] = 1;

		local_ref[2][5][0] = -1+2*step;
		local_ref[2][5][1] = -1+step;
		local_ref[2][5][2] = 1;

		local_ref[2][6][0] = -1+2*step;
		local_ref[2][6][1] = -1+2*step;
		local_ref[2][6][2] = 1;

		local_ref[2][7][0] = -1+step;
		local_ref[2][7][1] = -1+2*step;
		local_ref[2][7][2] = 1;

		//element 3
		local_ref[3][0][0] = -1;
		local_ref[3][0][1] = -1;
		local_ref[3][0][2] = -1;

		local_ref[3][1][0] =  1;
		local_ref[3][1][1] = -1;
		local_ref[3][1][2] = -1;

		local_ref[3][2][0] =  1;
		local_ref[3][2][1] = -1+step;
		local_ref[3][2][2] = -1+step;

		local_ref[3][3][0] = -1;
		local_ref[3][3][1] = -1+step;
		local_ref[3][3][2] = -1+step;

		local_ref[3][4][0] = -1+step;
		local_ref[3][4][1] = -1;
		local_ref[3][4][2] = -1+step;

		local_ref[3][5][0] = -1+2*step;
		local_ref[3][5][1] = -1;
		local_ref[3][5][2] = -1+step;

		local_ref[3][6][0] = -1+2*step;
		local_ref[3][6][1] = -1+1*step;
		local_ref[3][6][2] = -1+2*step;

		local_ref[3][7][0] = -1+step;
		local_ref[3][7][1] = -1+1*step;
		local_ref[3][7][2] = -1+2*step;

		//element 4
		local_ref[4][0][0] = -1;
		local_ref[4][0][1] = -1+2*step;
		local_ref[4][0][2] = -1+step;

		local_ref[4][1][0] =  1;
		local_ref[4][1][1] = -1+2*step;
		local_ref[4][1][2] = -1+step;

		local_ref[4][2][0] =  1;
		local_ref[4][2][1] =  1;
		local_ref[4][2][2] = -1;

		local_ref[4][3][0] = -1;
		local_ref[4][3][1] =  1;
		local_ref[4][3][2] = -1;

		local_ref[4][4][0] = -1+step;
		local_ref[4][4][1] = -1+2*step;
		local_ref[4][4][2] = -1+2*step;

		local_ref[4][5][0] = -1+2*step;
		local_ref[4][5][1] = -1+2*step;
		local_ref[4][5][2] = -1+2*step;

		local_ref[4][6][0] = -1+2*step;
		local_ref[4][6][1] = 1;
		local_ref[4][6][2] = -1+step;

		local_ref[4][7][0] = -1+step;
		local_ref[4][7][1] = 1;
		local_ref[4][7][2] = -1+step;

		//element 5
		local_ref[5][0][0] = -1;
		local_ref[5][0][1] = -1;
		local_ref[5][0][2] = -1;

		local_ref[5][1][0] = -1+step;
		local_ref[5][1][1] = -1;
		local_ref[5][1][2] = -1+step;

		local_ref[5][2][0] = -1+step;
		local_ref[5][2][1] = -1+step;
		local_ref[5][2][2] = -1+2*step;

		local_ref[5][3][0] = -1;
		local_ref[5][3][1] = -1+step;
		local_ref[5][3][2] = -1+step;

		local_ref[5][4][0] = -1;
		local_ref[5][4][1] = -1;
		local_ref[5][4][2] = 1;

		local_ref[5][5][0] = -1+step;
		local_ref[5][5][1] = -1;
		local_ref[5][5][2] = 1;

		local_ref[5][6][0] = -1+step;
		local_ref[5][6][1] = -1+step;
		local_ref[5][6][2] = 1;

		local_ref[5][7][0] = -1;
		local_ref[5][7][1] = -1+step;
		local_ref[5][7][2] = 1;

		//element 6
		local_ref[6][0][0] = -1+step;
		local_ref[6][0][1] = -1;
		local_ref[6][0][2] = -1+1*step;

		local_ref[6][1][0] = -1+2*step;
		local_ref[6][1][1] = -1;
		local_ref[6][1][2] = -1+step;

		local_ref[6][2][0] = -1+2*step;
		local_ref[6][2][1] = -1+step;
		local_ref[6][2][2] = -1+2*step;

		local_ref[6][3][0] = -1+step;
		local_ref[6][3][1] = -1+step;
		local_ref[6][3][2] = -1+2*step;

		local_ref[6][4][0] = -1+step;
		local_ref[6][4][1] = -1;
		local_ref[6][4][2] = 1;

		local_ref[6][5][0] = -1+2*step;
		local_ref[6][5][1] = -1;
		local_ref[6][5][2] = 1;

		local_ref[6][6][0] = -1+2*step;
		local_ref[6][6][1] = -1+1*step;
		local_ref[6][6][2] = 1;

		local_ref[6][7][0] = -1+step;
		local_ref[6][7][1] = -1+1*step;
		local_ref[6][7][2] = 1;

		//element 7
		local_ref[7][0][0] = -1+2*step;
		local_ref[7][0][1] = -1;
		local_ref[7][0][2] = -1+step;

		local_ref[7][1][0] = 1;
		local_ref[7][1][1] = -1;
		local_ref[7][1][2] = -1;

		local_ref[7][2][0] =  1;
		local_ref[7][2][1] = -1+step;
		local_ref[7][2][2] = -1+1*step;

		local_ref[7][3][0] = -1+2*step;
		local_ref[7][3][1] = -1+1*step;
		local_ref[7][3][2] = -1+2*step;

		local_ref[7][4][0] = -1+2*step;
		local_ref[7][4][1] = -1;
		local_ref[7][4][2] = 1;

		local_ref[7][5][0] = 1;
		local_ref[7][5][1] = -1;
		local_ref[7][5][2] = 1;

		local_ref[7][6][0] = 1;
		local_ref[7][6][1] = -1+1*step;
		local_ref[7][6][2] = 1;

		local_ref[7][7][0] = -1+2*step;
		local_ref[7][7][1] = -1+1*step;
		local_ref[7][7][2] = 1;

		//element 8
		local_ref[8][0][0] = -1+2*step;
		local_ref[8][0][1] = -1+1*step;
		local_ref[8][0][2] = -1+2*step;

		local_ref[8][1][0] = 1;
		local_ref[8][1][1] = -1+1*step;
		local_ref[8][1][2] = -1+1*step;

		local_ref[8][2][0] =  1;
		local_ref[8][2][1] = -1+2*step;
		local_ref[8][2][2] = -1+1*step;

		local_ref[8][3][0] = -1+2*step;
		local_ref[8][3][1] = -1+2*step;
		local_ref[8][3][2] = -1+2*step;

		local_ref[8][4][0] = -1+2*step;
		local_ref[8][4][1] = -1+1*step;
		local_ref[8][4][2] = 1;

		local_ref[8][5][0] = 1;
		local_ref[8][5][1] = -1+1*step;
		local_ref[8][5][2] = 1;

		local_ref[8][6][0] = 1;
		local_ref[8][6][1] = -1+2*step;
		local_ref[8][6][2] = 1;

		local_ref[8][7][0] = -1+2*step;
		local_ref[8][7][1] = -1+2*step;
		local_ref[8][7][2] = 1;

		//element 9
		local_ref[9][0][0] = -1+2*step;
		local_ref[9][0][1] = -1+2*step;
		local_ref[9][0][2] = -1+2*step;

		local_ref[9][1][0] = 1;
		local_ref[9][1][1] = -1+2*step;
		local_ref[9][1][2] = -1+1*step;

		local_ref[9][2][0] =  1;
		local_ref[9][2][1] =  1;
		local_ref[9][2][2] = -1;

		local_ref[9][3][0] = -1+2*step;
		local_ref[9][3][1] =  1;
		local_ref[9][3][2] = -1+1*step;

		local_ref[9][4][0] = -1+2*step;
		local_ref[9][4][1] = -1+2*step;
		local_ref[9][4][2] = 1;

		local_ref[9][5][0] = 1;
		local_ref[9][5][1] = -1+2*step;
		local_ref[9][5][2] = 1;

		local_ref[9][6][0] = 1;
		local_ref[9][6][1] = 1;
		local_ref[9][6][2] = 1;

		local_ref[9][7][0] = -1+2*step;
		local_ref[9][7][1] = 1;
		local_ref[9][7][2] = 1;

		//element 10
		local_ref[10][0][0] = -1+1*step;
		local_ref[10][0][1] = -1+2*step;
		local_ref[10][0][2] = -1+2*step;

		local_ref[10][1][0] = -1+2*step;
		local_ref[10][1][1] = -1+2*step;
		local_ref[10][1][2] = -1+2*step;

		local_ref[10][2][0] =  -1+2*step;
		local_ref[10][2][1] =  1;
		local_ref[10][2][2] = -1+step;

		local_ref[10][3][0] = -1+1*step;
		local_ref[10][3][1] =  1;
		local_ref[10][3][2] = -1+1*step;

		local_ref[10][4][0] = -1+1*step;
		local_ref[10][4][1] = -1+2*step;
		local_ref[10][4][2] = 1;

		local_ref[10][5][0] = -1+2*step;
		local_ref[10][5][1] = -1+2*step;
		local_ref[10][5][2] = 1;

		local_ref[10][6][0] = -1+2*step;
		local_ref[10][6][1] = 1;
		local_ref[10][6][2] = 1;

		local_ref[10][7][0] = -1+1*step;
		local_ref[10][7][1] = 1;
		local_ref[10][7][2] = 1;

		//element 11
		local_ref[11][0][0] = -1;
		local_ref[11][0][1] = -1+2*step;
		local_ref[11][0][2] = -1+1*step;

		local_ref[11][1][0] = -1+1*step;
		local_ref[11][1][1] = -1+2*step;
		local_ref[11][1][2] = -1+2*step;

		local_ref[11][2][0] =  -1+1*step;
		local_ref[11][2][1] =  1;
		local_ref[11][2][2] = -1+step;

		local_ref[11][3][0] = -1;
		local_ref[11][3][1] =  1;
		local_ref[11][3][2] = -1;

		local_ref[11][4][0] = -1;
		local_ref[11][4][1] = -1+2*step;
		local_ref[11][4][2] = 1;

		local_ref[11][5][0] = -1+1*step;
		local_ref[11][5][1] = -1+2*step;
		local_ref[11][5][2] = 1;

		local_ref[11][6][0] = -1+1*step;
		local_ref[11][6][1] = 1;
		local_ref[11][6][2] = 1;

		local_ref[11][7][0] = -1;
		local_ref[11][7][1] = 1;
		local_ref[11][7][2] = 1;

		//element 12
		local_ref[12][0][0] = -1;
		local_ref[12][0][1] = -1+1*step;
		local_ref[12][0][2] = -1+1*step;

		local_ref[12][1][0] = -1+1*step;
		local_ref[12][1][1] = -1+1*step;
		local_ref[12][1][2] = -1+2*step;

		local_ref[12][2][0] =  -1+1*step;
		local_ref[12][2][1] =  -1+2*step;
		local_ref[12][2][2] = -1+2*step;

		local_ref[12][3][0] = -1;
		local_ref[12][3][1] =  -1+2*step;
		local_ref[12][3][2] = -1+step;

		local_ref[12][4][0] = -1;
		local_ref[12][4][1] = -1+1*step;
		local_ref[12][4][2] = 1;

		local_ref[12][5][0] = -1+1*step;
		local_ref[12][5][1] = -1+1*step;
		local_ref[12][5][2] = 1;

		local_ref[12][6][0] = -1+1*step;
		local_ref[12][6][1] = -1+2*step;
		local_ref[12][6][2] = 1;

		local_ref[12][7][0] = -1;
		local_ref[12][7][1] = -1+2*step;
		local_ref[12][7][2] = 1;

	}

	//define the rotation of the reference element
	int rot[3];
	int sym[3];
	int id_node[8];
	vector<int> ord;

	if(elem->pad==49 || elem->pad==55 || elem->pad==56 || elem->pad==57 || elem->pad==58
			|| elem->pad==79 || elem->pad==80 || elem->pad==81 || elem->pad==82){
		//edge 0 1 2 3
		rot[0] = 0;
		rot[1] = 0;
		rot[2] = 0;

		sym[0] = 0;
		sym[1] = 0;
		sym[2] = 1;

		ord = RotateHex(rot,sym);

		for(int node_id=0;node_id<8;node_id++){
			id_node[node_id] = elem->nodes[ord[node_id]].id;
		}
	}else if(elem->pad==50 || elem->pad==59 || elem->pad==60 || elem->pad==61 || elem->pad==62
			|| elem->pad==83 || elem->pad==84 || elem->pad==85 || elem->pad==86){
		//edge 8 9 10 11
		rot[0] = 0;
		rot[1] = 0;
		rot[2] = 0;

		sym[0] = 0;
		sym[1] = 0;
		sym[2] = 0;

		ord = RotateHex(rot,sym);

		for(int node_id=0;node_id<8;node_id++){
			id_node[node_id] = elem->nodes[ord[node_id]].id;
		}
	} else if(elem->pad==51 || elem->pad==63 || elem->pad==64 || elem->pad==65 || elem->pad==66
			|| elem->pad==87 || elem->pad==88 || elem->pad==89 || elem->pad==90){
		//edge 0 4 5 8
		rot[0] = -1;
		rot[1] = 0;
		rot[2] = 0;

		sym[0] = 0;
		sym[1] = 0;
		sym[2] = 0;

		ord = RotateHex(rot,sym);

		for(int node_id=0;node_id<8;node_id++){
			id_node[node_id] = elem->nodes[ord[node_id]].id;
		}
	}else if(elem->pad==52 || elem->pad==67 || elem->pad==68 || elem->pad==69 || elem->pad==70
			|| elem->pad==91 || elem->pad==92 || elem->pad==93 || elem->pad==94){
		//edge 1 5 6 9
		rot[0] = 0;
		rot[1] = -1;
		rot[2] = 0;

		sym[0] = 0;
		sym[1] = 0;
		sym[2] = 0;

		ord = RotateHex(rot,sym);

		for(int node_id=0;node_id<8;node_id++){
			id_node[node_id] = elem->nodes[ord[node_id]].id;
		}
	}else if(elem->pad==53 || elem->pad==71 || elem->pad==72 || elem->pad==73 || elem->pad==74
			|| elem->pad==95 || elem->pad==96 || elem->pad==97 || elem->pad==98){
		//edge 2 6 7 10
		rot[0] = 1;
		rot[1] = 0;
		rot[2] = 0;

		sym[0] = 0;
		sym[1] = 0;
		sym[2] = 0;

		ord = RotateHex(rot,sym);

		for(int node_id=0;node_id<8;node_id++){
			id_node[node_id] = elem->nodes[ord[node_id]].id;
		}
	}else if(elem->pad==54 || elem->pad==75 || elem->pad==76 || elem->pad==77 || elem->pad==78
			|| elem->pad==99 || elem->pad==100 || elem->pad==101 || elem->pad==102){
		//edge 3 4 7 11
		rot[0] = 0;
		rot[1] = 1;
		rot[2] = 0;

		sym[0] = 0;
		sym[1] = 0;
		sym[2] = 0;

		ord = RotateHex(rot,sym);

		for(int node_id=0;node_id<8;node_id++){
			id_node[node_id] = elem->nodes[ord[node_id]].id;
		}
	}else{
		printf ("Error in template 5\n");
		exit (EXIT_FAILURE);
	}



	for(int i=0;i<13;i++){
		int conn_p[8];
		GtsPoint* point[8]={NULL};

		double cord_in_x[8],cord_in_y[8],cord_in_z[8];
		//add the nodes in the coord vector
		for (int ii = 0; ii < 8; ii++){
			cord_in_x[ii]=coords[3*id_node[ii]] ;
			cord_in_y[ii]=coords[3*id_node[ii]+1] ;
			cord_in_z[ii]=coords[3*id_node[ii]+2] ;
			////fprintf(mesh->fdbg,"coord in: %f, %f, %f, in the node: %d\n",cord_in_x[ii],cord_in_y[ii],cord_in_z[ii],elem->nodes[ii].id);
		}

		//add the new nodes in the
		for(int ii=0;ii<8;ii++){
			conn_p[ii] = 0;

			if((local_ref[i][ii][0]==1 || local_ref[i][ii][0]==-1) &&
					(local_ref[i][ii][1]==1 || local_ref[i][ii][1]==-1) &&
					(local_ref[i][ii][2]==1 || local_ref[i][ii][2]==-1)){

				if(local_ref[i][ii][0]==-1 && local_ref[i][ii][1]==-1 && local_ref[i][ii][2]==-1){
					conn_p[ii] = id_node[0];
				}else if(local_ref[i][ii][0]==1 && local_ref[i][ii][1]==-1 && local_ref[i][ii][2]==-1){
					conn_p[ii] = id_node[1];
				}else if(local_ref[i][ii][0]==1 && local_ref[i][ii][1]==1 && local_ref[i][ii][2]==-1){
					conn_p[ii] = id_node[2];
				}else if(local_ref[i][ii][0]==-1 && local_ref[i][ii][1]==1 && local_ref[i][ii][2]==-1){
					conn_p[ii] = id_node[3];
				}else if(local_ref[i][ii][0]==-1 && local_ref[i][ii][1]==-1 && local_ref[i][ii][2]==1){
					conn_p[ii] = id_node[4];
				}else if(local_ref[i][ii][0]==1 && local_ref[i][ii][1]==-1 && local_ref[i][ii][2]==1){
					conn_p[ii] = id_node[5];
				}else if(local_ref[i][ii][0]==1 && local_ref[i][ii][1]==1 && local_ref[i][ii][2]==1){
					conn_p[ii] = id_node[6];
				}else if(local_ref[i][ii][0]==-1 && local_ref[i][ii][1]==1 && local_ref[i][ii][2]==1){
					conn_p[ii] = id_node[7];
				}
				//fprintf(mesh->fdbg,"coord out: %f, %f, %f, in the node: %d\n",coords[3*conn_p[ii]],coords[3*conn_p[ii]+1],coords[3*conn_p[ii]+2],conn_p[ii]);
			}else{
				cord_in_ref[0] = local_ref[i][ii][0];
				cord_in_ref[1] = local_ref[i][ii][1];
				cord_in_ref[2] = local_ref[i][ii][2];

				point[ii] = LinearMapHex(cord_in_ref, cord_in_x,cord_in_y,cord_in_z);
				double var[3];
				var[0] = point[ii]->x;
				var[1] = point[ii]->y;
				var[2] = point[ii]->z;
				conn_p[ii] = AddPoint( mesh, hash_nodes, point[ii] , coords);
				//fprintf(mesh->fdbg,"coord out: %f, %f, %f, in the node: %d\n",var[0],var[1],var[2],conn_p[ii]);

			}
		}

		int aux[8];
		for(int ino=0;ino<8;ino++){
			aux[ino] = conn_p[ino];
		}
		for(int ino=0;ino<8;ino++){
			conn_p[ord[ino]] = aux[ino];
		}

		//add the new elements and make the connectivity
		if(i==0){
			octant_t *elem1 = (octant_t*) sc_array_index(&mesh->elements, elements_ids[iel]);

			elem1->nodes[0].id = conn_p[0];
			elem1->nodes[1].id = conn_p[1];
			elem1->nodes[2].id = conn_p[2];
			elem1->nodes[3].id = conn_p[3];

			elem1->nodes[4].id = conn_p[4];
			elem1->nodes[5].id = conn_p[5];
			elem1->nodes[6].id = conn_p[6];
			elem1->nodes[7].id = conn_p[7];

			CopyPropEl(mesh,id,elem1);

		}else{
			octant_t* elem2 = (octant_t*) sc_array_push(&mesh->elements);

			elem2->nodes[0].id = conn_p[0];
			elem2->nodes[1].id = conn_p[1];
			elem2->nodes[2].id = conn_p[2];
			elem2->nodes[3].id = conn_p[3];

			elem2->nodes[4].id = conn_p[4];
			elem2->nodes[5].id = conn_p[5];
			elem2->nodes[6].id = conn_p[6];
			elem2->nodes[7].id = conn_p[7];

			CopyPropEl(mesh,id,elem2);

		}

	}
	sc_array_reset(&toto);
}

void ApplyTemplate6(hexa_tree_t* mesh, std::vector<double>& coords, std::vector<int>& elements_ids, int iel, sc_hash_array_t* hash_nodes, double step){

	sc_array_t toto;
	sc_array_init(&toto, sizeof(octant_t));

	octant_t* elemOrig = (octant_t*) sc_array_index(&mesh->elements, elements_ids[iel]);
	octant_t * elem = (octant_t*) sc_array_push(&toto);
	hexa_element_copy(elemOrig,elem);

	int id = elements_ids[iel];

	double cord_in_ref[3];
	cord_in_ref[0] = 0;
	cord_in_ref[1] = 0;
	cord_in_ref[2] = 0;

	//reference element edge 0
	double local_ref[10][8][3];
	if(true){
		//element 0
		local_ref[0][0][0] = -1;
		local_ref[0][0][1] = -1;
		local_ref[0][0][2] = -1;

		local_ref[0][1][0] = -1+step;
		local_ref[0][1][1] = -1;
		local_ref[0][1][2] = -1;

		local_ref[0][2][0] = -1+2*3*step/double(4);
		local_ref[0][2][1] = -1+step;
		local_ref[0][2][2] = -1;

		local_ref[0][3][0] = -1+3*step/double(4);
		local_ref[0][3][1] = -1+step;
		local_ref[0][3][2] = -1;

		local_ref[0][4][0] = -1;
		local_ref[0][4][1] = -1;
		local_ref[0][4][2] = 1;

		local_ref[0][5][0] = -1+step;
		local_ref[0][5][1] = -1;
		local_ref[0][5][2] = 1;

		local_ref[0][6][0] = -1+2*3*step/double(4);
		local_ref[0][6][1] = -1+step;
		local_ref[0][6][2] = 1;

		local_ref[0][7][0] = -1+3*step/double(4);
		local_ref[0][7][1] = -1+step;
		local_ref[0][7][2] = 1;

		//element 1
		local_ref[1][0][0] = -1+1*step;
		local_ref[1][0][1] = -1;
		local_ref[1][0][2] = -1;

		local_ref[1][1][0] = -1+2*step;
		local_ref[1][1][1] = -1;
		local_ref[1][1][2] = -1;

		local_ref[1][2][0] = -1+3*3*step/double(4);
		local_ref[1][2][1] = -1+step;
		local_ref[1][2][2] = -1;

		local_ref[1][3][0] = -1+2*3*step/double(4);
		local_ref[1][3][1] = -1+step;
		local_ref[1][3][2] = -1;

		local_ref[1][4][0] = -1+1*step;
		local_ref[1][4][1] = -1;
		local_ref[1][4][2] = 1;

		local_ref[1][5][0] = -1+2*step;
		local_ref[1][5][1] = -1;
		local_ref[1][5][2] = 1;

		local_ref[1][6][0] = -1+3*3*step/double(4);
		local_ref[1][6][1] = -1+step;
		local_ref[1][6][2] = 1;

		local_ref[1][7][0] = -1+2*3*step/double(4);
		local_ref[1][7][1] = -1+step;
		local_ref[1][7][2] = 1;


		//element 2
		local_ref[2][0][0] = -1+2*step;
		local_ref[2][0][1] = -1;
		local_ref[2][0][2] = -1;

		local_ref[2][1][0] =  1;
		local_ref[2][1][1] = -1;
		local_ref[2][1][2] = -1;

		local_ref[2][2][0] =  1;
		local_ref[2][2][1] =  -1+step;
		local_ref[2][2][2] = -1;

		local_ref[2][3][0] = -1+3*3*step/double(4);
		local_ref[2][3][1] = -1+step;
		local_ref[2][3][2] = -1;

		local_ref[2][4][0] = -1+2*step;
		local_ref[2][4][1] = -1;
		local_ref[2][4][2] = 1;

		local_ref[2][5][0] =  1;
		local_ref[2][5][1] = -1;
		local_ref[2][5][2] =  1;

		local_ref[2][6][0] = 1;
		local_ref[2][6][1] = -1+step;
		local_ref[2][6][2] = 1;

		local_ref[2][7][0] = -1+3*3*step/double(4);
		local_ref[2][7][1] = -1+step;
		local_ref[2][7][2] = 1;

		//element 3
		local_ref[3][0][0] = -1+1*3*step/double(4);
		local_ref[3][0][1] = -1+step;
		local_ref[3][0][2] = -1;

		local_ref[3][1][0] = -1+2*3*step/double(4);
		local_ref[3][1][1] = -1+step;
		local_ref[3][1][2] = -1;

		local_ref[3][2][0] = -1+2*3*step/double(4);
		local_ref[3][2][1] = -1+2*step;
		local_ref[3][2][2] = -1;

		local_ref[3][3][0] = -1+1*3*step/double(4);
		local_ref[3][3][1] = -1+2*step;
		local_ref[3][3][2] = -1;

		local_ref[3][4][0] = -1+1*3*step/double(4);
		local_ref[3][4][1] = -1+step;
		local_ref[3][4][2] = 1;

		local_ref[3][5][0] = -1+2*3*step/double(4);
		local_ref[3][5][1] = -1+step;
		local_ref[3][5][2] = 1;

		local_ref[3][6][0] = -1+2*3*step/double(4);
		local_ref[3][6][1] = -1+2*step;
		local_ref[3][6][2] = 1;

		local_ref[3][7][0] = -1+1*3*step/double(4);
		local_ref[3][7][1] = -1+2*step;
		local_ref[3][7][2] = 1;

		//element 4
		local_ref[4][0][0] = -1+2*3*step/double(4);
		local_ref[4][0][1] = -1+step;
		local_ref[4][0][2] = -1;

		local_ref[4][1][0] = -1+3*3*step/double(4);
		local_ref[4][1][1] = -1+step;
		local_ref[4][1][2] = -1;

		local_ref[4][2][0] = -1+3*3*step/double(4);
		local_ref[4][2][1] = -1+2*step;
		local_ref[4][2][2] = -1;

		local_ref[4][3][0] = -1+2*3*step/double(4);
		local_ref[4][3][1] = -1+2*step;
		local_ref[4][3][2] = -1;

		local_ref[4][4][0] = -1+2*3*step/double(4);
		local_ref[4][4][1] = -1+step;
		local_ref[4][4][2] = 1;

		local_ref[4][5][0] = -1+3*3*step/double(4);
		local_ref[4][5][1] = -1+step;
		local_ref[4][5][2] = 1;

		local_ref[4][6][0] = -1+3*3*step/double(4);
		local_ref[4][6][1] = -1+2*step;
		local_ref[4][6][2] = 1;

		local_ref[4][7][0] = -1+2*3*step/double(4);
		local_ref[4][7][1] = -1+2*step;
		local_ref[4][7][2] = 1;

		//element 5
		local_ref[5][0][0] = -1+3*3*step/double(4);
		local_ref[5][0][1] = -1+step;
		local_ref[5][0][2] = -1;

		local_ref[5][1][0] = 1;
		local_ref[5][1][1] = -1+step;
		local_ref[5][1][2] = -1;

		local_ref[5][2][0] = 1;
		local_ref[5][2][1] = -1+2*step;
		local_ref[5][2][2] = -1;

		local_ref[5][3][0] = -1+3*3*step/double(4);
		local_ref[5][3][1] = -1+2*step;
		local_ref[5][3][2] = -1;

		local_ref[5][4][0] = -1+3*3*step/double(4);;
		local_ref[5][4][1] = -1+step;
		local_ref[5][4][2] = 1;

		local_ref[5][5][0] = 1;
		local_ref[5][5][1] = -1+step;
		local_ref[5][5][2] = 1;

		local_ref[5][6][0] = 1;
		local_ref[5][6][1] = -1+2*step;
		local_ref[5][6][2] = 1;

		local_ref[5][7][0] = -1+3*3*step/double(4);
		local_ref[5][7][1] = -1+2*step;
		local_ref[5][7][2] = 1;

		//element 6
		local_ref[6][0][0] = -1+1*3*step/double(4);
		local_ref[6][0][1] = -1+2*step;
		local_ref[6][0][2] = -1;

		local_ref[6][1][0] = -1+2*3*step/double(4);
		local_ref[6][1][1] = -1+2*step;
		local_ref[6][1][2] = -1;

		local_ref[6][2][0] = -1+step;
		local_ref[6][2][1] = 1;
		local_ref[6][2][2] = -1;

		local_ref[6][3][0] = -1;
		local_ref[6][3][1] = 1;
		local_ref[6][3][2] = -1;

		local_ref[6][4][0] = -1+1*3*step/double(4);
		local_ref[6][4][1] = -1+2*step;
		local_ref[6][4][2] = 1;

		local_ref[6][5][0] = -1+2*3*step/double(4);
		local_ref[6][5][1] = -1+2*step;
		local_ref[6][5][2] = 1;

		local_ref[6][6][0] = -1+step;
		local_ref[6][6][1] = 1;
		local_ref[6][6][2] = 1;

		local_ref[6][7][0] = -1;
		local_ref[6][7][1] = 1;
		local_ref[6][7][2] = 1;

		//element 7
		local_ref[7][0][0] = -1+2*3*step/double(4);
		local_ref[7][0][1] = -1+2*step;
		local_ref[7][0][2] = -1;

		local_ref[7][1][0] = -1+3*3*step/double(4);
		local_ref[7][1][1] = -1+2*step;
		local_ref[7][1][2] = -1;

		local_ref[7][2][0] = -1+2*step;
		local_ref[7][2][1] = 1;
		local_ref[7][2][2] = -1;

		local_ref[7][3][0] = -1+step;
		local_ref[7][3][1] = 1;
		local_ref[7][3][2] = -1;

		local_ref[7][4][0] = -1+2*3*step/double(4);
		local_ref[7][4][1] = -1+2*step;
		local_ref[7][4][2] = 1;

		local_ref[7][5][0] = -1+3*3*step/double(4);
		local_ref[7][5][1] = -1+2*step;
		local_ref[7][5][2] = 1;

		local_ref[7][6][0] = -1+2*step;
		local_ref[7][6][1] = 1;
		local_ref[7][6][2] = 1;

		local_ref[7][7][0] = -1+step;
		local_ref[7][7][1] = 1;
		local_ref[7][7][2] = 1;


		//element 8
		local_ref[8][0][0] = -1+3*3*step/double(4);
		local_ref[8][0][1] = -1+2*step;
		local_ref[8][0][2] = -1;

		local_ref[8][1][0] = 1;
		local_ref[8][1][1] = -1+2*step;
		local_ref[8][1][2] = -1;

		local_ref[8][2][0] = 1;
		local_ref[8][2][1] = 1;
		local_ref[8][2][2] = -1;

		local_ref[8][3][0] = -1+2*step;
		local_ref[8][3][1] = 1;
		local_ref[8][3][2] = -1;

		local_ref[8][4][0] = -1+3*3*step/double(4);
		local_ref[8][4][1] = -1+2*step;
		local_ref[8][4][2] = 1;

		local_ref[8][5][0] = 1;
		local_ref[8][5][1] = -1+2*step;
		local_ref[8][5][2] = 1;

		local_ref[8][6][0] = 1;
		local_ref[8][6][1] = 1;
		local_ref[8][6][2] = 1;

		local_ref[8][7][0] = -1+2*step;
		local_ref[8][7][1] = 1;
		local_ref[8][7][2] = 1;

		//element 9
		local_ref[9][0][0] = -1;
		local_ref[9][0][1] = -1;
		local_ref[9][0][2] = -1;

		local_ref[9][1][0] = -1+1*3*step/double(4);
		local_ref[9][1][1] = -1+step;
		local_ref[9][1][2] = -1;

		local_ref[9][2][0] = -1+1*3*step/double(4);
		local_ref[9][2][1] = -1+2*step;
		local_ref[9][2][2] = -1;

		local_ref[9][3][0] = -1;
		local_ref[9][3][1] = 1;
		local_ref[9][3][2] = -1;

		local_ref[9][4][0] = -1;
		local_ref[9][4][1] = -1;
		local_ref[9][4][2] = 1;

		local_ref[9][5][0] = -1+1*3*step/double(4);
		local_ref[9][5][1] = -1+step;
		local_ref[9][5][2] = 1;

		local_ref[9][6][0] = -1+1*3*step/double(4);
		local_ref[9][6][1] = -1+2*step;
		local_ref[9][6][2] = 1;

		local_ref[9][7][0] = -1;
		local_ref[9][7][1] = 1;
		local_ref[9][7][2] = 1;
	}

	//define the rotation of the reference element
	int rot[3];
	int sym[3];
	int id_node[8];
	vector<int> ord;

	if(elem->pad==103){
		//edge 0 1 3 8 9 11
		rot[0] = 0;
		rot[1] = 0;
		rot[2] = 1;

		sym[0] = 0;
		sym[1] = 0;
		sym[2] = 0;

		ord = RotateHex(rot,sym);

		for(int node_id=0;node_id<8;node_id++){
			id_node[node_id] = elem->nodes[ord[node_id]].id;
		}

	}else if(elem->pad==104){
		//edge 0 1 2 8 9 10
		rot[0] = 0;
		rot[1] = 0;
		rot[2] = 0;

		sym[0] = 0;
		sym[1] = 0;
		sym[2] = 0;

		ord = RotateHex(rot,sym);

		for(int node_id=0;node_id<8;node_id++){
			id_node[node_id] = elem->nodes[ord[node_id]].id;
		}

	}else if(elem->pad==105){
		//edge 1 2 3 9 10 11
		rot[0] = 0;
		rot[1] = 0;
		rot[2] = -1;

		sym[0] = 0;
		sym[1] = 0;
		sym[2] = 0;

		ord = RotateHex(rot,sym);

		for(int node_id=0;node_id<8;node_id++){
			id_node[node_id] = elem->nodes[ord[node_id]].id;
		}

	}else if(elem->pad==106){
		//edge 0 2 3 8 10 11
		rot[0] = 0;
		rot[1] = 0;
		rot[2] = 0;

		sym[0] = 1;
		sym[1] = 0;
		sym[2] = 0;

		ord = RotateHex(rot,sym);

		for(int node_id=0;node_id<8;node_id++){
			id_node[node_id] = elem->nodes[ord[node_id]].id;
		}

	}else if(elem->pad==107){
		//edge 0 5 8 2 6 10
		rot[0] = 1;
		rot[1] = 0;
		rot[2] = 0;

		sym[0] = 0;
		sym[1] = 0;
		sym[2] = 0;

		ord = RotateHex(rot,sym);

		for(int node_id=0;node_id<8;node_id++){
			id_node[node_id] = elem->nodes[ord[node_id]].id;
		}

	}else if(elem->pad==108){
		//edge 5 8 4 6 10 7
		rot[0] = 1;
		rot[1] = 0;
		rot[2] = 1;

		sym[0] = 0;
		sym[1] = 0;
		sym[2] = 0;

		ord = RotateHex(rot,sym);

		for(int node_id=0;node_id<8;node_id++){
			id_node[node_id] = elem->nodes[ord[node_id]].id;
		}

	}else if(elem->pad==109){
		//edge 0 4 8 2 7 10
		rot[0] = 1;
		rot[1] = 0;
		rot[2] = 0;

		sym[0] = 1;
		sym[1] = 0;
		sym[2] = 0;

		ord = RotateHex(rot,sym);

		for(int node_id=0;node_id<8;node_id++){
			id_node[node_id] = elem->nodes[ord[node_id]].id;
		}

	}else if(elem->pad==110){
		//edge 5 0 4 6 2 7
		rot[0] = 1;
		rot[1] = 0;
		rot[2] = -1;

		sym[0] = 0;
		sym[1] = 0;
		sym[2] = 0;

		ord = RotateHex(rot,sym);

		for(int node_id=0;node_id<8;node_id++){
			id_node[node_id] = elem->nodes[ord[node_id]].id;
		}

	}else if(elem->pad==111){
		//edge 4 3 7 5 1 6
		rot[0] = 0;
		rot[1] = 1;
		rot[2] = 0;

		sym[0] = 1;
		sym[1] = 0;
		sym[2] = 0;

		ord = RotateHex(rot,sym);

		for(int node_id=0;node_id<8;node_id++){
			id_node[node_id] = elem->nodes[ord[node_id]].id;
		}

	}else if(elem->pad==112){
		//edge 1 5 9 3 4 11
		rot[0] = 1;
		rot[1] = 1;
		rot[2] = 0;

		sym[0] = 1;
		sym[1] = 0;
		sym[2] = 0;

		ord = RotateHex(rot,sym);

		for(int node_id=0;node_id<8;node_id++){
			id_node[node_id] = elem->nodes[ord[node_id]].id;
		}

	}else if(elem->pad==113){
		//edge 5 9 6 4 11 7
		rot[0] = 0;
		rot[1] = 1;
		rot[2] = 0;

		sym[0] = 0;
		sym[1] = 0;
		sym[2] = 0;

		ord = RotateHex(rot,sym);

		for(int node_id=0;node_id<8;node_id++){
			id_node[node_id] = elem->nodes[ord[node_id]].id;
		}

	}else if(elem->pad==114){
		//edge 1 6 9 3 7 11
		rot[0] = 1;
		rot[1] = 1;
		rot[2] = 0;

		sym[0] = 0;
		sym[1] = 0;
		sym[2] = 0;

		ord = RotateHex(rot,sym);

		for(int node_id=0;node_id<8;node_id++){
			id_node[node_id] = elem->nodes[ord[node_id]].id;
		}

	}else{
		printf ("Error in template 6\n");
		exit (EXIT_FAILURE);
	}


	for(int i=0;i<10;i++){
		int conn_p[8];
		GtsPoint* point[8]={NULL};

		double cord_in_x[8],cord_in_y[8],cord_in_z[8];
		//add the nodes in the coord vector
		for (int ii = 0; ii < 8; ii++){
			cord_in_x[ii]=coords[3*id_node[ii]] ;
			cord_in_y[ii]=coords[3*id_node[ii]+1] ;
			cord_in_z[ii]=coords[3*id_node[ii]+2] ;
			//fprintf(mesh->fdbg,"coord in: %f, %f, %f, in the node: %d\n",cord_in_x[ii],cord_in_y[ii],cord_in_z[ii],elem->nodes[ii].id);
		}

		//add the new nodes in the
		for(int ii=0;ii<8;ii++){
			conn_p[ii] = 0;

			if((local_ref[i][ii][0]==1 || local_ref[i][ii][0]==-1) &&
					(local_ref[i][ii][1]==1 || local_ref[i][ii][1]==-1) &&
					(local_ref[i][ii][2]==1 || local_ref[i][ii][2]==-1)){

				if(local_ref[i][ii][0]==-1 && local_ref[i][ii][1]==-1 && local_ref[i][ii][2]==-1){
					conn_p[ii] = id_node[0];
				}else if(local_ref[i][ii][0]==1 && local_ref[i][ii][1]==-1 && local_ref[i][ii][2]==-1){
					conn_p[ii] = id_node[1];
				}else if(local_ref[i][ii][0]==1 && local_ref[i][ii][1]==1 && local_ref[i][ii][2]==-1){
					conn_p[ii] = id_node[2];
				}else if(local_ref[i][ii][0]==-1 && local_ref[i][ii][1]==1 && local_ref[i][ii][2]==-1){
					conn_p[ii] = id_node[3];
				}else if(local_ref[i][ii][0]==-1 && local_ref[i][ii][1]==-1 && local_ref[i][ii][2]==1){
					conn_p[ii] = id_node[4];
				}else if(local_ref[i][ii][0]==1 && local_ref[i][ii][1]==-1 && local_ref[i][ii][2]==1){
					conn_p[ii] = id_node[5];
				}else if(local_ref[i][ii][0]==1 && local_ref[i][ii][1]==1 && local_ref[i][ii][2]==1){
					conn_p[ii] = id_node[6];
				}else if(local_ref[i][ii][0]==-1 && local_ref[i][ii][1]==1 && local_ref[i][ii][2]==1){
					conn_p[ii] = id_node[7];
				}
				//fprintf(mesh->fdbg,"coord out: %f, %f, %f, in the node: %d\n",coords[3*conn_p[ii]],coords[3*conn_p[ii]+1],coords[3*conn_p[ii]+2],conn_p[ii]);
			}else{
				cord_in_ref[0] = local_ref[i][ii][0];
				cord_in_ref[1] = local_ref[i][ii][1];
				cord_in_ref[2] = local_ref[i][ii][2];

				point[ii] = LinearMapHex(cord_in_ref, cord_in_x,cord_in_y,cord_in_z);
				double var[3];
				var[0] = point[ii]->x;
				var[1] = point[ii]->y;
				var[2] = point[ii]->z;
				conn_p[ii] = AddPoint( mesh, hash_nodes, point[ii] , coords);
				//fprintf(mesh->fdbg,"coord out: %f, %f, %f, in the node: %d\n",var[0],var[1],var[2],conn_p[ii]);

			}
		}

		int aux[8];
		for(int ino=0;ino<8;ino++){
			aux[ino] = conn_p[ino];
		}
		for(int ino=0;ino<8;ino++){
			conn_p[ord[ino]] = aux[ino];
		}

		//add the new elements and make the connectivity
		if(i==0){
			octant_t *elem1 = (octant_t*) sc_array_index(&mesh->elements, elements_ids[iel]);

			elem1->nodes[0].id = conn_p[0];
			elem1->nodes[1].id = conn_p[1];
			elem1->nodes[2].id = conn_p[2];
			elem1->nodes[3].id = conn_p[3];

			elem1->nodes[4].id = conn_p[4];
			elem1->nodes[5].id = conn_p[5];
			elem1->nodes[6].id = conn_p[6];
			elem1->nodes[7].id = conn_p[7];

			CopyPropEl(mesh,id,elem1);

		}else{
			octant_t* elem2 = (octant_t*) sc_array_push(&mesh->elements);

			elem2->nodes[0].id = conn_p[0];
			elem2->nodes[1].id = conn_p[1];
			elem2->nodes[2].id = conn_p[2];
			elem2->nodes[3].id = conn_p[3];

			elem2->nodes[4].id = conn_p[4];
			elem2->nodes[5].id = conn_p[5];
			elem2->nodes[6].id = conn_p[6];
			elem2->nodes[7].id = conn_p[7];

			CopyPropEl(mesh,id,elem2);

		}
	}
	sc_array_reset(&toto);
}

void ApplyTemplate7(hexa_tree_t* mesh, std::vector<double>& coords, std::vector<int>& elements_ids, int iel, sc_hash_array_t* hash_nodes, double step){

	sc_array_t toto;
	sc_array_init(&toto, sizeof(octant_t));

	octant_t* elemOrig = (octant_t*) sc_array_index(&mesh->elements, elements_ids[iel]);
	octant_t * elem = (octant_t*) sc_array_push(&toto);
	hexa_element_copy(elemOrig,elem);

	int id = elements_ids[iel];

	double cord_in_ref[3];
	cord_in_ref[0] = 0;
	cord_in_ref[1] = 0;
	cord_in_ref[2] = 0;

	//reference element edge 0
	double local_ref[7][8][3];
	if(true){
		//element 0
		local_ref[0][0][0] = -1;
		local_ref[0][0][1] = -1;
		local_ref[0][0][2] = -1;

		local_ref[0][1][0] = -1+step;
		local_ref[0][1][1] = -1;
		local_ref[0][1][2] = -1;

		local_ref[0][2][0] = -1+step;
		local_ref[0][2][1] = -1+step;
		local_ref[0][2][2] = -1;

		local_ref[0][3][0] = -1;
		local_ref[0][3][1] = -1+step;
		local_ref[0][3][2] = -1;

		local_ref[0][4][0] = -1;
		local_ref[0][4][1] = -1;
		local_ref[0][4][2] =  -1+step;

		local_ref[0][5][0] = -1+step;
		local_ref[0][5][1] = -1;
		local_ref[0][5][2] = -1+step;

		local_ref[0][6][0] = -1+step;
		local_ref[0][6][1] = -1+step;
		local_ref[0][6][2] = -1+step;

		local_ref[0][7][0] = -1;
		local_ref[0][7][1] = -1+step;
		local_ref[0][7][2] =  -1+step;

		//element 1
		local_ref[1][0][0] = -1+step;
		local_ref[1][0][1] = -1;
		local_ref[1][0][2] = -1;

		local_ref[1][1][0] = -1+2*step;
		local_ref[1][1][1] = -1;
		local_ref[1][1][2] = -1;

		local_ref[1][2][0] = -1+2*step;
		local_ref[1][2][1] = -1+2*step;
		local_ref[1][2][2] = -1;

		local_ref[1][3][0] = -1+step;
		local_ref[1][3][1] = -1+step;
		local_ref[1][3][2] = -1;

		local_ref[1][4][0] = -1+step;
		local_ref[1][4][1] = -1;
		local_ref[1][4][2] = -1+step;

		local_ref[1][5][0] = -1+2*step;
		local_ref[1][5][1] = -1;
		local_ref[1][5][2] = -1+2*step;

		local_ref[1][6][0] = -1+2*step;
		local_ref[1][6][1] = -1+2*step;
		local_ref[1][6][2] = -1+2*step;;

		local_ref[1][7][0] = -1+step;
		local_ref[1][7][1] = -1+step;
		local_ref[1][7][2] = -1+step;


		//element 2
		local_ref[2][0][0] = -1+2*step;
		local_ref[2][0][1] = -1;
		local_ref[2][0][2] = -1;

		local_ref[2][1][0] = 1;
		local_ref[2][1][1] = -1;
		local_ref[2][1][2] = -1;

		local_ref[2][2][0] = 1;
		local_ref[2][2][1] = 1;
		local_ref[2][2][2] = -1;

		local_ref[2][3][0] = -1+2*step;
		local_ref[2][3][1] = -1+2*step;
		local_ref[2][3][2] = -1;

		local_ref[2][4][0] = -1+2*step;
		local_ref[2][4][1] = -1;
		local_ref[2][4][2] = -1+2*step;

		local_ref[2][5][0] =  1;
		local_ref[2][5][1] =  -1;
		local_ref[2][5][2] =  1;

		local_ref[2][6][0] = 1;
		local_ref[2][6][1] = 1;
		local_ref[2][6][2] = 1;

		local_ref[2][7][0] = -1+2*step;
		local_ref[2][7][1] = -1+2*step;
		local_ref[2][7][2] = -1+2*step;


		//element 3
		local_ref[3][0][0] = -1;
		local_ref[3][0][1] = -1+step;
		local_ref[3][0][2] = -1;

		local_ref[3][1][0] = -1+step;
		local_ref[3][1][1] = -1+step;
		local_ref[3][1][2] = -1;

		local_ref[3][2][0] = -1+2*step;
		local_ref[3][2][1] = -1+2*step;
		local_ref[3][2][2] = -1;

		local_ref[3][3][0] = -1;
		local_ref[3][3][1] = -1+2*step;
		local_ref[3][3][2] = -1;

		local_ref[3][4][0] = -1;
		local_ref[3][4][1] = -1+step;
		local_ref[3][4][2] = -1+step;

		local_ref[3][5][0] = -1+step;
		local_ref[3][5][1] = -1+step;
		local_ref[3][5][2] = -1+step;

		local_ref[3][6][0] = -1+2*step;
		local_ref[3][6][1] = -1+2*step;
		local_ref[3][6][2] = -1+2*step;

		local_ref[3][7][0] = -1;
		local_ref[3][7][1] = -1+2*step;
		local_ref[3][7][2] = -1+2*step;

		//element 4
		local_ref[4][0][0] = -1;
		local_ref[4][0][1] = -1;
		local_ref[4][0][2] = -1+1*step;

		local_ref[4][1][0] = -1+1*step;
		local_ref[4][1][1] = -1;
		local_ref[4][1][2] = -1+1*step;

		local_ref[4][2][0] = -1+1*step;
		local_ref[4][2][1] = -1+1*step;
		local_ref[4][2][2] = -1+1*step;

		local_ref[4][3][0] = -1;
		local_ref[4][3][1] = -1+1*step;
		local_ref[4][3][2] = -1+1*step;

		local_ref[4][4][0] = -1;
		local_ref[4][4][1] = -1;
		local_ref[4][4][2] = -1+2*step;

		local_ref[4][5][0] = -1+2*step;
		local_ref[4][5][1] = -1+0*step;
		local_ref[4][5][2] = -1+2*step;

		local_ref[4][6][0] = -1+2*step;
		local_ref[4][6][1] = -1+2*step;
		local_ref[4][6][2] = -1+2*step;

		local_ref[4][7][0] = -1+0*step;
		local_ref[4][7][1] = -1+2*step;
		local_ref[4][7][2] = -1+2*step;

		//element 5
		local_ref[5][0][0] = -1;
		local_ref[5][0][1] = -1;
		local_ref[5][0][2] = -1+2*step;

		local_ref[5][1][0] =  -1+2*step;
		local_ref[5][1][1] = -1;
		local_ref[5][1][2] = -1+2*step;

		local_ref[5][2][0] = -1+2*step;
		local_ref[5][2][1] = -1+2*step;
		local_ref[5][2][2] = -1+2*step;

		local_ref[5][3][0] = -1;
		local_ref[5][3][1] = -1+2*step;
		local_ref[5][3][2] = -1+2*step;

		local_ref[5][4][0] = -1;
		local_ref[5][4][1] = -1;
		local_ref[5][4][2] =  1;

		local_ref[5][5][0] = 1;
		local_ref[5][5][1] = -1;
		local_ref[5][5][2] = 1;

		local_ref[5][6][0] = 1;
		local_ref[5][6][1] = 1;
		local_ref[5][6][2] = 1;

		local_ref[5][7][0] = -1;
		local_ref[5][7][1] = 1;
		local_ref[5][7][2] = 1;

		//element 6
		local_ref[6][0][0] = -1;
		local_ref[6][0][1] = -1+2*step;
		local_ref[6][0][2] = -1;

		local_ref[6][1][0] = -1+2*step;
		local_ref[6][1][1] = -1+2*step;
		local_ref[6][1][2] = -1;

		local_ref[6][2][0] = 1;
		local_ref[6][2][1] = 1;
		local_ref[6][2][2] = -1;

		local_ref[6][3][0] = -1;
		local_ref[6][3][1] = 1;
		local_ref[6][3][2] = -1;

		local_ref[6][4][0] = -1;
		local_ref[6][4][1] = -1+2*step;
		local_ref[6][4][2] =  -1+2*step;;

		local_ref[6][5][0] = -1+2*step;;
		local_ref[6][5][1] = -1+2*step;
		local_ref[6][5][2] = -1+2*step;;

		local_ref[6][6][0] = 1;
		local_ref[6][6][1] = 1;
		local_ref[6][6][2] = 1;

		local_ref[6][7][0] = -1;
		local_ref[6][7][1] = 1;
		local_ref[6][7][2] = 1;
	}

	//define the rotation of the reference element
	int rot[3];
	int sym[3];
	int id_node[8];
	vector<int> ord;

	if(elem->pad==115){
		//edge 0 4 3
		rot[0] = 0;
		rot[1] = 0;
		rot[2] = 0;

		sym[0] = 0;
		sym[1] = 0;
		sym[2] = 0;

		ord = RotateHex(rot,sym);

		for(int node_id=0;node_id<8;node_id++){
			id_node[node_id] = elem->nodes[ord[node_id]].id;
		}

	}else if(elem->pad==116){
		//edge 0 1 5
		rot[0] = 0;
		rot[1] = 0;
		rot[2] = -1;

		sym[0] = 0;
		sym[1] = 0;
		sym[2] = 0;

		ord = RotateHex(rot,sym);

		for(int node_id=0;node_id<8;node_id++){
			id_node[node_id] = elem->nodes[ord[node_id]].id;
		}

	}else if(elem->pad==117){
		//edge 1 2 6
		rot[0] = 0;
		rot[1] = 0;
		rot[2] = -1;

		sym[0] = 1;
		sym[1] = 0;
		sym[2] = 0;

		ord = RotateHex(rot,sym);

		for(int node_id=0;node_id<8;node_id++){
			id_node[node_id] = elem->nodes[ord[node_id]].id;
		}

	}else if(elem->pad==118){
		//edge 2 3 7
		rot[0] = 0;
		rot[1] = 0;
		rot[2] = 1;

		sym[0] = 0;
		sym[1] = 0;
		sym[2] = 0;

		ord = RotateHex(rot,sym);

		for(int node_id=0;node_id<8;node_id++){
			id_node[node_id] = elem->nodes[ord[node_id]].id;
		}

	}else if(elem->pad==119){
		//edge 4 8 11
		rot[0] = 0;
		rot[1] = 0;
		rot[2] = 0;

		sym[0] = 0;
		sym[1] = 0;
		sym[2] = 1;

		ord = RotateHex(rot,sym);

		for(int node_id=0;node_id<8;node_id++){
			id_node[node_id] = elem->nodes[ord[node_id]].id;
		}

	}else if(elem->pad==120){
		//edge 8 9 5
		rot[0] = 0;
		rot[1] = 0;
		rot[2] = -1;

		sym[0] = 0;
		sym[1] = 0;
		sym[2] = 1;

		ord = RotateHex(rot,sym);

		for(int node_id=0;node_id<8;node_id++){
			id_node[node_id] = elem->nodes[ord[node_id]].id;
		}

	}else if(elem->pad==121){
		//edge 9 10 6
		rot[0] = 0;
		rot[1] = 0;
		rot[2] = -1;

		sym[0] = 1;
		sym[1] = 0;
		sym[2] = 1;

		ord = RotateHex(rot,sym);

		for(int node_id=0;node_id<8;node_id++){
			id_node[node_id] = elem->nodes[ord[node_id]].id;
		}

	}else if(elem->pad==122){
		//edge 11 10 7
		rot[0] = 0;
		rot[1] = 0;
		rot[2] = 1;

		sym[0] = 0;
		sym[1] = 0;
		sym[2] = 1;

		ord = RotateHex(rot,sym);

		for(int node_id=0;node_id<8;node_id++){
			id_node[node_id] = elem->nodes[ord[node_id]].id;
		}

	}else{
		printf ("Error in template 7\n");
		exit (EXIT_FAILURE);
	}


	for(int i=0;i<7;i++){
		int conn_p[8];
		GtsPoint* point[8]={NULL};

		double cord_in_x[8],cord_in_y[8],cord_in_z[8];
		//add the nodes in the coord vector
		for (int ii = 0; ii < 8; ii++){
			cord_in_x[ii]=coords[3*id_node[ii]] ;
			cord_in_y[ii]=coords[3*id_node[ii]+1] ;
			cord_in_z[ii]=coords[3*id_node[ii]+2] ;
			//fprintf(mesh->fdbg,"coord in: %f, %f, %f, in the node: %d\n",cord_in_x[ii],cord_in_y[ii],cord_in_z[ii],elem->nodes[ii].id);
		}

		//add the new nodes in the
		for(int ii=0;ii<8;ii++){
			conn_p[ii] = 0;

			if((local_ref[i][ii][0]==1 || local_ref[i][ii][0]==-1) &&
					(local_ref[i][ii][1]==1 || local_ref[i][ii][1]==-1) &&
					(local_ref[i][ii][2]==1 || local_ref[i][ii][2]==-1)){

				if(local_ref[i][ii][0]==-1 && local_ref[i][ii][1]==-1 && local_ref[i][ii][2]==-1){
					conn_p[ii] = id_node[0];
				}else if(local_ref[i][ii][0]==1 && local_ref[i][ii][1]==-1 && local_ref[i][ii][2]==-1){
					conn_p[ii] = id_node[1];
				}else if(local_ref[i][ii][0]==1 && local_ref[i][ii][1]==1 && local_ref[i][ii][2]==-1){
					conn_p[ii] = id_node[2];
				}else if(local_ref[i][ii][0]==-1 && local_ref[i][ii][1]==1 && local_ref[i][ii][2]==-1){
					conn_p[ii] = id_node[3];
				}else if(local_ref[i][ii][0]==-1 && local_ref[i][ii][1]==-1 && local_ref[i][ii][2]==1){
					conn_p[ii] = id_node[4];
				}else if(local_ref[i][ii][0]==1 && local_ref[i][ii][1]==-1 && local_ref[i][ii][2]==1){
					conn_p[ii] = id_node[5];
				}else if(local_ref[i][ii][0]==1 && local_ref[i][ii][1]==1 && local_ref[i][ii][2]==1){
					conn_p[ii] = id_node[6];
				}else if(local_ref[i][ii][0]==-1 && local_ref[i][ii][1]==1 && local_ref[i][ii][2]==1){
					conn_p[ii] = id_node[7];
				}
				//fprintf(mesh->fdbg,"coord out: %f, %f, %f, in the node: %d\n",coords[3*conn_p[ii]],coords[3*conn_p[ii]+1],coords[3*conn_p[ii]+2],conn_p[ii]);
			}else{
				cord_in_ref[0] = local_ref[i][ii][0];
				cord_in_ref[1] = local_ref[i][ii][1];
				cord_in_ref[2] = local_ref[i][ii][2];

				point[ii] = LinearMapHex(cord_in_ref, cord_in_x,cord_in_y,cord_in_z);
				double var[3];
				var[0] = point[ii]->x;
				var[1] = point[ii]->y;
				var[2] = point[ii]->z;
				conn_p[ii] = AddPoint( mesh, hash_nodes, point[ii] , coords);
				//fprintf(mesh->fdbg,"coord out: %f, %f, %f, in the node: %d\n",var[0],var[1],var[2],conn_p[ii]);

			}
		}

		int aux[8];
		for(int ino=0;ino<8;ino++){
			aux[ino] = conn_p[ino];
		}
		for(int ino=0;ino<8;ino++){
			conn_p[ord[ino]] = aux[ino];
		}

		//add the new elements and make the connectivity
		if(i==0){
			octant_t *elem1 = (octant_t*) sc_array_index(&mesh->elements, elements_ids[iel]);

			elem1->nodes[0].id = conn_p[0];
			elem1->nodes[1].id = conn_p[1];
			elem1->nodes[2].id = conn_p[2];
			elem1->nodes[3].id = conn_p[3];

			elem1->nodes[4].id = conn_p[4];
			elem1->nodes[5].id = conn_p[5];
			elem1->nodes[6].id = conn_p[6];
			elem1->nodes[7].id = conn_p[7];

			CopyPropEl(mesh,id,elem1);

		}else{
			octant_t* elem2 = (octant_t*) sc_array_push(&mesh->elements);

			elem2->nodes[0].id = conn_p[0];
			elem2->nodes[1].id = conn_p[1];
			elem2->nodes[2].id = conn_p[2];
			elem2->nodes[3].id = conn_p[3];

			elem2->nodes[4].id = conn_p[4];
			elem2->nodes[5].id = conn_p[5];
			elem2->nodes[6].id = conn_p[6];
			elem2->nodes[7].id = conn_p[7];

			CopyPropEl(mesh,id,elem2);

		}
	}
	sc_array_reset(&toto);
}

void ApplyTemplate8(hexa_tree_t* mesh, std::vector<double>& coords, std::vector<int>& elements_ids, int iel, sc_hash_array_t* hash_nodes, double step){

	sc_array_t toto;
	sc_array_init(&toto, sizeof(octant_t));

	octant_t* elemOrig = (octant_t*) sc_array_index(&mesh->elements, elements_ids[iel]);
	octant_t * elem = (octant_t*) sc_array_push(&toto);
	hexa_element_copy(elemOrig,elem);

	int id = elements_ids[iel];

	double cord_in_ref[3];
	cord_in_ref[0] = 0;
	cord_in_ref[1] = 0;
	cord_in_ref[2] = 0;

	//reference element edge 0
	double local_ref[9][8][3];
	if(true){
		//element 0
		local_ref[0][0][0] = -1;
		local_ref[0][0][1] = -1;
		local_ref[0][0][2] = -1;

		local_ref[0][1][0] = -1+step;
		local_ref[0][1][1] = -1;
		local_ref[0][1][2] = -1;

		local_ref[0][2][0] = -1+step;
		local_ref[0][2][1] = -1+step;
		local_ref[0][2][2] = -1;

		local_ref[0][3][0] = -1;
		local_ref[0][3][1] = -1+step;
		local_ref[0][3][2] = -1;

		local_ref[0][4][0] = -1;
		local_ref[0][4][1] = -1;
		local_ref[0][4][2] = 1;

		local_ref[0][5][0] = -1+step;
		local_ref[0][5][1] = -1;
		local_ref[0][5][2] = 1;

		local_ref[0][6][0] = -1+step;
		local_ref[0][6][1] = -1+step;
		local_ref[0][6][2] = 1;

		local_ref[0][7][0] = -1;
		local_ref[0][7][1] = -1+step;
		local_ref[0][7][2] = 1;

		//element 1
		local_ref[1][0][0] = -1+step;
		local_ref[1][0][1] = -1;
		local_ref[1][0][2] = -1;

		local_ref[1][1][0] = -1+2*step;
		local_ref[1][1][1] = -1;
		local_ref[1][1][2] = -1;

		local_ref[1][2][0] = -1+2*step;
		local_ref[1][2][1] = -1+step;
		local_ref[1][2][2] = -1;

		local_ref[1][3][0] = -1+step;
		local_ref[1][3][1] = -1+step;
		local_ref[1][3][2] = -1;

		local_ref[1][4][0] = -1+step;
		local_ref[1][4][1] = -1;
		local_ref[1][4][2] = 1;

		local_ref[1][5][0] = -1+2*step;
		local_ref[1][5][1] = -1;
		local_ref[1][5][2] = 1;

		local_ref[1][6][0] = -1+2*step;
		local_ref[1][6][1] = -1+step;
		local_ref[1][6][2] = 1;

		local_ref[1][7][0] = -1+step;
		local_ref[1][7][1] = -1+step;
		local_ref[1][7][2] = 1;

		//element 2
		local_ref[2][0][0] = -1+2*step;
		local_ref[2][0][1] = -1;
		local_ref[2][0][2] = -1;

		local_ref[2][1][0] = 1;
		local_ref[2][1][1] = -1;
		local_ref[2][1][2] = -1;

		local_ref[2][2][0] = 1;
		local_ref[2][2][1] = -1+step;
		local_ref[2][2][2] = -1;

		local_ref[2][3][0] = -1+2*step;
		local_ref[2][3][1] = -1+step;
		local_ref[2][3][2] = -1;

		local_ref[2][4][0] = -1+2*step;
		local_ref[2][4][1] = -1;
		local_ref[2][4][2] = 1;

		local_ref[2][5][0] = 1;
		local_ref[2][5][1] = -1;
		local_ref[2][5][2] = 1;

		local_ref[2][6][0] = 1;
		local_ref[2][6][1] = -1+step;
		local_ref[2][6][2] = 1;

		local_ref[2][7][0] = -1+2*step;
		local_ref[2][7][1] = -1+step;
		local_ref[2][7][2] = 1;

		//element 3
		local_ref[3][0][0] = -1;
		local_ref[3][0][1] = -1+step;
		local_ref[3][0][2] = -1;

		local_ref[3][1][0] = -1+step;
		local_ref[3][1][1] = -1+step;
		local_ref[3][1][2] = -1;

		local_ref[3][2][0] = -1+step;
		local_ref[3][2][1] = -1+2*step;
		local_ref[3][2][2] = -1;

		local_ref[3][3][0] = -1;
		local_ref[3][3][1] = -1+2*step;
		local_ref[3][3][2] = -1;

		local_ref[3][4][0] = -1;
		local_ref[3][4][1] = -1+step;
		local_ref[3][4][2] = 1;

		local_ref[3][5][0] = -1+step;
		local_ref[3][5][1] = -1+step;
		local_ref[3][5][2] = 1;

		local_ref[3][6][0] = -1+step;
		local_ref[3][6][1] = -1+2*step;
		local_ref[3][6][2] = 1;

		local_ref[3][7][0] = -1;
		local_ref[3][7][1] = -1+2*step;
		local_ref[3][7][2] = 1;

		//element 4
		local_ref[4][0][0] = -1+step;
		local_ref[4][0][1] = -1+step;
		local_ref[4][0][2] = -1;

		local_ref[4][1][0] = -1+2*step;
		local_ref[4][1][1] = -1+step;
		local_ref[4][1][2] = -1;

		local_ref[4][2][0] = -1+2*step;
		local_ref[4][2][1] = -1+2*step;
		local_ref[4][2][2] = -1;

		local_ref[4][3][0] = -1+step;
		local_ref[4][3][1] = -1+2*step;
		local_ref[4][3][2] = -1;

		local_ref[4][4][0] = -1+step;
		local_ref[4][4][1] = -1+step;
		local_ref[4][4][2] = 1;

		local_ref[4][5][0] = -1+2*step;
		local_ref[4][5][1] = -1+step;
		local_ref[4][5][2] = 1;

		local_ref[4][6][0] = -1+2*step;
		local_ref[4][6][1] = -1+2*step;
		local_ref[4][6][2] = 1;

		local_ref[4][7][0] = -1+step;
		local_ref[4][7][1] = -1+2*step;
		local_ref[4][7][2] = 1;

		//element 5
		local_ref[5][0][0] = -1+2*step;
		local_ref[5][0][1] = -1+step;
		local_ref[5][0][2] = -1;

		local_ref[5][1][0] = 1;
		local_ref[5][1][1] = -1+step;
		local_ref[5][1][2] = -1;

		local_ref[5][2][0] = 1;
		local_ref[5][2][1] = -1+2*step;
		local_ref[5][2][2] = -1;

		local_ref[5][3][0] = -1+2*step;
		local_ref[5][3][1] = -1+2*step;
		local_ref[5][3][2] = -1;

		local_ref[5][4][0] = -1+2*step;
		local_ref[5][4][1] = -1+step;
		local_ref[5][4][2] = 1;

		local_ref[5][5][0] = 1;
		local_ref[5][5][1] = -1+step;
		local_ref[5][5][2] = 1;

		local_ref[5][6][0] = 1;
		local_ref[5][6][1] = -1+2*step;
		local_ref[5][6][2] = 1;

		local_ref[5][7][0] = -1+2*step;
		local_ref[5][7][1] = -1+2*step;
		local_ref[5][7][2] = 1;

		//element 6
		local_ref[6][0][0] = -1;
		local_ref[6][0][1] = -1+2*step;
		local_ref[6][0][2] = -1;

		local_ref[6][1][0] = -1+step;
		local_ref[6][1][1] = -1+2*step;
		local_ref[6][1][2] = -1;

		local_ref[6][2][0] = -1+step;
		local_ref[6][2][1] = 1;
		local_ref[6][2][2] = -1;

		local_ref[6][3][0] = -1;
		local_ref[6][3][1] = 1;
		local_ref[6][3][2] = -1;

		local_ref[6][4][0] = -1;
		local_ref[6][4][1] = -1+2*step;
		local_ref[6][4][2] = 1;

		local_ref[6][5][0] = -1+step;
		local_ref[6][5][1] = -1+2*step;
		local_ref[6][5][2] = 1;

		local_ref[6][6][0] = -1+step;
		local_ref[6][6][1] = 1;
		local_ref[6][6][2] = 1;

		local_ref[6][7][0] = -1;
		local_ref[6][7][1] = 1;
		local_ref[6][7][2] = 1;

		//element 7
		local_ref[7][0][0] = -1+step;
		local_ref[7][0][1] = -1+2*step;
		local_ref[7][0][2] = -1;

		local_ref[7][1][0] = -1+2*step;
		local_ref[7][1][1] = -1+2*step;
		local_ref[7][1][2] = -1;

		local_ref[7][2][0] = -1+2*step;
		local_ref[7][2][1] = 1;
		local_ref[7][2][2] = -1;

		local_ref[7][3][0] = -1+step;
		local_ref[7][3][1] = 1;
		local_ref[7][3][2] = -1;

		local_ref[7][4][0] = -1+step;
		local_ref[7][4][1] = -1+2*step;
		local_ref[7][4][2] = 1;

		local_ref[7][5][0] = -1+2*step;
		local_ref[7][5][1] = -1+2*step;
		local_ref[7][5][2] = 1;

		local_ref[7][6][0] = -1+2*step;
		local_ref[7][6][1] = 1;
		local_ref[7][6][2] = 1;

		local_ref[7][7][0] = -1+step;
		local_ref[7][7][1] = 1;
		local_ref[7][7][2] = 1;

		//element 8
		local_ref[8][0][0] = -1+2*step;
		local_ref[8][0][1] = -1+2*step;
		local_ref[8][0][2] = -1;

		local_ref[8][1][0] = 1;
		local_ref[8][1][1] = -1+2*step;
		local_ref[8][1][2] = -1;

		local_ref[8][2][0] = 1;
		local_ref[8][2][1] = 1;
		local_ref[8][2][2] = -1;

		local_ref[8][3][0] = -1+2*step;
		local_ref[8][3][1] = 1;
		local_ref[8][3][2] = -1;

		local_ref[8][4][0] = -1+2*step;
		local_ref[8][4][1] = -1+2*step;
		local_ref[8][4][2] = 1;

		local_ref[8][5][0] = 1;
		local_ref[8][5][1] = -1+2*step;
		local_ref[8][5][2] = 1;

		local_ref[8][6][0] = 1;
		local_ref[8][6][1] = 1;
		local_ref[8][6][2] = 1;

		local_ref[8][7][0] = -1+2*step;
		local_ref[8][7][1] = 1;
		local_ref[8][7][2] = 1;
	}

	//define the rotation of the reference element
	int rot[3];
	int sym[3];
	int id_node[8];
	vector<int> ord;

	if(elem->pad==123){
		//edge 0 1 2 3 8 9 10 11
		rot[0] = 0;
		rot[1] = 0;
		rot[2] = 0;

		sym[0] = 0;
		sym[1] = 0;
		sym[2] = 0;

		ord = RotateHex(rot,sym);

		for(int node_id=0;node_id<8;node_id++){
			id_node[node_id] = elem->nodes[ord[node_id]].id;
		}

	}else if(elem->pad==124){
		//edge 0 5 4 8 2 6 10 7
		rot[0] = 1;
		rot[1] = 0;
		rot[2] = 0;

		sym[0] = 0;
		sym[1] = 0;
		sym[2] = 0;

		ord = RotateHex(rot,sym);

		for(int node_id=0;node_id<8;node_id++){
			id_node[node_id] = elem->nodes[ord[node_id]].id;
		}

	}else if(elem->pad==125){
		//edge 1 6 8 5 3 7 11 4
		rot[0] = 0;
		rot[1] = 1;
		rot[2] = 0;

		sym[0] = 0;
		sym[1] = 0;
		sym[2] = 0;

		ord = RotateHex(rot,sym);

		for(int node_id=0;node_id<8;node_id++){
			id_node[node_id] = elem->nodes[ord[node_id]].id;
		}

	}else{
		printf ("Error in template 8\n");
		exit (EXIT_FAILURE);
	}


	for(int i=0;i<9;i++){
		int conn_p[8];
		GtsPoint* point[8]={NULL};

		double cord_in_x[8],cord_in_y[8],cord_in_z[8];
		//add the nodes in the coord vector
		for (int ii = 0; ii < 8; ii++){
			cord_in_x[ii]=coords[3*id_node[ii]] ;
			cord_in_y[ii]=coords[3*id_node[ii]+1] ;
			cord_in_z[ii]=coords[3*id_node[ii]+2] ;
			//fprintf(mesh->fdbg,"coord in: %f, %f, %f, in the node: %d\n",cord_in_x[ii],cord_in_y[ii],cord_in_z[ii],elem->nodes[ii].id);
		}

		//add the new nodes in the
		for(int ii=0;ii<8;ii++){
			conn_p[ii] = 0;

			if((local_ref[i][ii][0]==1 || local_ref[i][ii][0]==-1) &&
					(local_ref[i][ii][1]==1 || local_ref[i][ii][1]==-1) &&
					(local_ref[i][ii][2]==1 || local_ref[i][ii][2]==-1)){

				if(local_ref[i][ii][0]==-1 && local_ref[i][ii][1]==-1 && local_ref[i][ii][2]==-1){
					conn_p[ii] = id_node[0];
				}else if(local_ref[i][ii][0]==1 && local_ref[i][ii][1]==-1 && local_ref[i][ii][2]==-1){
					conn_p[ii] = id_node[1];
				}else if(local_ref[i][ii][0]==1 && local_ref[i][ii][1]==1 && local_ref[i][ii][2]==-1){
					conn_p[ii] = id_node[2];
				}else if(local_ref[i][ii][0]==-1 && local_ref[i][ii][1]==1 && local_ref[i][ii][2]==-1){
					conn_p[ii] = id_node[3];
				}else if(local_ref[i][ii][0]==-1 && local_ref[i][ii][1]==-1 && local_ref[i][ii][2]==1){
					conn_p[ii] = id_node[4];
				}else if(local_ref[i][ii][0]==1 && local_ref[i][ii][1]==-1 && local_ref[i][ii][2]==1){
					conn_p[ii] = id_node[5];
				}else if(local_ref[i][ii][0]==1 && local_ref[i][ii][1]==1 && local_ref[i][ii][2]==1){
					conn_p[ii] = id_node[6];
				}else if(local_ref[i][ii][0]==-1 && local_ref[i][ii][1]==1 && local_ref[i][ii][2]==1){
					conn_p[ii] = id_node[7];
				}
				//fprintf(mesh->fdbg,"coord out: %f, %f, %f, in the node: %d\n",coords[3*conn_p[ii]],coords[3*conn_p[ii]+1],coords[3*conn_p[ii]+2],conn_p[ii]);
			}else{
				cord_in_ref[0] = local_ref[i][ii][0];
				cord_in_ref[1] = local_ref[i][ii][1];
				cord_in_ref[2] = local_ref[i][ii][2];

				point[ii] = LinearMapHex(cord_in_ref, cord_in_x,cord_in_y,cord_in_z);
				double var[3];
				var[0] = point[ii]->x;
				var[1] = point[ii]->y;
				var[2] = point[ii]->z;
				conn_p[ii] = AddPoint( mesh, hash_nodes, point[ii] , coords);
				//fprintf(mesh->fdbg,"coord out: %f, %f, %f, in the node: %d\n",var[0],var[1],var[2],conn_p[ii]);

			}
		}

		int aux[8];
		for(int ino=0;ino<8;ino++){
			aux[ino] = conn_p[ino];
		}
		for(int ino=0;ino<8;ino++){
			conn_p[ord[ino]] = aux[ino];
		}

		//add the new elements and make the connectivity
		if(i==0){
			octant_t *elem1 = (octant_t*) sc_array_index(&mesh->elements, elements_ids[iel]);

			elem1->nodes[0].id = conn_p[0];
			elem1->nodes[1].id = conn_p[1];
			elem1->nodes[2].id = conn_p[2];
			elem1->nodes[3].id = conn_p[3];

			elem1->nodes[4].id = conn_p[4];
			elem1->nodes[5].id = conn_p[5];
			elem1->nodes[6].id = conn_p[6];
			elem1->nodes[7].id = conn_p[7];

			CopyPropEl(mesh,id,elem1);

		}else{
			octant_t* elem2 = (octant_t*) sc_array_push(&mesh->elements);

			elem2->nodes[0].id = conn_p[0];
			elem2->nodes[1].id = conn_p[1];
			elem2->nodes[2].id = conn_p[2];
			elem2->nodes[3].id = conn_p[3];

			elem2->nodes[4].id = conn_p[4];
			elem2->nodes[5].id = conn_p[5];
			elem2->nodes[6].id = conn_p[6];
			elem2->nodes[7].id = conn_p[7];

			CopyPropEl(mesh,id,elem2);

		}
	}
	sc_array_reset(&toto);
}

void ApplyTemplate9(hexa_tree_t* mesh, std::vector<double>& coords, std::vector<int>& elements_ids, int iel, sc_hash_array_t* hash_nodes, double step){

	sc_array_t toto;
	sc_array_init(&toto, sizeof(octant_t));

	octant_t* elemOrig = (octant_t*) sc_array_index(&mesh->elements, elements_ids[iel]);
	octant_t * elem = (octant_t*) sc_array_push(&toto);
	hexa_element_copy(elemOrig,elem);

	int id = elements_ids[iel];

	double cord_in_ref[3];
	cord_in_ref[0] = 0;
	cord_in_ref[1] = 0;
	cord_in_ref[2] = 0;

	//reference element edge 0
	double local_ref[31][8][3];
	if(true){
		//element 0
		local_ref[0][0][0] = -1;
		local_ref[0][0][1] = -1;
		local_ref[0][0][2] = -1;

		local_ref[0][1][0] = 1;
		local_ref[0][1][1] = -1;
		local_ref[0][1][2] = -1;

		local_ref[0][2][0] = 1;
		local_ref[0][2][1] = 1;
		local_ref[0][2][2] = -1;

		local_ref[0][3][0] = -1;
		local_ref[0][3][1] = 1;
		local_ref[0][3][2] = -1;

		local_ref[0][4][0] = -1;
		local_ref[0][4][1] = -1+step;
		local_ref[0][4][2] = -1+1*3*step/double(4);

		local_ref[0][5][0] = 1;
		local_ref[0][5][1] = -1+step;
		local_ref[0][5][2] = -1+1*3*step/double(4);

		local_ref[0][6][0] = 1;
		local_ref[0][6][1] = -1+2*step;
		local_ref[0][6][2] = -1+1*3*step/double(4);

		local_ref[0][7][0] = -1;
		local_ref[0][7][1] = -1+2*step;
		local_ref[0][7][2] = -1+1*3*step/double(4);

		//element 1
		local_ref[1][0][0] = -1;
		local_ref[1][0][1] = -1+step;
		local_ref[1][0][2] = -1+1*3*step/double(4);

		local_ref[1][1][0] = 1;
		local_ref[1][1][1] = -1+step;
		local_ref[1][1][2] = -1+1*3*step/double(4);

		local_ref[1][2][0] = 1;
		local_ref[1][2][1] = -1+2*step;
		local_ref[1][2][2] = -1+1*3*step/double(4);

		local_ref[1][3][0] = -1;
		local_ref[1][3][1] = -1+2*step;
		local_ref[1][3][2] = -1+1*3*step/double(4);

		local_ref[1][4][0] = -1+step;
		local_ref[1][4][1] = -1+step;
		local_ref[1][4][2] = -1+3*3*step/double(8);//-1+2*3*step/double(4);

		local_ref[1][5][0] = -1+2*step;
		local_ref[1][5][1] = -1+step;
		local_ref[1][5][2] = -1+3*3*step/double(8);//-1+2*3*step/double(4);

		local_ref[1][6][0] = -1+2*step;
		local_ref[1][6][1] = -1+2*step;
		local_ref[1][6][2] = -1+3*3*step/double(8);//-1+2*3*step/double(4);

		local_ref[1][7][0] = -1+step;
		local_ref[1][7][1] = -1+2*step;
		local_ref[1][7][2] = -1+3*3*step/double(8);//-1+2*3*step/double(4);


		//element 2
		local_ref[2][0][0] = -1;
		local_ref[2][0][1] = -1;
		local_ref[2][0][2] = -1;

		local_ref[2][1][0] =  1;
		local_ref[2][1][1] = -1;
		local_ref[2][1][2] = -1;

		local_ref[2][2][0] =  1;
		local_ref[2][2][1] = -1+step;
		local_ref[2][2][2] = -1+1*3*step/double(4);

		local_ref[2][3][0] = -1;
		local_ref[2][3][1] = -1+step;
		local_ref[2][3][2] = -1+1*3*step/double(4);

		local_ref[2][4][0] = -1+step;
		local_ref[2][4][1] = -1;
		local_ref[2][4][2] = -1+1*3*step/double(4);

		local_ref[2][5][0] = -1+2*step;
		local_ref[2][5][1] = -1;
		local_ref[2][5][2] = -1+1*3*step/double(4);

		local_ref[2][6][0] = -1+2*step;
		local_ref[2][6][1] = -1+1*step;
		local_ref[2][6][2] = -1+3*3*step/double(8);//-1+2*3*step/double(4);

		local_ref[2][7][0] = -1+step;
		local_ref[2][7][1] = -1+1*step;
		local_ref[2][7][2] = -1+3*3*step/double(8);//-1+2*3*step/double(4);

		//element 3
		local_ref[3][0][0] = -1;
		local_ref[3][0][1] = -1+2*step;
		local_ref[3][0][2] = -1+1*3*step/double(4);

		local_ref[3][1][0] =  1;
		local_ref[3][1][1] = -1+2*step;
		local_ref[3][1][2] = -1+1*3*step/double(4);

		local_ref[3][2][0] =  1;
		local_ref[3][2][1] =  1;
		local_ref[3][2][2] = -1;

		local_ref[3][3][0] = -1;
		local_ref[3][3][1] =  1;
		local_ref[3][3][2] = -1;

		local_ref[3][4][0] = -1+step;
		local_ref[3][4][1] = -1+2*step;
		local_ref[3][4][2] = -1+3*3*step/double(8);//-1+2*3*step/double(4);

		local_ref[3][5][0] = -1+2*step;
		local_ref[3][5][1] = -1+2*step;
		local_ref[3][5][2] = -1+3*3*step/double(8);//-1+2*3*step/double(4);

		local_ref[3][6][0] = -1+2*step;
		local_ref[3][6][1] = 1;
		local_ref[3][6][2] = -1+1*3*step/double(4);

		local_ref[3][7][0] = -1+step;
		local_ref[3][7][1] = 1;
		local_ref[3][7][2] = -1+1*3*step/double(4);

		//element 4
		local_ref[4][0][0] = -1+step;
		local_ref[4][0][1] = -1+step;
		local_ref[4][0][2] = -1+3*3*step/double(8);

		local_ref[4][1][0] = -1+2*step;
		local_ref[4][1][1] = -1+step;
		local_ref[4][1][2] = -1+3*3*step/double(8);

		local_ref[4][2][0] = -1+2*step;
		local_ref[4][2][1] = -1+2*step;
		local_ref[4][2][2] = -1+3*3*step/double(8);

		local_ref[4][3][0] = -1+step;
		local_ref[4][3][1] = -1+2*step;
		local_ref[4][3][2] = -1+3*3*step/double(8);

		local_ref[4][4][0] = -1+step;
		local_ref[4][4][1] = -1+step;
		local_ref[4][4][2] = -1+2*3*step/double(4);

		local_ref[4][5][0] = -1+2*step;
		local_ref[4][5][1] = -1+step;
		local_ref[4][5][2] = -1+2*3*step/double(4);

		local_ref[4][6][0] = -1+2*step;
		local_ref[4][6][1] = -1+2*step;
		local_ref[4][6][2] = -1+2*3*step/double(4);

		local_ref[4][7][0] = -1+step;
		local_ref[4][7][1] = -1+2*step;
		local_ref[4][7][2] = -1+2*3*step/double(4);

		//element 5
		local_ref[5][0][0] = -1+step;
		local_ref[5][0][1] = -1+step;
		local_ref[5][0][2] = -1+2*3*step/double(4);

		local_ref[5][1][0] = -1+2*step;
		local_ref[5][1][1] = -1+step;
		local_ref[5][1][2] = -1+2*3*step/double(4);

		local_ref[5][2][0] = -1+2*step;
		local_ref[5][2][1] = -1+2*step;
		local_ref[5][2][2] = -1+2*3*step/double(4);

		local_ref[5][3][0] = -1+step;
		local_ref[5][3][1] = -1+2*step;
		local_ref[5][3][2] = -1+2*3*step/double(4);

		local_ref[5][4][0] = -1+step;
		local_ref[5][4][1] = -1+step;
		local_ref[5][4][2] = -1+3*3*step/double(4);

		local_ref[5][5][0] = -1+2*step;
		local_ref[5][5][1] = -1+step;
		local_ref[5][5][2] = -1+3*3*step/double(4);

		local_ref[5][6][0] = -1+2*step;
		local_ref[5][6][1] = -1+2*step;
		local_ref[5][6][2] = -1+3*3*step/double(4);

		local_ref[5][7][0] = -1+step;
		local_ref[5][7][1] = -1+2*step;
		local_ref[5][7][2] = -1+3*3*step/double(4);

		//element 6
		local_ref[6][0][0] = -1+step;
		local_ref[6][0][1] = -1+step;
		local_ref[6][0][2] = -1+3*3*step/double(4);

		local_ref[6][1][0] = -1+2*step;
		local_ref[6][1][1] = -1+step;
		local_ref[6][1][2] = -1+3*3*step/double(4);

		local_ref[6][2][0] = -1+2*step;
		local_ref[6][2][1] = -1+2*step;
		local_ref[6][2][2] = -1+3*3*step/double(4);

		local_ref[6][3][0] = -1+step;
		local_ref[6][3][1] = -1+2*step;
		local_ref[6][3][2] = -1+3*3*step/double(4);

		local_ref[6][4][0] = -1+step;
		local_ref[6][4][1] = -1+step;
		local_ref[6][4][2] = 1;

		local_ref[6][5][0] = -1+2*step;
		local_ref[6][5][1] = -1+step;
		local_ref[6][5][2] = 1;

		local_ref[6][6][0] = -1+2*step;
		local_ref[6][6][1] = -1+2*step;
		local_ref[6][6][2] = 1;

		local_ref[6][7][0] = -1+step;
		local_ref[6][7][1] = -1+2*step;
		local_ref[6][7][2] = 1;

		//element 7
		local_ref[7][0][0] = -1;
		local_ref[7][0][1] = -1;
		local_ref[7][0][2] = -1;

		local_ref[7][1][0] = -1+step;
		local_ref[7][1][1] = -1;
		local_ref[7][1][2] = -1+1*3*step/double(4);

		local_ref[7][2][0] = -1+step;
		local_ref[7][2][1] = -1+step;
		local_ref[7][2][2] = -1+3*3*step/double(8);

		local_ref[7][3][0] = -1;
		local_ref[7][3][1] = -1+step;
		local_ref[7][3][2] = -1+1*3*step/double(4);

		local_ref[7][4][0] = -1;
		local_ref[7][4][1] = -1;
		local_ref[7][4][2] = -1+step;//-1+2*3*step/double(4);

		local_ref[7][5][0] = -1+step;
		local_ref[7][5][1] = -1;
		local_ref[7][5][2] = -1+2*3*step/double(4);

		local_ref[7][6][0] = -1+step;
		local_ref[7][6][1] = -1+step;
		local_ref[7][6][2] = -1+2*3*step/double(4);

		local_ref[7][7][0] = -1;
		local_ref[7][7][1] = -1+step;
		local_ref[7][7][2] = -1+2*3*step/double(4);


		//element 8
		local_ref[8][0][0] = -1+step;
		local_ref[8][0][1] = -1;
		local_ref[8][0][2] = -1+1*3*step/double(4);;

		local_ref[8][1][0] = -1+2*step;
		local_ref[8][1][1] = -1;
		local_ref[8][1][2] = -1+1*3*step/double(4);;

		local_ref[8][2][0] = -1+2*step;
		local_ref[8][2][1] = -1+step;
		local_ref[8][2][2] = -1+3*3*step/double(8);

		local_ref[8][3][0] = -1+step;
		local_ref[8][3][1] = -1+step;
		local_ref[8][3][2] = -1+3*3*step/double(8);

		local_ref[8][4][0] = -1+step;
		local_ref[8][4][1] = -1;
		local_ref[8][4][2] = -1+2*3*step/double(4);

		local_ref[8][5][0] = -1+2*step;
		local_ref[8][5][1] = -1;
		local_ref[8][5][2] = -1+2*3*step/double(4);

		local_ref[8][6][0] = -1+2*step;
		local_ref[8][6][1] = -1+1*step;
		local_ref[8][6][2] = -1+2*3*step/double(4);

		local_ref[8][7][0] = -1+step;
		local_ref[8][7][1] = -1+1*step;
		local_ref[8][7][2] = -1+2*3*step/double(4);

		//element 9
		local_ref[9][0][0] = -1+2*step;
		local_ref[9][0][1] = -1;
		local_ref[9][0][2] = -1+1*3*step/double(4);

		local_ref[9][1][0] = 1;
		local_ref[9][1][1] = -1;
		local_ref[9][1][2] = -1;

		local_ref[9][2][0] =  1;
		local_ref[9][2][1] = -1+step;
		local_ref[9][2][2] = -1+1*3*step/double(4);

		local_ref[9][3][0] = -1+2*step;
		local_ref[9][3][1] = -1+1*step;
		local_ref[9][3][2] = -1+3*3*step/double(8);

		local_ref[9][4][0] = -1+2*step;
		local_ref[9][4][1] = -1;
		local_ref[9][4][2] = -1+2*3*step/double(4);

		local_ref[9][5][0] = 1;
		local_ref[9][5][1] = -1;
		local_ref[9][5][2] = -1+step;//-1+2*3*step/double(4);

		local_ref[9][6][0] = 1;
		local_ref[9][6][1] = -1+1*step;
		local_ref[9][6][2] = -1+2*3*step/double(4);

		local_ref[9][7][0] = -1+2*step;
		local_ref[9][7][1] = -1+1*step;
		local_ref[9][7][2] = -1+2*3*step/double(4);

		//element 10
		local_ref[10][0][0] = -1+2*step;
		local_ref[10][0][1] = -1+1*step;
		local_ref[10][0][2] = -1+3*3*step/double(8);

		local_ref[10][1][0] = 1;
		local_ref[10][1][1] = -1+1*step;
		local_ref[10][1][2] = -1+1*3*step/double(4);

		local_ref[10][2][0] =  1;
		local_ref[10][2][1] = -1+2*step;
		local_ref[10][2][2] = -1+1*3*step/double(4);

		local_ref[10][3][0] = -1+2*step;
		local_ref[10][3][1] = -1+2*step;
		local_ref[10][3][2] = -1+3*3*step/double(8);

		local_ref[10][4][0] = -1+2*step;
		local_ref[10][4][1] = -1+1*step;
		local_ref[10][4][2] = -1+2*3*step/double(4);

		local_ref[10][5][0] = 1;
		local_ref[10][5][1] = -1+1*step;
		local_ref[10][5][2] = -1+2*3*step/double(4);

		local_ref[10][6][0] = 1;
		local_ref[10][6][1] = -1+2*step;
		local_ref[10][6][2] = -1+2*3*step/double(4);

		local_ref[10][7][0] = -1+2*step;
		local_ref[10][7][1] = -1+2*step;
		local_ref[10][7][2] = -1+2*3*step/double(4);

		//element 11
		local_ref[11][0][0] = -1+2*step;
		local_ref[11][0][1] = -1+2*step;
		local_ref[11][0][2] = -1+3*3*step/double(8);

		local_ref[11][1][0] = 1;
		local_ref[11][1][1] = -1+2*step;
		local_ref[11][1][2] = -1+1*3*step/double(4);

		local_ref[11][2][0] =  1;
		local_ref[11][2][1] =  1;
		local_ref[11][2][2] = -1;

		local_ref[11][3][0] = -1+2*step;
		local_ref[11][3][1] =  1;
		local_ref[11][3][2] = -1+1*3*step/double(4);

		local_ref[11][4][0] = -1+2*step;
		local_ref[11][4][1] = -1+2*step;
		local_ref[11][4][2] = -1+2*3*step/double(4);

		local_ref[11][5][0] = 1;
		local_ref[11][5][1] = -1+2*step;
		local_ref[11][5][2] = -1+2*3*step/double(4);

		local_ref[11][6][0] = 1;
		local_ref[11][6][1] = 1;
		local_ref[11][6][2] = -1+step;//-1+2*3*step/double(4);

		local_ref[11][7][0] = -1+2*step;
		local_ref[11][7][1] = 1;
		local_ref[11][7][2] = -1+2*3*step/double(4);

		//element 12
		local_ref[12][0][0] = -1+1*step;
		local_ref[12][0][1] = -1+2*step;
		local_ref[12][0][2] = -1+3*3*step/double(8);

		local_ref[12][1][0] = -1+2*step;
		local_ref[12][1][1] = -1+2*step;
		local_ref[12][1][2] = -1+3*3*step/double(8);

		local_ref[12][2][0] =  -1+2*step;
		local_ref[12][2][1] =  1;
		local_ref[12][2][2] = -1+1*3*step/double(4);

		local_ref[12][3][0] = -1+1*step;
		local_ref[12][3][1] =  1;
		local_ref[12][3][2] = -1+1*3*step/double(4);

		local_ref[12][4][0] = -1+1*step;
		local_ref[12][4][1] = -1+2*step;
		local_ref[12][4][2] = -1+2*3*step/double(4);

		local_ref[12][5][0] = -1+2*step;
		local_ref[12][5][1] = -1+2*step;
		local_ref[12][5][2] = -1+2*3*step/double(4);

		local_ref[12][6][0] = -1+2*step;
		local_ref[12][6][1] = 1;
		local_ref[12][6][2] = -1+2*3*step/double(4);

		local_ref[12][7][0] = -1+1*step;
		local_ref[12][7][1] = 1;
		local_ref[12][7][2] = -1+2*3*step/double(4);

		//element 13
		local_ref[13][0][0] = -1;
		local_ref[13][0][1] = -1+2*step;
		local_ref[13][0][2] = -1+1*3*step/double(4);

		local_ref[13][1][0] = -1+1*step;
		local_ref[13][1][1] = -1+2*step;
		local_ref[13][1][2] = -1+3*3*step/double(8);

		local_ref[13][2][0] =  -1+1*step;
		local_ref[13][2][1] =  1;
		local_ref[13][2][2] = -1+1*3*step/double(4);

		local_ref[13][3][0] = -1;
		local_ref[13][3][1] =  1;
		local_ref[13][3][2] = -1;

		local_ref[13][4][0] = -1;
		local_ref[13][4][1] = -1+2*step;
		local_ref[13][4][2] = -1+2*3*step/double(4);

		local_ref[13][5][0] = -1+1*step;
		local_ref[13][5][1] = -1+2*step;
		local_ref[13][5][2] = -1+2*3*step/double(4);

		local_ref[13][6][0] = -1+1*step;
		local_ref[13][6][1] = 1;
		local_ref[13][6][2] = -1+2*3*step/double(4);

		local_ref[13][7][0] = -1;
		local_ref[13][7][1] = 1;
		local_ref[13][7][2] = -1+step;//-1+2*3*step/double(4);

		//element 14
		local_ref[14][0][0] = -1;
		local_ref[14][0][1] = -1+1*step;
		local_ref[14][0][2] = -1+1*3*step/double(4);

		local_ref[14][1][0] = -1+1*step;
		local_ref[14][1][1] = -1+1*step;
		local_ref[14][1][2] = -1+3*3*step/double(8);

		local_ref[14][2][0] =  -1+1*step;
		local_ref[14][2][1] =  -1+2*step;
		local_ref[14][2][2] = -1+3*3*step/double(8);

		local_ref[14][3][0] = -1;
		local_ref[14][3][1] =  -1+2*step;
		local_ref[14][3][2] = -1+1*3*step/double(4);

		local_ref[14][4][0] = -1;
		local_ref[14][4][1] = -1+1*step;
		local_ref[14][4][2] = -1+2*3*step/double(4);

		local_ref[14][5][0] = -1+1*step;
		local_ref[14][5][1] = -1+1*step;
		local_ref[14][5][2] = -1+2*3*step/double(4);

		local_ref[14][6][0] = -1+1*step;
		local_ref[14][6][1] = -1+2*step;
		local_ref[14][6][2] = -1+2*3*step/double(4);

		local_ref[14][7][0] = -1;
		local_ref[14][7][1] = -1+2*step;
		local_ref[14][7][2] = -1+2*3*step/double(4);

		//element 15
		local_ref[15][0][0] = -1;
		local_ref[15][0][1] = -1;
		local_ref[15][0][2] = -1+step;

		local_ref[15][1][0] = -1+step;
		local_ref[15][1][1] = -1;
		local_ref[15][1][2] = -1+2*3*step/double(4);

		local_ref[15][2][0] = -1+step;
		local_ref[15][2][1] = -1+step;
		local_ref[15][2][2] = -1+2*3*step/double(4);

		local_ref[15][3][0] = -1;
		local_ref[15][3][1] = -1+step;
		local_ref[15][3][2] = -1+2*3*step/double(4);

		local_ref[15][4][0] = -1;
		local_ref[15][4][1] = -1;
		local_ref[15][4][2] = -1+2*step;//-1+2*3*step/double(4);

		local_ref[15][5][0] = -1+step;
		local_ref[15][5][1] = -1;
		local_ref[15][5][2] = -1+3*3*step/double(4);

		local_ref[15][6][0] = -1+step;
		local_ref[15][6][1] = -1+step;
		local_ref[15][6][2] = -1+3*3*step/double(4);

		local_ref[15][7][0] = -1;
		local_ref[15][7][1] = -1+step;
		local_ref[15][7][2] = -1+3*3*step/double(4);


		//element 16
		local_ref[16][0][0] = -1+step;
		local_ref[16][0][1] = -1;
		local_ref[16][0][2] = -1+2*3*step/double(4);

		local_ref[16][1][0] = -1+2*step;
		local_ref[16][1][1] = -1;
		local_ref[16][1][2] = -1+2*3*step/double(4);

		local_ref[16][2][0] = -1+2*step;
		local_ref[16][2][1] = -1+step;
		local_ref[16][2][2] = -1+2*3*step/double(4);

		local_ref[16][3][0] = -1+step;
		local_ref[16][3][1] = -1+step;
		local_ref[16][3][2] = -1+2*3*step/double(4);

		local_ref[16][4][0] = -1+step;
		local_ref[16][4][1] = -1;
		local_ref[16][4][2] = -1+3*3*step/double(4);

		local_ref[16][5][0] = -1+2*step;
		local_ref[16][5][1] = -1;
		local_ref[16][5][2] = -1+3*3*step/double(4);

		local_ref[16][6][0] = -1+2*step;
		local_ref[16][6][1] = -1+1*step;
		local_ref[16][6][2] = -1+3*3*step/double(4);

		local_ref[16][7][0] = -1+step;
		local_ref[16][7][1] = -1+1*step;
		local_ref[16][7][2] = -1+3*3*step/double(4);

		//element 17
		local_ref[17][0][0] = -1+2*step;
		local_ref[17][0][1] = -1;
		local_ref[17][0][2] = -1+2*3*step/double(4);

		local_ref[17][1][0] = 1;
		local_ref[17][1][1] = -1;
		local_ref[17][1][2] = -1+step;

		local_ref[17][2][0] =  1;
		local_ref[17][2][1] = -1+step;
		local_ref[17][2][2] = -1+2*3*step/double(4);

		local_ref[17][3][0] = -1+2*step;
		local_ref[17][3][1] = -1+1*step;
		local_ref[17][3][2] = -1+2*3*step/double(4);

		local_ref[17][4][0] = -1+2*step;
		local_ref[17][4][1] = -1;
		local_ref[17][4][2] = -1+3*3*step/double(4);

		local_ref[17][5][0] = 1;
		local_ref[17][5][1] = -1;
		local_ref[17][5][2] = -1+2*step;//-1+2*3*step/double(4);

		local_ref[17][6][0] = 1;
		local_ref[17][6][1] = -1+1*step;
		local_ref[17][6][2] = -1+3*3*step/double(4);

		local_ref[17][7][0] = -1+2*step;
		local_ref[17][7][1] = -1+1*step;
		local_ref[17][7][2] = -1+3*3*step/double(4);

		//element 18
		local_ref[18][0][0] = -1+2*step;
		local_ref[18][0][1] = -1+1*step;
		local_ref[18][0][2] = -1+2*3*step/double(4);

		local_ref[18][1][0] = 1;
		local_ref[18][1][1] = -1+1*step;
		local_ref[18][1][2] = -1+2*3*step/double(4);

		local_ref[18][2][0] =  1;
		local_ref[18][2][1] = -1+2*step;
		local_ref[18][2][2] = -1+2*3*step/double(4);

		local_ref[18][3][0] = -1+2*step;
		local_ref[18][3][1] = -1+2*step;
		local_ref[18][3][2] = -1+2*3*step/double(4);

		local_ref[18][4][0] = -1+2*step;
		local_ref[18][4][1] = -1+1*step;
		local_ref[18][4][2] = -1+3*3*step/double(4);

		local_ref[18][5][0] = 1;
		local_ref[18][5][1] = -1+1*step;
		local_ref[18][5][2] = -1+3*3*step/double(4);

		local_ref[18][6][0] = 1;
		local_ref[18][6][1] = -1+2*step;
		local_ref[18][6][2] = -1+3*3*step/double(4);

		local_ref[18][7][0] = -1+2*step;
		local_ref[18][7][1] = -1+2*step;
		local_ref[18][7][2] = -1+3*3*step/double(4);

		//element 19
		local_ref[19][0][0] = -1+2*step;
		local_ref[19][0][1] = -1+2*step;
		local_ref[19][0][2] = -1+2*3*step/double(4);

		local_ref[19][1][0] = 1;
		local_ref[19][1][1] = -1+2*step;
		local_ref[19][1][2] = -1+2*3*step/double(4);

		local_ref[19][2][0] =  1;
		local_ref[19][2][1] =  1;
		local_ref[19][2][2] = -1+step;

		local_ref[19][3][0] = -1+2*step;
		local_ref[19][3][1] =  1;
		local_ref[19][3][2] = -1+2*3*step/double(4);

		local_ref[19][4][0] = -1+2*step;
		local_ref[19][4][1] = -1+2*step;
		local_ref[19][4][2] = -1+3*3*step/double(4);

		local_ref[19][5][0] = 1;
		local_ref[19][5][1] = -1+2*step;
		local_ref[19][5][2] = -1+3*3*step/double(4);

		local_ref[19][6][0] = 1;
		local_ref[19][6][1] = 1;
		local_ref[19][6][2] = -1+2*step;

		local_ref[19][7][0] = -1+2*step;
		local_ref[19][7][1] = 1;
		local_ref[19][7][2] = -1+3*3*step/double(4);

		//element 20
		local_ref[20][0][0] = -1+1*step;
		local_ref[20][0][1] = -1+2*step;
		local_ref[20][0][2] = -1+2*3*step/double(4);

		local_ref[20][1][0] = -1+2*step;
		local_ref[20][1][1] = -1+2*step;
		local_ref[20][1][2] = -1+2*3*step/double(4);

		local_ref[20][2][0] =  -1+2*step;
		local_ref[20][2][1] =  1;
		local_ref[20][2][2] = -1+2*3*step/double(4);

		local_ref[20][3][0] = -1+1*step;
		local_ref[20][3][1] =  1;
		local_ref[20][3][2] = -1+2*3*step/double(4);

		local_ref[20][4][0] = -1+1*step;
		local_ref[20][4][1] = -1+2*step;
		local_ref[20][4][2] = -1+3*3*step/double(4);

		local_ref[20][5][0] = -1+2*step;
		local_ref[20][5][1] = -1+2*step;
		local_ref[20][5][2] = -1+3*3*step/double(4);

		local_ref[20][6][0] = -1+2*step;
		local_ref[20][6][1] = 1;
		local_ref[20][6][2] = -1+3*3*step/double(4);

		local_ref[20][7][0] = -1+1*step;
		local_ref[20][7][1] = 1;
		local_ref[20][7][2] = -1+3*3*step/double(4);

		//element 21
		local_ref[21][0][0] = -1;
		local_ref[21][0][1] = -1+2*step;
		local_ref[21][0][2] = -1+2*3*step/double(4);

		local_ref[21][1][0] = -1+1*step;
		local_ref[21][1][1] = -1+2*step;
		local_ref[21][1][2] = -1+2*3*step/double(4);

		local_ref[21][2][0] =  -1+1*step;
		local_ref[21][2][1] =  1;
		local_ref[21][2][2] = -1+2*3*step/double(4);

		local_ref[21][3][0] = -1;
		local_ref[21][3][1] =  1;
		local_ref[21][3][2] = -1+step;

		local_ref[21][4][0] = -1;
		local_ref[21][4][1] = -1+2*step;
		local_ref[21][4][2] = -1+3*3*step/double(4);

		local_ref[21][5][0] = -1+1*step;
		local_ref[21][5][1] = -1+2*step;
		local_ref[21][5][2] = -1+3*3*step/double(4);

		local_ref[21][6][0] = -1+1*step;
		local_ref[21][6][1] = 1;
		local_ref[21][6][2] = -1+3*3*step/double(4);

		local_ref[21][7][0] = -1;
		local_ref[21][7][1] = 1;
		local_ref[21][7][2] = -1+2*step;//-1+2*3*step/double(4);

		//element 22
		local_ref[22][0][0] = -1;
		local_ref[22][0][1] = -1+1*step;
		local_ref[22][0][2] = -1+2*3*step/double(4);

		local_ref[22][1][0] = -1+1*step;
		local_ref[22][1][1] = -1+1*step;
		local_ref[22][1][2] = -1+2*3*step/double(4);

		local_ref[22][2][0] =  -1+1*step;
		local_ref[22][2][1] =  -1+2*step;
		local_ref[22][2][2] = -1+2*3*step/double(4);

		local_ref[22][3][0] = -1;
		local_ref[22][3][1] =  -1+2*step;
		local_ref[22][3][2] = -1+2*3*step/double(4);

		local_ref[22][4][0] = -1;
		local_ref[22][4][1] = -1+1*step;
		local_ref[22][4][2] = -1+3*3*step/double(4);

		local_ref[22][5][0] = -1+1*step;
		local_ref[22][5][1] = -1+1*step;
		local_ref[22][5][2] = -1+3*3*step/double(4);

		local_ref[22][6][0] = -1+1*step;
		local_ref[22][6][1] = -1+2*step;
		local_ref[22][6][2] = -1+3*3*step/double(4);

		local_ref[22][7][0] = -1;
		local_ref[22][7][1] = -1+2*step;
		local_ref[22][7][2] = -1+3*3*step/double(4);

		//element 23
		local_ref[23][0][0] = -1;
		local_ref[23][0][1] = -1;
		local_ref[23][0][2] = -1+2*step;

		local_ref[23][1][0] = -1+step;
		local_ref[23][1][1] = -1;
		local_ref[23][1][2] = -1+3*3*step/double(4);

		local_ref[23][2][0] = -1+step;
		local_ref[23][2][1] = -1+step;
		local_ref[23][2][2] = -1+3*3*step/double(4);

		local_ref[23][3][0] = -1;
		local_ref[23][3][1] = -1+step;
		local_ref[23][3][2] = -1+3*3*step/double(4);

		local_ref[23][4][0] = -1;
		local_ref[23][4][1] = -1;
		local_ref[23][4][2] = 1;

		local_ref[23][5][0] = -1+step;
		local_ref[23][5][1] = -1;
		local_ref[23][5][2] = 1;

		local_ref[23][6][0] = -1+step;
		local_ref[23][6][1] = -1+step;
		local_ref[23][6][2] = 1;

		local_ref[23][7][0] = -1;
		local_ref[23][7][1] = -1+step;
		local_ref[23][7][2] = 1;


		//element 24
		local_ref[24][0][0] = -1+step;
		local_ref[24][0][1] = -1;
		local_ref[24][0][2] = -1+3*3*step/double(4);

		local_ref[24][1][0] = -1+2*step;
		local_ref[24][1][1] = -1;
		local_ref[24][1][2] = -1+3*3*step/double(4);

		local_ref[24][2][0] = -1+2*step;
		local_ref[24][2][1] = -1+step;
		local_ref[24][2][2] = -1+3*3*step/double(4);

		local_ref[24][3][0] = -1+step;
		local_ref[24][3][1] = -1+step;
		local_ref[24][3][2] = -1+3*3*step/double(4);

		local_ref[24][4][0] = -1+step;
		local_ref[24][4][1] = -1;
		local_ref[24][4][2] = 1;

		local_ref[24][5][0] = -1+2*step;
		local_ref[24][5][1] = -1;
		local_ref[24][5][2] = 1;

		local_ref[24][6][0] = -1+2*step;
		local_ref[24][6][1] = -1+1*step;
		local_ref[24][6][2] = 1;

		local_ref[24][7][0] = -1+step;
		local_ref[24][7][1] = -1+1*step;
		local_ref[24][7][2] = 1;

		//element 25
		local_ref[25][0][0] = -1+2*step;
		local_ref[25][0][1] = -1;
		local_ref[25][0][2] = -1+3*3*step/double(4);

		local_ref[25][1][0] = 1;
		local_ref[25][1][1] = -1;
		local_ref[25][1][2] = -1+2*step;

		local_ref[25][2][0] =  1;
		local_ref[25][2][1] = -1+step;
		local_ref[25][2][2] = -1+3*3*step/double(4);

		local_ref[25][3][0] = -1+2*step;
		local_ref[25][3][1] = -1+1*step;
		local_ref[25][3][2] = -1+3*3*step/double(4);

		local_ref[25][4][0] = -1+2*step;
		local_ref[25][4][1] = -1;
		local_ref[25][4][2] = 1;

		local_ref[25][5][0] = 1;
		local_ref[25][5][1] = -1;
		local_ref[25][5][2] = 1;

		local_ref[25][6][0] = 1;
		local_ref[25][6][1] = -1+1*step;
		local_ref[25][6][2] = 1;

		local_ref[25][7][0] = -1+2*step;
		local_ref[25][7][1] = -1+1*step;
		local_ref[25][7][2] = 1;

		//element 26
		local_ref[26][0][0] = -1+2*step;
		local_ref[26][0][1] = -1+1*step;
		local_ref[26][0][2] = -1+3*3*step/double(4);

		local_ref[26][1][0] = 1;
		local_ref[26][1][1] = -1+1*step;
		local_ref[26][1][2] = -1+3*3*step/double(4);

		local_ref[26][2][0] =  1;
		local_ref[26][2][1] = -1+2*step;
		local_ref[26][2][2] = -1+3*3*step/double(4);

		local_ref[26][3][0] = -1+2*step;
		local_ref[26][3][1] = -1+2*step;
		local_ref[26][3][2] = -1+3*3*step/double(4);

		local_ref[26][4][0] = -1+2*step;
		local_ref[26][4][1] = -1+1*step;
		local_ref[26][4][2] = 1;

		local_ref[26][5][0] = 1;
		local_ref[26][5][1] = -1+1*step;
		local_ref[26][5][2] = 1;

		local_ref[26][6][0] = 1;
		local_ref[26][6][1] = -1+2*step;
		local_ref[26][6][2] = 1;

		local_ref[26][7][0] = -1+2*step;
		local_ref[26][7][1] = -1+2*step;
		local_ref[26][7][2] = 1;

		//element 27
		local_ref[27][0][0] = -1+2*step;
		local_ref[27][0][1] = -1+2*step;
		local_ref[27][0][2] = -1+3*3*step/double(4);

		local_ref[27][1][0] = 1;
		local_ref[27][1][1] = -1+2*step;
		local_ref[27][1][2] = -1+3*3*step/double(4);

		local_ref[27][2][0] =  1;
		local_ref[27][2][1] =  1;
		local_ref[27][2][2] = -1+2*step;

		local_ref[27][3][0] = -1+2*step;
		local_ref[27][3][1] =  1;
		local_ref[27][3][2] = -1+3*3*step/double(4);

		local_ref[27][4][0] = -1+2*step;
		local_ref[27][4][1] = -1+2*step;
		local_ref[27][4][2] = 1;

		local_ref[27][5][0] = 1;
		local_ref[27][5][1] = -1+2*step;
		local_ref[27][5][2] = 1;

		local_ref[27][6][0] = 1;
		local_ref[27][6][1] = 1;
		local_ref[27][6][2] = 1;

		local_ref[27][7][0] = -1+2*step;
		local_ref[27][7][1] = 1;
		local_ref[27][7][2] = 1;

		//element 28
		local_ref[28][0][0] = -1+1*step;
		local_ref[28][0][1] = -1+2*step;
		local_ref[28][0][2] = -1+3*3*step/double(4);

		local_ref[28][1][0] = -1+2*step;
		local_ref[28][1][1] = -1+2*step;
		local_ref[28][1][2] = -1+3*3*step/double(4);

		local_ref[28][2][0] =  -1+2*step;
		local_ref[28][2][1] =  1;
		local_ref[28][2][2] = -1+3*3*step/double(4);

		local_ref[28][3][0] = -1+1*step;
		local_ref[28][3][1] =  1;
		local_ref[28][3][2] = -1+3*3*step/double(4);

		local_ref[28][4][0] = -1+1*step;
		local_ref[28][4][1] = -1+2*step;
		local_ref[28][4][2] = 1;

		local_ref[28][5][0] = -1+2*step;
		local_ref[28][5][1] = -1+2*step;
		local_ref[28][5][2] = 1;

		local_ref[28][6][0] = -1+2*step;
		local_ref[28][6][1] = 1;
		local_ref[28][6][2] = 1;

		local_ref[28][7][0] = -1+1*step;
		local_ref[28][7][1] = 1;
		local_ref[28][7][2] = 1;

		//element 29
		local_ref[29][0][0] = -1;
		local_ref[29][0][1] = -1+2*step;
		local_ref[29][0][2] = -1+3*3*step/double(4);

		local_ref[29][1][0] = -1+1*step;
		local_ref[29][1][1] = -1+2*step;
		local_ref[29][1][2] = -1+3*3*step/double(4);

		local_ref[29][2][0] =  -1+1*step;
		local_ref[29][2][1] =  1;
		local_ref[29][2][2] = -1+3*3*step/double(4);

		local_ref[29][3][0] = -1;
		local_ref[29][3][1] =  1;
		local_ref[29][3][2] = -1+2*step;

		local_ref[29][4][0] = -1;
		local_ref[29][4][1] = -1+2*step;
		local_ref[29][4][2] = 1;

		local_ref[29][5][0] = -1+1*step;
		local_ref[29][5][1] = -1+2*step;
		local_ref[29][5][2] = 1;

		local_ref[29][6][0] = -1+1*step;
		local_ref[29][6][1] = 1;
		local_ref[29][6][2] = 1;

		local_ref[29][7][0] = -1;
		local_ref[29][7][1] = 1;
		local_ref[29][7][2] = 1;

		//element 30
		local_ref[30][0][0] = -1;
		local_ref[30][0][1] = -1+1*step;
		local_ref[30][0][2] = -1+3*3*step/double(4);

		local_ref[30][1][0] = -1+1*step;
		local_ref[30][1][1] = -1+1*step;
		local_ref[30][1][2] = -1+3*3*step/double(4);

		local_ref[30][2][0] =  -1+1*step;
		local_ref[30][2][1] =  -1+2*step;
		local_ref[30][2][2] = -1+3*3*step/double(4);

		local_ref[30][3][0] = -1;
		local_ref[30][3][1] =  -1+2*step;
		local_ref[30][3][2] = -1+3*3*step/double(4);

		local_ref[30][4][0] = -1;
		local_ref[30][4][1] = -1+1*step;
		local_ref[30][4][2] = 1;

		local_ref[30][5][0] = -1+1*step;
		local_ref[30][5][1] = -1+1*step;
		local_ref[30][5][2] = 1;

		local_ref[30][6][0] = -1+1*step;
		local_ref[30][6][1] = -1+2*step;
		local_ref[30][6][2] = 1;

		local_ref[30][7][0] = -1;
		local_ref[30][7][1] = -1+2*step;
		local_ref[30][7][2] = 1;

	}

	//define the rotation of the reference element
	int rot[3];
	int sym[3];
	int id_node[8];
	vector<int> ord;

	if(elem->pad==126){
		//edge 0 1 2 3 4 5 6 7
		rot[0] = 0;
		rot[1] = 0;
		rot[2] = 0;

		sym[0] = 0;
		sym[1] = 0;
		sym[2] = 1;

		ord = RotateHex(rot,sym);

		for(int node_id=0;node_id<8;node_id++){
			id_node[node_id] = elem->nodes[ord[node_id]].id;
		}
	}else if(elem->pad==127){
		//edge 5 6 7 8 9 10 11
		rot[0] = 0;
		rot[1] = 0;
		rot[2] = 0;

		sym[0] = 0;
		sym[1] = 0;
		sym[2] = 0;

		ord = RotateHex(rot,sym);

		for(int node_id=0;node_id<8;node_id++){
			id_node[node_id] = elem->nodes[ord[node_id]].id;
		}
	}else if(elem->pad==128){
		//edge 0 5 8 4 1 9 3 11

		rot[0] = -1;
		rot[1] = 0;
		rot[2] = 0;

		sym[0] = 0;
		sym[1] = 0;
		sym[2] = 0;

		ord = RotateHex(rot,sym);

		for(int node_id=0;node_id<8;node_id++){
			id_node[node_id] = elem->nodes[ord[node_id]].id;
		}
	}else if(elem->pad==129){
		//edge 2 6 10 7 1 9 3 11
		rot[0] = 1;
		rot[1] = 0;
		rot[2] = 0;

		sym[0] = 0;
		sym[1] = 0;
		sym[2] = 0;

		ord = RotateHex(rot,sym);

		for(int node_id=0;node_id<8;node_id++){
			id_node[node_id] = elem->nodes[ord[node_id]].id;
		}
	}else if(elem->pad==130){
		//edge 1 5 9 6 0 8 10 2
		rot[0] = 0;
		rot[1] = -1;
		rot[2] = 0;

		sym[0] = 0;
		sym[1] = 0;
		sym[2] = 0;

		ord = RotateHex(rot,sym);

		for(int node_id=0;node_id<8;node_id++){
			id_node[node_id] = elem->nodes[ord[node_id]].id;
		}
	}else if(elem->pad==131){
		//edge 3 7 11 4 0 8 10 2
		rot[0] = 0;
		rot[1] = 1;
		rot[2] = 0;

		sym[0] = 0;
		sym[1] = 0;
		sym[2] = 0;

		ord = RotateHex(rot,sym);

		for(int node_id=0;node_id<8;node_id++){
			id_node[node_id] = elem->nodes[ord[node_id]].id;
		}
	}else{
		printf ("Error in template 9\n");
		exit (EXIT_FAILURE);
	}


	for(int i=0;i<31;i++){
		int conn_p[8];
		GtsPoint* point[8]={NULL};

		double cord_in_x[8],cord_in_y[8],cord_in_z[8];
		//add the nodes in the coord vector
		for (int ii = 0; ii < 8; ii++){
			cord_in_x[ii]=coords[3*id_node[ii]] ;
			cord_in_y[ii]=coords[3*id_node[ii]+1] ;
			cord_in_z[ii]=coords[3*id_node[ii]+2] ;
			//fprintf(mesh->fdbg,"coord in: %f, %f, %f, in the node: %d\n",cord_in_x[ii],cord_in_y[ii],cord_in_z[ii],elem->nodes[ii].id);
		}

		//add the new nodes in the
		for(int ii=0;ii<8;ii++){
			conn_p[ii] = 0;

			if((local_ref[i][ii][0]==1 || local_ref[i][ii][0]==-1) &&
					(local_ref[i][ii][1]==1 || local_ref[i][ii][1]==-1) &&
					(local_ref[i][ii][2]==1 || local_ref[i][ii][2]==-1)){

				if(local_ref[i][ii][0]==-1 && local_ref[i][ii][1]==-1 && local_ref[i][ii][2]==-1){
					conn_p[ii] = id_node[0];
				}else if(local_ref[i][ii][0]==1 && local_ref[i][ii][1]==-1 && local_ref[i][ii][2]==-1){
					conn_p[ii] = id_node[1];
				}else if(local_ref[i][ii][0]==1 && local_ref[i][ii][1]==1 && local_ref[i][ii][2]==-1){
					conn_p[ii] = id_node[2];
				}else if(local_ref[i][ii][0]==-1 && local_ref[i][ii][1]==1 && local_ref[i][ii][2]==-1){
					conn_p[ii] = id_node[3];
				}else if(local_ref[i][ii][0]==-1 && local_ref[i][ii][1]==-1 && local_ref[i][ii][2]==1){
					conn_p[ii] = id_node[4];
				}else if(local_ref[i][ii][0]==1 && local_ref[i][ii][1]==-1 && local_ref[i][ii][2]==1){
					conn_p[ii] = id_node[5];
				}else if(local_ref[i][ii][0]==1 && local_ref[i][ii][1]==1 && local_ref[i][ii][2]==1){
					conn_p[ii] = id_node[6];
				}else if(local_ref[i][ii][0]==-1 && local_ref[i][ii][1]==1 && local_ref[i][ii][2]==1){
					conn_p[ii] = id_node[7];
				}
				//fprintf(mesh->fdbg,"coord out: %f, %f, %f, in the node: %d\n",coords[3*conn_p[ii]],coords[3*conn_p[ii]+1],coords[3*conn_p[ii]+2],conn_p[ii]);
			}else{
				cord_in_ref[0] = local_ref[i][ii][0];
				cord_in_ref[1] = local_ref[i][ii][1];
				cord_in_ref[2] = local_ref[i][ii][2];

				point[ii] = LinearMapHex(cord_in_ref, cord_in_x,cord_in_y,cord_in_z);
				double var[3];
				var[0] = point[ii]->x;
				var[1] = point[ii]->y;
				var[2] = point[ii]->z;
				conn_p[ii] = AddPoint( mesh, hash_nodes, point[ii] , coords);
				//fprintf(mesh->fdbg,"coord out: %f, %f, %f, in the node: %d\n",var[0],var[1],var[2],conn_p[ii]);

			}
		}

		int aux[8];
		for(int ino=0;ino<8;ino++){
			aux[ino] = conn_p[ino];
		}
		for(int ino=0;ino<8;ino++){
			conn_p[ord[ino]] = aux[ino];
		}

		//add the new elements and make the connectivity
		if(i==0){
			octant_t *elem1 = (octant_t*) sc_array_index(&mesh->elements, elements_ids[iel]);

			elem1->nodes[0].id = conn_p[0];
			elem1->nodes[1].id = conn_p[1];
			elem1->nodes[2].id = conn_p[2];
			elem1->nodes[3].id = conn_p[3];

			elem1->nodes[4].id = conn_p[4];
			elem1->nodes[5].id = conn_p[5];
			elem1->nodes[6].id = conn_p[6];
			elem1->nodes[7].id = conn_p[7];

			CopyPropEl(mesh,id,elem1);

		}else{
			octant_t* elem2 = (octant_t*) sc_array_push(&mesh->elements);

			elem2->nodes[0].id = conn_p[0];
			elem2->nodes[1].id = conn_p[1];
			elem2->nodes[2].id = conn_p[2];
			elem2->nodes[3].id = conn_p[3];

			elem2->nodes[4].id = conn_p[4];
			elem2->nodes[5].id = conn_p[5];
			elem2->nodes[6].id = conn_p[6];
			elem2->nodes[7].id = conn_p[7];

			CopyPropEl(mesh,id,elem2);

		}

	}
	sc_array_reset(&toto);
}

void ApplyTemplate10(hexa_tree_t* mesh, std::vector<double>& coords, std::vector<int>& elements_ids, int iel, sc_hash_array_t* hash_nodes, double step){

	sc_array_t toto;
	sc_array_init(&toto, sizeof(octant_t));

	octant_t* elemOrig = (octant_t*) sc_array_index(&mesh->elements, elements_ids[iel]);
	octant_t * elem = (octant_t*) sc_array_push(&toto);
	hexa_element_copy(elemOrig,elem);

	int id = elements_ids[iel];

	double cord_in_ref[3];
	cord_in_ref[0] = 0;
	cord_in_ref[1] = 0;
	cord_in_ref[2] = 0;

	//reference element edge 0
	double local_ref[17][8][3];
	if(true){
		//element 0
		local_ref[0][0][0] = -1;
		local_ref[0][0][1] = -1;
		local_ref[0][0][2] = -1;

		local_ref[0][1][0] = -1+step;
		local_ref[0][1][1] = -1;
		local_ref[0][1][2] = -1;

		local_ref[0][2][0] = -1+step;
		local_ref[0][2][1] = -1+1*3*step/double(4);
		local_ref[0][2][2] = -1;

		local_ref[0][3][0] = -1;
		local_ref[0][3][1] = -1+step;
		local_ref[0][3][2] = -1;

		local_ref[0][4][0] = -1;
		local_ref[0][4][1] = -1;
		local_ref[0][4][2] = -1+step;

		local_ref[0][5][0] = -1+step;
		local_ref[0][5][1] = -1;
		local_ref[0][5][2] = -1+1*3*step/double(4);

		local_ref[0][6][0] = -1+step;
		local_ref[0][6][1] = -1+1*3*step/double(4);
		local_ref[0][6][2] = -1+1*3*step/double(4);

		local_ref[0][7][0] = -1;
		local_ref[0][7][1] = -1+step;
		local_ref[0][7][2] = -1+step;

		//element 1
		local_ref[1][0][0] = -1+1*step;
		local_ref[1][0][1] = -1;
		local_ref[1][0][2] = -1;

		local_ref[1][1][0] = -1+2*step;
		local_ref[1][1][1] = -1;
		local_ref[1][1][2] = -1;

		local_ref[1][2][0] = -1+2*step;
		local_ref[1][2][1] = -1+1*3*step/double(4);
		local_ref[1][2][2] = -1;

		local_ref[1][3][0] = -1+1*step;
		local_ref[1][3][1] = -1+1*3*step/double(4);
		local_ref[1][3][2] = -1;

		local_ref[1][4][0] = -1+1*step;
		local_ref[1][4][1] = -1;
		local_ref[1][4][2] = -1+1*3*step/double(4);

		local_ref[1][5][0] = -1+2*step;
		local_ref[1][5][1] = -1;
		local_ref[1][5][2] = -1+1*3*step/double(4);

		local_ref[1][6][0] = -1+2*step;
		local_ref[1][6][1] = -1+1*3*step/double(4);
		local_ref[1][6][2] = -1+1*3*step/double(4);

		local_ref[1][7][0] = -1+1*step;
		local_ref[1][7][1] = -1+1*3*step/double(4);
		local_ref[1][7][2] = -1+1*3*step/double(4);


		//element 2
		local_ref[2][0][0] = -1+2*step;
		local_ref[2][0][1] = -1;
		local_ref[2][0][2] = -1;

		local_ref[2][1][0] =  1;
		local_ref[2][1][1] = -1;
		local_ref[2][1][2] = -1;

		local_ref[2][2][0] =  1;
		local_ref[2][2][1] = -1+step;
		local_ref[2][2][2] = -1;

		local_ref[2][3][0] = -1+2*step;
		local_ref[2][3][1] = -1+1*3*step/double(4);
		local_ref[2][3][2] = -1;

		local_ref[2][4][0] = -1+2*step;
		local_ref[2][4][1] = -1;
		local_ref[2][4][2] = -1+1*3*step/double(4);

		local_ref[2][5][0] =  1;
		local_ref[2][5][1] = -1;
		local_ref[2][5][2] =  -1+step;

		local_ref[2][6][0] = 1;
		local_ref[2][6][1] = -1+step;
		local_ref[2][6][2] = -1+step;

		local_ref[2][7][0] = -1+2*step;
		local_ref[2][7][1] = -1+1*3*step/double(4);
		local_ref[2][7][2] = -1+1*3*step/double(4);

		//element 3
		local_ref[3][0][0] = -1;
		local_ref[3][0][1] = -1;
		local_ref[3][0][2] = -1+step;

		local_ref[3][1][0] = -1+step;
		local_ref[3][1][1] = -1;
		local_ref[3][1][2] = -1+1*3*step/double(4);

		local_ref[3][2][0] = -1+step;
		local_ref[3][2][1] = -1+1*3*step/double(4);
		local_ref[3][2][2] = -1+1*3*step/double(4);

		local_ref[3][3][0] = -1;
		local_ref[3][3][1] = -1+step;
		local_ref[3][3][2] = -1+step;

		local_ref[3][4][0] = -1;
		local_ref[3][4][1] = -1;
		local_ref[3][4][2] = -1+2*step;

		local_ref[3][5][0] = -1+step;
		local_ref[3][5][1] = -1;
		local_ref[3][5][2] = -1+2*3*step/double(4);

		local_ref[3][6][0] = -1+step;
		local_ref[3][6][1] = -1+2*3*step/double(4);
		local_ref[3][6][2] = -1+2*3*step/double(4);

		local_ref[3][7][0] = -1;
		local_ref[3][7][1] = -1+2*step;
		local_ref[3][7][2] = -1+2*step;

		//element 4
		local_ref[4][0][0] = -1+step;
		local_ref[4][0][1] = -1;
		local_ref[4][0][2] = -1+1*3*step/double(4);

		local_ref[4][1][0] = -1+2*step;
		local_ref[4][1][1] = -1;
		local_ref[4][1][2] = -1+1*3*step/double(4);

		local_ref[4][2][0] = -1+2*step;
		local_ref[4][2][1] = -1+1*3*step/double(4);
		local_ref[4][2][2] = -1+1*3*step/double(4);

		local_ref[4][3][0] = -1+step;
		local_ref[4][3][1] = -1+1*3*step/double(4);
		local_ref[4][3][2] = -1+1*3*step/double(4);

		local_ref[4][4][0] = -1+step;
		local_ref[4][4][1] = -1;
		local_ref[4][4][2] = -1+2*3*step/double(4);

		local_ref[4][5][0] = -1+2*step;
		local_ref[4][5][1] = -1;
		local_ref[4][5][2] = -1+2*3*step/double(4);

		local_ref[4][6][0] = -1+2*step;
		local_ref[4][6][1] = -1+2*3*step/double(4);
		local_ref[4][6][2] = -1+2*3*step/double(4);

		local_ref[4][7][0] = -1+step;
		local_ref[4][7][1] = -1+2*3*step/double(4);
		local_ref[4][7][2] = -1+2*3*step/double(4);

		//element 5
		local_ref[5][0][0] = -1+2*step;
		local_ref[5][0][1] = -1;
		local_ref[5][0][2] = -1+1*3*step/double(4);

		local_ref[5][1][0] =  1;
		local_ref[5][1][1] = -1;
		local_ref[5][1][2] = -1+step;

		local_ref[5][2][0] = 1;
		local_ref[5][2][1] = -1+step;
		local_ref[5][2][2] = -1+step;

		local_ref[5][3][0] = -1+2*step;
		local_ref[5][3][1] = -1+1*3*step/double(4);
		local_ref[5][3][2] = -1+1*3*step/double(4);

		local_ref[5][4][0] = -1+2*step;
		local_ref[5][4][1] = -1;
		local_ref[5][4][2] = -1+2*3*step/double(4);

		local_ref[5][5][0] = 1;
		local_ref[5][5][1] = -1;
		local_ref[5][5][2] = -1+2*step;

		local_ref[5][6][0] = 1;
		local_ref[5][6][1] = -1+2*step;
		local_ref[5][6][2] = -1+2*step;

		local_ref[5][7][0] = -1+2*step;
		local_ref[5][7][1] = -1+2*3*step/double(4);
		local_ref[5][7][2] = -1+2*3*step/double(4);

		//element 6
		local_ref[6][0][0] = -1;
		local_ref[6][0][1] = -1;
		local_ref[6][0][2] = -1+2*step;

		local_ref[6][1][0] = -1+step;
		local_ref[6][1][1] = -1;
		local_ref[6][1][2] = -1+2*3*step/double(4);

		local_ref[6][2][0] = -1+step;
		local_ref[6][2][1] = -1+2*3*step/double(4);
		local_ref[6][2][2] = -1+2*3*step/double(4);

		local_ref[6][3][0] = -1;
		local_ref[6][3][1] = -1+2*step;
		local_ref[6][3][2] = -1+2*step;

		local_ref[6][4][0] = -1;
		local_ref[6][4][1] = -1;
		local_ref[6][4][2] = 1;

		local_ref[6][5][0] = -1+step;
		local_ref[6][5][1] = -1;
		local_ref[6][5][2] = -1+3*3*step/double(4);

		local_ref[6][6][0] = -1+step;
		local_ref[6][6][1] = -1+3*3*step/double(4);
		local_ref[6][6][2] = -1+3*3*step/double(4);

		local_ref[6][7][0] = -1;
		local_ref[6][7][1] = 1;
		local_ref[6][7][2] = 1;

		//element 7
		local_ref[7][0][0] = -1+step;
		local_ref[7][0][1] = -1;
		local_ref[7][0][2] = -1+2*3*step/double(4);

		local_ref[7][1][0] = -1+2*step;
		local_ref[7][1][1] = -1;
		local_ref[7][1][2] = -1+2*3*step/double(4);

		local_ref[7][2][0] = -1+2*step;
		local_ref[7][2][1] = -1+2*3*step/double(4);
		local_ref[7][2][2] = -1+2*3*step/double(4);

		local_ref[7][3][0] = -1+step;
		local_ref[7][3][1] = -1+2*3*step/double(4);
		local_ref[7][3][2] = -1+2*3*step/double(4);

		local_ref[7][4][0] = -1+step;
		local_ref[7][4][1] = -1;
		local_ref[7][4][2] = -1+3*3*step/double(4);

		local_ref[7][5][0] = -1+2*step;
		local_ref[7][5][1] = -1;
		local_ref[7][5][2] = -1+3*3*step/double(4);

		local_ref[7][6][0] = -1+2*step;
		local_ref[7][6][1] = -1+3*3*step/double(4);
		local_ref[7][6][2] = -1+3*3*step/double(4);

		local_ref[7][7][0] = -1+step;
		local_ref[7][7][1] = -1+3*3*step/double(4);
		local_ref[7][7][2] = -1+3*3*step/double(4);

		//element 8
		local_ref[8][0][0] = -1+2*step;
		local_ref[8][0][1] = -1;
		local_ref[8][0][2] = -1+2*3*step/double(4);

		local_ref[8][1][0] =  1;
		local_ref[8][1][1] = -1;
		local_ref[8][1][2] = -1+2*step;

		local_ref[8][2][0] = 1;
		local_ref[8][2][1] = -1+2*step;
		local_ref[8][2][2] = -1+2*step;

		local_ref[8][3][0] = -1+2*step;
		local_ref[8][3][1] = -1+2*3*step/double(4);
		local_ref[8][3][2] = -1+2*3*step/double(4);

		local_ref[8][4][0] = -1+2*step;
		local_ref[8][4][1] = -1;
		local_ref[8][4][2] = -1+3*3*step/double(4);

		local_ref[8][5][0] = 1;
		local_ref[8][5][1] = -1;
		local_ref[8][5][2] = 1;

		local_ref[8][6][0] = 1;
		local_ref[8][6][1] = 1;
		local_ref[8][6][2] = 1;

		local_ref[8][7][0] = -1+2*step;
		local_ref[8][7][1] = -1+3*3*step/double(4);
		local_ref[8][7][2] = -1+3*3*step/double(4);

		//element 9
		local_ref[9][0][0] = -1;
		local_ref[9][0][1] = -1+step;
		local_ref[9][0][2] = -1;

		local_ref[9][1][0] = -1+step;
		local_ref[9][1][1] = -1+1*3*step/double(4);
		local_ref[9][1][2] = -1;

		local_ref[9][2][0] = -1+step;
		local_ref[9][2][1] = -1+2*3*step/double(4);
		local_ref[9][2][2] = -1;

		local_ref[9][3][0] = -1;
		local_ref[9][3][1] = -1+2*step;;
		local_ref[9][3][2] = -1;

		local_ref[9][4][0] = -1;
		local_ref[9][4][1] = -1+step;
		local_ref[9][4][2] = -1+step;

		local_ref[9][5][0] = -1+step;
		local_ref[9][5][1] = -1+1*3*step/double(4);
		local_ref[9][5][2] = -1+1*3*step/double(4);

		local_ref[9][6][0] = -1+step;
		local_ref[9][6][1] = -1+2*3*step/double(4);
		local_ref[9][6][2] = -1+2*3*step/double(4);

		local_ref[9][7][0] = -1;
		local_ref[9][7][1] = -1+2*step;
		local_ref[9][7][2] = -1+2*step;

		//element 10
		local_ref[10][0][0] = -1+step;;
		local_ref[10][0][1] = -1+1*3*step/double(4);
		local_ref[10][0][2] = -1;

		local_ref[10][1][0] = -1+2*step;
		local_ref[10][1][1] = -1+1*3*step/double(4);
		local_ref[10][1][2] = -1;

		local_ref[10][2][0] = -1+2*step;
		local_ref[10][2][1] = -1+2*3*step/double(4);
		local_ref[10][2][2] = -1;

		local_ref[10][3][0] = -1+step;
		local_ref[10][3][1] = -1+2*3*step/double(4);
		local_ref[10][3][2] = -1;

		local_ref[10][4][0] = -1+step;
		local_ref[10][4][1] = -1+1*3*step/double(4);
		local_ref[10][4][2] = -1+1*3*step/double(4);

		local_ref[10][5][0] = -1+2*step;
		local_ref[10][5][1] = -1+1*3*step/double(4);
		local_ref[10][5][2] = -1+1*3*step/double(4);

		local_ref[10][6][0] = -1+2*step;
		local_ref[10][6][1] = -1+2*3*step/double(4);
		local_ref[10][6][2] = -1+2*3*step/double(4);

		local_ref[10][7][0] = -1+step;
		local_ref[10][7][1] = -1+2*3*step/double(4);
		local_ref[10][7][2] = -1+2*3*step/double(4);

		//element 11
		local_ref[11][0][0] = -1+2*step;
		local_ref[11][0][1] = -1+1*3*step/double(4);
		local_ref[11][0][2] = -1;

		local_ref[11][1][0] =  1;
		local_ref[11][1][1] = -1+step;
		local_ref[11][1][2] = -1;

		local_ref[11][2][0] = 1;
		local_ref[11][2][1] = -1+2*step;
		local_ref[11][2][2] = -1;

		local_ref[11][3][0] = -1+2*step;
		local_ref[11][3][1] = -1+2*3*step/double(4);
		local_ref[11][3][2] = -1;

		local_ref[11][4][0] =  -1+2*step;
		local_ref[11][4][1] = -1+1*3*step/double(4);
		local_ref[11][4][2] = -1+1*3*step/double(4);

		local_ref[11][5][0] =  1;
		local_ref[11][5][1] = -1+step;
		local_ref[11][5][2] = -1+step;

		local_ref[11][6][0] = 1;
		local_ref[11][6][1] = -1+2*step;
		local_ref[11][6][2] = -1+2*step;

		local_ref[11][7][0] = -1+2*step;
		local_ref[11][7][1] = -1+2*3*step/double(4);
		local_ref[11][7][2] = -1+2*3*step/double(4);

		//element 12
		local_ref[12][0][0] = -1;
		local_ref[12][0][1] = -1+2*step;
		local_ref[12][0][2] = -1;

		local_ref[12][1][0] = -1+step;
		local_ref[12][1][1] = -1+2*3*step/double(4);
		local_ref[12][1][2] = -1;

		local_ref[12][2][0] = -1+step;
		local_ref[12][2][1] = -1+3*3*step/double(4);
		local_ref[12][2][2] = -1;

		local_ref[12][3][0] = -1;
		local_ref[12][3][1] = 1;
		local_ref[12][3][2] = -1;

		local_ref[12][4][0] = -1;
		local_ref[12][4][1] = -1+2*step;
		local_ref[12][4][2] = -1+2*step;

		local_ref[12][5][0] = -1+step;
		local_ref[12][5][1] = -1+2*3*step/double(4);
		local_ref[12][5][2] = -1+2*3*step/double(4);

		local_ref[12][6][0] = -1+step;
		local_ref[12][6][1] = -1+3*3*step/double(4);
		local_ref[12][6][2] = -1+3*3*step/double(4);

		local_ref[12][7][0] = -1;
		local_ref[12][7][1] = 1;
		local_ref[12][7][2] = 1;

		//element 13
		local_ref[13][0][0] = -1+step;;
		local_ref[13][0][1] = -1+2*3*step/double(4);
		local_ref[13][0][2] = -1;

		local_ref[13][1][0] = -1+2*step;
		local_ref[13][1][1] = -1+2*3*step/double(4);
		local_ref[13][1][2] = -1;

		local_ref[13][2][0] = -1+2*step;
		local_ref[13][2][1] = -1+3*3*step/double(4);
		local_ref[13][2][2] = -1;

		local_ref[13][3][0] = -1+step;
		local_ref[13][3][1] = -1+3*3*step/double(4);
		local_ref[13][3][2] = -1;

		local_ref[13][4][0] = -1+step;
		local_ref[13][4][1] = -1+2*3*step/double(4);
		local_ref[13][4][2] = -1+2*3*step/double(4);

		local_ref[13][5][0] = -1+2*step;
		local_ref[13][5][1] = -1+2*3*step/double(4);
		local_ref[13][5][2] = -1+2*3*step/double(4);

		local_ref[13][6][0] = -1+2*step;
		local_ref[13][6][1] = -1+3*3*step/double(4);
		local_ref[13][6][2] = -1+3*3*step/double(4);

		local_ref[13][7][0] = -1+step;
		local_ref[13][7][1] = -1+3*3*step/double(4);
		local_ref[13][7][2] = -1+3*3*step/double(4);

		//element 14
		local_ref[14][0][0] = -1+2*step;
		local_ref[14][0][1] = -1+2*3*step/double(4);
		local_ref[14][0][2] = -1;

		local_ref[14][1][0] =  1;
		local_ref[14][1][1] = -1+2*step;
		local_ref[14][1][2] = -1;

		local_ref[14][2][0] = 1;
		local_ref[14][2][1] = 1;
		local_ref[14][2][2] = -1;

		local_ref[14][3][0] = -1+2*step;
		local_ref[14][3][1] = -1+3*3*step/double(4);
		local_ref[14][3][2] = -1;

		local_ref[14][4][0] =  -1+2*step;
		local_ref[14][4][1] = -1+2*3*step/double(4);
		local_ref[14][4][2] = -1+2*3*step/double(4);

		local_ref[14][5][0] =  1;
		local_ref[14][5][1] = -1+2*step;
		local_ref[14][5][2] = -1+2*step;

		local_ref[14][6][0] = 1;
		local_ref[14][6][1] = 1;
		local_ref[14][6][2] = 1;

		local_ref[14][7][0] = -1+2*step;
		local_ref[14][7][1] = -1+3*3*step/double(4);
		local_ref[14][7][2] = -1+3*3*step/double(4);

		//element 15
		local_ref[15][0][0] = -1+step;;
		local_ref[15][0][1] = -1;
		local_ref[15][0][2] = -1+3*3*step/double(4);

		local_ref[15][1][0] = -1+2*step;
		local_ref[15][1][1] = -1;
		local_ref[15][1][2] = -1+3*3*step/double(4);

		local_ref[15][2][0] = -1+2*step;
		local_ref[15][2][1] = -1+3*3*step/double(4);
		local_ref[15][2][2] = -1+3*3*step/double(4);

		local_ref[15][3][0] = -1+step;
		local_ref[15][3][1] = -1+3*3*step/double(4);
		local_ref[15][3][2] = -1+3*3*step/double(4);

		local_ref[15][4][0] = -1;
		local_ref[15][4][1] = -1;
		local_ref[15][4][2] = 1;

		local_ref[15][5][0] = 1;
		local_ref[15][5][1] = -1;
		local_ref[15][5][2] = 1;

		local_ref[15][6][0] = 1;
		local_ref[15][6][1] = 1;
		local_ref[15][6][2] = 1;

		local_ref[15][7][0] = -1;
		local_ref[15][7][1] = 1;
		local_ref[15][7][2] = 1;

		//element 16
		local_ref[16][0][0] = -1+1*step;
		local_ref[16][0][1] = -1+3*3*step/double(4);
		local_ref[16][0][2] = -1+3*3*step/double(4);

		local_ref[16][1][0] =  -1+2*step;
		local_ref[16][1][1] = -1+3*3*step/double(4);
		local_ref[16][1][2] = -1+3*3*step/double(4);

		local_ref[16][2][0] = 1;
		local_ref[16][2][1] = 1;
		local_ref[16][2][2] = 1;

		local_ref[16][3][0] = -1;
		local_ref[16][3][1] = 1;
		local_ref[16][3][2] = 1;

		local_ref[16][4][0] = -1+1*step;
		local_ref[16][4][1] = -1+3*3*step/double(4);
		local_ref[16][4][2] = -1;

		local_ref[16][5][0] =  -1+2*step;
		local_ref[16][5][1] = -1+3*3*step/double(4);
		local_ref[16][5][2] = -1;

		local_ref[16][6][0] = 1;
		local_ref[16][6][1] = 1;
		local_ref[16][6][2] = -1;

		local_ref[16][7][0] = -1;
		local_ref[16][7][1] = 1;
		local_ref[16][7][2] = -1;
	}

	//define the rotation of the reference element
	int rot[3];
	int sym[3];
	int id_node[8];
	vector<int> ord;

	if(elem->pad==132){
		//edge 3 0 1 4 5
		rot[0] = 0;
		rot[1] = 0;
		rot[2] = 0;

		sym[0] = 0;
		sym[1] = 0;
		sym[2] = 0;

		ord = RotateHex(rot,sym);

		for(int node_id=0;node_id<8;node_id++){
			id_node[node_id] = elem->nodes[ord[node_id]].id;
		}

	}else if(elem->pad==133){
		//edge 0 1 2 5 6
		rot[0] = 0;
		rot[1] = 0;
		rot[2] = -1;

		sym[0] = 0;
		sym[1] = 0;
		sym[2] = 0;

		ord = RotateHex(rot,sym);

		for(int node_id=0;node_id<8;node_id++){
			id_node[node_id] = elem->nodes[ord[node_id]].id;
		}

	}else if(elem->pad==134){
		//edge 1 2 3 6 7
		rot[0] = 0;
		rot[1] = 0;
		rot[2] = 0;

		sym[0] = 0;
		sym[1] = 1;
		sym[2] = 0;

		ord = RotateHex(rot,sym);

		for(int node_id=0;node_id<8;node_id++){
			id_node[node_id] = elem->nodes[ord[node_id]].id;
		}

	}else if(elem->pad==135){
		//edge 0 2 3 4 7
		rot[0] = 0;
		rot[1] = 0;
		rot[2] = 1;

		sym[0] = 0;
		sym[1] = 0;
		sym[2] = 0;

		ord = RotateHex(rot,sym);

		for(int node_id=0;node_id<8;node_id++){
			id_node[node_id] = elem->nodes[ord[node_id]].id;
		}

	}else if(elem->pad==136){
		//edge 0 3 4 8 11
		rot[0] = 0;
		rot[1] = -1;
		rot[2] = 0;

		sym[0] = 0;
		sym[1] = 0;
		sym[2] = 0;

		ord = RotateHex(rot,sym);

		for(int node_id=0;node_id<8;node_id++){
			id_node[node_id] = elem->nodes[ord[node_id]].id;
		}

	}else if(elem->pad==137){
		//edge 0 1 5 8 9
		rot[0] = 0;
		rot[1] = 1;
		rot[2] = 0;

		sym[0] = 0;
		sym[1] = 0;
		sym[2] = 0;

		ord = RotateHex(rot,sym);

		for(int node_id=0;node_id<8;node_id++){
			id_node[node_id] = elem->nodes[ord[node_id]].id;
		}

	}else if(elem->pad==138){
		//edge 1 2 6 9 10
		rot[0] = 0;
		rot[1] = 1;
		rot[2] = 0;

		sym[0] = 0;
		sym[1] = 1;
		sym[2] = 0;

		ord = RotateHex(rot,sym);

		for(int node_id=0;node_id<8;node_id++){
			id_node[node_id] = elem->nodes[ord[node_id]].id;
		}

	}else if(elem->pad==139){
		//edge 2 3 7 10 11
		rot[0] = 0;
		rot[1] = -1;
		rot[2] = 0;

		sym[0] = 0;
		sym[1] = 1;
		sym[2] = 0;

		ord = RotateHex(rot,sym);

		for(int node_id=0;node_id<8;node_id++){
			id_node[node_id] = elem->nodes[ord[node_id]].id;
		}

	}else if(elem->pad==140){
		//edge 8 9 11 4 5
		rot[0] = 0;
		rot[1] = 0;
		rot[2] = 0;

		sym[0] = 0;
		sym[1] = 0;
		sym[2] = 1;

		ord = RotateHex(rot,sym);

		for(int node_id=0;node_id<8;node_id++){
			id_node[node_id] = elem->nodes[ord[node_id]].id;
		}

	}else if(elem->pad==141){
		//edge 8 9 10 5 6
		rot[0] = 0;
		rot[1] = 0;
		rot[2] = -1;

		sym[0] = 0;
		sym[1] = 0;
		sym[2] = 1;

		ord = RotateHex(rot,sym);

		for(int node_id=0;node_id<8;node_id++){
			id_node[node_id] = elem->nodes[ord[node_id]].id;
		}

	}else if(elem->pad==142){
		//edge 9 10 11 4 7
		rot[0] = 0;
		rot[1] = 0;
		rot[2] = 0;

		sym[0] = 0;
		sym[1] = 1;
		sym[2] = 1;

		ord = RotateHex(rot,sym);

		for(int node_id=0;node_id<8;node_id++){
			id_node[node_id] = elem->nodes[ord[node_id]].id;
		}

	}else if(elem->pad==143){
		//edge 8 10 11 4 7
		rot[0] = 0;
		rot[1] = 0;
		rot[2] = 1;

		sym[0] = 0;
		sym[1] = 0;
		sym[2] = 1;

		ord = RotateHex(rot,sym);

		for(int node_id=0;node_id<8;node_id++){
			id_node[node_id] = elem->nodes[ord[node_id]].id;
		}

	}else{
		printf ("Error in template 10\n");
		exit (EXIT_FAILURE);
	}


	for(int i=0;i<17;i++){
		int conn_p[8];
		GtsPoint* point[8]={NULL};

		double cord_in_x[8],cord_in_y[8],cord_in_z[8];
		//add the nodes in the coord vector
		for (int ii = 0; ii < 8; ii++){
			cord_in_x[ii]=coords[3*id_node[ii]] ;
			cord_in_y[ii]=coords[3*id_node[ii]+1] ;
			cord_in_z[ii]=coords[3*id_node[ii]+2] ;
			//fprintf(mesh->fdbg,"coord in: %f, %f, %f, in the node: %d\n",cord_in_x[ii],cord_in_y[ii],cord_in_z[ii],elem->nodes[ii].id);
		}

		//add the new nodes in the
		for(int ii=0;ii<8;ii++){
			conn_p[ii] = 0;

			if((local_ref[i][ii][0]==1 || local_ref[i][ii][0]==-1) &&
					(local_ref[i][ii][1]==1 || local_ref[i][ii][1]==-1) &&
					(local_ref[i][ii][2]==1 || local_ref[i][ii][2]==-1)){

				if(local_ref[i][ii][0]==-1 && local_ref[i][ii][1]==-1 && local_ref[i][ii][2]==-1){
					conn_p[ii] = id_node[0];
				}else if(local_ref[i][ii][0]==1 && local_ref[i][ii][1]==-1 && local_ref[i][ii][2]==-1){
					conn_p[ii] = id_node[1];
				}else if(local_ref[i][ii][0]==1 && local_ref[i][ii][1]==1 && local_ref[i][ii][2]==-1){
					conn_p[ii] = id_node[2];
				}else if(local_ref[i][ii][0]==-1 && local_ref[i][ii][1]==1 && local_ref[i][ii][2]==-1){
					conn_p[ii] = id_node[3];
				}else if(local_ref[i][ii][0]==-1 && local_ref[i][ii][1]==-1 && local_ref[i][ii][2]==1){
					conn_p[ii] = id_node[4];
				}else if(local_ref[i][ii][0]==1 && local_ref[i][ii][1]==-1 && local_ref[i][ii][2]==1){
					conn_p[ii] = id_node[5];
				}else if(local_ref[i][ii][0]==1 && local_ref[i][ii][1]==1 && local_ref[i][ii][2]==1){
					conn_p[ii] = id_node[6];
				}else if(local_ref[i][ii][0]==-1 && local_ref[i][ii][1]==1 && local_ref[i][ii][2]==1){
					conn_p[ii] = id_node[7];
				}
				//fprintf(mesh->fdbg,"coord out: %f, %f, %f, in the node: %d\n",coords[3*conn_p[ii]],coords[3*conn_p[ii]+1],coords[3*conn_p[ii]+2],conn_p[ii]);
			}else{
				cord_in_ref[0] = local_ref[i][ii][0];
				cord_in_ref[1] = local_ref[i][ii][1];
				cord_in_ref[2] = local_ref[i][ii][2];

				point[ii] = LinearMapHex(cord_in_ref, cord_in_x,cord_in_y,cord_in_z);
				double var[3];
				var[0] = point[ii]->x;
				var[1] = point[ii]->y;
				var[2] = point[ii]->z;
				conn_p[ii] = AddPoint( mesh, hash_nodes, point[ii] , coords);
				//fprintf(mesh->fdbg,"coord out: %f, %f, %f, in the node: %d\n",var[0],var[1],var[2],conn_p[ii]);

			}
		}

		int aux[8];
		for(int ino=0;ino<8;ino++){
			aux[ino] = conn_p[ino];
		}
		for(int ino=0;ino<8;ino++){
			conn_p[ord[ino]] = aux[ino];
		}

		//add the new elements and make the connectivity
		if(i==0){
			octant_t *elem1 = (octant_t*) sc_array_index(&mesh->elements, elements_ids[iel]);

			elem1->nodes[0].id = conn_p[0];
			elem1->nodes[1].id = conn_p[1];
			elem1->nodes[2].id = conn_p[2];
			elem1->nodes[3].id = conn_p[3];

			elem1->nodes[4].id = conn_p[4];
			elem1->nodes[5].id = conn_p[5];
			elem1->nodes[6].id = conn_p[6];
			elem1->nodes[7].id = conn_p[7];

			CopyPropEl(mesh,id,elem1);

		}else{
			octant_t* elem2 = (octant_t*) sc_array_push(&mesh->elements);

			elem2->nodes[0].id = conn_p[0];
			elem2->nodes[1].id = conn_p[1];
			elem2->nodes[2].id = conn_p[2];
			elem2->nodes[3].id = conn_p[3];

			elem2->nodes[4].id = conn_p[4];
			elem2->nodes[5].id = conn_p[5];
			elem2->nodes[6].id = conn_p[6];
			elem2->nodes[7].id = conn_p[7];

			CopyPropEl(mesh,id,elem2);

		}
	}
	sc_array_reset(&toto);
}

void ApplyTemplate11(hexa_tree_t* mesh, std::vector<double>& coords, std::vector<int>& elements_ids, int iel, sc_hash_array_t* hash_nodes, double step){

	sc_array_t toto;
	sc_array_init(&toto, sizeof(octant_t));

	octant_t* elemOrig = (octant_t*) sc_array_index(&mesh->elements, elements_ids[iel]);
	octant_t * elem = (octant_t*) sc_array_push(&toto);
	hexa_element_copy(elemOrig,elem);

	int id = elements_ids[iel];
	double cord_in_ref[3];
	cord_in_ref[0] = 0;
	cord_in_ref[1] = 0;
	cord_in_ref[2] = 0;

	//reference element edge 0
	double local_ref[27][8][3];
	if(true){
		//element 0
		local_ref[0][0][0] = -1+0*step;
		local_ref[0][0][1] = -1+0*step;
		local_ref[0][0][2] = -1+0*step;

		local_ref[0][1][0] = -1+1*step;
		local_ref[0][1][1] = -1+0*step;
		local_ref[0][1][2] = -1+0*step;

		local_ref[0][2][0] = -1+1*step;
		local_ref[0][2][1] = -1+1*step;
		local_ref[0][2][2] = -1+0*step;

		local_ref[0][3][0] = -1+0*step;
		local_ref[0][3][1] = -1+1*step;
		local_ref[0][3][2] = -1+0*step;

		local_ref[0][4][0] = -1+0*step;
		local_ref[0][4][1] = -1+0*step;
		local_ref[0][4][2] = -1+1*step;

		local_ref[0][5][0] = -1+1*step;
		local_ref[0][5][1] = -1+0*step;
		local_ref[0][5][2] = -1+1*step;

		local_ref[0][6][0] = -1+1*step;
		local_ref[0][6][1] = -1+1*step;
		local_ref[0][6][2] = -1+1*step;

		local_ref[0][7][0] = -1+0*step;
		local_ref[0][7][1] = -1+1*step;
		local_ref[0][7][2] = -1+1*step;

		//element 1
		local_ref[1][0][0] = -1+1*step;
		local_ref[1][0][1] = -1+0*step;
		local_ref[1][0][2] = -1+0*step;

		local_ref[1][1][0] = -1+2*step;
		local_ref[1][1][1] = -1+0*step;
		local_ref[1][1][2] = -1+0*step;

		local_ref[1][2][0] = -1+2*step;
		local_ref[1][2][1] = -1+1*step;
		local_ref[1][2][2] = -1+0*step;

		local_ref[1][3][0] = -1+1*step;
		local_ref[1][3][1] = -1+1*step;
		local_ref[1][3][2] = -1+0*step;

		local_ref[1][4][0] = -1+1*step;
		local_ref[1][4][1] = -1+0*step;
		local_ref[1][4][2] = -1+1*step;

		local_ref[1][5][0] = -1+2*step;
		local_ref[1][5][1] = -1+0*step;
		local_ref[1][5][2] = -1+1*step;

		local_ref[1][6][0] = -1+2*step;
		local_ref[1][6][1] = -1+1*step;
		local_ref[1][6][2] = -1+1*step;

		local_ref[1][7][0] = -1+1*step;
		local_ref[1][7][1] = -1+1*step;
		local_ref[1][7][2] = -1+1*step;

		//element 2
		local_ref[2][0][0] = -1+2*step;
		local_ref[2][0][1] = -1+0*step;
		local_ref[2][0][2] = -1+0*step;

		local_ref[2][1][0] = -1+3*step;
		local_ref[2][1][1] = -1+0*step;
		local_ref[2][1][2] = -1+0*step;

		local_ref[2][2][0] = -1+3*step;
		local_ref[2][2][1] = -1+1*step;
		local_ref[2][2][2] = -1+0*step;

		local_ref[2][3][0] = -1+2*step;
		local_ref[2][3][1] = -1+1*step;
		local_ref[2][3][2] = -1+0*step;

		local_ref[2][4][0] = -1+2*step;
		local_ref[2][4][1] = -1+0*step;
		local_ref[2][4][2] = -1+1*step;

		local_ref[2][5][0] = -1+3*step;
		local_ref[2][5][1] = -1+0*step;
		local_ref[2][5][2] = -1+1*step;

		local_ref[2][6][0] = -1+3*step;
		local_ref[2][6][1] = -1+1*step;
		local_ref[2][6][2] = -1+1*step;

		local_ref[2][7][0] = -1+2*step;
		local_ref[2][7][1] = -1+1*step;
		local_ref[2][7][2] = -1+1*step;

		//element 3
		local_ref[3][0][0] = -1+0*step;
		local_ref[3][0][1] = -1+1*step;
		local_ref[3][0][2] = -1+0*step;

		local_ref[3][1][0] = -1+1*step;
		local_ref[3][1][1] = -1+1*step;
		local_ref[3][1][2] = -1+0*step;

		local_ref[3][2][0] = -1+1*step;
		local_ref[3][2][1] = -1+2*step;
		local_ref[3][2][2] = -1+0*step;

		local_ref[3][3][0] = -1+0*step;
		local_ref[3][3][1] = -1+2*step;
		local_ref[3][3][2] = -1+0*step;

		local_ref[3][4][0] = -1+0*step;
		local_ref[3][4][1] = -1+1*step;
		local_ref[3][4][2] = -1+1*step;

		local_ref[3][5][0] = -1+1*step;
		local_ref[3][5][1] = -1+1*step;
		local_ref[3][5][2] = -1+1*step;

		local_ref[3][6][0] = -1+1*step;
		local_ref[3][6][1] = -1+2*step;
		local_ref[3][6][2] = -1+1*step;

		local_ref[3][7][0] = -1+0*step;
		local_ref[3][7][1] = -1+2*step;
		local_ref[3][7][2] = -1+1*step;

		//element 4
		local_ref[4][0][0] = -1+1*step;
		local_ref[4][0][1] = -1+1*step;
		local_ref[4][0][2] = -1+0*step;

		local_ref[4][1][0] = -1+2*step;
		local_ref[4][1][1] = -1+1*step;
		local_ref[4][1][2] = -1+0*step;

		local_ref[4][2][0] = -1+2*step;
		local_ref[4][2][1] = -1+2*step;
		local_ref[4][2][2] = -1+0*step;

		local_ref[4][3][0] = -1+1*step;
		local_ref[4][3][1] = -1+2*step;
		local_ref[4][3][2] = -1+0*step;

		local_ref[4][4][0] = -1+1*step;
		local_ref[4][4][1] = -1+1*step;
		local_ref[4][4][2] = -1+1*step;

		local_ref[4][5][0] = -1+2*step;
		local_ref[4][5][1] = -1+1*step;
		local_ref[4][5][2] = -1+1*step;

		local_ref[4][6][0] = -1+2*step;
		local_ref[4][6][1] = -1+2*step;
		local_ref[4][6][2] = -1+1*step;

		local_ref[4][7][0] = -1+1*step;
		local_ref[4][7][1] = -1+2*step;
		local_ref[4][7][2] = -1+1*step;

		//element 5
		local_ref[5][0][0] = -1+2*step;
		local_ref[5][0][1] = -1+1*step;
		local_ref[5][0][2] = -1+0*step;

		local_ref[5][1][0] = -1+3*step;
		local_ref[5][1][1] = -1+1*step;
		local_ref[5][1][2] = -1+0*step;

		local_ref[5][2][0] = -1+3*step;
		local_ref[5][2][1] = -1+2*step;
		local_ref[5][2][2] = -1+0*step;

		local_ref[5][3][0] = -1+2*step;
		local_ref[5][3][1] = -1+2*step;
		local_ref[5][3][2] = -1+0*step;

		local_ref[5][4][0] = -1+2*step;
		local_ref[5][4][1] = -1+1*step;
		local_ref[5][4][2] = -1+1*step;

		local_ref[5][5][0] = -1+3*step;
		local_ref[5][5][1] = -1+1*step;
		local_ref[5][5][2] = -1+1*step;

		local_ref[5][6][0] = -1+3*step;
		local_ref[5][6][1] = -1+2*step;
		local_ref[5][6][2] = -1+1*step;

		local_ref[5][7][0] = -1+2*step;
		local_ref[5][7][1] = -1+2*step;
		local_ref[5][7][2] = -1+1*step;

		//element 6
		local_ref[6][0][0] = -1+0*step;
		local_ref[6][0][1] = -1+2*step;
		local_ref[6][0][2] = -1+0*step;

		local_ref[6][1][0] = -1+1*step;
		local_ref[6][1][1] = -1+2*step;
		local_ref[6][1][2] = -1+0*step;

		local_ref[6][2][0] = -1+1*step;
		local_ref[6][2][1] = -1+3*step;
		local_ref[6][2][2] = -1+0*step;

		local_ref[6][3][0] = -1+0*step;
		local_ref[6][3][1] = -1+3*step;
		local_ref[6][3][2] = -1+0*step;

		local_ref[6][4][0] = -1+0*step;
		local_ref[6][4][1] = -1+2*step;
		local_ref[6][4][2] = -1+1*step;

		local_ref[6][5][0] = -1+1*step;
		local_ref[6][5][1] = -1+2*step;
		local_ref[6][5][2] = -1+1*step;

		local_ref[6][6][0] = -1+1*step;
		local_ref[6][6][1] = -1+3*step;
		local_ref[6][6][2] = -1+1*step;

		local_ref[6][7][0] = -1+0*step;
		local_ref[6][7][1] = -1+3*step;
		local_ref[6][7][2] = -1+1*step;

		//element 7
		local_ref[7][0][0] = -1+1*step;
		local_ref[7][0][1] = -1+2*step;
		local_ref[7][0][2] = -1+0*step;

		local_ref[7][1][0] = -1+2*step;
		local_ref[7][1][1] = -1+2*step;
		local_ref[7][1][2] = -1+0*step;

		local_ref[7][2][0] = -1+2*step;
		local_ref[7][2][1] = -1+3*step;
		local_ref[7][2][2] = -1+0*step;

		local_ref[7][3][0] = -1+1*step;
		local_ref[7][3][1] = -1+3*step;
		local_ref[7][3][2] = -1+0*step;

		local_ref[7][4][0] = -1+1*step;
		local_ref[7][4][1] = -1+2*step;
		local_ref[7][4][2] = -1+1*step;

		local_ref[7][5][0] = -1+2*step;
		local_ref[7][5][1] = -1+2*step;
		local_ref[7][5][2] = -1+1*step;

		local_ref[7][6][0] = -1+2*step;
		local_ref[7][6][1] = -1+3*step;
		local_ref[7][6][2] = -1+1*step;

		local_ref[7][7][0] = -1+1*step;
		local_ref[7][7][1] = -1+3*step;
		local_ref[7][7][2] = -1+1*step;

		//element 8
		local_ref[8][0][0] = -1+2*step;
		local_ref[8][0][1] = -1+2*step;
		local_ref[8][0][2] = -1+0*step;

		local_ref[8][1][0] = -1+3*step;
		local_ref[8][1][1] = -1+2*step;
		local_ref[8][1][2] = -1+0*step;

		local_ref[8][2][0] = -1+3*step;
		local_ref[8][2][1] = -1+3*step;
		local_ref[8][2][2] = -1+0*step;

		local_ref[8][3][0] = -1+2*step;
		local_ref[8][3][1] = -1+3*step;
		local_ref[8][3][2] = -1+0*step;

		local_ref[8][4][0] = -1+2*step;
		local_ref[8][4][1] = -1+2*step;
		local_ref[8][4][2] = -1+1*step;

		local_ref[8][5][0] = -1+3*step;
		local_ref[8][5][1] = -1+2*step;
		local_ref[8][5][2] = -1+1*step;

		local_ref[8][6][0] = -1+3*step;
		local_ref[8][6][1] = -1+3*step;
		local_ref[8][6][2] = -1+1*step;

		local_ref[8][7][0] = -1+2*step;
		local_ref[8][7][1] = -1+3*step;
		local_ref[8][7][2] = -1+1*step;

		//element 9
		local_ref[9][0][0] = -1+0*step;
		local_ref[9][0][1] = -1+0*step;
		local_ref[9][0][2] = -1+1*step;

		local_ref[9][1][0] = -1+1*step;
		local_ref[9][1][1] = -1+0*step;
		local_ref[9][1][2] = -1+1*step;

		local_ref[9][2][0] = -1+1*step;
		local_ref[9][2][1] = -1+1*step;
		local_ref[9][2][2] = -1+1*step;

		local_ref[9][3][0] = -1+0*step;
		local_ref[9][3][1] = -1+1*step;
		local_ref[9][3][2] = -1+1*step;

		local_ref[9][4][0] = -1+0*step;
		local_ref[9][4][1] = -1+0*step;
		local_ref[9][4][2] = -1+2*step;

		local_ref[9][5][0] = -1+1*step;
		local_ref[9][5][1] = -1+0*step;
		local_ref[9][5][2] = -1+2*step;

		local_ref[9][6][0] = -1+1*step;
		local_ref[9][6][1] = -1+1*step;
		local_ref[9][6][2] = -1+2*step;

		local_ref[9][7][0] = -1+0*step;
		local_ref[9][7][1] = -1+1*step;
		local_ref[9][7][2] = -1+2*step;

		//element 10
		local_ref[10][0][0] = -1+1*step;
		local_ref[10][0][1] = -1+0*step;
		local_ref[10][0][2] = -1+1*step;

		local_ref[10][1][0] = -1+2*step;
		local_ref[10][1][1] = -1+0*step;
		local_ref[10][1][2] = -1+1*step;

		local_ref[10][2][0] = -1+2*step;
		local_ref[10][2][1] = -1+1*step;
		local_ref[10][2][2] = -1+1*step;

		local_ref[10][3][0] = -1+1*step;
		local_ref[10][3][1] = -1+1*step;
		local_ref[10][3][2] = -1+1*step;

		local_ref[10][4][0] = -1+1*step;
		local_ref[10][4][1] = -1+0*step;
		local_ref[10][4][2] = -1+2*step;

		local_ref[10][5][0] = -1+2*step;
		local_ref[10][5][1] = -1+0*step;
		local_ref[10][5][2] = -1+2*step;

		local_ref[10][6][0] = -1+2*step;
		local_ref[10][6][1] = -1+1*step;
		local_ref[10][6][2] = -1+2*step;

		local_ref[10][7][0] = -1+1*step;
		local_ref[10][7][1] = -1+1*step;
		local_ref[10][7][2] = -1+2*step;

		//element 11
		local_ref[11][0][0] = -1+2*step;
		local_ref[11][0][1] = -1+0*step;
		local_ref[11][0][2] = -1+1*step;

		local_ref[11][1][0] = -1+3*step;
		local_ref[11][1][1] = -1+0*step;
		local_ref[11][1][2] = -1+1*step;

		local_ref[11][2][0] = -1+3*step;
		local_ref[11][2][1] = -1+1*step;
		local_ref[11][2][2] = -1+1*step;

		local_ref[11][3][0] = -1+2*step;
		local_ref[11][3][1] = -1+1*step;
		local_ref[11][3][2] = -1+1*step;

		local_ref[11][4][0] = -1+2*step;
		local_ref[11][4][1] = -1+0*step;
		local_ref[11][4][2] = -1+2*step;

		local_ref[11][5][0] = -1+3*step;
		local_ref[11][5][1] = -1+0*step;
		local_ref[11][5][2] = -1+2*step;

		local_ref[11][6][0] = -1+3*step;
		local_ref[11][6][1] = -1+1*step;
		local_ref[11][6][2] = -1+2*step;

		local_ref[11][7][0] = -1+2*step;
		local_ref[11][7][1] = -1+1*step;
		local_ref[11][7][2] = -1+2*step;

		//element 12
		local_ref[12][0][0] = -1+0*step;
		local_ref[12][0][1] = -1+1*step;
		local_ref[12][0][2] = -1+1*step;

		local_ref[12][1][0] = -1+1*step;
		local_ref[12][1][1] = -1+1*step;
		local_ref[12][1][2] = -1+1*step;

		local_ref[12][2][0] = -1+1*step;
		local_ref[12][2][1] = -1+2*step;
		local_ref[12][2][2] = -1+1*step;

		local_ref[12][3][0] = -1+0*step;
		local_ref[12][3][1] = -1+2*step;
		local_ref[12][3][2] = -1+1*step;

		local_ref[12][4][0] = -1+0*step;
		local_ref[12][4][1] = -1+1*step;
		local_ref[12][4][2] = -1+2*step;

		local_ref[12][5][0] = -1+1*step;
		local_ref[12][5][1] = -1+1*step;
		local_ref[12][5][2] = -1+2*step;

		local_ref[12][6][0] = -1+1*step;
		local_ref[12][6][1] = -1+2*step;
		local_ref[12][6][2] = -1+2*step;

		local_ref[12][7][0] = -1+0*step;
		local_ref[12][7][1] = -1+2*step;
		local_ref[12][7][2] = -1+2*step;

		//element 13
		local_ref[13][0][0] = -1+1*step;
		local_ref[13][0][1] = -1+1*step;
		local_ref[13][0][2] = -1+1*step;

		local_ref[13][1][0] = -1+2*step;
		local_ref[13][1][1] = -1+1*step;
		local_ref[13][1][2] = -1+1*step;

		local_ref[13][2][0] = -1+2*step;
		local_ref[13][2][1] = -1+2*step;
		local_ref[13][2][2] = -1+1*step;

		local_ref[13][3][0] = -1+1*step;
		local_ref[13][3][1] = -1+2*step;
		local_ref[13][3][2] = -1+1*step;

		local_ref[13][4][0] = -1+1*step;
		local_ref[13][4][1] = -1+1*step;
		local_ref[13][4][2] = -1+2*step;

		local_ref[13][5][0] = -1+2*step;
		local_ref[13][5][1] = -1+1*step;
		local_ref[13][5][2] = -1+2*step;

		local_ref[13][6][0] = -1+2*step;
		local_ref[13][6][1] = -1+2*step;
		local_ref[13][6][2] = -1+2*step;

		local_ref[13][7][0] = -1+1*step;
		local_ref[13][7][1] = -1+2*step;
		local_ref[13][7][2] = -1+2*step;

		//element 14
		local_ref[14][0][0] = -1+2*step;
		local_ref[14][0][1] = -1+1*step;
		local_ref[14][0][2] = -1+1*step;

		local_ref[14][1][0] = -1+3*step;
		local_ref[14][1][1] = -1+1*step;
		local_ref[14][1][2] = -1+1*step;

		local_ref[14][2][0] = -1+3*step;
		local_ref[14][2][1] = -1+2*step;
		local_ref[14][2][2] = -1+1*step;

		local_ref[14][3][0] = -1+2*step;
		local_ref[14][3][1] = -1+2*step;
		local_ref[14][3][2] = -1+1*step;

		local_ref[14][4][0] = -1+2*step;
		local_ref[14][4][1] = -1+1*step;
		local_ref[14][4][2] = -1+2*step;

		local_ref[14][5][0] = -1+3*step;
		local_ref[14][5][1] = -1+1*step;
		local_ref[14][5][2] = -1+2*step;

		local_ref[14][6][0] = -1+3*step;
		local_ref[14][6][1] = -1+2*step;
		local_ref[14][6][2] = -1+2*step;

		local_ref[14][7][0] = -1+2*step;
		local_ref[14][7][1] = -1+2*step;
		local_ref[14][7][2] = -1+2*step;

		//element 15
		local_ref[15][0][0] = -1+0*step;
		local_ref[15][0][1] = -1+2*step;
		local_ref[15][0][2] = -1+1*step;

		local_ref[15][1][0] = -1+1*step;
		local_ref[15][1][1] = -1+2*step;
		local_ref[15][1][2] = -1+1*step;

		local_ref[15][2][0] = -1+1*step;
		local_ref[15][2][1] = -1+3*step;
		local_ref[15][2][2] = -1+1*step;

		local_ref[15][3][0] = -1+0*step;
		local_ref[15][3][1] = -1+3*step;
		local_ref[15][3][2] = -1+1*step;

		local_ref[15][4][0] = -1+0*step;
		local_ref[15][4][1] = -1+2*step;
		local_ref[15][4][2] = -1+2*step;

		local_ref[15][5][0] = -1+1*step;
		local_ref[15][5][1] = -1+2*step;
		local_ref[15][5][2] = -1+2*step;

		local_ref[15][6][0] = -1+1*step;
		local_ref[15][6][1] = -1+3*step;
		local_ref[15][6][2] = -1+2*step;

		local_ref[15][7][0] = -1+0*step;
		local_ref[15][7][1] = -1+3*step;
		local_ref[15][7][2] = -1+2*step;

		//element 16
		local_ref[16][0][0] = -1+1*step;
		local_ref[16][0][1] = -1+2*step;
		local_ref[16][0][2] = -1+1*step;

		local_ref[16][1][0] = -1+2*step;
		local_ref[16][1][1] = -1+2*step;
		local_ref[16][1][2] = -1+1*step;

		local_ref[16][2][0] = -1+2*step;
		local_ref[16][2][1] = -1+3*step;
		local_ref[16][2][2] = -1+1*step;

		local_ref[16][3][0] = -1+1*step;
		local_ref[16][3][1] = -1+3*step;
		local_ref[16][3][2] = -1+1*step;

		local_ref[16][4][0] = -1+1*step;
		local_ref[16][4][1] = -1+2*step;
		local_ref[16][4][2] = -1+2*step;

		local_ref[16][5][0] = -1+2*step;
		local_ref[16][5][1] = -1+2*step;
		local_ref[16][5][2] = -1+2*step;

		local_ref[16][6][0] = -1+2*step;
		local_ref[16][6][1] = -1+3*step;
		local_ref[16][6][2] = -1+2*step;

		local_ref[16][7][0] = -1+1*step;
		local_ref[16][7][1] = -1+3*step;
		local_ref[16][7][2] = -1+2*step;

		//element 17
		local_ref[17][0][0] = -1+2*step;
		local_ref[17][0][1] = -1+2*step;
		local_ref[17][0][2] = -1+1*step;

		local_ref[17][1][0] = -1+3*step;
		local_ref[17][1][1] = -1+2*step;
		local_ref[17][1][2] = -1+1*step;

		local_ref[17][2][0] = -1+3*step;
		local_ref[17][2][1] = -1+3*step;
		local_ref[17][2][2] = -1+1*step;

		local_ref[17][3][0] = -1+2*step;
		local_ref[17][3][1] = -1+3*step;
		local_ref[17][3][2] = -1+1*step;

		local_ref[17][4][0] = -1+2*step;
		local_ref[17][4][1] = -1+2*step;
		local_ref[17][4][2] = -1+2*step;

		local_ref[17][5][0] = -1+3*step;
		local_ref[17][5][1] = -1+2*step;
		local_ref[17][5][2] = -1+2*step;

		local_ref[17][6][0] = -1+3*step;
		local_ref[17][6][1] = -1+3*step;
		local_ref[17][6][2] = -1+2*step;

		local_ref[17][7][0] = -1+2*step;
		local_ref[17][7][1] = -1+3*step;
		local_ref[17][7][2] = -1+2*step;

		//element 18
		local_ref[18][0][0] = -1+0*step;
		local_ref[18][0][1] = -1+0*step;
		local_ref[18][0][2] = -1+2*step;

		local_ref[18][1][0] = -1+1*step;
		local_ref[18][1][1] = -1+0*step;
		local_ref[18][1][2] = -1+2*step;

		local_ref[18][2][0] = -1+1*step;
		local_ref[18][2][1] = -1+1*step;
		local_ref[18][2][2] = -1+2*step;

		local_ref[18][3][0] = -1+0*step;
		local_ref[18][3][1] = -1+1*step;
		local_ref[18][3][2] = -1+2*step;

		local_ref[18][4][0] = -1+0*step;
		local_ref[18][4][1] = -1+0*step;
		local_ref[18][4][2] = -1+3*step;

		local_ref[18][5][0] = -1+1*step;
		local_ref[18][5][1] = -1+0*step;
		local_ref[18][5][2] = -1+3*step;

		local_ref[18][6][0] = -1+1*step;
		local_ref[18][6][1] = -1+1*step;
		local_ref[18][6][2] = -1+3*step;

		local_ref[18][7][0] = -1+0*step;
		local_ref[18][7][1] = -1+1*step;
		local_ref[18][7][2] = -1+3*step;

		//element 19
		local_ref[19][0][0] = -1+1*step;
		local_ref[19][0][1] = -1+0*step;
		local_ref[19][0][2] = -1+2*step;

		local_ref[19][1][0] = -1+2*step;
		local_ref[19][1][1] = -1+0*step;
		local_ref[19][1][2] = -1+2*step;

		local_ref[19][2][0] = -1+2*step;
		local_ref[19][2][1] = -1+1*step;
		local_ref[19][2][2] = -1+2*step;

		local_ref[19][3][0] = -1+1*step;
		local_ref[19][3][1] = -1+1*step;
		local_ref[19][3][2] = -1+2*step;

		local_ref[19][4][0] = -1+1*step;
		local_ref[19][4][1] = -1+0*step;
		local_ref[19][4][2] = -1+3*step;

		local_ref[19][5][0] = -1+2*step;
		local_ref[19][5][1] = -1+0*step;
		local_ref[19][5][2] = -1+3*step;

		local_ref[19][6][0] = -1+2*step;
		local_ref[19][6][1] = -1+1*step;
		local_ref[19][6][2] = -1+3*step;

		local_ref[19][7][0] = -1+1*step;
		local_ref[19][7][1] = -1+1*step;
		local_ref[19][7][2] = -1+3*step;

		//element 20
		local_ref[20][0][0] = -1+2*step;
		local_ref[20][0][1] = -1+0*step;
		local_ref[20][0][2] = -1+2*step;

		local_ref[20][1][0] = -1+3*step;
		local_ref[20][1][1] = -1+0*step;
		local_ref[20][1][2] = -1+2*step;

		local_ref[20][2][0] = -1+3*step;
		local_ref[20][2][1] = -1+1*step;
		local_ref[20][2][2] = -1+2*step;

		local_ref[20][3][0] = -1+2*step;
		local_ref[20][3][1] = -1+1*step;
		local_ref[20][3][2] = -1+2*step;

		local_ref[20][4][0] = -1+2*step;
		local_ref[20][4][1] = -1+0*step;
		local_ref[20][4][2] = -1+3*step;

		local_ref[20][5][0] = -1+3*step;
		local_ref[20][5][1] = -1+0*step;
		local_ref[20][5][2] = -1+3*step;

		local_ref[20][6][0] = -1+3*step;
		local_ref[20][6][1] = -1+1*step;
		local_ref[20][6][2] = -1+3*step;

		local_ref[20][7][0] = -1+2*step;
		local_ref[20][7][1] = -1+1*step;
		local_ref[20][7][2] = -1+3*step;

		//element 21
		local_ref[21][0][0] = -1+0*step;
		local_ref[21][0][1] = -1+1*step;
		local_ref[21][0][2] = -1+2*step;

		local_ref[21][1][0] = -1+1*step;
		local_ref[21][1][1] = -1+1*step;
		local_ref[21][1][2] = -1+2*step;

		local_ref[21][2][0] = -1+1*step;
		local_ref[21][2][1] = -1+2*step;
		local_ref[21][2][2] = -1+2*step;

		local_ref[21][3][0] = -1+0*step;
		local_ref[21][3][1] = -1+2*step;
		local_ref[21][3][2] = -1+2*step;

		local_ref[21][4][0] = -1+0*step;
		local_ref[21][4][1] = -1+1*step;
		local_ref[21][4][2] = -1+3*step;

		local_ref[21][5][0] = -1+1*step;
		local_ref[21][5][1] = -1+1*step;
		local_ref[21][5][2] = -1+3*step;

		local_ref[21][6][0] = -1+1*step;
		local_ref[21][6][1] = -1+2*step;
		local_ref[21][6][2] = -1+3*step;

		local_ref[21][7][0] = -1+0*step;
		local_ref[21][7][1] = -1+2*step;
		local_ref[21][7][2] = -1+3*step;

		//element 22
		local_ref[22][0][0] = -1+1*step;
		local_ref[22][0][1] = -1+1*step;
		local_ref[22][0][2] = -1+2*step;

		local_ref[22][1][0] = -1+2*step;
		local_ref[22][1][1] = -1+1*step;
		local_ref[22][1][2] = -1+2*step;

		local_ref[22][2][0] = -1+2*step;
		local_ref[22][2][1] = -1+2*step;
		local_ref[22][2][2] = -1+2*step;

		local_ref[22][3][0] = -1+1*step;
		local_ref[22][3][1] = -1+2*step;
		local_ref[22][3][2] = -1+2*step;

		local_ref[22][4][0] = -1+1*step;
		local_ref[22][4][1] = -1+1*step;
		local_ref[22][4][2] = -1+3*step;

		local_ref[22][5][0] = -1+2*step;
		local_ref[22][5][1] = -1+1*step;
		local_ref[22][5][2] = -1+3*step;

		local_ref[22][6][0] = -1+2*step;
		local_ref[22][6][1] = -1+2*step;
		local_ref[22][6][2] = -1+3*step;

		local_ref[22][7][0] = -1+1*step;
		local_ref[22][7][1] = -1+2*step;
		local_ref[22][7][2] = -1+3*step;

		//element 23
		local_ref[23][0][0] = -1+2*step;
		local_ref[23][0][1] = -1+1*step;
		local_ref[23][0][2] = -1+2*step;

		local_ref[23][1][0] = -1+3*step;
		local_ref[23][1][1] = -1+1*step;
		local_ref[23][1][2] = -1+2*step;

		local_ref[23][2][0] = -1+3*step;
		local_ref[23][2][1] = -1+2*step;
		local_ref[23][2][2] = -1+2*step;

		local_ref[23][3][0] = -1+2*step;
		local_ref[23][3][1] = -1+2*step;
		local_ref[23][3][2] = -1+2*step;

		local_ref[23][4][0] = -1+2*step;
		local_ref[23][4][1] = -1+1*step;
		local_ref[23][4][2] = -1+3*step;

		local_ref[23][5][0] = -1+3*step;
		local_ref[23][5][1] = -1+1*step;
		local_ref[23][5][2] = -1+3*step;

		local_ref[23][6][0] = -1+3*step;
		local_ref[23][6][1] = -1+2*step;
		local_ref[23][6][2] = -1+3*step;

		local_ref[23][7][0] = -1+2*step;
		local_ref[23][7][1] = -1+2*step;
		local_ref[23][7][2] = -1+3*step;

		//element 24
		local_ref[24][0][0] = -1+0*step;
		local_ref[24][0][1] = -1+2*step;
		local_ref[24][0][2] = -1+2*step;

		local_ref[24][1][0] = -1+1*step;
		local_ref[24][1][1] = -1+2*step;
		local_ref[24][1][2] = -1+2*step;

		local_ref[24][2][0] = -1+1*step;
		local_ref[24][2][1] = -1+3*step;
		local_ref[24][2][2] = -1+2*step;

		local_ref[24][3][0] = -1+0*step;
		local_ref[24][3][1] = -1+3*step;
		local_ref[24][3][2] = -1+2*step;

		local_ref[24][4][0] = -1+0*step;
		local_ref[24][4][1] = -1+2*step;
		local_ref[24][4][2] = -1+3*step;

		local_ref[24][5][0] = -1+1*step;
		local_ref[24][5][1] = -1+2*step;
		local_ref[24][5][2] = -1+3*step;

		local_ref[24][6][0] = -1+1*step;
		local_ref[24][6][1] = -1+3*step;
		local_ref[24][6][2] = -1+3*step;

		local_ref[24][7][0] = -1+0*step;
		local_ref[24][7][1] = -1+3*step;
		local_ref[24][7][2] = -1+3*step;

		//element 25
		local_ref[25][0][0] = -1+1*step;
		local_ref[25][0][1] = -1+2*step;
		local_ref[25][0][2] = -1+2*step;

		local_ref[25][1][0] = -1+2*step;
		local_ref[25][1][1] = -1+2*step;
		local_ref[25][1][2] = -1+2*step;

		local_ref[25][2][0] = -1+2*step;
		local_ref[25][2][1] = -1+3*step;
		local_ref[25][2][2] = -1+2*step;

		local_ref[25][3][0] = -1+1*step;
		local_ref[25][3][1] = -1+3*step;
		local_ref[25][3][2] = -1+2*step;

		local_ref[25][4][0] = -1+1*step;
		local_ref[25][4][1] = -1+2*step;
		local_ref[25][4][2] = -1+3*step;

		local_ref[25][5][0] = -1+2*step;
		local_ref[25][5][1] = -1+2*step;
		local_ref[25][5][2] = -1+3*step;

		local_ref[25][6][0] = -1+2*step;
		local_ref[25][6][1] = -1+3*step;
		local_ref[25][6][2] = -1+3*step;

		local_ref[25][7][0] = -1+1*step;
		local_ref[25][7][1] = -1+3*step;
		local_ref[25][7][2] = -1+3*step;

		//element 26
		local_ref[26][0][0] = -1+2*step;
		local_ref[26][0][1] = -1+2*step;
		local_ref[26][0][2] = -1+2*step;

		local_ref[26][1][0] = -1+3*step;
		local_ref[26][1][1] = -1+2*step;
		local_ref[26][1][2] = -1+2*step;

		local_ref[26][2][0] = -1+3*step;
		local_ref[26][2][1] = -1+3*step;
		local_ref[26][2][2] = -1+2*step;

		local_ref[26][3][0] = -1+2*step;
		local_ref[26][3][1] = -1+3*step;
		local_ref[26][3][2] = -1+2*step;

		local_ref[26][4][0] = -1+2*step;
		local_ref[26][4][1] = -1+2*step;
		local_ref[26][4][2] = -1+3*step;

		local_ref[26][5][0] = -1+3*step;
		local_ref[26][5][1] = -1+2*step;
		local_ref[26][5][2] = -1+3*step;

		local_ref[26][6][0] = -1+3*step;
		local_ref[26][6][1] = -1+3*step;
		local_ref[26][6][2] = -1+3*step;

		local_ref[26][7][0] = -1+2*step;
		local_ref[26][7][1] = -1+3*step;
		local_ref[26][7][2] = -1+3*step;
	}

	int id_node[8];
	for(int ino=0;ino<8;ino++){
		id_node[ino] = elem->nodes[ino].id;
	}

	for(int i=0;i<27;i++){
		int conn_p[8];
		GtsPoint* point[8]={NULL};

		double cord_in_x[8],cord_in_y[8],cord_in_z[8];
		//add the nodes in the coord vector
		for (int ii = 0; ii < 8; ii++){
			cord_in_x[ii]=coords[3*id_node[ii]] ;
			cord_in_y[ii]=coords[3*id_node[ii]+1] ;
			cord_in_z[ii]=coords[3*id_node[ii]+2] ;
			//fprintf(mesh->fdbg,"coord in: %f, %f, %f, in the node: %d\n",cord_in_x[ii],cord_in_y[ii],cord_in_z[ii],elem->nodes[ii].id);
		}

		//add the new nodes in the
		for(int ii=0;ii<8;ii++){
			conn_p[ii] = 0;

			if((local_ref[i][ii][0]==1 || local_ref[i][ii][0]==-1) &&
					(local_ref[i][ii][1]==1 || local_ref[i][ii][1]==-1) &&
					(local_ref[i][ii][2]==1 || local_ref[i][ii][2]==-1)){

				if(local_ref[i][ii][0]==-1 && local_ref[i][ii][1]==-1 && local_ref[i][ii][2]==-1){
					conn_p[ii] = id_node[0];
				}else if(local_ref[i][ii][0]==1 && local_ref[i][ii][1]==-1 && local_ref[i][ii][2]==-1){
					conn_p[ii] = id_node[1];
				}else if(local_ref[i][ii][0]==1 && local_ref[i][ii][1]==1 && local_ref[i][ii][2]==-1){
					conn_p[ii] = id_node[2];
				}else if(local_ref[i][ii][0]==-1 && local_ref[i][ii][1]==1 && local_ref[i][ii][2]==-1){
					conn_p[ii] = id_node[3];
				}else if(local_ref[i][ii][0]==-1 && local_ref[i][ii][1]==-1 && local_ref[i][ii][2]==1){
					conn_p[ii] = id_node[4];
				}else if(local_ref[i][ii][0]==1 && local_ref[i][ii][1]==-1 && local_ref[i][ii][2]==1){
					conn_p[ii] = id_node[5];
				}else if(local_ref[i][ii][0]==1 && local_ref[i][ii][1]==1 && local_ref[i][ii][2]==1){
					conn_p[ii] = id_node[6];
				}else if(local_ref[i][ii][0]==-1 && local_ref[i][ii][1]==1 && local_ref[i][ii][2]==1){
					conn_p[ii] = id_node[7];
				}
				//fprintf(mesh->fdbg,"coord out: %f, %f, %f, in the node: %d\n",coords[3*conn_p[ii]],coords[3*conn_p[ii]+1],coords[3*conn_p[ii]+2],conn_p[ii]);
			}else{
				cord_in_ref[0] = local_ref[i][ii][0];
				cord_in_ref[1] = local_ref[i][ii][1];
				cord_in_ref[2] = local_ref[i][ii][2];

				point[ii] = LinearMapHex(cord_in_ref, cord_in_x,cord_in_y,cord_in_z);
				double var[3];
				var[0] = point[ii]->x;
				var[1] = point[ii]->y;
				var[2] = point[ii]->z;
				conn_p[ii] = AddPoint( mesh, hash_nodes, point[ii] , coords);
				//fprintf(mesh->fdbg,"coord out: %f, %f, %f, in the node: %d\n",var[0],var[1],var[2],conn_p[ii]);

			}
		}

		//add the new elements and make the connectivity
		if(i==0){
			octant_t *elem1 = (octant_t*) sc_array_index(&mesh->elements, elements_ids[iel]);

			elem1->nodes[0].id = conn_p[0];
			elem1->nodes[1].id = conn_p[1];
			elem1->nodes[2].id = conn_p[2];
			elem1->nodes[3].id = conn_p[3];

			elem1->nodes[4].id = conn_p[4];
			elem1->nodes[5].id = conn_p[5];
			elem1->nodes[6].id = conn_p[6];
			elem1->nodes[7].id = conn_p[7];

			CopyPropEl(mesh,id,elem1);

		}else{
			octant_t* elem2 = (octant_t*) sc_array_push(&mesh->elements);

			elem2->nodes[0].id = conn_p[0];
			elem2->nodes[1].id = conn_p[1];
			elem2->nodes[2].id = conn_p[2];
			elem2->nodes[3].id = conn_p[3];

			elem2->nodes[4].id = conn_p[4];
			elem2->nodes[5].id = conn_p[5];
			elem2->nodes[6].id = conn_p[6];
			elem2->nodes[7].id = conn_p[7];

			CopyPropEl(mesh,id,elem2);

		}
	}
	sc_array_reset(&toto);









	/*
	int id = elements_ids[iel];
	int conn_p[64];
	GtsPoint* point[64]={NULL};

	double cord_in_x[8],cord_in_y[8],cord_in_z[8];
	double cord_in_ref[3];
	cord_in_ref[0] = -1;
	cord_in_ref[1] = -1;
	cord_in_ref[2] = -1;

	//add the nodes in the coord vector
	for (int i = 0; i < 8; i++){
		cord_in_x[i]=coords[3*elem->nodes[i].id] ;
		cord_in_y[i]=coords[3*elem->nodes[i].id+1] ;
		cord_in_z[i]=coords[3*elem->nodes[i].id+2] ;
		//fprintf(mesh->fdbg,"coord in: %f, %f, %f, in the node: %d\n",cord_in_x[i],cord_in_y[i],cord_in_z[i],i);
	}

	for (int i = 0; i < 4; ++i) {
		cord_in_ref[1] = -1;
		for (int ii = 0; ii < 4; ++ii) {
			cord_in_ref[2] = -1;
			for (int iii = 0; iii < 4; ++iii) {

				//fprintf(mesh->fdbg,"coord ref: %f, %f, %f\n",cord_in_ref[0],cord_in_ref[1],cord_in_ref[2]);

				if((i==0 || i==3) && (ii==0 || ii==3) && (iii==0 || iii==3) ){
					if(i==0 && ii==0 && iii==0){
						conn_p[i*16+ii*4+iii] = elem->nodes[0].id;
					}else if(i==3 && ii==0 && iii==0){
						conn_p[i*16+ii*4+iii] = elem->nodes[1].id;
					}else if(i==0 && ii==3 && iii==0){
						conn_p[i*16+ii*4+iii] = elem->nodes[3].id;
					}else if(i==3 && ii==3 && iii==0){
						conn_p[i*16+ii*4+iii] = elem->nodes[2].id;
					}else if(i==0 && ii==0 && iii==3){
						conn_p[i*16+ii*4+iii] = elem->nodes[4].id;
					}else if(i==3 && ii==0 && iii==3){
						conn_p[i*16+ii*4+iii] = elem->nodes[5].id;
					}else if(i==0 && ii==3 && iii==3){
						conn_p[i*16+ii*4+iii] = elem->nodes[7].id;
					}else if(i==3 && ii==3 && iii==3){
						conn_p[i*16+ii*4+iii] = elem->nodes[6].id;
					}
				}else{
					point[i*16+ii*4+iii] = LinearMapHex(cord_in_ref, cord_in_x,cord_in_y,cord_in_z);
					double var[3];
					var[0] = point[i*16+ii*4+iii]->x;
					var[1] = point[i*16+ii*4+iii]->y;
					var[2] = point[i*16+ii*4+iii]->z;
					conn_p[i*16+ii*4+iii] = AddPoint( mesh, hash_nodes, point[i*16+ii*4+iii] , coords);
				}

				//fprintf(mesh->fdbg,"id do no: %d\n",conn_p[i*16+ii*4+iii]);
				//double xxx = coords[3*conn_p[i*16+ii*4+iii]];
				//double yyy = coords[3*conn_p[i*16+ii*4+iii]+1];
				//double zzz = coords[3*conn_p[i*16+ii*4+iii]+2];
				//fprintf(mesh->fdbg,"no vetor Coords x: %f, y:%f, z:%f\n", xxx, yyy, zzz );
				cord_in_ref[2] = cord_in_ref[2] + step;
			}
			cord_in_ref[1] = cord_in_ref[1] + step;
		}
		cord_in_ref[0] = cord_in_ref[0] + step;
	}

	//add the elements in the octants
	for (int i = 0; i < 3; ++i) {
		for (int ii = 0; ii < 3; ++ii) {
			for (int iii = 0; iii < 3; ++iii) {

				if(i==0 && ii==0 && iii==0 ){

					octant_t *elem1 = (octant_t*) sc_array_index(&mesh->elements, elements_ids[iel]);

					elem1->nodes[0].id = conn_p[i*16+ii*4+iii];
					elem1->nodes[1].id = conn_p[(i+1)*16+ii*4+iii];
					elem1->nodes[2].id = conn_p[(i+1)*16+(ii+1)*4+iii];
					elem1->nodes[3].id = conn_p[i*16+(ii+1)*4+iii];

					elem1->nodes[4].id = conn_p[i*16+ii*4+iii+1];
					elem1->nodes[5].id = conn_p[(i+1)*16+ii*4+iii+1];
					elem1->nodes[6].id = conn_p[(i+1)*16+(ii+1)*4+iii+1];
					elem1->nodes[7].id = conn_p[i*16+(ii+1)*4+iii+1];

					CopyPropEl(mesh,id,elem1);
				} else{

					octant_t* elem2 = (octant_t*) sc_array_push(&mesh->elements);

					elem2->nodes[0].id = conn_p[i*16+ii*4+iii];
					elem2->nodes[1].id = conn_p[(i+1)*16+ii*4+iii];
					elem2->nodes[2].id = conn_p[(i+1)*16+(ii+1)*4+iii];
					elem2->nodes[3].id = conn_p[i*16+(ii+1)*4+iii];

					elem2->nodes[4].id = conn_p[i*16+ii*4+iii+1];
					elem2->nodes[5].id = conn_p[(i+1)*16+ii*4+iii+1];
					elem2->nodes[6].id = conn_p[(i+1)*16+(ii+1)*4+iii+1];
					elem2->nodes[7].id = conn_p[i*16+(ii+1)*4+iii+1];

					CopyPropEl(mesh,id,elem2);
				}
			}
		}
	}
	sc_array_reset(&toto);
	 */
}




void ApplyOctreeTemplate(hexa_tree_t* mesh, std::vector<double>& coords) {

	bool clamped = true;
	sc_hash_array_t* hash_nodes = sc_hash_array_new(sizeof(node_t), edge_hash_fn, edge_equal_fn, &clamped);

	for(int ino = 0; ino<mesh->nodes.elem_count; ino++){
		size_t position;
		node_t *r;
		node_t key;
		octant_node_t* node = (octant_node_t*) sc_array_index (&mesh->nodes, ino);
		key.coord[0] = coords[3*node->id];
		key.coord[1] = coords[3*node->id+1];
		key.coord[2] = coords[3*node->id+2];
		key.node_id = node->id;

		r = (node_t*) sc_hash_array_insert_unique(hash_nodes, &key, &position);
	}

	std::vector<int> elements_ids;
	//for(int ioc = 0; ioc <mesh->oct.elem_count; ioc++){
	//	octree_t * oc = (octree_t*) sc_array_index (&mesh->oct, ioc);
	//	for(int iel = 0; iel<8; iel++) elements_ids.push_back(oc->id[iel]);
	//}
	for(int iel = 0; iel <mesh->elements.elem_count; iel++){
		octant_t * elem = (octant_t*) sc_array_index (&mesh->elements, iel);
		if(elem->pad!=0) elements_ids.push_back(elem->id);
	}

	double xs;
	double xe;
	double ys;
	double ye;

	for(int n = 0;n<mesh->nodes.elem_count;n++){
		octant_node_t* node = (octant_node_t*) sc_array_index (&mesh->nodes, n);
		if(node->y == mesh->y_start){ys=coords[3*node->id+1];}
		if(node->y == mesh->y_end)  {ye=coords[3*node->id+1];}
		if(node->x == mesh->x_start)  {xs=coords[3*node->id];}
		if(node->x == mesh->x_end)    {xe=coords[3*node->id];}
	}

	for (int iel = 0; iel < elements_ids.size(); ++iel) {

		double step = double(2)/double(3);

		octant_t *elem = (octant_t*) sc_array_index(&mesh->elements, elements_ids[iel]);
		int id = elements_ids[iel];

		//fprintf(mesh->fdbg,"Element: %d\n", elements_ids[iel]);

		//printf("Element: %d, pad: %d, temp: %d, level: %d\n",elements_ids[iel],elem->pad,elem->tem,elem->level);

		//template 1
		if(elem->tem==1){
			ApplyTemplate1(mesh, coords, elements_ids,iel, hash_nodes, step);
		}

		//template 2
		else if(elem->tem==2){
			ApplyTemplate2(mesh, coords, elements_ids,iel, hash_nodes, step);
		}

		//template 3
		else if(elem->tem==3){
			ApplyTemplate3(mesh, coords, elements_ids,iel, hash_nodes, step);
		}

		//template 4
		else if(elem->tem==4){
			ApplyTemplate4(mesh, coords, elements_ids,iel, hash_nodes, step);
		}

		//template 5
		else if(elem->tem==5){
			ApplyTemplate5(mesh, coords, elements_ids,iel, hash_nodes, step);
		}

		//template 6
		else if(elem->tem==6){
			ApplyTemplate6(mesh, coords, elements_ids,iel, hash_nodes, step);
		}

		//template 7
		else if(elem->tem==7){
			ApplyTemplate7(mesh, coords, elements_ids,iel, hash_nodes, step);
		}

		//template 8
		else if(elem->tem==8){
			ApplyTemplate8(mesh, coords, elements_ids,iel, hash_nodes, step);
		}

		//template 9
		else if(elem->tem==9){
			ApplyTemplate9(mesh, coords, elements_ids,iel, hash_nodes, step);
		}

		//template 10
		else if(elem->tem==10){
			ApplyTemplate10(mesh, coords, elements_ids,iel, hash_nodes, step);
		}

		//template 11
		else if(elem->tem==11){
			ApplyTemplate11(mesh, coords, elements_ids,iel, hash_nodes, step);
		}

	}

	//update the vectors
	mesh->local_n_elements = mesh->elements.elem_count;
	mesh->local_n_nodes = mesh->nodes.elem_count;
	MPI_Allreduce(&mesh->local_n_elements, &mesh->total_n_elements, 1, MPI_LONG, MPI_SUM, MPI_COMM_WORLD);
	MPI_Allreduce(&mesh->local_n_nodes, &mesh->total_n_nodes, 1, MPI_LONG, MPI_SUM, MPI_COMM_WORLD);

	if (mesh->mpi_rank == 0 && false) {
		printf("Total number of elements: %lld\n", mesh->local_n_elements);
		printf("Total number of nodes: %lld\n", mesh->local_n_nodes);
	}

	//redefine the connectivity data for write the vtk file
	mesh->part_nodes = NULL;
	mesh->part_nodes = (int*) malloc (mesh->local_n_nodes*sizeof(int));
	for(int i =0; i < mesh->local_n_nodes; i++)
		mesh->part_nodes[i] = mesh->mpi_rank;

	sc_hash_array_t* shared_nodes    = (sc_hash_array_t *)sc_hash_array_new(sizeof (shared_node_t), node_shared_hash_fn, node_shared_equal_fn, &clamped);
	//insert the shared nodes in the hash_array
	for(int i = 0; i < mesh->nodes.elem_count; i++)
	{
		octant_node_t* node = (octant_node_t*) sc_array_index (&mesh->nodes, i);
		if(node->y == mesh->y_start)
		{
			if(node->x == mesh->x_start)
				hexa_insert_shared_node(shared_nodes,node,mesh->neighbors[0]);
			if(node->x == mesh->x_end)
				hexa_insert_shared_node(shared_nodes,node,mesh->neighbors[2]);

			hexa_insert_shared_node(shared_nodes,node,mesh->neighbors[1]);
			continue;
		}

		if(node->y == mesh->y_end) {
			if(node->x == mesh->x_start)
				hexa_insert_shared_node(shared_nodes,node,mesh->neighbors[6]);
			if(node->x == mesh->x_end)
				hexa_insert_shared_node(shared_nodes,node,mesh->neighbors[8]);

			hexa_insert_shared_node(shared_nodes,node,mesh->neighbors[7]);
			continue;
		}

		if(node->x == mesh->x_start)
			hexa_insert_shared_node(shared_nodes,node,mesh->neighbors[3]);
		if(node->x == mesh->x_end)
			hexa_insert_shared_node(shared_nodes,node,mesh->neighbors[5]);
	}

	double fy = (ye-ys)/1000;
	double fx = (xe-xs)/1000;

	ys += fy;
	xs += fx;
	ye -= fy;
	xe -= fx;

	for (int iel = 0; iel < mesh->elements.elem_count; ++iel) {
		octant_t *elem = (octant_t*) sc_array_index(&mesh->elements, iel);
		elem->id = iel;
		if(false){
			fprintf(mesh->fdbg,"El:%d\n",elem->id);
			for(int i = 0; i < 8; i++){
				if(ys>=coords[3*elem->nodes[i].id+1]){
					if(xs>=coords[3*elem->nodes[i].id]){
						hexa_insert_shared_node(shared_nodes,&elem->nodes[i],mesh->neighbors[0]);
					}
					if(xe<=coords[3*elem->nodes[i].id]){
						hexa_insert_shared_node(shared_nodes,&elem->nodes[i],mesh->neighbors[2]);
					}
					hexa_insert_shared_node(shared_nodes,&elem->nodes[i],mesh->neighbors[1]);
					//continue;
				}
				if(ye<=coords[3*elem->nodes[i].id+1]){
					if(xs>=coords[3*elem->nodes[i].id]){
						hexa_insert_shared_node(shared_nodes,&elem->nodes[i],mesh->neighbors[6]);
					}
					if(xe<=coords[3*elem->nodes[i].id]){
						hexa_insert_shared_node(shared_nodes,&elem->nodes[i],mesh->neighbors[8]);
					}
					hexa_insert_shared_node(shared_nodes,&elem->nodes[i],mesh->neighbors[7]);
					//continue;
				}

				if(xs>=coords[3*elem->nodes[i].id]){
					hexa_insert_shared_node(shared_nodes,&elem->nodes[i],mesh->neighbors[3]);
				}
				if(xe<=coords[3*elem->nodes[i].id]){
					hexa_insert_shared_node(shared_nodes,&elem->nodes[i],mesh->neighbors[5]);
				}
			}
		}
	}


#ifdef HEXA_DEBUG_
	if(1){
		fprintf(mesh->fdbg, "Shared Nodes: \n");
		fprintf(mesh->fdbg, "Total: %d\n",shared_nodes->a.elem_count);
		for(int i = 0; i < shared_nodes->a.elem_count; ++i)
		{
			shared_node_t* sn = (shared_node_t*) sc_array_index(&shared_nodes->a,i);
			fprintf(mesh->fdbg, "(%d): %d %d %d\n", sn->id, sn->x, sn->y, sn->z);
			fprintf(mesh->fdbg, "     shared with processors: ");
			for(int j = 0; j < sn->listSz; j++){
				fprintf(mesh->fdbg, "%d ", sn->rankList[j]);
			}
			fprintf(mesh->fdbg, "\n");
		}
	}
#endif

	////////////////
	//extract the share nodes from shared_nodes
	sc_hash_array_rip (shared_nodes, &mesh->shared_nodes);
	sc_array_sort(&mesh->shared_nodes,node_comp);


	int local[3];
	size_t              position;

	local[0] = mesh->local_n_nodes    = mesh->nodes.elem_count;
	local[1] = mesh->local_n_elements = mesh->elements.elem_count;

	// node map
	int not_my_nodes    = 0;
	int my_own_nodes    = 0;
	mesh->global_id     = (int64_t*)malloc(sizeof(int64_t)*mesh->local_n_nodes);
	memset(mesh->global_id,-2,mesh->local_n_nodes*sizeof(int64_t));
	sc_hash_array_t* SendTo   = (sc_hash_array_t *) sc_hash_array_new(sizeof(message_t), processors_hash_fn, processors_equal_fn, &clamped);
	sc_hash_array_t* RecvFrom = (sc_hash_array_t *) sc_hash_array_new(sizeof(message_t), processors_hash_fn, processors_equal_fn, &clamped);

	for(int i = 0; i < mesh->shared_nodes.elem_count; ++i)
	{
		shared_node_t* sn = (shared_node_t*) sc_array_index(&mesh->shared_nodes,i);
		for(int j = 0; j < sn->listSz; j++)
		{
			if(sn->rankList[j] < mesh->mpi_rank) {
				message_t* m = (message_t*)sc_hash_array_insert_unique(SendTo,&sn->rankList[j],&position);
				if(m!=NULL)
				{
					m->rank  = sn->rankList[j];
					sc_array_init(&m->idxs, sizeof(uint32_t));
					uint32_t* p = (uint32_t*) sc_array_push(&m->idxs);
					*p = sn->id;
				} else
				{
					message_t* m = (message_t*)sc_array_index(&SendTo->a, position);
					uint32_t* p = (uint32_t*) sc_array_push(&m->idxs);
					*p = sn->id;
				}
				mesh->global_id[sn->id] = -3;
			}
			else if (sn->rankList[j] > mesh->mpi_rank)
			{
				message_t *m = (message_t*)sc_hash_array_insert_unique(RecvFrom,&sn->rankList[j],&position);
				if(m!=NULL)
				{
					m->rank  = sn->rankList[j];
					sc_array_init(&m->idxs, sizeof(uint32_t));
					uint32_t* p = (uint32_t*) sc_array_push(&m->idxs);
					*p = sn->id;
				} else
				{
					message_t* m = (message_t*)sc_array_index(&RecvFrom->a, position);
					uint32_t* p = (uint32_t*) sc_array_push(&m->idxs);
					*p = sn->id;
				}
				mesh->global_id[sn->id] = -1;
			}
		}

		if(mesh->global_id[sn->id] == -1) not_my_nodes++;
		if(mesh->global_id[sn->id] == -3) my_own_nodes++;

	}

	local[0] -= not_my_nodes;

	sc_hash_array_rip(RecvFrom, &mesh->comm_map.RecvFrom );
	sc_hash_array_rip(SendTo  , &mesh->comm_map.SendTo   );

#ifdef HEXA_DEBUG_
	if(0){
		//nodes
		fprintf(mesh->fdbg,"Nodes:\n");
		fprintf(mesh->fdbg, "Recv from %ld processors\n", mesh->comm_map.RecvFrom.elem_count);
		for(int i=0; i < mesh->comm_map.RecvFrom.elem_count; i++)
		{
			message_t* m = (message_t*) sc_array_index(&mesh->comm_map.RecvFrom, i);
			fprintf(mesh->fdbg, "  \n Recv %ld nodes from %d\n", m->idxs.elem_count, m->rank);
			for(int j =0; j < m->idxs.elem_count; j++)
			{
				int32_t *id = (int32_t*) sc_array_index(&m->idxs,j);
				fprintf(mesh->fdbg, "%d ", *id);
				if( (j+1) % 5 == 0 )fprintf(mesh->fdbg, "\n");
			}

		}
		fprintf(mesh->fdbg,"\n");
		fprintf(mesh->fdbg, "Send to %ld processors\n", mesh->comm_map.SendTo.elem_count);
		for(int i=0; i < mesh->comm_map.SendTo.elem_count; i++)
		{
			message_t* m = (message_t*) sc_array_index(&mesh->comm_map.SendTo, i);
			fprintf(mesh->fdbg, "\n Sending %ld nodes from %d\n", m->idxs.elem_count, m->rank);
			for(int j =0; j < m->idxs.elem_count; j++)
			{
				int32_t *id = (int32_t*) sc_array_index(&m->idxs,j);
				fprintf(mesh->fdbg, "%d ", *id);
				if( (j+1) % 5 == 0 )fprintf(mesh->fdbg, "\n");
			}

		}
	}
#endif

	//size for the nodes message
	mesh->comm_map.max_recvbuf_size = 0;
	for(int i = 0; i < mesh->comm_map.RecvFrom.elem_count; i++)
	{
		message_t* m = (message_t*) sc_array_index(&mesh->comm_map.RecvFrom, i);
		mesh->comm_map.max_recvbuf_size +=  m->idxs.elem_count;
	}

	mesh->comm_map.max_sendbuf_size = 0;
	for(int i = 0; i < mesh->comm_map.SendTo.elem_count; i++)
	{
		message_t* m = (message_t*) sc_array_index(&mesh->comm_map.SendTo, i);
		mesh->comm_map.max_sendbuf_size +=  m->idxs.elem_count;
	}

	mesh->comm_map.nrequests = mesh->comm_map.RecvFrom.elem_count +
			mesh->comm_map.SendTo.elem_count;

	int offset = 0;
	MPI_Scan(&local[0], &offset, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

	offset = offset-local[0];
	int count = offset;
	for(int i=0; i < mesh->local_n_nodes; ++i)
		if(mesh->global_id[i] != -1) mesh->global_id[i] = count++;

	int          n_requests;

	MPI_Request *requests;
	MPI_Status  *statuses;
	long long    *recvbuf;
	long long   *sendbuf;

	n_requests = mesh->comm_map.nrequests;
	recvbuf    = (long long*)malloc(mesh->comm_map.max_recvbuf_size*sizeof(long long));
	sendbuf    = (long long*)malloc(mesh->comm_map.max_sendbuf_size*sizeof(long long));

	requests = (MPI_Request*) malloc (n_requests*sizeof(MPI_Request));
	statuses = (MPI_Status*)  malloc (n_requests*sizeof(MPI_Status));
	int c = 0;

	offset = 0;

	// post all non-blocking receives
	for(int i = 0; i < mesh->comm_map.RecvFrom.elem_count; ++i) {
		message_t *m = (message_t*) sc_array_index(&mesh->comm_map.RecvFrom, i);
		MPI_Irecv(&recvbuf[offset], m->idxs.elem_count, MPI_LONG_LONG, m->rank,0,MPI_COMM_WORLD, &requests[c]);
		offset += m->idxs.elem_count;
		c++;
	}

	assert(offset == mesh->comm_map.max_recvbuf_size);

	offset = 0;
	for(int i = 0; i < mesh->comm_map.SendTo.elem_count; ++i) {
		message_t *m = (message_t*) sc_array_index(&mesh->comm_map.SendTo, i);
		for(int j = 0; j < m->idxs.elem_count; ++j)
		{
			int32_t *id = (int32_t*) sc_array_index(&m->idxs,j);
			sendbuf[offset+j] = (long long) mesh->global_id[*id];
		}
		MPI_Isend(&sendbuf[offset], m->idxs.elem_count, MPI_LONG_LONG, m->rank,0,MPI_COMM_WORLD, &requests[c]);
		offset += m->idxs.elem_count;
		c++;
	}
	assert(offset == mesh->comm_map.max_sendbuf_size);

	assert(c == n_requests);

	MPI_Waitall(n_requests,requests,statuses);

	offset = 0;
	for(int i = 0; i < mesh->comm_map.RecvFrom.elem_count; ++i) {
		message_t *m = (message_t*) sc_array_index(&mesh->comm_map.RecvFrom, i);
		for(int j = 0; j < m->idxs.elem_count; ++j)
		{
			int32_t *id = (int32_t*) sc_array_index(&m->idxs,j);
			mesh->global_id[*id] = recvbuf[offset+j];
		}
		offset += m->idxs.elem_count;
	}

	free(recvbuf);
	free(sendbuf);
	free(requests);
	free(statuses);

	mesh->part_nodes = (int*) malloc (mesh->local_n_nodes*sizeof(int));
	for(int i =0; i < mesh->local_n_nodes; i++)
		mesh->part_nodes[i] = mesh->mpi_rank;

	for(int i=0; i < mesh->comm_map.RecvFrom.elem_count; i++)
	{
		message_t* m = (message_t*) sc_array_index(&mesh->comm_map.RecvFrom, i);
		for(int j =0; j < m->idxs.elem_count; j++)
		{
			int32_t *id = (int32_t*) sc_array_index(&m->idxs,j);
			mesh->part_nodes[*id] = m->rank;
		}
	}

}


void PillowingInterface(hexa_tree_t* mesh, std::vector<double>& coords, std::vector<int>& nodes_b_mat){
	int elem_old = mesh->elements.elem_count;
	int nodes_old = mesh->nodes.elem_count;
	bool clamped = true;

	auto start = std::chrono::steady_clock::now( );

	//criando a hash de nos para evitar nos duplicados no pillowing
	sc_hash_array_t*   hash_nodes  = (sc_hash_array_t *)sc_hash_array_new(sizeof(node_t), edge_hash_fn, edge_equal_fn, &clamped);
	//hash nodes nodes_b_mat
	sc_hash_array_t*   hash_b_mat  = (sc_hash_array_t *)sc_hash_array_new(sizeof(node_t), edge_hash_fn, edge_equal_fn, &clamped);
	//vertex hash
	sc_hash_array_t*	  vertex_hash  = (sc_hash_array_t *)sc_hash_array_new(sizeof (octant_vertex_t), vertex_hash_id, vertex_equal_id, &clamped);
	//edge hash
	sc_hash_array_t* hash_edge_ref = sc_hash_array_new(sizeof (octant_edge_t), id_hash, id_equal, &clamped);

	printf("     Building hashs\n");
	BuildHash(mesh, coords, nodes_b_mat, hash_nodes, hash_b_mat, vertex_hash, hash_edge_ref);
	auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>( std::chrono::steady_clock::now( ) - start );

	std::cout << "Time in hash "<< elapsed.count() <<" millisecond(s)."<< std::endl;

	start = std::chrono::steady_clock::now( );
	printf("     Check octree template\n");
	CheckOctreeTemplate(mesh, hash_edge_ref);
	elapsed = std::chrono::duration_cast<std::chrono::milliseconds>( std::chrono::steady_clock::now( ) - start );

	std::cout << "Time check "<< elapsed.count() <<" millisecond(s)."<< std::endl;

	start = std::chrono::steady_clock::now( );
	printf("     Apply octree template\n");
	ApplyOctreeTemplate(mesh, coords);
	elapsed = std::chrono::duration_cast<std::chrono::milliseconds>( std::chrono::steady_clock::now( ) - start );

	std::cout << "Time apply "<< elapsed.count() <<" millisecond(s)."<< std::endl;

	FILE* bla;
	char treta[80];
	sprintf(treta,"edges_%04d_%04d.txt", mesh->mpi_size, mesh->mpi_rank);
	bla = fopen(treta,"w");
	for (int iel = 0; iel < mesh->elements.elem_count; ++iel) {
		octant_t *elem = (octant_t*) sc_array_index(&mesh->elements, iel);
		for (int edge = 0; edge < 12; ++edge) {
			fprintf(bla,"El:%d, edge:%d, coord[0]:%d, coord[1]:%d\n",iel,elem->edge[edge].id,elem->edge[edge].coord[0],elem->edge[edge].coord[1]);
		}
	}
	fclose(bla);

	//update the vectors
	mesh->local_n_elements = mesh->elements.elem_count;
	mesh->local_n_nodes = mesh->nodes.elem_count;
	MPI_Allreduce(&mesh->local_n_elements, &mesh->total_n_elements, 1, MPI_LONG, MPI_SUM, MPI_COMM_WORLD);
	MPI_Allreduce(&mesh->local_n_nodes, &mesh->total_n_nodes, 1, MPI_LONG, MPI_SUM, MPI_COMM_WORLD);
	printf("     %d nodes were created in the pillowing process\n",mesh->nodes.elem_count - nodes_old);
	printf("     %d elements were created in the pillowing process\n",mesh->elements.elem_count - elem_old);

	free(mesh->part_nodes);
	mesh->part_nodes = (int*) malloc (mesh->local_n_nodes*sizeof(int));
	for (int ino = 0; ino < mesh->local_n_nodes; ino++) {
		mesh->part_nodes[ino]=0;
	}

	for(int ino = 0; ino < nodes_b_mat.size(); ino++){
		mesh->part_nodes[nodes_b_mat[ino]] = 1;
	}
}
