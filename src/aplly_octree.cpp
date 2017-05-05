
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

#include "hexa.h"
#include "hilbert.h"
#include "refinement.h"




/*
typedef struct {
	bitmask_t coord[2];
	int node_id;
} node_in_edge_t;

unsigned edge_hash_fn(const void *v, const void *u) {
	const node_in_edge_t *q = (const node_in_edge_t*) v;
	uint32_t a, b, c;

	a = (uint32_t) q->coord[0];
	b = (uint32_t) q->coord[1];
	c = (uint32_t) 0;
	sc_hash_mix(a, b, c);
	sc_hash_final(a, b, c);
	return (unsigned) c;
}

int edge_equal_fn(const void *v, const void *u, const void *w) {
	const node_in_edge_t *e1 = (const node_in_edge_t*) v;
	const node_in_edge_t *e2 = (const node_in_edge_t*) u;

	return (unsigned) ((e1->coord[0] == e2->coord[0]) && (e2->coord[1] == e2->coord[1]));

}

int AddPointOnEdge(int* nodes, sc_hash_array_t* hash, int &npoints, GtsPoint *p, std::vector<double> &coords) {
	size_t position;
	node_in_edge_t *r;
	node_in_edge_t key;
	key.coord[0] = nodes[0];
	key.coord[1] = nodes[1];

	r = (node_in_edge_t*) sc_hash_array_insert_unique(hash, &key, &position);
	if (r != NULL) {
		r->coord[0] = key.coord[0];
		r->coord[1] = key.coord[1];
		r->node_id = npoints;
		npoints++;
		coords.push_back(p->x);
		coords.push_back(p->y);
		coords.push_back(p->z);
		return r->node_id;
	} else {
		r = (node_in_edge_t*) sc_array_index(&hash->a, position);
		return r->node_id;
	}
}

int EdgeCreateSurfMap[24][2] = {
		{0, 1}, // Edge 0
		{1, 2}, // Edge 1
		{2, 3}, //      2
		{3, 0},
		{0, 4},
		{1, 5},
		{2, 6},
		{3, 7},
		{4, 5},
		{5, 6},
		{6, 7},
		{7, 4}
};
 */

typedef struct {
	bitmask_t coord[3];
	int node_id;
} node_t;

unsigned edge_hash_fn(const void *v, const void *u) {
	const node_t *q = (const node_t*) v;
	uint32_t a, b, c;

	a = (uint32_t) q->coord[0];
	b = (uint32_t) q->coord[1];
	c = (uint32_t) 0;
	sc_hash_mix(a, b, c);
	a += (uint32_t) q->coord[3];
	sc_hash_final(a, b, c);
	return (unsigned) c;
}

/*
 * static unsigned p8est_mesh_indep_hash_fn (const void *v, const void *u)
{
  const p4est_locidx_t *indep = (p4est_locidx_t *) v;
  uint32_t            a, b, c;

  a = (uint32_t) indep[0];
  b = (uint32_t) indep[1];
  c = (uint32_t) indep[2];
  sc_hash_mix (a, b, c);
  a += (uint32_t) indep[3];
  sc_hash_final (a, b, c);

  return (unsigned) c;
}*/

int edge_equal_fn(const void *v, const void *u, const void *w) {
	const node_t *e1 = (const node_t*) v;
	const node_t *e2 = (const node_t*) u;

	return (unsigned) ((e1->coord[0] == e2->coord[0]) &&
			(e1->coord[1] == e2->coord[1]) &&
			(e1->coord[2] == e2->coord[2]));

}

int AddPoint(double* nodes, sc_hash_array_t* hash, int &npoints, GtsPoint *p, std::vector<double> &coords) {
	size_t position;
	node_t *r;
	node_t key;
	key.coord[0] = nodes[0];
	key.coord[1] = nodes[1];
	key.coord[2] = nodes[2];

	r = (node_t*) sc_hash_array_insert_unique(hash, &key, &position);
	if (r != NULL) {
		r->coord[0] = key.coord[0];
		r->coord[1] = key.coord[1];
		r->coord[2] = key.coord[2];
		r->node_id = npoints;
		npoints++;
		coords.push_back(p->x);
		coords.push_back(p->y);
		coords.push_back(p->z);
		return r->node_id;
	} else {
		r = (node_t*) sc_array_index(&hash->a, position);
		return r->node_id;
	}
}

GtsPoint* LinearMapHex(const double* cord_in_ref, const double* cord_in_x, const double* cord_in_y, const double* cord_in_z){

	double N[8];
	GtsPoint* point;
	double out[3];


	N[0] = (1-cord_in_ref[0])*(1-cord_in_ref[1])*(1-cord_in_ref[2])/double(8);
	N[1] = (1+cord_in_ref[0])*(1-cord_in_ref[1])*(1-cord_in_ref[2])/double(8);
	N[2] = (1+cord_in_ref[0])*(1+cord_in_ref[1])*(1-cord_in_ref[2])/double(8);
	N[3] = (1-cord_in_ref[0])*(1+cord_in_ref[1])*(1-cord_in_ref[2])/double(8);

	N[4] = (1-cord_in_ref[0])*(1-cord_in_ref[1])*(1+cord_in_ref[2])/double(8);
	N[5] = (1+cord_in_ref[0])*(1-cord_in_ref[1])*(1+cord_in_ref[2])/double(8);
	N[6] = (1+cord_in_ref[0])*(1+cord_in_ref[1])*(1+cord_in_ref[2])/double(8);
	N[7] = (1-cord_in_ref[0])*(1+cord_in_ref[1])*(1+cord_in_ref[2])/double(8);

	out[0] = 0;
	out[1] = 0;
	out[2] = 0;

	for(int i=0;i<8;i++){
		out[0] = N[i]*cord_in_x[i] + out[0] ;
		out[1] = N[i]*cord_in_y[i] + out[1];
		out[2] = N[i]*cord_in_z[i] + out[2];
	}

	point = gts_point_new(gts_point_class(),out[0],out[1],out[2]);

	return point;
}

vector<int> RotateHex(int* rot, int* sym){

	vector<int> order;

	for(int k = 0;k<8;k++){
		order.push_back(k);
	}

	if(rot[0]==-1){
		//order = {3,2,6,7,0,1,5,4};
		order.clear();
		order.push_back(3);
		order.push_back(2);
		order.push_back(6);
		order.push_back(7);

		order.push_back(0);
		order.push_back(1);
		order.push_back(5);
		order.push_back(4);

	}
	if(rot[0]==1){

		//order = {4,0,7,3,5,1,2,6};
		order.clear();

		order.push_back(4);
		order.push_back(5);
		order.push_back(1);
		order.push_back(0);

		order.push_back(7);
		order.push_back(6);
		order.push_back(2);
		order.push_back(3);

	}
	if(rot[1]==1){

		//order = {1,5,6,2,0,4,7,3};
		order.clear();

		order.push_back(1);
		order.push_back(5);
		order.push_back(6);
		order.push_back(2);

		order.push_back(0);
		order.push_back(4);
		order.push_back(7);
		order.push_back(3);


	}
	if(rot[1]==-1){

		//order = {4,5,1,0,7,6,2,3};
		order.clear();

		order.push_back(4);
		order.push_back(0);
		order.push_back(3);
		order.push_back(7);

		order.push_back(5);
		order.push_back(1);
		order.push_back(2);
		order.push_back(6);

	}
	if(rot[2]==1){

		//order = {3,0,1,2,7,4,5,6};
		order.clear();

		order.push_back(3);
		order.push_back(0);
		order.push_back(1);
		order.push_back(2);

		order.push_back(7);
		order.push_back(4);
		order.push_back(5);
		order.push_back(6);

	}
	if(rot[2]==-1){

		//order = {1,2,3,0,5,6,7,4};
		order.clear();

		order.push_back(1);
		order.push_back(2);
		order.push_back(3);
		order.push_back(0);

		order.push_back(5);
		order.push_back(6);
		order.push_back(7);
		order.push_back(4);

	}

	int aux[8];
	for(int i=0;i<8;i++){
		aux[i]=order[i];
	}

	if(sym[0]==1){
		order.clear();
		//order = {aux[1],aux[0],aux[3],aux[2],aux[5],aux[4],aux[7],aux[6]};
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
		order.clear();
		//order = {aux[3],aux[2],aux[6],aux[7],aux[0],aux[1],aux[5],aux[6]};
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
		order.clear();
		//order = {aux[4],aux[5],aux[6],aux[7],aux[0],aux[1],aux[2],aux[3]};
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

void CopyPropEl(octant_t *elem, octant_t *elem1){
	elem1->level = elem->level;
	elem1->ref = elem->ref+1;
	elem1->tem = elem->tem;
	elem1->pad = elem->pad;
	elem1->n_mat = elem->n_mat;
	elem1->pml_id = elem->pml_id;


	//elem1->x=elem->x;
	//elem1->y=elem->y;
	//elem1->z=elem->z;
}

void ApplyOctreeTemplate(hexa_tree_t* mesh, std::vector<double>& coords, std::vector<int>& elements_ids) {

	bool clamped = true;
	FILE * fdbg;
	char filename[80];
	sprintf(filename, "Nodeid_deb_%04d.txt", mesh->mpi_rank);
	fdbg = fopen(filename, "w");

	for (int iel = 0; iel < elements_ids.size(); ++iel) {

		double step = double(2)/double(3);

		octant_t *elem = (octant_t*) sc_array_index(&mesh->elements, elements_ids[iel]);

		sc_hash_array_t* hash_nodes = sc_hash_array_new(sizeof(node_t), edge_hash_fn, edge_equal_fn, &clamped);

		fprintf(fdbg,"Element: %d\n", elements_ids[iel]);

		if(elem->tem==1){

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

				/*
				id_node[0] = elem->nodes[0].id;
				id_node[1] = elem->nodes[1].id;
				id_node[2] = elem->nodes[2].id;
				id_node[3] = elem->nodes[3].id;

				id_node[4] = elem->nodes[4].id;
				id_node[5] = elem->nodes[5].id;
				id_node[6] = elem->nodes[6].id;
				id_node[7] = elem->nodes[7].id;
				 */
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

				id_node[0] = elem->nodes[1].id;
				id_node[1] = elem->nodes[2].id;
				id_node[2] = elem->nodes[3].id;
				id_node[3] = elem->nodes[0].id;

				id_node[4] = elem->nodes[5].id;
				id_node[5] = elem->nodes[6].id;
				id_node[6] = elem->nodes[7].id;
				id_node[7] = elem->nodes[4].id;

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

				id_node[0] = elem->nodes[2].id;
				id_node[1] = elem->nodes[3].id;
				id_node[2] = elem->nodes[0].id;
				id_node[3] = elem->nodes[1].id;

				id_node[4] = elem->nodes[6].id;
				id_node[5] = elem->nodes[7].id;
				id_node[6] = elem->nodes[4].id;
				id_node[7] = elem->nodes[5].id;

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

				id_node[0] = elem->nodes[3].id;
				id_node[1] = elem->nodes[0].id;
				id_node[2] = elem->nodes[1].id;
				id_node[3] = elem->nodes[2].id;

				id_node[4] = elem->nodes[7].id;
				id_node[5] = elem->nodes[4].id;
				id_node[6] = elem->nodes[5].id;
				id_node[7] = elem->nodes[6].id;

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

				id_node[0] = elem->nodes[0].id;
				id_node[1] = elem->nodes[4].id;
				id_node[2] = elem->nodes[7].id;
				id_node[3] = elem->nodes[3].id;

				id_node[4] = elem->nodes[1].id;
				id_node[5] = elem->nodes[2].id;
				id_node[6] = elem->nodes[5].id;
				id_node[7] = elem->nodes[6].id;

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

				id_node[0] = elem->nodes[1].id;
				id_node[1] = elem->nodes[5].id;
				id_node[2] = elem->nodes[6].id;
				id_node[3] = elem->nodes[2].id;

				id_node[4] = elem->nodes[0].id;
				id_node[5] = elem->nodes[4].id;
				id_node[6] = elem->nodes[7].id;
				id_node[7] = elem->nodes[3].id;

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

				id_node[0] = elem->nodes[2].id;
				id_node[1] = elem->nodes[6].id;
				id_node[2] = elem->nodes[7].id;
				id_node[3] = elem->nodes[3].id;

				id_node[4] = elem->nodes[1].id;
				id_node[5] = elem->nodes[5].id;
				id_node[6] = elem->nodes[4].id;
				id_node[7] = elem->nodes[0].id;

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

				id_node[0] = elem->nodes[3].id;
				id_node[1] = elem->nodes[7].id;
				id_node[2] = elem->nodes[6].id;
				id_node[3] = elem->nodes[2].id;

				id_node[4] = elem->nodes[0].id;
				id_node[5] = elem->nodes[4].id;
				id_node[6] = elem->nodes[5].id;
				id_node[7] = elem->nodes[1].id;

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

				id_node[0] = elem->nodes[4].id;
				id_node[1] = elem->nodes[5].id;
				id_node[2] = elem->nodes[6].id;
				id_node[3] = elem->nodes[7].id;

				id_node[4] = elem->nodes[0].id;
				id_node[5] = elem->nodes[1].id;
				id_node[6] = elem->nodes[2].id;
				id_node[7] = elem->nodes[3].id;

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

				id_node[0] = elem->nodes[5].id;
				id_node[1] = elem->nodes[6].id;
				id_node[2] = elem->nodes[7].id;
				id_node[3] = elem->nodes[4].id;

				id_node[4] = elem->nodes[1].id;
				id_node[5] = elem->nodes[2].id;
				id_node[6] = elem->nodes[3].id;
				id_node[7] = elem->nodes[0].id;

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

				id_node[0] = elem->nodes[6].id;
				id_node[1] = elem->nodes[7].id;
				id_node[2] = elem->nodes[4].id;
				id_node[3] = elem->nodes[5].id;

				id_node[4] = elem->nodes[2].id;
				id_node[5] = elem->nodes[3].id;
				id_node[6] = elem->nodes[0].id;
				id_node[7] = elem->nodes[1].id;

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

				id_node[0] = elem->nodes[7].id;
				id_node[1] = elem->nodes[4].id;
				id_node[2] = elem->nodes[5].id;
				id_node[3] = elem->nodes[6].id;

				id_node[4] = elem->nodes[3].id;
				id_node[5] = elem->nodes[0].id;
				id_node[6] = elem->nodes[1].id;
				id_node[7] = elem->nodes[2].id;

				ord = RotateHex(rot,sym);

				for(int node_id=0;node_id<8;node_id++){
					id_node[node_id] = elem->nodes[ord[node_id]].id;
				}

			}else{

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
					fprintf(fdbg,"coord in: %f, %f, %f, in the node: %d\n",cord_in_x[ii],cord_in_y[ii],cord_in_z[ii],elem->nodes[ii].id);
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
						fprintf(fdbg,"coord out: %f, %f, %f, in the node: %d\n",coords[3*conn_p[ii]],coords[3*conn_p[ii]+1],coords[3*conn_p[ii]+2],conn_p[ii]);
					}else{
						cord_in_ref[0] = local_ref[i][ii][0];
						cord_in_ref[1] = local_ref[i][ii][1];
						cord_in_ref[2] = local_ref[i][ii][2];

						point[ii] = LinearMapHex(cord_in_ref, cord_in_x,cord_in_y,cord_in_z);
						double var[3];
						var[0] = point[ii]->x;
						var[1] = point[ii]->y;
						var[2] = point[ii]->z;
						conn_p[ii] = AddPoint(var, hash_nodes, mesh->local_n_nodes, point[ii] , coords);
						fprintf(fdbg,"coord out: %f, %f, %f, in the node: %d\n",var[0],var[1],var[2],conn_p[ii]);

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

					CopyPropEl(elem,elem1);

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

					CopyPropEl(elem,elem2);

				}
			}
		}

		if(elem->tem==2){

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
			elem->pad=22;
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

				id_node[0] = elem->nodes[1].id;
				id_node[1] = elem->nodes[2].id;
				id_node[2] = elem->nodes[3].id;
				id_node[3] = elem->nodes[0].id;

				id_node[4] = elem->nodes[5].id;
				id_node[5] = elem->nodes[6].id;
				id_node[6] = elem->nodes[7].id;
				id_node[7] = elem->nodes[4].id;

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

				id_node[0] = elem->nodes[2].id;
				id_node[1] = elem->nodes[3].id;
				id_node[2] = elem->nodes[0].id;
				id_node[3] = elem->nodes[1].id;

				id_node[4] = elem->nodes[6].id;
				id_node[5] = elem->nodes[7].id;
				id_node[6] = elem->nodes[4].id;
				id_node[7] = elem->nodes[5].id;

				ord = RotateHex(rot,sym);

				for(int node_id=0;node_id<8;node_id++){
					id_node[node_id] = elem->nodes[ord[node_id]].id;
				}

			}else{

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
					fprintf(fdbg,"coord in: %f, %f, %f, in the node: %d\n",cord_in_x[ii],cord_in_y[ii],cord_in_z[ii],elem->nodes[ii].id);
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
						fprintf(fdbg,"coord out: %f, %f, %f, in the node: %d\n",coords[3*conn_p[ii]],coords[3*conn_p[ii]+1],coords[3*conn_p[ii]+2],conn_p[ii]);
					}else{
						cord_in_ref[0] = local_ref[i][ii][0];
						cord_in_ref[1] = local_ref[i][ii][1];
						cord_in_ref[2] = local_ref[i][ii][2];

						point[ii] = LinearMapHex(cord_in_ref, cord_in_x,cord_in_y,cord_in_z);
						double var[3];
						var[0] = point[ii]->x;
						var[1] = point[ii]->y;
						var[2] = point[ii]->z;
						conn_p[ii] = AddPoint(var, hash_nodes, mesh->local_n_nodes, point[ii] , coords);
						fprintf(fdbg,"coord out: %f, %f, %f, in the node: %d\n",var[0],var[1],var[2],conn_p[ii]);

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

					CopyPropEl(elem,elem1);

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

					CopyPropEl(elem,elem2);

				}
			}
		}


		if(elem->tem==5){

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

			if(elem->pad==50 || elem->pad==59 || elem->pad==60 || elem->pad==61 || elem->pad==62
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
			}else if(elem->pad==49 || elem->pad==55 || elem->pad==56 || elem->pad==57 || elem->pad==58
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
			}else if(elem->pad==51 || elem->pad==63 || elem->pad==64 || elem->pad==65 || elem->pad==66){
				//edge 0 4 5 8
				//code with a bug here, and the identifcation process missed something
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
					|| elem->pad==87 || elem->pad==88 || elem->pad==89 || elem->pad==90){
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
					|| elem->pad==91 || elem->pad==92 || elem->pad==93 || elem->pad==94){
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
					|| elem->pad==95 || elem->pad==96 || elem->pad==97 || elem->pad==98){
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
					fprintf(fdbg,"coord in: %f, %f, %f, in the node: %d\n",cord_in_x[ii],cord_in_y[ii],cord_in_z[ii],elem->nodes[ii].id);
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
						fprintf(fdbg,"coord out: %f, %f, %f, in the node: %d\n",coords[3*conn_p[ii]],coords[3*conn_p[ii]+1],coords[3*conn_p[ii]+2],conn_p[ii]);
					}else{
						cord_in_ref[0] = local_ref[i][ii][0];
						cord_in_ref[1] = local_ref[i][ii][1];
						cord_in_ref[2] = local_ref[i][ii][2];

						point[ii] = LinearMapHex(cord_in_ref, cord_in_x,cord_in_y,cord_in_z);
						double var[3];
						var[0] = point[ii]->x;
						var[1] = point[ii]->y;
						var[2] = point[ii]->z;
						conn_p[ii] = AddPoint(var, hash_nodes, mesh->local_n_nodes, point[ii] , coords);
						fprintf(fdbg,"coord out: %f, %f, %f, in the node: %d\n",var[0],var[1],var[2],conn_p[ii]);

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

					CopyPropEl(elem,elem1);

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

					CopyPropEl(elem,elem2);

				}

			}

		}



		//template 11
		if(elem->tem==11){

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
				fprintf(fdbg,"coord in: %f, %f, %f, in the node: %d\n",cord_in_x[i],cord_in_y[i],cord_in_z[i],i);
			}

			for (int i = 0; i < 4; ++i) {
				cord_in_ref[1] = -1;
				for (int ii = 0; ii < 4; ++ii) {
					cord_in_ref[2] = -1;
					for (int iii = 0; iii < 4; ++iii) {

						fprintf(fdbg,"coord ref: %f, %f, %f\n",cord_in_ref[0],cord_in_ref[1],cord_in_ref[2]);

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
							conn_p[i*16+ii*4+iii] = AddPoint(var, hash_nodes, mesh->local_n_nodes, point[i*16+ii*4+iii] , coords);
						}

						fprintf(fdbg,"id do no: %d\n",conn_p[i*16+ii*4+iii]);
						double xxx = coords[3*conn_p[i*16+ii*4+iii]];
						double yyy = coords[3*conn_p[i*16+ii*4+iii]+1];
						double zzz = coords[3*conn_p[i*16+ii*4+iii]+2];
						fprintf(fdbg,"no vetor Coords x: %f, y:%f, z:%f\n", xxx, yyy, zzz );
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

							//segfautl with the function
							//CopyPropEl(elem,elem1);
							//elem1->level = elem->level;
							//elem1->ref = elem->ref+1;
							//elem1->tem = elem->tem;
							//elem1->pad = elem->pad;
							//elem1->n_mat = elem->n_mat;
							//elem1->pml_id = elem->pml_id;
							elem1->pad=140;
							elem1->tem=11;
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

							//segfautl with the function
							//CopyPropEl(elem,elem2);
							//elem2->level = elem->level;
							//elem2->ref = elem->ref+1;
							//elem2->tem = elem->tem;
							//elem2->pad = elem->pad;
							//elem2->n_mat = elem->n_mat;
							//elem2->pml_id = elem->pml_id;
							elem2->pad=140;
							elem2->tem=11;
						}
					}
				}
			}
		}






	}
	fclose(fdbg);


	//update the vectors
	mesh->local_n_elements = mesh->elements.elem_count;
	MPI_Allreduce(&mesh->local_n_elements, &mesh->total_n_elements, 1, MPI_LONG, MPI_SUM, MPI_COMM_WORLD);
	MPI_Allreduce(&mesh->local_n_nodes, &mesh->total_n_nodes, 1, MPI_LONG, MPI_SUM, MPI_COMM_WORLD);
	mesh->nodes.elem_count =  mesh->local_n_nodes;

	if (mesh->mpi_rank == 0) {
		printf("Total number of elements: %lld\n", mesh->local_n_elements);
		printf("Total number of nodes: %lld\n", mesh->local_n_nodes);
	}

}


// Old Code... will be removed...
void SideEL(hexa_tree_t *mesh,int x, int y, int z ,std::vector<int>& element_ids_local){

	sc_array_t *elements = &mesh->elements;

	for(int iel = 0; iel < elements->elem_count; ++iel) {

		octant_t *elem = (octant_t*) sc_array_index(&mesh->elements, iel);

		if(z==0){
			if((elem->x==x-1)&&(elem->y==y-1)&&(elem->z==z+1)){
				element_ids_local.push_back(iel);
			}else if((elem->x==x-1)&&(elem->y==y)&&(elem->z==z+1)) {
				element_ids_local.push_back(iel);
			}else if((elem->x==x-1)&&(elem->y==y+1)&&(elem->z==z+1)) {
				element_ids_local.push_back(iel);
			}else if((elem->x==x)&&(elem->y==y-1)&&(elem->z==z+1)) {
				element_ids_local.push_back(iel);
			}else if((elem->x==x)&&(elem->y==y)&&(elem->z==z+1)) {
				element_ids_local.push_back(iel);
			}else if((elem->x==x)&&(elem->y==y+1)&&(elem->z==z+1)) {
				element_ids_local.push_back(iel);
			}else if((elem->x==x+1)&&(elem->y==y-1)&&(elem->z==z+1)) {
				element_ids_local.push_back(iel);
			}else if((elem->x==x+1)&&(elem->y==y)&&(elem->z==z+1)) {
				element_ids_local.push_back(iel);
			}else if((elem->x==x+1)&&(elem->y==y+1)&&(elem->z==z+1)) {
				element_ids_local.push_back(iel);
			}else if((elem->x==x-1)&&(elem->y==y-1)&&(elem->z==z)){
				element_ids_local.push_back(iel);
			}else if((elem->x==x-1)&&(elem->y==y)&&(elem->z==z)) {
				element_ids_local.push_back(iel);
			}else if((elem->x==x-1)&&(elem->y==y+1)&&(elem->z==z)) {
				element_ids_local.push_back(iel);
			}else if((elem->x==x)&&(elem->y==y-1)&&(elem->z==z)) {
				element_ids_local.push_back(iel);
			}else if((elem->x==x)&&(elem->y==y)&&(elem->z==z)) {
				element_ids_local.push_back(iel);
			}else if((elem->x==x)&&(elem->y==y+1)&&(elem->z==z)) {
				element_ids_local.push_back(iel);
			}else if((elem->x==x+1)&&(elem->y==y-1)&&(elem->z==z)) {
				element_ids_local.push_back(iel);
			}else if((elem->x==x+1)&&(elem->y==y)&&(elem->z==z)) {
				element_ids_local.push_back(iel);
			}else if((elem->x==x+1)&&(elem->y==y+1)&&(elem->z==z)) {
				element_ids_local.push_back(iel);
			}
		}else{
			if((elem->x==x-1)&&(elem->y==y-1)&&(elem->z==z+1)){
				element_ids_local.push_back(iel);
			}else if((elem->x==x-1)&&(elem->y==y)&&(elem->z==z+1)) {
				element_ids_local.push_back(iel);
			}else if((elem->x==x-1)&&(elem->y==y+1)&&(elem->z==z+1)) {
				element_ids_local.push_back(iel);
			}else if((elem->x==x)&&(elem->y==y-1)&&(elem->z==z+1)) {
				element_ids_local.push_back(iel);
			}else if((elem->x==x)&&(elem->y==y)&&(elem->z==z+1)) {
				element_ids_local.push_back(iel);
			}else if((elem->x==x)&&(elem->y==y+1)&&(elem->z==z+1)) {
				element_ids_local.push_back(iel);
			}else if((elem->x==x+1)&&(elem->y==y-1)&&(elem->z==z+1)) {
				element_ids_local.push_back(iel);
			}else if((elem->x==x+1)&&(elem->y==y)&&(elem->z==z+1)) {
				element_ids_local.push_back(iel);
			}else if((elem->x==x+1)&&(elem->y==y+1)&&(elem->z==z+1)) {
				element_ids_local.push_back(iel);
			}else if((elem->x==x-1)&&(elem->y==y-1)&&(elem->z==z)){
				element_ids_local.push_back(iel);
			}else if((elem->x==x-1)&&(elem->y==y)&&(elem->z==z)) {
				element_ids_local.push_back(iel);
			}else if((elem->x==x-1)&&(elem->y==y+1)&&(elem->z==z)) {
				element_ids_local.push_back(iel);
			}else if((elem->x==x)&&(elem->y==y-1)&&(elem->z==z)) {
				element_ids_local.push_back(iel);
			}else if((elem->x==x)&&(elem->y==y)&&(elem->z==z)) {
				element_ids_local.push_back(iel);
			}else if((elem->x==x)&&(elem->y==y+1)&&(elem->z==z)) {
				element_ids_local.push_back(iel);
			}else if((elem->x==x+1)&&(elem->y==y-1)&&(elem->z==z)) {
				element_ids_local.push_back(iel);
			}else if((elem->x==x+1)&&(elem->y==y)&&(elem->z==z)) {
				element_ids_local.push_back(iel);
			}else if((elem->x==x+1)&&(elem->y==y+1)&&(elem->z==z)) {
				element_ids_local.push_back(iel);
			}else if((elem->x==x-1)&&(elem->y==y-1)&&(elem->z==z-1)){
				element_ids_local.push_back(iel);
			}else if((elem->x==x-1)&&(elem->y==y)&&(elem->z==z-1)) {
				element_ids_local.push_back(iel);
			}else if((elem->x==x-1)&&(elem->y==y+1)&&(elem->z==z-1)) {
				element_ids_local.push_back(iel);
			}else if((elem->x==x)&&(elem->y==y-1)&&(elem->z==z-1)) {
				element_ids_local.push_back(iel);
			}else if((elem->x==x)&&(elem->y==y)&&(elem->z==z-1)) {
				element_ids_local.push_back(iel);
			}else if((elem->x==x)&&(elem->y==y+1)&&(elem->z==z-1)) {
				element_ids_local.push_back(iel);
			}else if((elem->x==x+1)&&(elem->y==y-1)&&(elem->z==z-1)) {
				element_ids_local.push_back(iel);
			}else if((elem->x==x+1)&&(elem->y==y)&&(elem->z==z-1)) {
				element_ids_local.push_back(iel);
			}else if((elem->x==x+1)&&(elem->y==y+1)&&(elem->z==z-1)) {
				element_ids_local.push_back(iel);
			}
		}
	}
}

void template2Rand(hexa_tree_t *mesh, std::vector<double>& coords,octant_t * elem, double v){

	double dx = coords[3*elem->nodes[1].id  ] - coords[3*elem->nodes[0].id  ] ;
	double dy = coords[3*elem->nodes[2].id+1] - coords[3*elem->nodes[0].id+1] ;
	double dz = coords[3*elem->nodes[0].id+2] - coords[3*elem->nodes[4].id+2] ;

	dx = fabs(dx);
	dy = fabs(dy);
	dz = fabs(dz);

	for(int i = 0;i<8;i++){
		int v1 = (-50 + rand() % 100)*0.02;
		int v2 = (-50 + rand() % 100)*0.02;
		int v3 = (-50 + rand() % 100)*0.02;
		int node_change = 3*elem->nodes[i].id;
		coords[node_change  ] = coords[node_change  ] + v*v1*dx;
		coords[node_change+1] = coords[node_change+1] + v*v2*dy;
		if(elem->z>0){
			coords[node_change+2] = coords[node_change+2] + v*v3*dz;
		}else if(elem->z==0 && i>=4){
			coords[node_change+2] = coords[node_change+2] + v*v3*dz;
		}
	}
}

void ChangeTemplate(hexa_tree_t* mesh, std::vector<double>& coords, std::vector<int>& elements_ids) {
	std::vector<int> element_ids_local;
	double coords_orgi[24];
	double v_1 = 1;
	double v_2 = 1;
	int iter_max = 200;
	int ref=0;

#if 0
	double tt=1.02;


	FILE * fdbg;
	fdbg = fopen("Change_template.txt", "w");

	int nao_sei = 0;
	int n_iter=0;
	int var_aux[27];
	int var_aux_1[27];
	int el_4 = 0;

	srand (time(NULL));

	sc_array_t *elements = &mesh->elements;

	for (int iel = 0; iel < elements_ids.size(); ++iel) {
		octant_t *elem = (octant_t*) sc_array_index(&mesh->elements, elements_ids[iel]);
		v_1 = 1;
		//v_2 = 1;

		if(elem->pad==4){

			//printf("El:%d\n",elements_ids[iel]);
			element_ids_local.clear();
			SideEL(mesh, elem->x, elem->y, elem->z ,element_ids_local);

			fprintf(fdbg,"El: %d Case type 4\n",elements_ids[iel]);
			fprintf(fdbg,"El:");
			for(int co = 0; co<element_ids_local.size();co++){
				fprintf(fdbg,"%d ",element_ids_local[co]);
				if(element_ids_local[co]==elements_ids[iel]){ref=co;}
			}
			fprintf(fdbg,"\n");

			fprintf(fdbg,"REF:%d \n",ref);

			fprintf(fdbg,"Pad Original: ");
			for(int c = 0; c < element_ids_local.size(); c++){
				octant_t *h = (octant_t*) sc_array_index(&mesh->elements, element_ids_local[c]);
				var_aux[c]=h->pad;
				fprintf(fdbg,"%d ",var_aux[c]);
			}
			fprintf(fdbg,"\n");

			for(int n_nodes = 0; n_nodes<8;n_nodes++){
				octant_t *elem = (octant_t*) sc_array_index(&mesh->elements, elements_ids[iel]);
				int node = 3*elem->nodes[n_nodes].id;
				coords_orgi[3*n_nodes  ] = coords[node  ];
				coords_orgi[3*n_nodes+1] = coords[node+1];
				coords_orgi[3*n_nodes+2] = coords[node+2];
			}

			n_iter = 0;
			bool flag1 = true;
			bool flag2 = true;

			while(n_iter<iter_max && flag1){

				// move one node... try to change to template 1
#if 0
				template2Rand(mesh, coords, elem, v_1);
				CheckTemplate(mesh, coords, element_ids_local,false);
				octant_t *elem = (octant_t*) sc_array_index(&mesh->elements, elements_ids[iel]);

				fprintf(fdbg,"Pad local:    ");
				for(int co = 0; co < element_ids_local.size(); co++){
					octant_t *h = (octant_t*) sc_array_index(&mesh->elements, element_ids_local[co]);
					var_aux_1[co]=h->pad;
					fprintf(fdbg,"%d ",h->pad);
				}
				fprintf(fdbg,"\n");

				if(elem->pad!=4 && elem->pad!=-10){
					flag2=false;
					for(int co = 0; co < element_ids_local.size(); co++){
						if(element_ids_local[co]!=elements_ids[iel]){
							if(var_aux[co]==var_aux_1[co]){
								flag2=false;
							}else{
								if(var_aux_1[co]==4){
									flag2=true;
									break;
								}else if(var_aux_1[co]==-10){
									flag2=true;
									break;
								}
							}
						}
					}
				}

				if(!flag2){
					flag1=false;
					flag2=false;
					printf("To livreeeee rand El:%d, n_iter:%d\n",elements_ids[iel],n_iter);
					break;
				}

				//if(elem->pad==4){
				//	CuttedEdges (mesh, coords, elements_ids, iel, &point[12], &edge_list[12]);

				GtsSegment * segments[12]={0};
				int Edge2GNode[12][2];
				GtsPoint * point[12] = {NULL};
				bool edge_list[12] = {false};

				for (int edge = 0; edge < 12; ++edge) {
					point[edge] = NULL;

					int node1 = elem->nodes[EdgeVerticesMap[edge][0]].id;
					int node2 = elem->nodes[EdgeVerticesMap[edge][1]].id;
					Edge2GNode[edge][0] = node1 <= node2 ? node1 : node2;
					Edge2GNode[edge][1] = node1 >= node2 ? node1 : node2;

					GtsVertex *v1 = gts_vertex_new(gts_vertex_class(), coords[node1 * 3], coords[node1 * 3 + 1], coords[node1 * 3 + 2]);
					GtsVertex *v2 = gts_vertex_new(gts_vertex_class(), coords[node2 * 3], coords[node2 * 3 + 1], coords[node2 * 3 + 2]);

					segments[edge] = gts_segment_new(gts_segment_class(), v1, v2);

					GtsBBox *sb = gts_bbox_segment(gts_bbox_class(), segments[edge]);
					GSList* list = gts_bb_tree_overlap(mesh->gdata.bbt, sb);

					edge_list[edge] = false;

					if (list == NULL) continue;
					while (list) {
						GtsBBox *b = GTS_BBOX(list->data);
						point[edge] = SegmentTriangleIntersection(segments[edge], GTS_TRIANGLE(b->bounded));
						if (point[edge]) {
							edge_list[edge] = true;
							break;
						}
						list = list->next;
					}
				}

				double d_c1=0;
				double d_c2=0;

				if(edge_list[0]&&edge_list[1]){
					int node_change = 3*elem->nodes[1].id;

					GtsPoint *p0 = gts_point_new(gts_point_class(),
							coords[node_change],
							coords[node_change + 1],
							coords[node_change + 2]);

					d_c1 = gts_point_distance(point[0], p0);
					d_c2 = gts_point_distance(point[1], p0);

					coords[node_change  ] = coords[node_change  ] - d_c1*tt;
					coords[node_change+1] = coords[node_change+1] + d_c2*tt;
					//move 1

				}else if(edge_list[1]&&edge_list[2]){
					int node_change = 3*elem->nodes[2].id;

					GtsPoint *p0 = gts_point_new(gts_point_class(),
							coords[node_change],
							coords[node_change + 1],
							coords[node_change + 2]);

					d_c1 = gts_point_distance(point[2], p0);
					d_c2 = gts_point_distance(point[1], p0);

					coords[node_change  ] = coords[node_change  ] - d_c1*tt;
					coords[node_change+1] = coords[node_change+1] - d_c2*tt;
					//move 2

				}else if(edge_list[2]&&edge_list[3]){
					int node_change = 3*elem->nodes[3].id;

					GtsPoint *p0 = gts_point_new(gts_point_class(),
							coords[node_change],
							coords[node_change + 1],
							coords[node_change + 2]);

					d_c1 = gts_point_distance(point[2], p0);
					d_c2 = gts_point_distance(point[3], p0);

					coords[node_change  ] = coords[node_change  ] + d_c1*tt;
					coords[node_change+1] = coords[node_change+1] - d_c2*tt;
					//move 3

				}else if(edge_list[3]&&edge_list[0]){
					int node_change = 3*elem->nodes[0].id;

					GtsPoint *p0 = gts_point_new(gts_point_class(),
							coords[node_change],
							coords[node_change + 1],
							coords[node_change + 2]);

					d_c1 = gts_point_distance(point[0], p0);
					d_c2 = gts_point_distance(point[3], p0);

					coords[node_change  ] = coords[node_change  ] + d_c1*tt;
					coords[node_change+1] = coords[node_change+1] + d_c2*tt;
					//move 0

				}else if(edge_list[8]&&edge_list[9]){
					int node_change = 3*elem->nodes[5].id;

					GtsPoint *p0 = gts_point_new(gts_point_class(),
							coords[node_change],
							coords[node_change + 1],
							coords[node_change + 2]);

					d_c1 = gts_point_distance(point[8], p0);
					d_c2 = gts_point_distance(point[9], p0);

					coords[node_change  ] = coords[node_change  ] - d_c1*tt;
					coords[node_change+1] = coords[node_change+1] + d_c2*tt;
					//move 5

				}else if(edge_list[9]&&edge_list[10]){
					int node_change = 3*elem->nodes[6].id;

					GtsPoint *p0 = gts_point_new(gts_point_class(),
							coords[node_change],
							coords[node_change + 1],
							coords[node_change + 2]);

					d_c1 = gts_point_distance(point[10], p0);
					d_c2 = gts_point_distance(point[9], p0);

					coords[node_change  ] = coords[node_change  ] - d_c1*tt;
					coords[node_change+1] = coords[node_change+1] - d_c2*tt;
					//move 6

				}else if(edge_list[10]&&edge_list[11]){
					int node_change = 3*elem->nodes[7].id;

					GtsPoint *p0 = gts_point_new(gts_point_class(),
							coords[node_change],
							coords[node_change + 1],
							coords[node_change + 2]);

					d_c1 = gts_point_distance(point[10], p0);
					d_c2 = gts_point_distance(point[11], p0);

					coords[node_change  ] = coords[node_change  ] + d_c1*tt;
					coords[node_change+1] = coords[node_change+1] - d_c2*tt;
					//move 7

				}else if(edge_list[11]&&edge_list[8]){
					int node_change = 3*elem->nodes[4].id;

					GtsPoint *p0 = gts_point_new(gts_point_class(),
							coords[node_change],
							coords[node_change + 1],
							coords[node_change + 2]);

					d_c1 = gts_point_distance(point[8], p0);
					d_c2 = gts_point_distance(point[11], p0);

					coords[node_change  ] = coords[node_change  ] + d_c1*tt;
					coords[node_change+1] = coords[node_change+1] + d_c2*tt;
					//move 4

				}else if( edge_list[0]&&edge_list[4] ){
					int node_change = 3*elem->nodes[0].id;

					GtsPoint *p0 = gts_point_new(gts_point_class(),
							coords[node_change],
							coords[node_change + 1],
							coords[node_change + 2]);

					d_c1 = gts_point_distance(point[0], p0);
					d_c2 = gts_point_distance(point[4], p0);

					coords[node_change  ] = coords[node_change  ] + d_c1*tt;
					coords[node_change+2] = coords[node_change+2] - d_c2*tt;
					//move 0

				}else if( edge_list[3]&&edge_list[4] ){
					int node_change = 3*elem->nodes[0].id;

					GtsPoint *p0 = gts_point_new(gts_point_class(),
							coords[node_change],
							coords[node_change + 1],
							coords[node_change + 2]);

					d_c1 = gts_point_distance(point[3], p0);
					d_c2 = gts_point_distance(point[4], p0);

					coords[node_change+1] = coords[node_change+1] + d_c1*tt;
					coords[node_change+2] = coords[node_change+2] - d_c2*tt;
					//move 0

				}else if(edge_list[0]&&edge_list[5]){
					int node_change = 3*elem->nodes[1].id;

					GtsPoint *p0 = gts_point_new(gts_point_class(),
							coords[node_change],
							coords[node_change + 1],
							coords[node_change + 2]);

					d_c1 = gts_point_distance(point[0], p0);
					d_c2 = gts_point_distance(point[5], p0);

					coords[node_change  ] = coords[node_change  ] - d_c1*tt;
					coords[node_change+2] = coords[node_change+2] - d_c2*tt;
					//move 1

				}else if(edge_list[1]&&edge_list[5]){
					int node_change = 3*elem->nodes[1].id;

					GtsPoint *p0 = gts_point_new(gts_point_class(),
							coords[node_change],
							coords[node_change + 1],
							coords[node_change + 2]);

					d_c1 = gts_point_distance(point[1], p0);
					d_c2 = gts_point_distance(point[5], p0);

					coords[node_change+1] = coords[node_change+1] + d_c1*tt;
					coords[node_change+2] = coords[node_change+2] - d_c2*tt;
					//move 1

				}else if(edge_list[1]&&edge_list[6]){
					int node_change = 3*elem->nodes[2].id;

					GtsPoint *p0 = gts_point_new(gts_point_class(),
							coords[node_change],
							coords[node_change + 1],
							coords[node_change + 2]);

					d_c1 = gts_point_distance(point[1], p0);
					d_c2 = gts_point_distance(point[6], p0);

					coords[node_change+1] = coords[node_change+1] - d_c1*tt;
					coords[node_change+2] = coords[node_change+2] - d_c2*tt;
					//move 2

				}else if(edge_list[2]&&edge_list[6]){
					int node_change = 3*elem->nodes[2].id;

					GtsPoint *p0 = gts_point_new(gts_point_class(),
							coords[node_change],
							coords[node_change + 1],
							coords[node_change + 2]);

					d_c1 = gts_point_distance(point[2], p0);
					d_c2 = gts_point_distance(point[6], p0);

					coords[node_change  ] = coords[node_change  ] - d_c1*tt;
					coords[node_change+2] = coords[node_change+2] - d_c2*tt;
					//move 2

				}else if(edge_list[2]&&edge_list[7]){
					int node_change = 3*elem->nodes[3].id;

					GtsPoint *p0 = gts_point_new(gts_point_class(),
							coords[node_change],
							coords[node_change + 1],
							coords[node_change + 2]);

					d_c1 = gts_point_distance(point[2], p0);
					d_c2 = gts_point_distance(point[7], p0);

					coords[node_change  ] = coords[node_change  ] + d_c1*tt;
					coords[node_change+2] = coords[node_change+2] - d_c2*tt;
					//move 3

				}else if(edge_list[3]&&edge_list[7]){
					int node_change = 3*elem->nodes[3].id;

					GtsPoint *p0 = gts_point_new(gts_point_class(),
							coords[node_change],
							coords[node_change + 1],
							coords[node_change + 2]);

					d_c1 = gts_point_distance(point[3], p0);
					d_c2 = gts_point_distance(point[7], p0);

					coords[node_change+1] = coords[node_change+1] - d_c1*tt;
					coords[node_change+2] = coords[node_change+2] - d_c2*tt;
					//move 3

				}else if(edge_list[8]&&edge_list[4]){
					int node_change = 3*elem->nodes[4].id;

					GtsPoint *p0 = gts_point_new(gts_point_class(),
							coords[node_change],
							coords[node_change + 1],
							coords[node_change + 2]);

					d_c1 = gts_point_distance(point[8], p0);
					d_c2 = gts_point_distance(point[4], p0);

					coords[node_change  ] = coords[node_change  ] + d_c1*tt;
					coords[node_change+2] = coords[node_change+2] + d_c2*tt;
					//move 4

				}else if(edge_list[11]&&edge_list[4]){
					int node_change = 3*elem->nodes[4].id;

					GtsPoint *p0 = gts_point_new(gts_point_class(),
							coords[node_change],
							coords[node_change + 1],
							coords[node_change + 2]);

					d_c1 = gts_point_distance(point[11], p0);
					d_c2 = gts_point_distance(point[4], p0);

					coords[node_change+1] = coords[node_change+1] + d_c1*tt;
					coords[node_change+2] = coords[node_change+2] + d_c2*tt;
					//move 4

				}else if(edge_list[8]&&edge_list[5]){
					int node_change = 3*elem->nodes[5].id;

					GtsPoint *p0 = gts_point_new(gts_point_class(),
							coords[node_change],
							coords[node_change + 1],
							coords[node_change + 2]);

					d_c1 = gts_point_distance(point[8], p0);
					d_c2 = gts_point_distance(point[5], p0);

					coords[node_change  ] = coords[node_change  ] - d_c1*tt;
					coords[node_change+2] = coords[node_change+2] + d_c2*tt;
					//move 5

				}else if(edge_list[9]&&edge_list[5]){
					int node_change = 3*elem->nodes[5].id;

					GtsPoint *p0 = gts_point_new(gts_point_class(),
							coords[node_change],
							coords[node_change + 1],
							coords[node_change + 2]);

					d_c1 = gts_point_distance(point[9], p0);
					d_c2 = gts_point_distance(point[5], p0);

					coords[node_change+1] = coords[node_change+1] + d_c1*tt;
					coords[node_change+2] = coords[node_change+2] + d_c2*tt;
					//move 5

				}else if(edge_list[9]&&edge_list[6]){
					int node_change = 3*elem->nodes[6].id;

					GtsPoint *p0 = gts_point_new(gts_point_class(),
							coords[node_change],
							coords[node_change + 1],
							coords[node_change + 2]);

					d_c1 = gts_point_distance(point[9], p0);
					d_c2 = gts_point_distance(point[6], p0);

					coords[node_change+1] = coords[node_change+1] - d_c1*tt;
					coords[node_change+2] = coords[node_change+2] + d_c2*tt;
					//move 6

				}else if(edge_list[10]&&edge_list[6]){
					int node_change = 3*elem->nodes[6].id;

					GtsPoint *p0 = gts_point_new(gts_point_class(),
							coords[node_change],
							coords[node_change + 1],
							coords[node_change + 2]);

					d_c1 = gts_point_distance(point[10], p0);
					d_c2 = gts_point_distance(point[6], p0);

					coords[node_change  ] = coords[node_change  ] - d_c1*tt;
					coords[node_change+2] = coords[node_change+2] + d_c2*tt;
					//move 6

				}else if(edge_list[10]&&edge_list[7]){
					int node_change = 3*elem->nodes[7].id;

					GtsPoint *p0 = gts_point_new(gts_point_class(),
							coords[node_change],
							coords[node_change + 1],
							coords[node_change + 2]);

					d_c1 = gts_point_distance(point[10], p0);
					d_c2 = gts_point_distance(point[7], p0);

					coords[node_change  ] = coords[node_change  ] + d_c1*tt;
					coords[node_change+2] = coords[node_change+2] + d_c2*tt;
					//move 7

				}else if(edge_list[11]&&edge_list[7]){
					int node_change = 3*elem->nodes[7].id;

					GtsPoint *p0 = gts_point_new(gts_point_class(),
							coords[node_change],
							coords[node_change + 1],
							coords[node_change + 2]);

					d_c1 = gts_point_distance(point[11], p0);
					d_c2 = gts_point_distance(point[7], p0);

					coords[node_change+1] = coords[node_change+1] - d_c1*tt;
					coords[node_change+2] = coords[node_change+2] + d_c2*tt;
					//move 7
				}

				CheckTemplate(mesh, coords, element_ids_local,false);
				elem = (octant_t*) sc_array_index(&mesh->elements, elements_ids[iel]);

				fprintf(fdbg,"Pad local:    ");
				for(int co = 0; co < element_ids_local.size(); co++){
					octant_t *h = (octant_t*) sc_array_index(&mesh->elements, element_ids_local[co]);
					var_aux_1[co]=h->pad;
					fprintf(fdbg,"%d ",h->pad);
				}
				fprintf(fdbg,"\n");


				if(elem->pad!=4 && elem->pad!=-10){
					flag2=false;
					for(int co = 0; co < element_ids_local.size(); co++){
						if(element_ids_local[co]!=elements_ids[iel]){
							if(var_aux[co]==var_aux_1[co]){
								flag2=false;
							}else{
								if(var_aux_1[co]==4){
									flag2=true;
									break;
								}else if(var_aux_1[co]==-10){
									flag2=true;
									break;
								}
							}
						}
					}
				}
				//}

				if(!flag2){
					flag1=false;
					flag2=false;
					printf("To livreeeee mov 1 El:%d, n_iter:%d\n",elements_ids[iel],n_iter);
				}

				if( flag2 ){
					for(int n_nodes = 0; n_nodes<8;n_nodes++){
						int node = 3*elem->nodes[n_nodes].id;
						coords[node]   = coords_orgi[3*n_nodes];
						coords[node+1] = coords_orgi[3*n_nodes+1];
						coords[node+2] = coords_orgi[3*n_nodes+2];
					}
				}
#endif

#if 0
				template2Rand(mesh, coords, elem, v_1);
				CheckTemplate(mesh, coords, element_ids_local,false);
				elem = (octant_t*) sc_array_index(&mesh->elements, elements_ids[iel]);

				fprintf(fdbg,"Pad local:    ");
				for(int co = 0; co < element_ids_local.size(); co++){
					octant_t *h = (octant_t*) sc_array_index(&mesh->elements, element_ids_local[co]);
					var_aux_1[co]=h->pad;
					fprintf(fdbg,"%d ",h->pad);
				}
				fprintf(fdbg,"\n");

				if(elem->pad!=4 && elem->pad!=-10){
					flag2=false;
					for(int co = 0; co < element_ids_local.size(); co++){
						if(element_ids_local[co]!=elements_ids[iel]){
							if(var_aux[co]==var_aux_1[co]){
								flag2=false;
							}else{
								if(var_aux_1[co]==4){
									flag2=true;
									break;
								}else if(var_aux_1[co]==-10){
									flag2=true;
									break;
								}
							}
						}
					}
				}

				if(!flag2){
					flag1=false;
					flag2=false;
					printf("To livreeeee rand El:%d, n_iter:%d\n",elements_ids[iel],n_iter);
					break;
				}

				//template4to1(mesh, coords, elements_ids, iel, &point[12], &edge_list[12],1.1);
				// move two nodes "or" one edge
				//if(elem->pad==4){
				//	CuttedEdges (mesh, coords, elements_ids, iel, &point[12], &edge_list[12]);

				GtsSegment * segments[12]={0};
				int Edge2GNode[12][2];
				GtsPoint * point[12] = {NULL};
				bool edge_list[12] = {false};

				for (int edge = 0; edge < 12; ++edge) {
					point[edge] = NULL;

					int node1 = elem->nodes[EdgeVerticesMap[edge][0]].id;
					int node2 = elem->nodes[EdgeVerticesMap[edge][1]].id;
					Edge2GNode[edge][0] = node1 <= node2 ? node1 : node2;
					Edge2GNode[edge][1] = node1 >= node2 ? node1 : node2;

					GtsVertex *v1 = gts_vertex_new(gts_vertex_class(), coords[node1 * 3], coords[node1 * 3 + 1], coords[node1 * 3 + 2]);
					GtsVertex *v2 = gts_vertex_new(gts_vertex_class(), coords[node2 * 3], coords[node2 * 3 + 1], coords[node2 * 3 + 2]);

					segments[edge] = 0;
					segments[edge] = gts_segment_new(gts_segment_class(), v1, v2);

					GtsBBox *sb = gts_bbox_segment(gts_bbox_class(), segments[edge]);
					GSList* list = gts_bb_tree_overlap(mesh->gdata.bbt, sb);

					edge_list[edge] = false;

					if (list == NULL) continue;
					while (list) {
						GtsBBox *b = GTS_BBOX(list->data);
						point[edge] = SegmentTriangleIntersection(segments[edge], GTS_TRIANGLE(b->bounded));
						if (point[edge]) {
							edge_list[edge] = true;
							break;
						}
						list = list->next;
					}
				}

				d_c1=0;
				d_c2=0;

				if(edge_list[0]&&edge_list[1]){
					int node_change = 3*elem->nodes[1].id;

					GtsPoint *p0 = gts_point_new(gts_point_class(),
							coords[node_change],
							coords[node_change + 1],
							coords[node_change + 2]);

					d_c1 = gts_point_distance(point[0], p0);
					d_c2 = gts_point_distance(point[1], p0);

					coords[node_change  ] = coords[node_change  ] - d_c1*tt;
					coords[node_change+1] = coords[node_change+1] + d_c2*tt;

					node_change = 3*elem->nodes[5].id;
					coords[node_change  ] = coords[node_change  ] - d_c1*tt;
					coords[node_change+1] = coords[node_change+1] + d_c2*tt;

					//move 1

				}else if(edge_list[1]&&edge_list[2]){
					int node_change = 3*elem->nodes[2].id;

					GtsPoint *p0 = gts_point_new(gts_point_class(),
							coords[node_change],
							coords[node_change + 1],
							coords[node_change + 2]);

					d_c1 = gts_point_distance(point[2], p0);
					d_c2 = gts_point_distance(point[1], p0);

					coords[node_change  ] = coords[node_change  ] - d_c1*tt;
					coords[node_change+1] = coords[node_change+1] - d_c2*tt;

					node_change = 3*elem->nodes[6].id;
					coords[node_change  ] = coords[node_change  ] - d_c1*tt;
					coords[node_change+1] = coords[node_change+1] - d_c2*tt;
					//move 2

				}else if(edge_list[2]&&edge_list[3]){
					int node_change = 3*elem->nodes[3].id;

					GtsPoint *p0 = gts_point_new(gts_point_class(),
							coords[node_change],
							coords[node_change + 1],
							coords[node_change + 2]);

					d_c1 = gts_point_distance(point[2], p0);
					d_c2 = gts_point_distance(point[3], p0);

					coords[node_change  ] = coords[node_change  ] + d_c1*tt;
					coords[node_change+1] = coords[node_change+1] - d_c2*tt;

					node_change = 3*elem->nodes[7].id;
					coords[node_change  ] = coords[node_change  ] + d_c1*tt;
					coords[node_change+1] = coords[node_change+1] - d_c2*tt;
					//move 3

				}else if(edge_list[3]&&edge_list[0]){
					int node_change = 3*elem->nodes[0].id;

					GtsPoint *p0 = gts_point_new(gts_point_class(),
							coords[node_change],
							coords[node_change + 1],
							coords[node_change + 2]);

					d_c1 = gts_point_distance(point[0], p0);
					d_c2 = gts_point_distance(point[3], p0);

					coords[node_change  ] = coords[node_change  ] + d_c1*tt;
					coords[node_change+1] = coords[node_change+1] + d_c2*tt;

					node_change = 3*elem->nodes[4].id;
					coords[node_change  ] = coords[node_change  ] + d_c1*tt;
					coords[node_change+1] = coords[node_change+1] + d_c2*tt;
					//move 0

				}else if(edge_list[8]&&edge_list[9]){
					int node_change = 3*elem->nodes[5].id;

					GtsPoint *p0 = gts_point_new(gts_point_class(),
							coords[node_change],
							coords[node_change + 1],
							coords[node_change + 2]);

					d_c1 = gts_point_distance(point[8], p0);
					d_c2 = gts_point_distance(point[9], p0);

					coords[node_change  ] = coords[node_change  ] - d_c1*tt;
					coords[node_change+1] = coords[node_change+1] + d_c2*tt;

					node_change = 3*elem->nodes[1].id;
					coords[node_change  ] = coords[node_change  ] - d_c1*tt;
					coords[node_change+1] = coords[node_change+1] + d_c2*tt;
					//move 5

				}else if(edge_list[9]&&edge_list[10]){
					int node_change = 3*elem->nodes[6].id;

					GtsPoint *p0 = gts_point_new(gts_point_class(),
							coords[node_change],
							coords[node_change + 1],
							coords[node_change + 2]);

					d_c1 = gts_point_distance(point[10], p0);
					d_c2 = gts_point_distance(point[9], p0);

					coords[node_change  ] = coords[node_change  ] - d_c1*tt;
					coords[node_change+1] = coords[node_change+1] - d_c2*tt;

					node_change = 3*elem->nodes[2].id;
					coords[node_change  ] = coords[node_change  ] - d_c1*tt;
					coords[node_change+1] = coords[node_change+1] - d_c2*tt;
					//move 6

				}else if(edge_list[10]&&edge_list[11]){
					int node_change = 3*elem->nodes[7].id;

					GtsPoint *p0 = gts_point_new(gts_point_class(),
							coords[node_change],
							coords[node_change + 1],
							coords[node_change + 2]);

					d_c1 = gts_point_distance(point[10], p0);
					d_c2 = gts_point_distance(point[11], p0);

					coords[node_change  ] = coords[node_change  ] + d_c1*tt;
					coords[node_change+1] = coords[node_change+1] - d_c2*tt;

					node_change = 3*elem->nodes[3].id;
					coords[node_change  ] = coords[node_change  ] + d_c1*tt;
					coords[node_change+1] = coords[node_change+1] - d_c2*tt;
					//move 7

				}else if(edge_list[11]&&edge_list[8]){
					int node_change = 3*elem->nodes[4].id;

					GtsPoint *p0 = gts_point_new(gts_point_class(),
							coords[node_change],
							coords[node_change + 1],
							coords[node_change + 2]);

					d_c1 = gts_point_distance(point[8], p0);
					d_c2 = gts_point_distance(point[11], p0);

					coords[node_change  ] = coords[node_change  ] + d_c1*tt;
					coords[node_change+1] = coords[node_change+1] + d_c2*tt;

					node_change = 3*elem->nodes[0].id;
					coords[node_change  ] = coords[node_change  ] + d_c1*tt;
					coords[node_change+1] = coords[node_change+1] + d_c2*tt;
					//move 4

				}else if( edge_list[0]&&edge_list[4] ){
					int node_change = 3*elem->nodes[0].id;

					GtsPoint *p0 = gts_point_new(gts_point_class(),
							coords[node_change],
							coords[node_change + 1],
							coords[node_change + 2]);

					d_c1 = gts_point_distance(point[0], p0);
					d_c2 = gts_point_distance(point[4], p0);

					coords[node_change  ] = coords[node_change  ] + d_c1*tt;
					coords[node_change+2] = coords[node_change+2] - d_c2*tt;

					node_change = 3*elem->nodes[3].id;
					coords[node_change  ] = coords[node_change  ] + d_c1*tt;
					coords[node_change+2] = coords[node_change+2] - d_c2*tt;
					//move 0

				}else if( edge_list[3]&&edge_list[4] ){
					int node_change = 3*elem->nodes[0].id;

					GtsPoint *p0 = gts_point_new(gts_point_class(),
							coords[node_change],
							coords[node_change + 1],
							coords[node_change + 2]);

					d_c1 = gts_point_distance(point[3], p0);
					d_c2 = gts_point_distance(point[4], p0);

					coords[node_change+1] = coords[node_change+1] + d_c1*tt;
					coords[node_change+2] = coords[node_change+2] - d_c2*tt;

					node_change = 3*elem->nodes[1].id;
					coords[node_change+1] = coords[node_change+1] + d_c1*tt;
					coords[node_change+2] = coords[node_change+2] - d_c2*tt;
					//move 0

				}else if(edge_list[0]&&edge_list[5]){
					int node_change = 3*elem->nodes[1].id;

					GtsPoint *p0 = gts_point_new(gts_point_class(),
							coords[node_change],
							coords[node_change + 1],
							coords[node_change + 2]);

					d_c1 = gts_point_distance(point[0], p0);
					d_c2 = gts_point_distance(point[5], p0);

					coords[node_change  ] = coords[node_change  ] - d_c1*tt;
					coords[node_change+2] = coords[node_change+2] - d_c2*tt;

					node_change = 3*elem->nodes[2].id;
					coords[node_change  ] = coords[node_change  ] - d_c1*tt;
					coords[node_change+2] = coords[node_change+2] - d_c2*tt;
					//move 1

				}else if(edge_list[1]&&edge_list[5]){
					int node_change = 3*elem->nodes[1].id;

					GtsPoint *p0 = gts_point_new(gts_point_class(),
							coords[node_change],
							coords[node_change + 1],
							coords[node_change + 2]);

					d_c1 = gts_point_distance(point[1], p0);
					d_c2 = gts_point_distance(point[5], p0);

					coords[node_change+1] = coords[node_change+1] + d_c1*tt;
					coords[node_change+2] = coords[node_change+2] - d_c2*tt;

					node_change = 3*elem->nodes[0].id;
					coords[node_change+1] = coords[node_change+1] + d_c1*tt;
					coords[node_change+2] = coords[node_change+2] - d_c2*tt;
					//move 1

				}else if(edge_list[1]&&edge_list[6]){
					int node_change = 3*elem->nodes[2].id;

					GtsPoint *p0 = gts_point_new(gts_point_class(),
							coords[node_change],
							coords[node_change + 1],
							coords[node_change + 2]);

					d_c1 = gts_point_distance(point[1], p0);
					d_c2 = gts_point_distance(point[6], p0);

					coords[node_change+1] = coords[node_change+1] - d_c1*tt;
					coords[node_change+2] = coords[node_change+2] - d_c2*tt;

					node_change = 3*elem->nodes[3].id;
					coords[node_change+1] = coords[node_change+1] - d_c1*tt;
					coords[node_change+2] = coords[node_change+2] - d_c2*tt;
					//move 2

				}else if(edge_list[2]&&edge_list[6]){
					int node_change = 3*elem->nodes[2].id;

					GtsPoint *p0 = gts_point_new(gts_point_class(),
							coords[node_change],
							coords[node_change + 1],
							coords[node_change + 2]);

					d_c1 = gts_point_distance(point[2], p0);
					d_c2 = gts_point_distance(point[6], p0);

					coords[node_change  ] = coords[node_change  ] - d_c1*tt;
					coords[node_change+2] = coords[node_change+2] - d_c2*tt;

					node_change = 3*elem->nodes[1].id;
					coords[node_change  ] = coords[node_change  ] - d_c1*tt;
					coords[node_change+2] = coords[node_change+2] - d_c2*tt;
					//move 2

				}else if(edge_list[2]&&edge_list[7]){
					int node_change = 3*elem->nodes[3].id;

					GtsPoint *p0 = gts_point_new(gts_point_class(),
							coords[node_change],
							coords[node_change + 1],
							coords[node_change + 2]);

					d_c1 = gts_point_distance(point[2], p0);
					d_c2 = gts_point_distance(point[7], p0);

					coords[node_change  ] = coords[node_change  ] + d_c1*tt;
					coords[node_change+2] = coords[node_change+2] - d_c2*tt;

					node_change = 3*elem->nodes[0].id;
					coords[node_change  ] = coords[node_change  ] + d_c1*tt;
					coords[node_change+2] = coords[node_change+2] - d_c2*tt;
					//move 3

				}else if(edge_list[3]&&edge_list[7]){
					int node_change = 3*elem->nodes[3].id;

					GtsPoint *p0 = gts_point_new(gts_point_class(),
							coords[node_change],
							coords[node_change + 1],
							coords[node_change + 2]);

					d_c1 = gts_point_distance(point[3], p0);
					d_c2 = gts_point_distance(point[7], p0);

					coords[node_change+1] = coords[node_change+1] - d_c1*tt;
					coords[node_change+2] = coords[node_change+2] - d_c2*tt;

					node_change = 3*elem->nodes[2].id;
					coords[node_change+1] = coords[node_change+1] - d_c1*tt;
					coords[node_change+2] = coords[node_change+2] - d_c2*tt;
					//move 3

				}else if(edge_list[8]&&edge_list[4]){
					int node_change = 3*elem->nodes[4].id;

					GtsPoint *p0 = gts_point_new(gts_point_class(),
							coords[node_change],
							coords[node_change + 1],
							coords[node_change + 2]);

					d_c1 = gts_point_distance(point[8], p0);
					d_c2 = gts_point_distance(point[4], p0);

					coords[node_change  ] = coords[node_change  ] + d_c1*tt;
					coords[node_change+2] = coords[node_change+2] + d_c2*tt;

					node_change = 3*elem->nodes[7].id;
					coords[node_change  ] = coords[node_change  ] + d_c1*tt;
					coords[node_change+2] = coords[node_change+2] + d_c2*tt;
					//move 4

				}else if(edge_list[11]&&edge_list[4]){
					int node_change = 3*elem->nodes[4].id;

					GtsPoint *p0 = gts_point_new(gts_point_class(),
							coords[node_change],
							coords[node_change + 1],
							coords[node_change + 2]);

					d_c1 = gts_point_distance(point[11], p0);
					d_c2 = gts_point_distance(point[4], p0);

					coords[node_change+1] = coords[node_change+1] + d_c1*tt;
					coords[node_change+2] = coords[node_change+2] + d_c2*tt;

					node_change = 3*elem->nodes[5].id;
					coords[node_change+1] = coords[node_change+1] + d_c1*tt;
					coords[node_change+2] = coords[node_change+2] + d_c2*tt;
					//move 4

				}else if(edge_list[8]&&edge_list[5]){
					int node_change = 3*elem->nodes[5].id;

					GtsPoint *p0 = gts_point_new(gts_point_class(),
							coords[node_change],
							coords[node_change + 1],
							coords[node_change + 2]);

					d_c1 = gts_point_distance(point[8], p0);
					d_c2 = gts_point_distance(point[5], p0);

					coords[node_change  ] = coords[node_change  ] - d_c1*tt;
					coords[node_change+2] = coords[node_change+2] + d_c2*tt;

					node_change = 3*elem->nodes[6].id;
					coords[node_change  ] = coords[node_change  ] - d_c1*tt;
					coords[node_change+2] = coords[node_change+2] + d_c2*tt;
					//move 5

				}else if(edge_list[9]&&edge_list[5]){
					int node_change = 3*elem->nodes[5].id;

					GtsPoint *p0 = gts_point_new(gts_point_class(),
							coords[node_change],
							coords[node_change + 1],
							coords[node_change + 2]);

					d_c1 = gts_point_distance(point[9], p0);
					d_c2 = gts_point_distance(point[5], p0);

					coords[node_change+1] = coords[node_change+1] + d_c1*tt;
					coords[node_change+2] = coords[node_change+2] + d_c2*tt;

					node_change = 3*elem->nodes[4].id;
					coords[node_change+1] = coords[node_change+1] + d_c1*tt;
					coords[node_change+2] = coords[node_change+2] + d_c2*tt;
					//move 5

				}else if(edge_list[9]&&edge_list[6]){
					int node_change = 3*elem->nodes[6].id;

					GtsPoint *p0 = gts_point_new(gts_point_class(),
							coords[node_change],
							coords[node_change + 1],
							coords[node_change + 2]);

					d_c1 = gts_point_distance(point[9], p0);
					d_c2 = gts_point_distance(point[6], p0);

					coords[node_change+1] = coords[node_change+1] - d_c1*tt;
					coords[node_change+2] = coords[node_change+2] + d_c2*tt;

					node_change = 3*elem->nodes[7].id;
					coords[node_change+1] = coords[node_change+1] - d_c1*tt;
					coords[node_change+2] = coords[node_change+2] + d_c2*tt;
					//move 6

				}else if(edge_list[10]&&edge_list[6]){
					int node_change = 3*elem->nodes[6].id;

					GtsPoint *p0 = gts_point_new(gts_point_class(),
							coords[node_change],
							coords[node_change + 1],
							coords[node_change + 2]);

					d_c1 = gts_point_distance(point[10], p0);
					d_c2 = gts_point_distance(point[6], p0);

					coords[node_change  ] = coords[node_change  ] - d_c1*tt;
					coords[node_change+2] = coords[node_change+2] + d_c2*tt;

					node_change = 3*elem->nodes[5].id;
					coords[node_change  ] = coords[node_change  ] - d_c1*tt;
					coords[node_change+2] = coords[node_change+2] + d_c2*tt;
					//move 6

				}else if(edge_list[10]&&edge_list[7]){
					int node_change = 3*elem->nodes[7].id;

					GtsPoint *p0 = gts_point_new(gts_point_class(),
							coords[node_change],
							coords[node_change + 1],
							coords[node_change + 2]);

					d_c1 = gts_point_distance(point[10], p0);
					d_c2 = gts_point_distance(point[7], p0);

					coords[node_change  ] = coords[node_change  ] + d_c1*tt;
					coords[node_change+2] = coords[node_change+2] + d_c2*tt;

					node_change = 3*elem->nodes[4].id;
					coords[node_change  ] = coords[node_change  ] + d_c1*tt;
					coords[node_change+2] = coords[node_change+2] + d_c2*tt;
					//move 7

				}else if(edge_list[11]&&edge_list[7]){
					int node_change = 3*elem->nodes[7].id;

					GtsPoint *p0 = gts_point_new(gts_point_class(),
							coords[node_change],
							coords[node_change + 1],
							coords[node_change + 2]);

					d_c1 = gts_point_distance(point[11], p0);
					d_c2 = gts_point_distance(point[7], p0);

					coords[node_change+1] = coords[node_change+1] - d_c1*tt;
					coords[node_change+2] = coords[node_change+2] + d_c2*tt;

					node_change = 3*elem->nodes[6].id;
					coords[node_change+1] = coords[node_change+1] - d_c1*tt;
					coords[node_change+2] = coords[node_change+2] + d_c2*tt;
					//move 7
				}

				CheckTemplate(mesh, coords, element_ids_local,false);
				elem = (octant_t*) sc_array_index(&mesh->elements, elements_ids[iel]);

				fprintf(fdbg,"Pad local:    ");
				for(int co = 0; co < element_ids_local.size(); co++){
					octant_t *h = (octant_t*) sc_array_index(&mesh->elements, element_ids_local[co]);
					var_aux_1[co]=h->pad;
					fprintf(fdbg,"%d ",h->pad);
				}
				fprintf(fdbg,"\n");


				if(elem->pad!=4 && elem->pad!=-10){
					flag2=false;
					for(int co = 0; co < element_ids_local.size(); co++){
						if(element_ids_local[co]!=elements_ids[iel]){
							if(var_aux[co]==var_aux_1[co]){
								flag2=false;
							}else{
								if(var_aux_1[co]==4){
									flag2=true;
									break;
								}else if(var_aux_1[co]==-10){
									flag2=true;
									break;
								}
							}
						}
					}
				}
				//}

				if(!flag2){
					flag1=false;
					flag2=false;
					printf("To livreeeee mov 2 El:%d, n_iter:%d\n",elements_ids[iel],n_iter);
				}

				if( flag2 ){
					for(int n_nodes = 0; n_nodes<8;n_nodes++){
						int node = 3*elem->nodes[n_nodes].id;
						coords[node]   = coords_orgi[3*n_nodes];
						coords[node+1] = coords_orgi[3*n_nodes+1];
						coords[node+2] = coords_orgi[3*n_nodes+2];
					}
				}
#endif

#if 0
				template2Rand(mesh, coords, elem, v_1);
				CheckTemplate(mesh, coords, element_ids_local,false);
				elem = (octant_t*) sc_array_index(&mesh->elements, elements_ids[iel]);

				fprintf(fdbg,"Pad local:    ");
				for(int co = 0; co < element_ids_local.size(); co++){
					octant_t *h = (octant_t*) sc_array_index(&mesh->elements, element_ids_local[co]);
					var_aux_1[co]=h->pad;
					fprintf(fdbg,"%d ",h->pad);
				}
				fprintf(fdbg,"\n");

				if(elem->pad!=4 && elem->pad!=-10){
					flag2=false;
					for(int co = 0; co < element_ids_local.size(); co++){
						if(element_ids_local[co]!=elements_ids[iel]){
							if(var_aux[co]==var_aux_1[co]){
								flag2=false;
							}else{
								if(var_aux_1[co]==4){
									flag2=true;
									break;
								}else if(var_aux_1[co]==-10){
									flag2=true;
									break;
								}
							}
						}
					}
				}

				if(!flag2){
					flag1=false;
					flag2=false;
					printf("To livreeeee rand El:%d, n_iter:%d\n",elements_ids[iel],n_iter);
					break;
				}

				//template4to2(mesh, coords, elements_ids, iel, &point[12], &edge_list[12],1.1);
				if(elem->pad==4){

					GtsSegment * segments[12]={0};
					int Edge2GNode[12][2];
					GtsPoint * point[12] = {NULL};
					bool edge_list[12] = {false};

					for (int edge = 0; edge < 12; ++edge) {
						point[edge] = NULL;

						int node1 = elem->nodes[EdgeVerticesMap[edge][0]].id;
						int node2 = elem->nodes[EdgeVerticesMap[edge][1]].id;
						Edge2GNode[edge][0] = node1 <= node2 ? node1 : node2;
						Edge2GNode[edge][1] = node1 >= node2 ? node1 : node2;

						GtsVertex *v1 = gts_vertex_new(gts_vertex_class(), coords[node1 * 3], coords[node1 * 3 + 1], coords[node1 * 3 + 2]);
						GtsVertex *v2 = gts_vertex_new(gts_vertex_class(), coords[node2 * 3], coords[node2 * 3 + 1], coords[node2 * 3 + 2]);

						segments[edge] = gts_segment_new(gts_segment_class(), v1, v2);

						GtsBBox *sb = gts_bbox_segment(gts_bbox_class(), segments[edge]);
						GSList* list = gts_bb_tree_overlap(mesh->gdata.bbt, sb);

						edge_list[edge] = false;

						if (list == NULL) continue;
						while (list) {
							GtsBBox *b = GTS_BBOX(list->data);
							point[edge] = SegmentTriangleIntersection(segments[edge], GTS_TRIANGLE(b->bounded));
							if (point[edge]) {
								edge_list[edge] = true;
								break;
							}
							list = list->next;
						}
					}

					double d_c1;
					double d_c2;

					if(edge_list[0]&&edge_list[1]){
						int node_change = 3*elem->nodes[1].id;

						GtsPoint *p0 = gts_point_new(gts_point_class(),
								coords[node_change],
								coords[node_change + 1],
								coords[node_change + 2]);

						d_c1 = gts_point_distance(point[0], p0);
						d_c2 = gts_point_distance(point[1], p0);

						if(d_c1<=d_c2){
							coords[node_change  ] = coords[node_change  ] - d_c1*tt;
						}else{
							coords[node_change+1] = coords[node_change+1] + d_c2*tt;
						}
						//move 1

					}else if(edge_list[1]&&edge_list[2]){
						int node_change = 3*elem->nodes[2].id;

						GtsPoint *p0 = gts_point_new(gts_point_class(),
								coords[node_change],
								coords[node_change + 1],
								coords[node_change + 2]);

						d_c1 = gts_point_distance(point[2], p0);
						d_c2 = gts_point_distance(point[1], p0);

						if(d_c1<=d_c2){
							coords[node_change  ] = coords[node_change  ] - d_c1*tt;
						}else{
							coords[node_change+1] = coords[node_change+1] - d_c2*tt;
						}
						//move 2

					}else if(edge_list[2]&&edge_list[3]){
						int node_change = 3*elem->nodes[3].id;

						GtsPoint *p0 = gts_point_new(gts_point_class(),
								coords[node_change],
								coords[node_change + 1],
								coords[node_change + 2]);

						d_c1 = gts_point_distance(point[2], p0);
						d_c2 = gts_point_distance(point[3], p0);

						if(d_c1<=d_c2){
							coords[node_change  ] = coords[node_change  ] + d_c1*tt;
						}else{
							coords[node_change+1] = coords[node_change+1] - d_c2*tt;
						}
						//move 3

					}else if(edge_list[3]&&edge_list[0]){
						int node_change = 3*elem->nodes[0].id;

						GtsPoint *p0 = gts_point_new(gts_point_class(),
								coords[node_change],
								coords[node_change + 1],
								coords[node_change + 2]);

						d_c1 = gts_point_distance(point[0], p0);
						d_c2 = gts_point_distance(point[3], p0);

						if(d_c1<=d_c2){
							coords[node_change  ] = coords[node_change  ] + d_c1*tt;
						}else{
							coords[node_change+1] = coords[node_change+1] + d_c2*tt;
						}
						//move 0

					}else if(edge_list[8]&&edge_list[9]){
						int node_change = 3*elem->nodes[5].id;

						GtsPoint *p0 = gts_point_new(gts_point_class(),
								coords[node_change],
								coords[node_change + 1],
								coords[node_change + 2]);

						d_c1 = gts_point_distance(point[8], p0);
						d_c2 = gts_point_distance(point[9], p0);

						if(d_c1<=d_c2){
							coords[node_change  ] = coords[node_change  ] - d_c1*tt;
						}else{
							coords[node_change+1] = coords[node_change+1] + d_c2*tt;
						}
						//move 5

					}else if(edge_list[9]&&edge_list[10]){
						int node_change = 3*elem->nodes[6].id;

						GtsPoint *p0 = gts_point_new(gts_point_class(),
								coords[node_change],
								coords[node_change + 1],
								coords[node_change + 2]);

						d_c1 = gts_point_distance(point[10], p0);
						d_c2 = gts_point_distance(point[9], p0);

						if(d_c1<=d_c2){
							coords[node_change  ] = coords[node_change  ] - d_c1*tt;
						}else{
							coords[node_change+1] = coords[node_change+1] - d_c2*tt;
						}
						//move 6

					}else if(edge_list[10]&&edge_list[11]){
						int node_change = 3*elem->nodes[7].id;

						GtsPoint *p0 = gts_point_new(gts_point_class(),
								coords[node_change],
								coords[node_change + 1],
								coords[node_change + 2]);

						d_c1 = gts_point_distance(point[10], p0);
						d_c2 = gts_point_distance(point[11], p0);

						if(d_c1<=d_c2){
							coords[node_change  ] = coords[node_change  ] + d_c1*tt;
						}else{
							coords[node_change+1] = coords[node_change+1] - d_c2*tt;
						}
						//move 7

					}else if(edge_list[11]&&edge_list[8]){
						int node_change = 3*elem->nodes[4].id;

						GtsPoint *p0 = gts_point_new(gts_point_class(),
								coords[node_change],
								coords[node_change + 1],
								coords[node_change + 2]);

						d_c1 = gts_point_distance(point[8], p0);
						d_c2 = gts_point_distance(point[11], p0);

						if(d_c1<=d_c2){
							coords[node_change  ] = coords[node_change  ] + d_c1*tt;
						}else{
							coords[node_change+1] = coords[node_change+1] + d_c2*tt;
						}
						//move 4

					}else if( edge_list[0]&&edge_list[4] ){
						int node_change = 3*elem->nodes[0].id;

						GtsPoint *p0 = gts_point_new(gts_point_class(),
								coords[node_change],
								coords[node_change + 1],
								coords[node_change + 2]);

						d_c1 = gts_point_distance(point[0], p0);
						d_c2 = gts_point_distance(point[4], p0);

						if(d_c1<=d_c2){
							coords[node_change  ] = coords[node_change  ] + d_c1*tt;
						}else{
							coords[node_change+2] = coords[node_change+2] - d_c2*tt;
						}
						//move 0

					}else if( edge_list[3]&&edge_list[4] ){
						int node_change = 3*elem->nodes[0].id;

						GtsPoint *p0 = gts_point_new(gts_point_class(),
								coords[node_change],
								coords[node_change + 1],
								coords[node_change + 2]);

						d_c1 = gts_point_distance(point[3], p0);
						d_c2 = gts_point_distance(point[4], p0);

						if(d_c1<=d_c2){
							coords[node_change+1] = coords[node_change+1] + d_c1*tt;
						}else{
							coords[node_change+2] = coords[node_change+2] - d_c2*tt;
						}
						//move 0

					}else if(edge_list[0]&&edge_list[5]){
						int node_change = 3*elem->nodes[1].id;

						GtsPoint *p0 = gts_point_new(gts_point_class(),
								coords[node_change],
								coords[node_change + 1],
								coords[node_change + 2]);

						d_c1 = gts_point_distance(point[0], p0);
						d_c2 = gts_point_distance(point[5], p0);

						if(d_c1<=d_c2){
							coords[node_change  ] = coords[node_change  ] - d_c1*tt;
						}else{
							coords[node_change+2] = coords[node_change+2] - d_c2*tt;
						}
						//move 1

					}else if(edge_list[1]&&edge_list[5]){
						int node_change = 3*elem->nodes[1].id;

						GtsPoint *p0 = gts_point_new(gts_point_class(),
								coords[node_change],
								coords[node_change + 1],
								coords[node_change + 2]);

						d_c1 = gts_point_distance(point[1], p0);
						d_c2 = gts_point_distance(point[5], p0);

						if(d_c1<=d_c2){
							coords[node_change+1] = coords[node_change+1] + d_c1*tt;
						}else{
							coords[node_change+2] = coords[node_change+2] - d_c2*tt;
						}
						//move 1

					}else if(edge_list[1]&&edge_list[6]){
						int node_change = 3*elem->nodes[2].id;

						GtsPoint *p0 = gts_point_new(gts_point_class(),
								coords[node_change],
								coords[node_change + 1],
								coords[node_change + 2]);

						d_c1 = gts_point_distance(point[1], p0);
						d_c2 = gts_point_distance(point[6], p0);

						if(d_c1<=d_c2){
							coords[node_change+1] = coords[node_change+1] - d_c1*tt;
						}else{
							coords[node_change+2] = coords[node_change+2] - d_c2*tt;
						}
						//move 2

					}else if(edge_list[2]&&edge_list[6]){
						int node_change = 3*elem->nodes[2].id;

						GtsPoint *p0 = gts_point_new(gts_point_class(),
								coords[node_change],
								coords[node_change + 1],
								coords[node_change + 2]);

						d_c1 = gts_point_distance(point[2], p0);
						d_c2 = gts_point_distance(point[6], p0);

						if(d_c1<=d_c2){
							coords[node_change  ] = coords[node_change  ] - d_c1*tt;
						}else{
							coords[node_change+2] = coords[node_change+2] - d_c2*tt;
						}
						//move 2

					}else if(edge_list[2]&&edge_list[7]){
						int node_change = 3*elem->nodes[3].id;

						GtsPoint *p0 = gts_point_new(gts_point_class(),
								coords[node_change],
								coords[node_change + 1],
								coords[node_change + 2]);

						d_c1 = gts_point_distance(point[2], p0);
						d_c2 = gts_point_distance(point[7], p0);

						if(d_c1<=d_c2){
							coords[node_change  ] = coords[node_change  ] + d_c1*tt;
						}else{
							coords[node_change+2] = coords[node_change+2] - d_c2*tt;
						}
						//move 3

					}else if(edge_list[3]&&edge_list[7]){
						int node_change = 3*elem->nodes[3].id;

						GtsPoint *p0 = gts_point_new(gts_point_class(),
								coords[node_change],
								coords[node_change + 1],
								coords[node_change + 2]);

						d_c1 = gts_point_distance(point[3], p0);
						d_c2 = gts_point_distance(point[7], p0);

						if(d_c1<=d_c2){
							coords[node_change+1] = coords[node_change+1] - d_c1*tt;
						}else{
							coords[node_change+2] = coords[node_change+2] - d_c2*tt;
						}
						//move 3

					}else if(edge_list[8]&&edge_list[4]){
						int node_change = 3*elem->nodes[4].id;

						GtsPoint *p0 = gts_point_new(gts_point_class(),
								coords[node_change],
								coords[node_change + 1],
								coords[node_change + 2]);

						d_c1 = gts_point_distance(point[8], p0);
						d_c2 = gts_point_distance(point[4], p0);

						if(d_c1<=d_c2){
							coords[node_change  ] = coords[node_change  ] + d_c1*tt;
						}else{
							coords[node_change+2] = coords[node_change+2] + d_c2*tt;
						}
						//move 4

					}else if(edge_list[11]&&edge_list[4]){
						int node_change = 3*elem->nodes[4].id;

						GtsPoint *p0 = gts_point_new(gts_point_class(),
								coords[node_change],
								coords[node_change + 1],
								coords[node_change + 2]);

						d_c1 = gts_point_distance(point[11], p0);
						d_c2 = gts_point_distance(point[4], p0);

						if(d_c1<=d_c2){
							coords[node_change+1] = coords[node_change+1] + d_c1*tt;
						}else{
							coords[node_change+2] = coords[node_change+2] + d_c2*tt;
						}
						//move 4

					}else if(edge_list[8]&&edge_list[5]){
						int node_change = 3*elem->nodes[5].id;

						GtsPoint *p0 = gts_point_new(gts_point_class(),
								coords[node_change],
								coords[node_change + 1],
								coords[node_change + 2]);

						d_c1 = gts_point_distance(point[8], p0);
						d_c2 = gts_point_distance(point[5], p0);

						if(d_c1<=d_c2){
							coords[node_change  ] = coords[node_change  ] - d_c1*tt;
						}else{
							coords[node_change+2] = coords[node_change+2] + d_c2*tt;
						}
						//move 5

					}else if(edge_list[9]&&edge_list[5]){
						int node_change = 3*elem->nodes[5].id;

						GtsPoint *p0 = gts_point_new(gts_point_class(),
								coords[node_change],
								coords[node_change + 1],
								coords[node_change + 2]);

						d_c1 = gts_point_distance(point[9], p0);
						d_c2 = gts_point_distance(point[5], p0);

						if(d_c1<=d_c2){
							coords[node_change+1] = coords[node_change+1] + d_c1*tt;
						}else{
							coords[node_change+2] = coords[node_change+2] + d_c2*tt;
						}
						//move 5

					}else if(edge_list[9]&&edge_list[6]){
						int node_change = 3*elem->nodes[6].id;

						GtsPoint *p0 = gts_point_new(gts_point_class(),
								coords[node_change],
								coords[node_change + 1],
								coords[node_change + 2]);

						d_c1 = gts_point_distance(point[9], p0);
						d_c2 = gts_point_distance(point[6], p0);

						if(d_c1<=d_c2){
							coords[node_change+1] = coords[node_change+1] - d_c1*tt;
						}else{
							coords[node_change+2] = coords[node_change+2] + d_c2*tt;
						}
						//move 6

					}else if(edge_list[10]&&edge_list[6]){
						int node_change = 3*elem->nodes[6].id;

						GtsPoint *p0 = gts_point_new(gts_point_class(),
								coords[node_change],
								coords[node_change + 1],
								coords[node_change + 2]);

						d_c1 = gts_point_distance(point[10], p0);
						d_c2 = gts_point_distance(point[6], p0);

						if(d_c1<=d_c2){
							coords[node_change  ] = coords[node_change  ] - d_c1*tt;
						}else{
							coords[node_change+2] = coords[node_change+2] + d_c2*tt;
						}
						//move 6

					}else if(edge_list[10]&&edge_list[7]){
						int node_change = 3*elem->nodes[7].id;

						GtsPoint *p0 = gts_point_new(gts_point_class(),
								coords[node_change],
								coords[node_change + 1],
								coords[node_change + 2]);

						d_c1 = gts_point_distance(point[10], p0);
						d_c2 = gts_point_distance(point[7], p0);

						if(d_c1<=d_c2){
							coords[node_change  ] = coords[node_change  ] + d_c1*tt;
						}else{
							coords[node_change+2] = coords[node_change+2] + d_c2*tt;
						}
						//move 7

					}else if(edge_list[11]&&edge_list[7]){
						int node_change = 3*elem->nodes[7].id;

						GtsPoint *p0 = gts_point_new(gts_point_class(),
								coords[node_change],
								coords[node_change + 1],
								coords[node_change + 2]);

						d_c1 = gts_point_distance(point[11], p0);
						d_c2 = gts_point_distance(point[7], p0);

						if(d_c1<=d_c2){
							coords[node_change+1] = coords[node_change+1] - d_c1*tt;
						}else{
							coords[node_change+2] = coords[node_change+2] + d_c2*tt;
						}
						//move 7
					}

					CheckTemplate(mesh, coords, element_ids_local,false);
					elem = (octant_t*) sc_array_index(&mesh->elements, elements_ids[iel]);

					fprintf(fdbg,"Pad local:    ");
					for(int co = 0; co < element_ids_local.size(); co++){
						octant_t *h = (octant_t*) sc_array_index(&mesh->elements, element_ids_local[co]);
						var_aux_1[co]=h->pad;
						fprintf(fdbg,"%d ",h->pad);
					}
					fprintf(fdbg,"\n");

					if(elem->pad!=4 && elem->pad!=-10){
						flag2=false;
						for(int co = 0; co < element_ids_local.size(); co++){
							if(element_ids_local[co]!=elements_ids[iel]){
								if(var_aux[co]==var_aux_1[co]){
									flag2=false;
								}else{
									if(var_aux_1[co]==4){
										flag2=true;
										break;
									}else if(var_aux_1[co]==-10){
										flag2=true;
										break;
									}
								}
							}
						}
					}
				}

				if(!flag2){
					flag1=false;
					flag2=false;
					printf("To livreeeee mov 3 El:%d, n_iter:%d\n",elements_ids[iel],n_iter);
				}

				if( flag2 ){
					for(int n_nodes = 0; n_nodes<8;n_nodes++){
						int node = 3*elem->nodes[n_nodes].id;
						coords[node]   = coords_orgi[3*n_nodes];
						coords[node+1] = coords_orgi[3*n_nodes+1];
						coords[node+2] = coords_orgi[3*n_nodes+2];
					}
				}
#endif

#if 0 //isso aqui não funciona...

				if (elem->z!=0){

					if(n_iter==0){
						printf("entrou aqui EL: %d\n",elements_ids[iel]);
					}

					elem = (octant_t*) sc_array_index(&mesh->elements, elements_ids[iel]);

					double dz = (coords[3*elem->nodes[4].id+ 3]-coords[3*elem->nodes[0].id+ 3]);
					dz = abs(dz);

					double nn = (50 - rand() % 100)*0.1;

					for(int co = 0; co<8; co++){
						int node_change = 3*elem->nodes[co].id;
						coords[node_change + 3] = coords[node_change + 3] + dz*nn;
					}

					CheckTemplate(mesh, coords, element_ids_local,false);
					elem = (octant_t*) sc_array_index(&mesh->elements, elements_ids[iel]);

					fprintf(fdbg,"Pad local:    ");
					for(int co = 0; co < element_ids_local.size(); co++){
						octant_t *h = (octant_t*) sc_array_index(&mesh->elements, element_ids_local[co]);
						var_aux_1[co]=h->pad;
						fprintf(fdbg,"%d ",h->pad);
					}
					fprintf(fdbg,"\n");

					if(elem->pad!=4 && elem->pad!=-10){
						flag2=false;
						for(int co = 0; co < element_ids_local.size(); co++){
							if(element_ids_local[co]!=elements_ids[iel]){
								if(var_aux[co]==var_aux_1[co]){
									flag2=false;
								}else{
									if(var_aux_1[co]==4){
										flag2=true;
										break;
									}else if(var_aux_1[co]==-10){
										flag2=true;
										break;
									}
								}
							}
						}
					}

					if(!flag2){
						flag1=false;
						flag2=false;
						printf("To livreeeee todos juntinhos El:%d, n_iter:%d\n",elements_ids[iel],n_iter);
						break;
					}
				}

				if( flag2 ){
					for(int n_nodes = 0; n_nodes<8;n_nodes++){
						int node = 3*elem->nodes[n_nodes].id;
						coords[node]   = coords_orgi[3*n_nodes];
						coords[node+1] = coords_orgi[3*n_nodes+1];
						coords[node+2] = coords_orgi[3*n_nodes+2];
					}
				}
#endif

#if 0

				template2Rand(mesh, coords, elem, v_1);
				CheckTemplate(mesh, coords, element_ids_local,false);
				elem = (octant_t*) sc_array_index(&mesh->elements, elements_ids[iel]);

				fprintf(fdbg,"Pad local:    ");
				for(int co = 0; co < element_ids_local.size(); co++){
					octant_t *h = (octant_t*) sc_array_index(&mesh->elements, element_ids_local[co]);
					var_aux_1[co]=h->pad;
					fprintf(fdbg,"%d ",h->pad);
				}
				fprintf(fdbg,"\n");

				if(elem->pad!=4 && elem->pad!=-10){
					flag2=false;
					for(int co = 0; co < element_ids_local.size(); co++){
						if(element_ids_local[co]!=elements_ids[iel]){
							if(var_aux[co]==var_aux_1[co]){
								flag2=false;
							}else{
								if(var_aux_1[co]==4){
									flag2=true;
									break;
								}else if(var_aux_1[co]==-10){
									flag2=true;
									break;
								}
							}
						}
					}
				}

				if(!flag2){
					flag1=false;
					flag2=false;
					printf("To livreeeee rand El:%d, n_iter:%d\n",elements_ids[iel],n_iter);
					break;
				}

				if(elem->pad==4){
					//	CuttedEdges (mesh, coords, elements_ids, iel, &point[12], &edge_list[12]);

					GtsSegment * segments[12]={0};
					int Edge2GNode[12][2];
					GtsPoint * point[12] = {NULL};
					bool edge_list[12] = {false};

					for (int edge = 0; edge < 12; ++edge) {
						point[edge] = NULL;

						int node1 = elem->nodes[EdgeVerticesMap[edge][0]].id;
						int node2 = elem->nodes[EdgeVerticesMap[edge][1]].id;
						Edge2GNode[edge][0] = node1 <= node2 ? node1 : node2;
						Edge2GNode[edge][1] = node1 >= node2 ? node1 : node2;

						GtsVertex *v1 = gts_vertex_new(gts_vertex_class(), coords[node1 * 3], coords[node1 * 3 + 1], coords[node1 * 3 + 2]);
						GtsVertex *v2 = gts_vertex_new(gts_vertex_class(), coords[node2 * 3], coords[node2 * 3 + 1], coords[node2 * 3 + 2]);

						segments[edge] = gts_segment_new(gts_segment_class(), v1, v2);

						GtsBBox *sb = gts_bbox_segment(gts_bbox_class(), segments[edge]);
						GSList* list = gts_bb_tree_overlap(mesh->gdata.bbt, sb);

						edge_list[edge] = false;

						if (list == NULL) continue;
						while (list) {
							GtsBBox *b = GTS_BBOX(list->data);
							point[edge] = SegmentTriangleIntersection(segments[edge], GTS_TRIANGLE(b->bounded));
							if (point[edge]) {
								edge_list[edge] = true;
								break;
							}
							list = list->next;
						}
					}

					double d_c1;
					double d_c2;

					double ref1 = 	coords[3*elem->nodes[0].id ] ;
					double ref2 = 	coords[3*elem->nodes[0].id + 1 ] ;
					double ref3 = coords[3*elem->nodes[0].id + 3 ] ;
					double dx = coords[3*elem->nodes[1].id ] - coords[3*elem->nodes[0].id ];
					double dy = coords[3*elem->nodes[2].id +1] - coords[3*elem->nodes[1].id+1 ];
					double dz = coords[3*elem->nodes[0].id +2] - coords[3*elem->nodes[4].id+2 ];

					if(edge_list[0]&&edge_list[1]){
						int node_change = 3*elem->nodes[1].id;

						GtsPoint *p0 = gts_point_new(gts_point_class(),
								coords[node_change],
								coords[node_change + 1],
								coords[node_change + 2]);

						d_c1 = gts_point_distance(point[0], p0);
						d_c2 = gts_point_distance(point[1], p0);

						double alpha = atan(d_c1/d_c2);

						for(int co = 0; co<8;co++){
							node_change = 3*elem->nodes[co].id;
							coords[node_change  ] = coords[node_change   ] - ref1 - dx*0.5;
							coords[node_change +1  ] =coords[node_change +1  ] -ref2-dy/2;
							coords[node_change  ] = coords[node_change  ]*cos(alpha)  - coords[node_change+1]*sin(alpha);
							coords[node_change +1  ] = coords[node_change]*sin(alpha) +coords[node_change +1  ]*cos(alpha);
							coords[node_change  ] = coords[node_change   ] + ref1 + dx*0.5;
							coords[node_change +1  ] =coords[node_change +1  ] + ref2 + dy/2;
						}
						//move 1

					}else if(edge_list[1]&&edge_list[2]){
						int node_change = 3*elem->nodes[2].id;

						GtsPoint *p0 = gts_point_new(gts_point_class(),
								coords[node_change],
								coords[node_change + 1],
								coords[node_change + 2]);

						d_c1 = gts_point_distance(point[2], p0);
						d_c2 = gts_point_distance(point[1], p0);

						double alpha = atan(d_c1/d_c2);

						for(int co = 0; co<8;co++){
							node_change = 3*elem->nodes[co].id;
							coords[node_change  ] = coords[node_change   ] - ref1 - dx*0.5;
							coords[node_change +1  ] =coords[node_change +1  ] -ref2-dy/2;
							coords[node_change  ] = coords[node_change  ]*cos(alpha)  - coords[node_change+1]*sin(alpha);
							coords[node_change +1  ] = coords[node_change  ]*sin(alpha) +coords[node_change +1  ]*cos(alpha);
							coords[node_change  ] = coords[node_change   ] + ref1 + dx*0.5;
							coords[node_change +1  ] =coords[node_change +1  ] + ref2 + dy/2;
						}
						//move 2

					}else if(edge_list[2]&&edge_list[3]){
						int node_change = 3*elem->nodes[3].id;

						GtsPoint *p0 = gts_point_new(gts_point_class(),
								coords[node_change],
								coords[node_change + 1],
								coords[node_change + 2]);

						d_c1 = gts_point_distance(point[2], p0);
						d_c2 = gts_point_distance(point[3], p0);

						double alpha = atan(d_c1/d_c2);

						for(int co = 0; co<8;co++){
							node_change = 3*elem->nodes[co].id;
							coords[node_change  ] = coords[node_change   ] - ref1 - dx*0.5;
							coords[node_change +1  ] =coords[node_change +1  ] -ref2-dy/2;
							coords[node_change  ] = coords[node_change  ]*cos(alpha)  - coords[node_change+1]*sin(alpha);
							coords[node_change +1  ] = coords[node_change  ]*sin(alpha) +coords[node_change +1  ]*cos(alpha);
							coords[node_change  ] = coords[node_change   ] + ref1 + dx*0.5;
							coords[node_change +1  ] =coords[node_change +1  ] + ref2 + dy/2;
						}
						//move 3

					}else if(edge_list[3]&&edge_list[0]){
						int node_change = 3*elem->nodes[0].id;

						GtsPoint *p0 = gts_point_new(gts_point_class(),
								coords[node_change],
								coords[node_change + 1],
								coords[node_change + 2]);

						d_c1 = gts_point_distance(point[0], p0);
						d_c2 = gts_point_distance(point[3], p0);

						double alpha = atan(d_c1/d_c2);

						for(int co = 0; co<8;co++){
							node_change = 3*elem->nodes[co].id;
							coords[node_change  ] = coords[node_change   ] - ref1 - dx*0.5;
							coords[node_change +1  ] =coords[node_change +1  ] -ref2-dy/2;
							coords[node_change  ] = coords[node_change  ]*cos(alpha)  - coords[node_change+1]*sin(alpha);
							coords[node_change +1  ] = coords[node_change  ]*sin(alpha) +coords[node_change +1  ]*cos(alpha);
							coords[node_change  ] = coords[node_change   ] + ref1 + dx*0.5;
							coords[node_change +1  ] =coords[node_change +1  ] + ref2 + dy/2;
						}
						//move 0

					}else if(edge_list[8]&&edge_list[9]){
						int node_change = 3*elem->nodes[5].id;

						GtsPoint *p0 = gts_point_new(gts_point_class(),
								coords[node_change],
								coords[node_change + 1],
								coords[node_change + 2]);

						d_c1 = gts_point_distance(point[8], p0);
						d_c2 = gts_point_distance(point[9], p0);

						double alpha = atan(d_c1/d_c2);

						for(int co = 0; co<8;co++){
							node_change = 3*elem->nodes[co].id;
							coords[node_change  ] = coords[node_change   ] - ref1 - dx*0.5;
							coords[node_change +1  ] =coords[node_change +1  ] -ref2-dy/2;
							coords[node_change  ] = coords[node_change  ]*cos(alpha)  - coords[node_change+1]*sin(alpha);
							coords[node_change +1  ] = coords[node_change]*sin(alpha) +coords[node_change +1  ]*cos(alpha);
							coords[node_change  ] = coords[node_change   ] + ref1 + dx*0.5;
							coords[node_change +1  ] =coords[node_change +1  ] + ref2 + dy/2;
						}
						//move 5

					}else if(edge_list[9]&&edge_list[10]){
						int node_change = 3*elem->nodes[6].id;

						GtsPoint *p0 = gts_point_new(gts_point_class(),
								coords[node_change],
								coords[node_change + 1],
								coords[node_change + 2]);

						d_c1 = gts_point_distance(point[10], p0);
						d_c2 = gts_point_distance(point[9], p0);

						double alpha = atan(d_c1/d_c2);

						for(int co = 0; co<8;co++){
							node_change = 3*elem->nodes[co].id;
							coords[node_change  ] = coords[node_change   ] - ref1 - dx*0.5;
							coords[node_change +1  ] =coords[node_change +1  ] -ref2-dy/2;
							coords[node_change  ] = coords[node_change  ]*cos(alpha)  - coords[node_change+1]*sin(alpha);
							coords[node_change +1  ] = coords[node_change ]*sin(alpha) +coords[node_change +1  ]*cos(alpha);
							coords[node_change  ] = coords[node_change   ] + ref1 + dx*0.5;
							coords[node_change +1  ] =coords[node_change +1  ] + ref2 + dy/2;
						}
						//move 6

					}else if(edge_list[10]&&edge_list[11]){
						int node_change = 3*elem->nodes[7].id;

						GtsPoint *p0 = gts_point_new(gts_point_class(),
								coords[node_change],
								coords[node_change + 1],
								coords[node_change + 2]);

						d_c1 = gts_point_distance(point[10], p0);
						d_c2 = gts_point_distance(point[11], p0);

						double alpha = atan(d_c1/d_c2);

						for(int co = 0; co<8;co++){
							node_change = 3*elem->nodes[co].id;
							coords[node_change  ] = coords[node_change   ] - ref1 - dx*0.5;
							coords[node_change +1  ] =coords[node_change +1  ] -ref2-dy/2;
							coords[node_change  ] = coords[node_change  ]*cos(alpha)  - coords[node_change+1]*sin(alpha);
							coords[node_change +1  ] = coords[node_change]*sin(alpha) +coords[node_change +1  ]*cos(alpha);
							coords[node_change  ] = coords[node_change   ] + ref1 + dx*0.5;
							coords[node_change +1  ] =coords[node_change +1  ] + ref2 + dy/2;
						}
						//move 7

					}else if(edge_list[11]&&edge_list[8]){
						int node_change = 3*elem->nodes[4].id;

						GtsPoint *p0 = gts_point_new(gts_point_class(),
								coords[node_change],
								coords[node_change + 1],
								coords[node_change + 2]);

						d_c1 = gts_point_distance(point[8], p0);
						d_c2 = gts_point_distance(point[11], p0);

						double alpha = atan(d_c1/d_c2);

						for(int co = 0; co<8;co++){
							node_change = 3*elem->nodes[co].id;
							coords[node_change  ] = coords[node_change   ] - ref1 - dx*0.5;
							coords[node_change +1  ] =coords[node_change +1  ] -ref2-dy/2;
							coords[node_change  ] = coords[node_change  ]*cos(alpha)  - coords[node_change+1]*sin(alpha);
							coords[node_change +1  ] = coords[node_change ]*sin(alpha) +coords[node_change +1  ]*cos(alpha);
							coords[node_change  ] = coords[node_change   ] + ref1 + dx*0.5;
							coords[node_change +1  ] =coords[node_change +1  ] + ref2 + dy/2;
						}
						//move 4

					}else if( edge_list[0]&&edge_list[4] ){
						int node_change = 3*elem->nodes[0].id;

						GtsPoint *p0 = gts_point_new(gts_point_class(),
								coords[node_change],
								coords[node_change + 1],
								coords[node_change + 2]);

						d_c1 = gts_point_distance(point[0], p0);
						d_c2 = gts_point_distance(point[4], p0);

						double alpha = atan(d_c1/d_c2);

						for(int co = 0; co<8;co++){
							node_change = 3*elem->nodes[co].id;
							coords[node_change  ] = coords[node_change   ] - ref1 - dx*0.5;
							coords[node_change +2  ] =coords[node_change +2  ] -ref3-dz/2;
							coords[node_change  ] = coords[node_change  ]*cos(alpha)  - coords[node_change+2]*sin(alpha);
							coords[node_change +2  ] = coords[node_change]*sin(alpha) +coords[node_change +2  ]*cos(alpha);
							coords[node_change  ] = coords[node_change   ] + ref1 + dx*0.5;
							coords[node_change +2  ] =coords[node_change +1  ] + ref3 + dz/2;
						}
						//move 0

					}else if( edge_list[3]&&edge_list[4] ){
						int node_change = 3*elem->nodes[0].id;

						GtsPoint *p0 = gts_point_new(gts_point_class(),
								coords[node_change],
								coords[node_change + 1],
								coords[node_change + 2]);

						d_c1 = gts_point_distance(point[3], p0);
						d_c2 = gts_point_distance(point[4], p0);

						double alpha = atan(d_c1/d_c2);

						for(int co = 0; co<8;co++){
							node_change = 3*elem->nodes[co].id;
							coords[node_change + 1] = coords[node_change + 1] - ref2 - dy*0.5;
							coords[node_change + 2] = coords[node_change + 2] - ref3 - dz*0.5;
							coords[node_change + 1] = coords[node_change + 1]*cos(alpha) - coords[node_change + 2]*sin(alpha);
							coords[node_change + 2] = coords[node_change + 1]*sin(alpha) + coords[node_change + 2]*cos(alpha);
							coords[node_change + 1] = coords[node_change + 1] + ref2 + dy*0.5;
							coords[node_change + 2] = coords[node_change + 2] + ref3 + dz*0.5;
						}
						//move 0

					}else if(edge_list[0]&&edge_list[5]){
						int node_change = 3*elem->nodes[1].id;

						GtsPoint *p0 = gts_point_new(gts_point_class(),
								coords[node_change],
								coords[node_change + 1],
								coords[node_change + 2]);

						d_c1 = gts_point_distance(point[0], p0);
						d_c2 = gts_point_distance(point[5], p0);

						double alpha = atan(d_c1/d_c2);

						for(int co = 0; co<8;co++){
							node_change = 3*elem->nodes[co].id;
							coords[node_change    ] = coords[node_change    ] - ref1 - dx*0.5;
							coords[node_change + 2] = coords[node_change + 2] - ref3 - dz*0.5;
							coords[node_change    ] = coords[node_change    ]*cos(alpha) - coords[node_change + 2]*sin(alpha);
							coords[node_change + 2] = coords[node_change    ]*sin(alpha) + coords[node_change + 2]*cos(alpha);
							coords[node_change    ] = coords[node_change    ] + ref1 + dx*0.5;
							coords[node_change + 2] = coords[node_change + 2] + ref3 + dz*0.5;
						}
						//move 1

					}else if(edge_list[1]&&edge_list[5]){
						int node_change = 3*elem->nodes[1].id;

						GtsPoint *p0 = gts_point_new(gts_point_class(),
								coords[node_change],
								coords[node_change + 1],
								coords[node_change + 2]);

						d_c1 = gts_point_distance(point[1], p0);
						d_c2 = gts_point_distance(point[5], p0);

						double alpha = atan(d_c1/d_c2);

						for(int co = 0; co<8;co++){
							node_change = 3*elem->nodes[co].id;
							coords[node_change + 1] = coords[node_change + 1] - ref2 - dy*0.5;
							coords[node_change + 2] = coords[node_change + 2] - ref3 - dz*0.5;
							coords[node_change + 1] = coords[node_change + 1]*cos(alpha) - coords[node_change + 2]*sin(alpha);
							coords[node_change + 2] = coords[node_change + 1]*sin(alpha) + coords[node_change + 2]*cos(alpha);
							coords[node_change + 1] = coords[node_change + 1] + ref2 + dy*0.5;
							coords[node_change + 2] = coords[node_change + 2] + ref3 + dz*0.5;
						}
						//move 1

					}else if(edge_list[1]&&edge_list[6]){
						int node_change = 3*elem->nodes[2].id;

						GtsPoint *p0 = gts_point_new(gts_point_class(),
								coords[node_change],
								coords[node_change + 1],
								coords[node_change + 2]);

						d_c1 = gts_point_distance(point[1], p0);
						d_c2 = gts_point_distance(point[6], p0);

						double alpha = atan(d_c1/d_c2);

						for(int co = 0; co<8;co++){
							node_change = 3*elem->nodes[co].id;
							coords[node_change + 1] = coords[node_change + 1] - ref2 - dy*0.5;
							coords[node_change + 2] = coords[node_change + 2] - ref3 - dz*0.5;
							coords[node_change + 1] = coords[node_change + 1]*cos(alpha) - coords[node_change + 2]*sin(alpha);
							coords[node_change + 2] = coords[node_change + 1]*sin(alpha) + coords[node_change + 2]*cos(alpha);
							coords[node_change + 1] = coords[node_change + 1] + ref2 + dy*0.5;
							coords[node_change + 2] = coords[node_change + 2] + ref3 + dz*0.5;
						}
						//move 2

					}else if(edge_list[2]&&edge_list[6]){
						int node_change = 3*elem->nodes[2].id;

						GtsPoint *p0 = gts_point_new(gts_point_class(),
								coords[node_change],
								coords[node_change + 1],
								coords[node_change + 2]);

						d_c1 = gts_point_distance(point[2], p0);
						d_c2 = gts_point_distance(point[6], p0);

						double alpha = atan(d_c1/d_c2);

						for(int co = 0; co<8;co++){
							node_change = 3*elem->nodes[co].id;
							coords[node_change    ] = coords[node_change    ] - ref1 - dx*0.5;
							coords[node_change + 2] = coords[node_change + 2] - ref3 - dz*0.5;
							coords[node_change    ] = coords[node_change    ]*cos(alpha) - coords[node_change + 2]*sin(alpha);
							coords[node_change + 2] = coords[node_change    ]*sin(alpha) + coords[node_change + 2]*cos(alpha);
							coords[node_change    ] = coords[node_change    ] + ref1 + dx*0.5;
							coords[node_change + 2] = coords[node_change + 2] + ref3 + dz*0.5;
						}
						//move 2

					}else if(edge_list[2]&&edge_list[7]){
						int node_change = 3*elem->nodes[3].id;

						GtsPoint *p0 = gts_point_new(gts_point_class(),
								coords[node_change],
								coords[node_change + 1],
								coords[node_change + 2]);

						d_c1 = gts_point_distance(point[2], p0);
						d_c2 = gts_point_distance(point[7], p0);

						double alpha = atan(d_c1/d_c2);

						for(int co = 0; co<8;co++){
							node_change = 3*elem->nodes[co].id;
							coords[node_change    ] = coords[node_change    ] - ref1 - dx*0.5;
							coords[node_change + 2] = coords[node_change + 2] - ref3 - dz*0.5;
							coords[node_change    ] = coords[node_change    ]*cos(alpha) - coords[node_change + 2]*sin(alpha);
							coords[node_change + 2] = coords[node_change    ]*sin(alpha) + coords[node_change + 2]*cos(alpha);
							coords[node_change    ] = coords[node_change    ] + ref1 + dx*0.5;
							coords[node_change + 2] = coords[node_change + 2] + ref3 + dz*0.5;
						}
						//move 3

					}else if(edge_list[3]&&edge_list[7]){
						int node_change = 3*elem->nodes[3].id;

						GtsPoint *p0 = gts_point_new(gts_point_class(),
								coords[node_change],
								coords[node_change + 1],
								coords[node_change + 2]);

						d_c1 = gts_point_distance(point[3], p0);
						d_c2 = gts_point_distance(point[7], p0);

						double alpha = atan(d_c1/d_c2);

						for(int co = 0; co<8;co++){
							node_change = 3*elem->nodes[co].id;
							coords[node_change + 1] = coords[node_change + 1] - ref2 - dy*0.5;
							coords[node_change + 2] = coords[node_change + 2] - ref3 - dz*0.5;
							coords[node_change + 1] = coords[node_change + 1]*cos(alpha) - coords[node_change + 2]*sin(alpha);
							coords[node_change + 2] = coords[node_change + 1]*sin(alpha) + coords[node_change + 2]*cos(alpha);
							coords[node_change + 1] = coords[node_change + 1] + ref2 + dy*0.5;
							coords[node_change + 2] = coords[node_change + 2] + ref3 + dz*0.5;
						}
						//move 3

					}else if(edge_list[8]&&edge_list[4]){
						int node_change = 3*elem->nodes[4].id;

						GtsPoint *p0 = gts_point_new(gts_point_class(),
								coords[node_change],
								coords[node_change + 1],
								coords[node_change + 2]);

						d_c1 = gts_point_distance(point[8], p0);
						d_c2 = gts_point_distance(point[4], p0);

						double alpha = atan(d_c1/d_c2);

						for(int co = 0; co<8;co++){
							node_change = 3*elem->nodes[co].id;
							coords[node_change    ] = coords[node_change    ] - ref1 - dx*0.5;
							coords[node_change + 2] = coords[node_change + 2] - ref3 - dz*0.5;
							coords[node_change    ] = coords[node_change    ]*cos(alpha) - coords[node_change + 2]*sin(alpha);
							coords[node_change + 2] = coords[node_change    ]*sin(alpha) + coords[node_change + 2]*cos(alpha);
							coords[node_change    ] = coords[node_change    ] + ref1 + dx*0.5;
							coords[node_change + 2] = coords[node_change + 2] + ref3 + dz*0.5;
						}
						//move 4

					}else if(edge_list[11]&&edge_list[4]){
						int node_change = 3*elem->nodes[4].id;

						GtsPoint *p0 = gts_point_new(gts_point_class(),
								coords[node_change],
								coords[node_change + 1],
								coords[node_change + 2]);

						d_c1 = gts_point_distance(point[11], p0);
						d_c2 = gts_point_distance(point[4], p0);

						double alpha = atan(d_c1/d_c2);

						for(int co = 0; co<8;co++){
							node_change = 3*elem->nodes[co].id;
							coords[node_change + 1] = coords[node_change + 1] - ref2 - dy*0.5;
							coords[node_change + 2] = coords[node_change + 2] - ref3 - dz*0.5;
							coords[node_change + 1] = coords[node_change + 1]*cos(alpha) - coords[node_change + 2]*sin(alpha);
							coords[node_change + 2] = coords[node_change + 1]*sin(alpha) + coords[node_change + 2]*cos(alpha);
							coords[node_change + 1] = coords[node_change + 1] + ref2 + dy*0.5;
							coords[node_change + 2] = coords[node_change + 2] + ref3 + dz*0.5;
						}
						//move 4

					}else if(edge_list[8]&&edge_list[5]){
						int node_change = 3*elem->nodes[5].id;

						GtsPoint *p0 = gts_point_new(gts_point_class(),
								coords[node_change],
								coords[node_change + 1],
								coords[node_change + 2]);

						d_c1 = gts_point_distance(point[8], p0);
						d_c2 = gts_point_distance(point[5], p0);

						double alpha = atan(d_c1/d_c2);

						for(int co = 0; co<8;co++){
							node_change = 3*elem->nodes[co].id;
							coords[node_change    ] = coords[node_change    ] - ref1 - dx*0.5;
							coords[node_change + 2] = coords[node_change + 2] - ref3 - dz*0.5;
							coords[node_change    ] = coords[node_change    ]*cos(alpha) - coords[node_change + 2]*sin(alpha);
							coords[node_change + 2] = coords[node_change    ]*sin(alpha) + coords[node_change + 2]*cos(alpha);
							coords[node_change    ] = coords[node_change    ] + ref1 + dx*0.5;
							coords[node_change + 2] = coords[node_change + 2] + ref3 + dz*0.5;
						}
						//move 5

					}else if(edge_list[9]&&edge_list[5]){
						int node_change = 3*elem->nodes[5].id;

						GtsPoint *p0 = gts_point_new(gts_point_class(),
								coords[node_change],
								coords[node_change + 1],
								coords[node_change + 2]);

						d_c1 = gts_point_distance(point[9], p0);
						d_c2 = gts_point_distance(point[5], p0);

						double alpha = atan(d_c1/d_c2);

						for(int co = 0; co<8;co++){
							node_change = 3*elem->nodes[co].id;
							coords[node_change + 1] = coords[node_change + 1] - ref2 - dy*0.5;
							coords[node_change + 2] = coords[node_change + 2] - ref3 - dz*0.5;
							coords[node_change + 1] = coords[node_change + 1]*cos(alpha) - coords[node_change + 2]*sin(alpha);
							coords[node_change + 2] = coords[node_change + 1]*sin(alpha) + coords[node_change + 2]*cos(alpha);
							coords[node_change + 1] = coords[node_change + 1] + ref2 + dy*0.5;
							coords[node_change + 2] = coords[node_change + 2] + ref3 + dz*0.5;
						}
						//move 5

					}else if(edge_list[9]&&edge_list[6]){
						int node_change = 3*elem->nodes[6].id;

						GtsPoint *p0 = gts_point_new(gts_point_class(),
								coords[node_change],
								coords[node_change + 1],
								coords[node_change + 2]);

						d_c1 = gts_point_distance(point[9], p0);
						d_c2 = gts_point_distance(point[6], p0);

						double alpha = atan(d_c1/d_c2);

						for(int co = 0; co<8;co++){
							node_change = 3*elem->nodes[co].id;
							coords[node_change + 1] = coords[node_change + 1] - ref2 - dy*0.5;
							coords[node_change + 2] = coords[node_change + 2] - ref3 - dz*0.5;
							coords[node_change + 1] = coords[node_change + 1]*cos(alpha) - coords[node_change + 2]*sin(alpha);
							coords[node_change + 2] = coords[node_change + 1]*sin(alpha) + coords[node_change + 2]*cos(alpha);
							coords[node_change + 1] = coords[node_change + 1] + ref2 + dy*0.5;
							coords[node_change + 2] = coords[node_change + 2] + ref3 + dz*0.5;
						}
						//move 6

					}else if(edge_list[10]&&edge_list[6]){
						int node_change = 3*elem->nodes[6].id;

						GtsPoint *p0 = gts_point_new(gts_point_class(),
								coords[node_change],
								coords[node_change + 1],
								coords[node_change + 2]);

						d_c1 = gts_point_distance(point[10], p0);
						d_c2 = gts_point_distance(point[6], p0);

						double alpha = atan(d_c1/d_c2);

						for(int co = 0; co<8;co++){
							node_change = 3*elem->nodes[co].id;
							coords[node_change    ] = coords[node_change    ] - ref1 - dx*0.5;
							coords[node_change + 2] = coords[node_change + 2] - ref3 - dz*0.5;
							coords[node_change    ] = coords[node_change    ]*cos(alpha) - coords[node_change + 2]*sin(alpha);
							coords[node_change + 2] = coords[node_change    ]*sin(alpha) + coords[node_change + 2]*cos(alpha);
							coords[node_change    ] = coords[node_change    ] + ref1 + dx*0.5;
							coords[node_change + 2] = coords[node_change + 2] + ref3 + dz*0.5;
						}
						//move 6

					}else if(edge_list[10]&&edge_list[7]){
						int node_change = 3*elem->nodes[7].id;

						GtsPoint *p0 = gts_point_new(gts_point_class(),
								coords[node_change],
								coords[node_change + 1],
								coords[node_change + 2]);

						d_c1 = gts_point_distance(point[10], p0);
						d_c2 = gts_point_distance(point[7], p0);

						double alpha = atan(d_c1/d_c2);

						for(int co = 0; co<8;co++){
							node_change = 3*elem->nodes[co].id;
							coords[node_change    ] = coords[node_change    ] - ref1 - dx*0.5;
							coords[node_change + 2] = coords[node_change + 2] - ref3 - dz*0.5;
							coords[node_change    ] = coords[node_change    ]*cos(alpha) - coords[node_change + 2]*sin(alpha);
							coords[node_change + 2] = coords[node_change    ]*sin(alpha) + coords[node_change + 2]*cos(alpha);
							coords[node_change    ] = coords[node_change    ] + ref1 + dx*0.5;
							coords[node_change + 2] = coords[node_change + 2] + ref3 + dz*0.5;
						}
						//move 7

					}else if(edge_list[11]&&edge_list[7]){
						int node_change = 3*elem->nodes[7].id;

						GtsPoint *p0 = gts_point_new(gts_point_class(),
								coords[node_change],
								coords[node_change + 1],
								coords[node_change + 2]);

						d_c1 = gts_point_distance(point[11], p0);
						d_c2 = gts_point_distance(point[7], p0);

						double alpha = atan(d_c1/d_c2);

						for(int co = 0; co<8;co++){
							node_change = 3*elem->nodes[co].id;
							coords[node_change + 1] = coords[node_change + 1] - ref2 - dy*0.5;
							coords[node_change + 2] = coords[node_change + 2] - ref3 - dz*0.5;
							coords[node_change + 1] = coords[node_change + 1]*cos(alpha) - coords[node_change + 2]*sin(alpha);
							coords[node_change + 2] = coords[node_change + 1]*sin(alpha) + coords[node_change + 2]*cos(alpha);
							coords[node_change + 1] = coords[node_change + 1] + ref2 + dy*0.5;
							coords[node_change + 2] = coords[node_change + 2] + ref3 + dz*0.5;
						}
						//move 7
					}


					CheckTemplate(mesh, coords, element_ids_local,false);
					elem = (octant_t*) sc_array_index(&mesh->elements, elements_ids[iel]);

					fprintf(fdbg,"Pad local:    ");
					for(int co = 0; co < element_ids_local.size(); co++){
						octant_t *h = (octant_t*) sc_array_index(&mesh->elements, element_ids_local[co]);
						var_aux_1[co]=h->pad;
						fprintf(fdbg,"%d ",h->pad);
					}
					fprintf(fdbg,"\n");

					if(elem->pad!=4 && elem->pad!=-10){
						flag2=false;
						for(int co = 0; co < element_ids_local.size(); co++){
							if(element_ids_local[co]!=elements_ids[iel]){
								if(var_aux[co]==var_aux_1[co]){
									flag2=false;
								}else{
									if(var_aux_1[co]==4){
										flag2=true;
										break;
									}else if(var_aux_1[co]==-10){
										flag2=true;
										break;
									}
								}
							}
						}
					}
				}

				if(!flag2){
					flag1=false;
					flag2=false;
					printf("To livreeeee rand tipo 4 rodei a bahiana El:%d, n_iter:%d\n",elements_ids[iel],n_iter);
					break;
				}

				if( flag2 ){
					for(int n_nodes = 0; n_nodes<8;n_nodes++){
						int node = 3*elem->nodes[n_nodes].id;
						coords[node]   = coords_orgi[3*n_nodes];
						coords[node+1] = coords_orgi[3*n_nodes+1];
						coords[node+2] = coords_orgi[3*n_nodes+2];
					}
				}


#endif

#if 0
				template2Rand(mesh, coords, elem, v_1);
				CheckTemplate(mesh, coords, element_ids_local,false);
				elem = (octant_t*) sc_array_index(&mesh->elements, elements_ids[iel]);

				fprintf(fdbg,"Pad local:    ");
				for(int co = 0; co < element_ids_local.size(); co++){
					octant_t *h = (octant_t*) sc_array_index(&mesh->elements, element_ids_local[co]);
					var_aux_1[co]=h->pad;
					fprintf(fdbg,"%d ",h->pad);
				}
				fprintf(fdbg,"\n");

				if(elem->pad!=4 && elem->pad!=-10){
					flag2=false;
					for(int co = 0; co < element_ids_local.size(); co++){
						if(element_ids_local[co]!=elements_ids[iel]){
							if(var_aux[co]==var_aux_1[co]){
								flag2=false;
							}else{
								if(var_aux_1[co]==4){
									flag2=true;
									break;
								}else if(var_aux_1[co]==-10){
									flag2=true;
									break;
								}
							}
						}
					}
				}

				if(!flag2){
					flag1=false;
					flag2=false;
					printf("To livreeeee rand El:%d, n_iter:%d\n",elements_ids[iel],n_iter);
					break;
				}

				if(elem->pad==4){
					//	CuttedEdges (mesh, coords, elements_ids, iel, &point[12], &edge_list[12]);

					GtsSegment * segments[12]={0};
					int Edge2GNode[12][2];
					GtsPoint * point[12] = {NULL};
					bool edge_list[12] = {false};

					for (int edge = 0; edge < 12; ++edge) {
						point[edge] = NULL;

						int node1 = elem->nodes[EdgeVerticesMap[edge][0]].id;
						int node2 = elem->nodes[EdgeVerticesMap[edge][1]].id;
						Edge2GNode[edge][0] = node1 <= node2 ? node1 : node2;
						Edge2GNode[edge][1] = node1 >= node2 ? node1 : node2;

						GtsVertex *v1 = gts_vertex_new(gts_vertex_class(), coords[node1 * 3], coords[node1 * 3 + 1], coords[node1 * 3 + 2]);
						GtsVertex *v2 = gts_vertex_new(gts_vertex_class(), coords[node2 * 3], coords[node2 * 3 + 1], coords[node2 * 3 + 2]);

						segments[edge] = gts_segment_new(gts_segment_class(), v1, v2);

						GtsBBox *sb = gts_bbox_segment(gts_bbox_class(), segments[edge]);
						GSList* list = gts_bb_tree_overlap(mesh->gdata.bbt, sb);

						edge_list[edge] = false;

						if (list == NULL) continue;
						while (list) {
							GtsBBox *b = GTS_BBOX(list->data);
							point[edge] = SegmentTriangleIntersection(segments[edge], GTS_TRIANGLE(b->bounded));
							if (point[edge]) {
								edge_list[edge] = true;
								break;
							}
							list = list->next;
						}
					}

					double d_c1;
					double d_c2;

					double ref1 = 	coords[3*elem->nodes[0].id ] ;
					double ref2 = 	coords[3*elem->nodes[0].id + 1 ] ;
					double ref3 = coords[3*elem->nodes[0].id + 3 ] ;
					double dx = coords[3*elem->nodes[1].id ] - coords[3*elem->nodes[0].id ];
					double dy = coords[3*elem->nodes[2].id +1] - coords[3*elem->nodes[1].id+1 ];
					double dz = coords[3*elem->nodes[0].id +2] - coords[3*elem->nodes[4].id+2 ];

					if(edge_list[0]&&edge_list[1]){
						int node_change = 3*elem->nodes[1].id;

						GtsPoint *p0 = gts_point_new(gts_point_class(),
								coords[node_change],
								coords[node_change + 1],
								coords[node_change + 2]);

						d_c1 = gts_point_distance(point[0], p0);
						d_c2 = gts_point_distance(point[1], p0);

						double alpha = atan(d_c2/d_c1);

						for(int co = 0; co<8;co++){
							node_change = 3*elem->nodes[co].id;
							coords[node_change  ] = coords[node_change   ] - ref1 - dx*0.5;
							coords[node_change +1  ] =coords[node_change +1  ] -ref2-dy/2;
							coords[node_change  ] = coords[node_change  ]*cos(alpha)  - coords[node_change+1]*sin(alpha);
							coords[node_change +1  ] = coords[node_change]*sin(alpha) +coords[node_change +1  ]*cos(alpha);
							coords[node_change  ] = coords[node_change   ] + ref1 + dx*0.5;
							coords[node_change +1  ] =coords[node_change +1  ] + ref2 + dy/2;
						}
						//move 1

					}else if(edge_list[1]&&edge_list[2]){
						int node_change = 3*elem->nodes[2].id;

						GtsPoint *p0 = gts_point_new(gts_point_class(),
								coords[node_change],
								coords[node_change + 1],
								coords[node_change + 2]);

						d_c1 = gts_point_distance(point[2], p0);
						d_c2 = gts_point_distance(point[1], p0);

						double alpha = atan(d_c2/d_c1);

						for(int co = 0; co<8;co++){
							node_change = 3*elem->nodes[co].id;
							coords[node_change  ] = coords[node_change   ] - ref1 - dx*0.5;
							coords[node_change +1  ] =coords[node_change +1  ] -ref2-dy/2;
							coords[node_change  ] = coords[node_change  ]*cos(alpha)  - coords[node_change+1]*sin(alpha);
							coords[node_change +1  ] = coords[node_change  ]*sin(alpha) +coords[node_change +1  ]*cos(alpha);
							coords[node_change  ] = coords[node_change   ] + ref1 + dx*0.5;
							coords[node_change +1  ] =coords[node_change +1  ] + ref2 + dy/2;
						}
						//move 2

					}else if(edge_list[2]&&edge_list[3]){
						int node_change = 3*elem->nodes[3].id;

						GtsPoint *p0 = gts_point_new(gts_point_class(),
								coords[node_change],
								coords[node_change + 1],
								coords[node_change + 2]);

						d_c1 = gts_point_distance(point[2], p0);
						d_c2 = gts_point_distance(point[3], p0);

						double alpha = atan(d_c2/d_c1);

						for(int co = 0; co<8;co++){
							node_change = 3*elem->nodes[co].id;
							coords[node_change  ] = coords[node_change   ] - ref1 - dx*0.5;
							coords[node_change +1  ] =coords[node_change +1  ] -ref2-dy/2;
							coords[node_change  ] = coords[node_change  ]*cos(alpha)  - coords[node_change+1]*sin(alpha);
							coords[node_change +1  ] = coords[node_change  ]*sin(alpha) +coords[node_change +1  ]*cos(alpha);
							coords[node_change  ] = coords[node_change   ] + ref1 + dx*0.5;
							coords[node_change +1  ] =coords[node_change +1  ] + ref2 + dy/2;
						}
						//move 3

					}else if(edge_list[3]&&edge_list[0]){
						int node_change = 3*elem->nodes[0].id;

						GtsPoint *p0 = gts_point_new(gts_point_class(),
								coords[node_change],
								coords[node_change + 1],
								coords[node_change + 2]);

						d_c1 = gts_point_distance(point[0], p0);
						d_c2 = gts_point_distance(point[3], p0);

						double alpha = atan(d_c2/d_c1);

						for(int co = 0; co<8;co++){
							node_change = 3*elem->nodes[co].id;
							coords[node_change  ] = coords[node_change   ] - ref1 - dx*0.5;
							coords[node_change +1  ] =coords[node_change +1  ] -ref2-dy/2;
							coords[node_change  ] = coords[node_change  ]*cos(alpha)  - coords[node_change+1]*sin(alpha);
							coords[node_change +1  ] = coords[node_change  ]*sin(alpha) +coords[node_change +1  ]*cos(alpha);
							coords[node_change  ] = coords[node_change   ] + ref1 + dx*0.5;
							coords[node_change +1  ] =coords[node_change +1  ] + ref2 + dy/2;
						}
						//move 0

					}else if(edge_list[8]&&edge_list[9]){
						int node_change = 3*elem->nodes[5].id;

						GtsPoint *p0 = gts_point_new(gts_point_class(),
								coords[node_change],
								coords[node_change + 1],
								coords[node_change + 2]);

						d_c1 = gts_point_distance(point[8], p0);
						d_c2 = gts_point_distance(point[9], p0);

						double alpha = atan(d_c2/d_c1);

						for(int co = 0; co<8;co++){
							node_change = 3*elem->nodes[co].id;
							coords[node_change  ] = coords[node_change   ] - ref1 - dx*0.5;
							coords[node_change +1  ] =coords[node_change +1  ] -ref2-dy/2;
							coords[node_change  ] = coords[node_change  ]*cos(alpha)  - coords[node_change+1]*sin(alpha);
							coords[node_change +1  ] = coords[node_change]*sin(alpha) +coords[node_change +1  ]*cos(alpha);
							coords[node_change  ] = coords[node_change   ] + ref1 + dx*0.5;
							coords[node_change +1  ] =coords[node_change +1  ] + ref2 + dy/2;
						}
						//move 5

					}else if(edge_list[9]&&edge_list[10]){
						int node_change = 3*elem->nodes[6].id;

						GtsPoint *p0 = gts_point_new(gts_point_class(),
								coords[node_change],
								coords[node_change + 1],
								coords[node_change + 2]);

						d_c1 = gts_point_distance(point[10], p0);
						d_c2 = gts_point_distance(point[9], p0);

						double alpha = atan(d_c2/d_c1);

						for(int co = 0; co<8;co++){
							node_change = 3*elem->nodes[co].id;
							coords[node_change  ] = coords[node_change   ] - ref1 - dx*0.5;
							coords[node_change +1  ] =coords[node_change +1  ] -ref2-dy/2;
							coords[node_change  ] = coords[node_change  ]*cos(alpha)  - coords[node_change+1]*sin(alpha);
							coords[node_change +1  ] = coords[node_change ]*sin(alpha) +coords[node_change +1  ]*cos(alpha);
							coords[node_change  ] = coords[node_change   ] + ref1 + dx*0.5;
							coords[node_change +1  ] =coords[node_change +1  ] + ref2 + dy/2;
						}
						//move 6

					}else if(edge_list[10]&&edge_list[11]){
						int node_change = 3*elem->nodes[7].id;

						GtsPoint *p0 = gts_point_new(gts_point_class(),
								coords[node_change],
								coords[node_change + 1],
								coords[node_change + 2]);

						d_c1 = gts_point_distance(point[10], p0);
						d_c2 = gts_point_distance(point[11], p0);

						double alpha = atan(d_c2/d_c1);

						for(int co = 0; co<8;co++){
							node_change = 3*elem->nodes[co].id;
							coords[node_change  ] = coords[node_change   ] - ref1 - dx*0.5;
							coords[node_change +1  ] =coords[node_change +1  ] -ref2-dy/2;
							coords[node_change  ] = coords[node_change  ]*cos(alpha)  - coords[node_change+1]*sin(alpha);
							coords[node_change +1  ] = coords[node_change]*sin(alpha) +coords[node_change +1  ]*cos(alpha);
							coords[node_change  ] = coords[node_change   ] + ref1 + dx*0.5;
							coords[node_change +1  ] =coords[node_change +1  ] + ref2 + dy/2;
						}
						//move 7

					}else if(edge_list[11]&&edge_list[8]){
						int node_change = 3*elem->nodes[4].id;

						GtsPoint *p0 = gts_point_new(gts_point_class(),
								coords[node_change],
								coords[node_change + 1],
								coords[node_change + 2]);

						d_c1 = gts_point_distance(point[8], p0);
						d_c2 = gts_point_distance(point[11], p0);

						double alpha = atan(d_c2/d_c1);

						for(int co = 0; co<8;co++){
							node_change = 3*elem->nodes[co].id;
							coords[node_change  ] = coords[node_change   ] - ref1 - dx*0.5;
							coords[node_change +1  ] =coords[node_change +1  ] -ref2-dy/2;
							coords[node_change  ] = coords[node_change  ]*cos(alpha)  - coords[node_change+1]*sin(alpha);
							coords[node_change +1  ] = coords[node_change ]*sin(alpha) +coords[node_change +1  ]*cos(alpha);
							coords[node_change  ] = coords[node_change   ] + ref1 + dx*0.5;
							coords[node_change +1  ] =coords[node_change +1  ] + ref2 + dy/2;
						}
						//move 4

					}else if( edge_list[0]&&edge_list[4] ){
						int node_change = 3*elem->nodes[0].id;

						GtsPoint *p0 = gts_point_new(gts_point_class(),
								coords[node_change],
								coords[node_change + 1],
								coords[node_change + 2]);

						d_c1 = gts_point_distance(point[0], p0);
						d_c2 = gts_point_distance(point[4], p0);

						double alpha = atan(d_c2/d_c1);

						for(int co = 0; co<8;co++){
							node_change = 3*elem->nodes[co].id;
							coords[node_change  ] = coords[node_change   ] - ref1 - dx*0.5;
							coords[node_change +2  ] =coords[node_change +2  ] -ref3-dz/2;
							coords[node_change  ] = coords[node_change  ]*cos(alpha)  - coords[node_change+2]*sin(alpha);
							coords[node_change +2  ] = coords[node_change]*sin(alpha) +coords[node_change +2  ]*cos(alpha);
							coords[node_change  ] = coords[node_change   ] + ref1 + dx*0.5;
							coords[node_change +2  ] =coords[node_change +1  ] + ref3 + dz/2;
						}
						//move 0

					}else if( edge_list[3]&&edge_list[4] ){
						int node_change = 3*elem->nodes[0].id;

						GtsPoint *p0 = gts_point_new(gts_point_class(),
								coords[node_change],
								coords[node_change + 1],
								coords[node_change + 2]);

						d_c1 = gts_point_distance(point[3], p0);
						d_c2 = gts_point_distance(point[4], p0);

						double alpha = atan(d_c2/d_c1);

						for(int co = 0; co<8;co++){
							node_change = 3*elem->nodes[co].id;
							coords[node_change + 1] = coords[node_change + 1] - ref2 - dy*0.5;
							coords[node_change + 2] = coords[node_change + 2] - ref3 - dz*0.5;
							coords[node_change + 1] = coords[node_change + 1]*cos(alpha) - coords[node_change + 2]*sin(alpha);
							coords[node_change + 2] = coords[node_change + 1]*sin(alpha) + coords[node_change + 2]*cos(alpha);
							coords[node_change + 1] = coords[node_change + 1] + ref2 + dy*0.5;
							coords[node_change + 2] = coords[node_change + 2] + ref3 + dz*0.5;
						}
						//move 0

					}else if(edge_list[0]&&edge_list[5]){
						int node_change = 3*elem->nodes[1].id;

						GtsPoint *p0 = gts_point_new(gts_point_class(),
								coords[node_change],
								coords[node_change + 1],
								coords[node_change + 2]);

						d_c1 = gts_point_distance(point[0], p0);
						d_c2 = gts_point_distance(point[5], p0);

						double alpha = atan(d_c2/d_c1);

						for(int co = 0; co<8;co++){
							node_change = 3*elem->nodes[co].id;
							coords[node_change    ] = coords[node_change    ] - ref1 - dx*0.5;
							coords[node_change + 2] = coords[node_change + 2] - ref3 - dz*0.5;
							coords[node_change    ] = coords[node_change    ]*cos(alpha) - coords[node_change + 2]*sin(alpha);
							coords[node_change + 2] = coords[node_change    ]*sin(alpha) + coords[node_change + 2]*cos(alpha);
							coords[node_change    ] = coords[node_change    ] + ref1 + dx*0.5;
							coords[node_change + 2] = coords[node_change + 2] + ref3 + dz*0.5;
						}
						//move 1

					}else if(edge_list[1]&&edge_list[5]){
						int node_change = 3*elem->nodes[1].id;

						GtsPoint *p0 = gts_point_new(gts_point_class(),
								coords[node_change],
								coords[node_change + 1],
								coords[node_change + 2]);

						d_c1 = gts_point_distance(point[1], p0);
						d_c2 = gts_point_distance(point[5], p0);

						double alpha = atan(d_c2/d_c1);

						for(int co = 0; co<8;co++){
							node_change = 3*elem->nodes[co].id;
							coords[node_change + 1] = coords[node_change + 1] - ref2 - dy*0.5;
							coords[node_change + 2] = coords[node_change + 2] - ref3 - dz*0.5;
							coords[node_change + 1] = coords[node_change + 1]*cos(alpha) - coords[node_change + 2]*sin(alpha);
							coords[node_change + 2] = coords[node_change + 1]*sin(alpha) + coords[node_change + 2]*cos(alpha);
							coords[node_change + 1] = coords[node_change + 1] + ref2 + dy*0.5;
							coords[node_change + 2] = coords[node_change + 2] + ref3 + dz*0.5;
						}
						//move 1

					}else if(edge_list[1]&&edge_list[6]){
						int node_change = 3*elem->nodes[2].id;

						GtsPoint *p0 = gts_point_new(gts_point_class(),
								coords[node_change],
								coords[node_change + 1],
								coords[node_change + 2]);

						d_c1 = gts_point_distance(point[1], p0);
						d_c2 = gts_point_distance(point[6], p0);

						double alpha = atan(d_c2/d_c1);

						for(int co = 0; co<8;co++){
							node_change = 3*elem->nodes[co].id;
							coords[node_change + 1] = coords[node_change + 1] - ref2 - dy*0.5;
							coords[node_change + 2] = coords[node_change + 2] - ref3 - dz*0.5;
							coords[node_change + 1] = coords[node_change + 1]*cos(alpha) - coords[node_change + 2]*sin(alpha);
							coords[node_change + 2] = coords[node_change + 1]*sin(alpha) + coords[node_change + 2]*cos(alpha);
							coords[node_change + 1] = coords[node_change + 1] + ref2 + dy*0.5;
							coords[node_change + 2] = coords[node_change + 2] + ref3 + dz*0.5;
						}
						//move 2

					}else if(edge_list[2]&&edge_list[6]){
						int node_change = 3*elem->nodes[2].id;

						GtsPoint *p0 = gts_point_new(gts_point_class(),
								coords[node_change],
								coords[node_change + 1],
								coords[node_change + 2]);

						d_c1 = gts_point_distance(point[2], p0);
						d_c2 = gts_point_distance(point[6], p0);

						double alpha = atan(d_c2/d_c1);

						for(int co = 0; co<8;co++){
							node_change = 3*elem->nodes[co].id;
							coords[node_change    ] = coords[node_change    ] - ref1 - dx*0.5;
							coords[node_change + 2] = coords[node_change + 2] - ref3 - dz*0.5;
							coords[node_change    ] = coords[node_change    ]*cos(alpha) - coords[node_change + 2]*sin(alpha);
							coords[node_change + 2] = coords[node_change    ]*sin(alpha) + coords[node_change + 2]*cos(alpha);
							coords[node_change    ] = coords[node_change    ] + ref1 + dx*0.5;
							coords[node_change + 2] = coords[node_change + 2] + ref3 + dz*0.5;
						}
						//move 2

					}else if(edge_list[2]&&edge_list[7]){
						int node_change = 3*elem->nodes[3].id;

						GtsPoint *p0 = gts_point_new(gts_point_class(),
								coords[node_change],
								coords[node_change + 1],
								coords[node_change + 2]);

						d_c1 = gts_point_distance(point[2], p0);
						d_c2 = gts_point_distance(point[7], p0);

						double alpha = atan(d_c2/d_c1);

						for(int co = 0; co<8;co++){
							node_change = 3*elem->nodes[co].id;
							coords[node_change    ] = coords[node_change    ] - ref1 - dx*0.5;
							coords[node_change + 2] = coords[node_change + 2] - ref3 - dz*0.5;
							coords[node_change    ] = coords[node_change    ]*cos(alpha) - coords[node_change + 2]*sin(alpha);
							coords[node_change + 2] = coords[node_change    ]*sin(alpha) + coords[node_change + 2]*cos(alpha);
							coords[node_change    ] = coords[node_change    ] + ref1 + dx*0.5;
							coords[node_change + 2] = coords[node_change + 2] + ref3 + dz*0.5;
						}
						//move 3

					}else if(edge_list[3]&&edge_list[7]){
						int node_change = 3*elem->nodes[3].id;

						GtsPoint *p0 = gts_point_new(gts_point_class(),
								coords[node_change],
								coords[node_change + 1],
								coords[node_change + 2]);

						d_c1 = gts_point_distance(point[3], p0);
						d_c2 = gts_point_distance(point[7], p0);

						double alpha = atan(d_c2/d_c1);

						for(int co = 0; co<8;co++){
							node_change = 3*elem->nodes[co].id;
							coords[node_change + 1] = coords[node_change + 1] - ref2 - dy*0.5;
							coords[node_change + 2] = coords[node_change + 2] - ref3 - dz*0.5;
							coords[node_change + 1] = coords[node_change + 1]*cos(alpha) - coords[node_change + 2]*sin(alpha);
							coords[node_change + 2] = coords[node_change + 1]*sin(alpha) + coords[node_change + 2]*cos(alpha);
							coords[node_change + 1] = coords[node_change + 1] + ref2 + dy*0.5;
							coords[node_change + 2] = coords[node_change + 2] + ref3 + dz*0.5;
						}
						//move 3

					}else if(edge_list[8]&&edge_list[4]){
						int node_change = 3*elem->nodes[4].id;

						GtsPoint *p0 = gts_point_new(gts_point_class(),
								coords[node_change],
								coords[node_change + 1],
								coords[node_change + 2]);

						d_c1 = gts_point_distance(point[8], p0);
						d_c2 = gts_point_distance(point[4], p0);

						double alpha = atan(d_c2/d_c1);

						for(int co = 0; co<8;co++){
							node_change = 3*elem->nodes[co].id;
							coords[node_change    ] = coords[node_change    ] - ref1 - dx*0.5;
							coords[node_change + 2] = coords[node_change + 2] - ref3 - dz*0.5;
							coords[node_change    ] = coords[node_change    ]*cos(alpha) - coords[node_change + 2]*sin(alpha);
							coords[node_change + 2] = coords[node_change    ]*sin(alpha) + coords[node_change + 2]*cos(alpha);
							coords[node_change    ] = coords[node_change    ] + ref1 + dx*0.5;
							coords[node_change + 2] = coords[node_change + 2] + ref3 + dz*0.5;
						}
						//move 4

					}else if(edge_list[11]&&edge_list[4]){
						int node_change = 3*elem->nodes[4].id;

						GtsPoint *p0 = gts_point_new(gts_point_class(),
								coords[node_change],
								coords[node_change + 1],
								coords[node_change + 2]);

						d_c1 = gts_point_distance(point[11], p0);
						d_c2 = gts_point_distance(point[4], p0);

						double alpha = atan(d_c2/d_c1);

						for(int co = 0; co<8;co++){
							node_change = 3*elem->nodes[co].id;
							coords[node_change + 1] = coords[node_change + 1] - ref2 - dy*0.5;
							coords[node_change + 2] = coords[node_change + 2] - ref3 - dz*0.5;
							coords[node_change + 1] = coords[node_change + 1]*cos(alpha) - coords[node_change + 2]*sin(alpha);
							coords[node_change + 2] = coords[node_change + 1]*sin(alpha) + coords[node_change + 2]*cos(alpha);
							coords[node_change + 1] = coords[node_change + 1] + ref2 + dy*0.5;
							coords[node_change + 2] = coords[node_change + 2] + ref3 + dz*0.5;
						}
						//move 4

					}else if(edge_list[8]&&edge_list[5]){
						int node_change = 3*elem->nodes[5].id;

						GtsPoint *p0 = gts_point_new(gts_point_class(),
								coords[node_change],
								coords[node_change + 1],
								coords[node_change + 2]);

						d_c1 = gts_point_distance(point[8], p0);
						d_c2 = gts_point_distance(point[5], p0);

						double alpha = atan(d_c2/d_c1);

						for(int co = 0; co<8;co++){
							node_change = 3*elem->nodes[co].id;
							coords[node_change    ] = coords[node_change    ] - ref1 - dx*0.5;
							coords[node_change + 2] = coords[node_change + 2] - ref3 - dz*0.5;
							coords[node_change    ] = coords[node_change    ]*cos(alpha) - coords[node_change + 2]*sin(alpha);
							coords[node_change + 2] = coords[node_change    ]*sin(alpha) + coords[node_change + 2]*cos(alpha);
							coords[node_change    ] = coords[node_change    ] + ref1 + dx*0.5;
							coords[node_change + 2] = coords[node_change + 2] + ref3 + dz*0.5;
						}
						//move 5

					}else if(edge_list[9]&&edge_list[5]){
						int node_change = 3*elem->nodes[5].id;

						GtsPoint *p0 = gts_point_new(gts_point_class(),
								coords[node_change],
								coords[node_change + 1],
								coords[node_change + 2]);

						d_c1 = gts_point_distance(point[9], p0);
						d_c2 = gts_point_distance(point[5], p0);

						double alpha = atan(d_c2/d_c1);

						for(int co = 0; co<8;co++){
							node_change = 3*elem->nodes[co].id;
							coords[node_change + 1] = coords[node_change + 1] - ref2 - dy*0.5;
							coords[node_change + 2] = coords[node_change + 2] - ref3 - dz*0.5;
							coords[node_change + 1] = coords[node_change + 1]*cos(alpha) - coords[node_change + 2]*sin(alpha);
							coords[node_change + 2] = coords[node_change + 1]*sin(alpha) + coords[node_change + 2]*cos(alpha);
							coords[node_change + 1] = coords[node_change + 1] + ref2 + dy*0.5;
							coords[node_change + 2] = coords[node_change + 2] + ref3 + dz*0.5;
						}
						//move 5

					}else if(edge_list[9]&&edge_list[6]){
						int node_change = 3*elem->nodes[6].id;

						GtsPoint *p0 = gts_point_new(gts_point_class(),
								coords[node_change],
								coords[node_change + 1],
								coords[node_change + 2]);

						d_c1 = gts_point_distance(point[9], p0);
						d_c2 = gts_point_distance(point[6], p0);

						double alpha = atan(d_c2/d_c1);

						for(int co = 0; co<8;co++){
							node_change = 3*elem->nodes[co].id;
							coords[node_change + 1] = coords[node_change + 1] - ref2 - dy*0.5;
							coords[node_change + 2] = coords[node_change + 2] - ref3 - dz*0.5;
							coords[node_change + 1] = coords[node_change + 1]*cos(alpha) - coords[node_change + 2]*sin(alpha);
							coords[node_change + 2] = coords[node_change + 1]*sin(alpha) + coords[node_change + 2]*cos(alpha);
							coords[node_change + 1] = coords[node_change + 1] + ref2 + dy*0.5;
							coords[node_change + 2] = coords[node_change + 2] + ref3 + dz*0.5;
						}
						//move 6

					}else if(edge_list[10]&&edge_list[6]){
						int node_change = 3*elem->nodes[6].id;

						GtsPoint *p0 = gts_point_new(gts_point_class(),
								coords[node_change],
								coords[node_change + 1],
								coords[node_change + 2]);

						d_c1 = gts_point_distance(point[10], p0);
						d_c2 = gts_point_distance(point[6], p0);

						double alpha = atan(d_c2/d_c1);

						for(int co = 0; co<8;co++){
							node_change = 3*elem->nodes[co].id;
							coords[node_change    ] = coords[node_change    ] - ref1 - dx*0.5;
							coords[node_change + 2] = coords[node_change + 2] - ref3 - dz*0.5;
							coords[node_change    ] = coords[node_change    ]*cos(alpha) - coords[node_change + 2]*sin(alpha);
							coords[node_change + 2] = coords[node_change    ]*sin(alpha) + coords[node_change + 2]*cos(alpha);
							coords[node_change    ] = coords[node_change    ] + ref1 + dx*0.5;
							coords[node_change + 2] = coords[node_change + 2] + ref3 + dz*0.5;
						}
						//move 6

					}else if(edge_list[10]&&edge_list[7]){
						int node_change = 3*elem->nodes[7].id;

						GtsPoint *p0 = gts_point_new(gts_point_class(),
								coords[node_change],
								coords[node_change + 1],
								coords[node_change + 2]);

						d_c1 = gts_point_distance(point[10], p0);
						d_c2 = gts_point_distance(point[7], p0);

						double alpha = atan(d_c2/d_c1);

						for(int co = 0; co<8;co++){
							node_change = 3*elem->nodes[co].id;
							coords[node_change    ] = coords[node_change    ] - ref1 - dx*0.5;
							coords[node_change + 2] = coords[node_change + 2] - ref3 - dz*0.5;
							coords[node_change    ] = coords[node_change    ]*cos(alpha) - coords[node_change + 2]*sin(alpha);
							coords[node_change + 2] = coords[node_change    ]*sin(alpha) + coords[node_change + 2]*cos(alpha);
							coords[node_change    ] = coords[node_change    ] + ref1 + dx*0.5;
							coords[node_change + 2] = coords[node_change + 2] + ref3 + dz*0.5;
						}
						//move 7

					}else if(edge_list[11]&&edge_list[7]){
						int node_change = 3*elem->nodes[7].id;

						GtsPoint *p0 = gts_point_new(gts_point_class(),
								coords[node_change],
								coords[node_change + 1],
								coords[node_change + 2]);

						d_c1 = gts_point_distance(point[11], p0);
						d_c2 = gts_point_distance(point[7], p0);

						double alpha = atan(d_c2/d_c1);

						for(int co = 0; co<8;co++){
							node_change = 3*elem->nodes[co].id;
							coords[node_change + 1] = coords[node_change + 1] - ref2 - dy*0.5;
							coords[node_change + 2] = coords[node_change + 2] - ref3 - dz*0.5;
							coords[node_change + 1] = coords[node_change + 1]*cos(alpha) - coords[node_change + 2]*sin(alpha);
							coords[node_change + 2] = coords[node_change + 1]*sin(alpha) + coords[node_change + 2]*cos(alpha);
							coords[node_change + 1] = coords[node_change + 1] + ref2 + dy*0.5;
							coords[node_change + 2] = coords[node_change + 2] + ref3 + dz*0.5;
						}
						//move 7
					}


					CheckTemplate(mesh, coords, element_ids_local,false);
					elem = (octant_t*) sc_array_index(&mesh->elements, elements_ids[iel]);

					fprintf(fdbg,"Pad local:    ");
					for(int co = 0; co < element_ids_local.size(); co++){
						octant_t *h = (octant_t*) sc_array_index(&mesh->elements, element_ids_local[co]);
						var_aux_1[co]=h->pad;
						fprintf(fdbg,"%d ",h->pad);
					}
					fprintf(fdbg,"\n");

					if(elem->pad!=4 && elem->pad!=-10){
						flag2=false;
						for(int co = 0; co < element_ids_local.size(); co++){
							if(element_ids_local[co]!=elements_ids[iel]){
								if(var_aux[co]==var_aux_1[co]){
									flag2=false;
								}else{
									if(var_aux_1[co]==4){
										flag2=true;
										break;
									}else if(var_aux_1[co]==-10){
										flag2=true;
										break;
									}
								}
							}
						}
					}
				}

				if(!flag2){
					flag1=false;
					flag2=false;
					printf("To livreeeee rand tipo 4 rodei a bahiana El:%d, n_iter:%d\n",elements_ids[iel],n_iter);
					break;
				}

				if( flag2 ){
					for(int n_nodes = 0; n_nodes<8;n_nodes++){
						int node = 3*elem->nodes[n_nodes].id;
						coords[node]   = coords_orgi[3*n_nodes];
						coords[node+1] = coords_orgi[3*n_nodes+1];
						coords[node+2] = coords_orgi[3*n_nodes+2];
					}
				}


#endif

#if 0
				template2Rand(mesh, coords, elem, v_1);
				CheckTemplate(mesh, coords, element_ids_local,false);
				elem = (octant_t*) sc_array_index(&mesh->elements, elements_ids[iel]);

				fprintf(fdbg,"Pad local:    ");
				for(int co = 0; co < element_ids_local.size(); co++){
					octant_t *h = (octant_t*) sc_array_index(&mesh->elements, element_ids_local[co]);
					var_aux_1[co]=h->pad;
					fprintf(fdbg,"%d ",h->pad);
				}
				fprintf(fdbg,"\n");

				if(elem->pad!=4 && elem->pad!=-10){
					flag2=false;
					for(int co = 0; co < element_ids_local.size(); co++){
						if(element_ids_local[co]!=elements_ids[iel]){
							if(var_aux[co]==var_aux_1[co]){
								flag2=false;
							}else{
								if(var_aux_1[co]==4){
									flag2=true;
									break;
								}else if(var_aux_1[co]==-10){
									flag2=true;
									break;
								}
							}
						}
					}
				}

				if(!flag2){
					flag1=false;
					flag2=false;
					printf("To livreeeee rand El:%d, n_iter:%d\n",elements_ids[iel],n_iter);
					break;
				}

				if(elem->pad==4){
					//	CuttedEdges (mesh, coords, elements_ids, iel, &point[12], &edge_list[12]);

					GtsSegment * segments[12]={0};
					int Edge2GNode[12][2];
					GtsPoint * point[12] = {NULL};
					bool edge_list[12] = {false};

					for (int edge = 0; edge < 12; ++edge) {
						point[edge] = NULL;

						int node1 = elem->nodes[EdgeVerticesMap[edge][0]].id;
						int node2 = elem->nodes[EdgeVerticesMap[edge][1]].id;
						Edge2GNode[edge][0] = node1 <= node2 ? node1 : node2;
						Edge2GNode[edge][1] = node1 >= node2 ? node1 : node2;

						GtsVertex *v1 = gts_vertex_new(gts_vertex_class(), coords[node1 * 3], coords[node1 * 3 + 1], coords[node1 * 3 + 2]);
						GtsVertex *v2 = gts_vertex_new(gts_vertex_class(), coords[node2 * 3], coords[node2 * 3 + 1], coords[node2 * 3 + 2]);

						segments[edge] = gts_segment_new(gts_segment_class(), v1, v2);

						GtsBBox *sb = gts_bbox_segment(gts_bbox_class(), segments[edge]);
						GSList* list = gts_bb_tree_overlap(mesh->gdata.bbt, sb);

						edge_list[edge] = false;

						if (list == NULL) continue;
						while (list) {
							GtsBBox *b = GTS_BBOX(list->data);
							point[edge] = SegmentTriangleIntersection(segments[edge], GTS_TRIANGLE(b->bounded));
							if (point[edge]) {
								edge_list[edge] = true;
								break;
							}
							list = list->next;
						}
					}

					double d_c1;
					double d_c2;

					double dx = coords[3*elem->nodes[1].id    ] - coords[3*elem->nodes[0].id    ];
					double dy = coords[3*elem->nodes[2].id + 1] - coords[3*elem->nodes[1].id + 1];
					double dz = coords[3*elem->nodes[0].id + 2] - coords[3*elem->nodes[4].id + 2];

					if(edge_list[0]&&edge_list[1]){
						int node_change = 3*elem->nodes[1].id;

						GtsPoint *p0 = gts_point_new(gts_point_class(),
								coords[node_change],
								coords[node_change + 1],
								coords[node_change + 2]);

						d_c1 = gts_point_distance(point[0], p0);
						d_c2 = gts_point_distance(point[1], p0);

						if(d_c1>d_c2){
							for(int co = 0; co<4;co++){
								node_change = 3*elem->nodes[FaceNodesMap[5][co]].id;
								if(dx*0.5>d_c1){
									coords[node_change  ] = coords[node_change  ]-(d_c1-dx*1.1);
								}else{
									coords[node_change  ] = coords[node_change  ]+(d_c1-dx*1.1);
								}
							}
						}else{
							for(int co = 0; co<4;co++){
								node_change = 3*elem->nodes[FaceNodesMap[5][co]].id;
								if(dy*0.5>d_c2){
									coords[node_change +1  ] = coords[node_change +1  ] - (d_c2-dy*1.1);
								}else{
									coords[node_change +1  ] = coords[node_change +1  ] + (d_c2-dy*1.1);
								}
							}
						}
						//move 1

					}else if(edge_list[1]&&edge_list[2]){
						int node_change = 3*elem->nodes[2].id;

						GtsPoint *p0 = gts_point_new(gts_point_class(),
								coords[node_change],
								coords[node_change + 1],
								coords[node_change + 2]);

						d_c1 = gts_point_distance(point[2], p0);
						d_c2 = gts_point_distance(point[1], p0);

						if(d_c1>d_c2){
							for(int co = 0; co<4;co++){
								node_change = 3*elem->nodes[FaceNodesMap[5][co]].id;
								if(dx*0.5>d_c1){
									coords[node_change  ] = coords[node_change  ]-(d_c1-dx*1.1);
								}else{
									coords[node_change  ] = coords[node_change  ]+(d_c1-dx*1.1);
								}
							}
						}else{
							for(int co = 0; co<4;co++){
								node_change = 3*elem->nodes[FaceNodesMap[5][co]].id;
								if(dy*0.5>d_c2){
									coords[node_change +1  ] = coords[node_change +1  ] - (d_c2-dy*1.1);
								}else{
									coords[node_change +1  ] = coords[node_change +1  ] + (d_c2-dy*1.1);
								}
							}
						}
						//move 2

					}else if(edge_list[2]&&edge_list[3]){
						int node_change = 3*elem->nodes[3].id;

						GtsPoint *p0 = gts_point_new(gts_point_class(),
								coords[node_change],
								coords[node_change + 1],
								coords[node_change + 2]);

						d_c1 = gts_point_distance(point[2], p0);
						d_c2 = gts_point_distance(point[3], p0);

						if(d_c1>d_c2){
							for(int co = 0; co<4;co++){
								node_change = 3*elem->nodes[FaceNodesMap[5][co]].id;
								if(dx*0.5>d_c1){
									coords[node_change  ] = coords[node_change  ]-(d_c1-dx*1.1);
								}else{
									coords[node_change  ] = coords[node_change  ]+(d_c1-dx*1.1);
								}
							}
						}else{
							for(int co = 0; co<4;co++){
								node_change = 3*elem->nodes[FaceNodesMap[5][co]].id;
								if(dy*0.5>d_c2){
									coords[node_change +1  ] = coords[node_change +1  ] - (d_c2-dy*1.1);
								}else{
									coords[node_change +1  ] = coords[node_change +1  ] + (d_c2-dy*1.1);
								}
							}
						}
						//move 3

					}else if(edge_list[3]&&edge_list[0]){
						int node_change = 3*elem->nodes[0].id;

						GtsPoint *p0 = gts_point_new(gts_point_class(),
								coords[node_change],
								coords[node_change + 1],
								coords[node_change + 2]);

						d_c1 = gts_point_distance(point[0], p0);
						d_c2 = gts_point_distance(point[3], p0);

						if(d_c1>d_c2){
							for(int co = 0; co<4;co++){
								node_change = 3*elem->nodes[FaceNodesMap[5][co]].id;
								if(dx*0.5>d_c1){
									coords[node_change  ] = coords[node_change  ]-(d_c1-dx*1.1);
								}else{
									coords[node_change  ] = coords[node_change  ]+(d_c1-dx*1.1);
								}
							}
						}else{
							for(int co = 0; co<4;co++){
								node_change = 3*elem->nodes[FaceNodesMap[5][co]].id;
								if(dy*0.5>d_c2){
									coords[node_change +1  ] = coords[node_change +1  ] - (d_c2-dy*1.1);
								}else{
									coords[node_change +1  ] = coords[node_change +1  ] + (d_c2-dy*1.1);
								}
							}
						}
						//move 0

					}else if(edge_list[8]&&edge_list[9]){
						int node_change = 3*elem->nodes[5].id;

						GtsPoint *p0 = gts_point_new(gts_point_class(),
								coords[node_change],
								coords[node_change + 1],
								coords[node_change + 2]);

						d_c1 = gts_point_distance(point[8], p0);
						d_c2 = gts_point_distance(point[9], p0);

						if(d_c1>d_c2){
							for(int co = 0; co<4;co++){
								node_change = 3*elem->nodes[FaceNodesMap[4][co]].id;
								if(dx*0.5>d_c1){
									coords[node_change  ] = coords[node_change  ]-(d_c1-dx*1.1);
								}else{
									coords[node_change  ] = coords[node_change  ]+(d_c1-dx*1.1);
								}
							}
						}else{
							for(int co = 0; co<4;co++){
								node_change = 3*elem->nodes[FaceNodesMap[4][co]].id;
								if(dy*0.5>d_c2){
									coords[node_change +1  ] = coords[node_change +1  ] - (d_c2-dy*1.1);
								}else{
									coords[node_change +1  ] = coords[node_change +1  ] + (d_c2-dy*1.1);
								}
							}
						}
						//move 5

					}else if(edge_list[9]&&edge_list[10]){
						int node_change = 3*elem->nodes[6].id;

						GtsPoint *p0 = gts_point_new(gts_point_class(),
								coords[node_change],
								coords[node_change + 1],
								coords[node_change + 2]);

						d_c1 = gts_point_distance(point[10], p0);
						d_c2 = gts_point_distance(point[9], p0);

						if(d_c1>d_c2){
							for(int co = 0; co<4;co++){
								node_change = 3*elem->nodes[FaceNodesMap[4][co]].id;
								if(dx*0.5>d_c1){
									coords[node_change  ] = coords[node_change  ]-(d_c1-dx*1.1);
								}else{
									coords[node_change  ] = coords[node_change  ]+(d_c1-dx*1.1);
								}
							}
						}else{
							for(int co = 0; co<4;co++){
								node_change = 3*elem->nodes[FaceNodesMap[4][co]].id;
								if(dy*0.5>d_c2){
									coords[node_change +1  ] = coords[node_change +1  ] - (d_c2-dy*1.1);
								}else{
									coords[node_change +1  ] = coords[node_change +1  ] + (d_c2-dy*1.1);
								}
							}
						}
						//move 6

					}else if(edge_list[10]&&edge_list[11]){
						int node_change = 3*elem->nodes[7].id;

						GtsPoint *p0 = gts_point_new(gts_point_class(),
								coords[node_change],
								coords[node_change + 1],
								coords[node_change + 2]);

						d_c1 = gts_point_distance(point[10], p0);
						d_c2 = gts_point_distance(point[11], p0);

						if(d_c1>d_c2){
							for(int co = 0; co<4;co++){
								node_change = 3*elem->nodes[FaceNodesMap[4][co]].id;
								if(dx*0.5>d_c1){
									coords[node_change  ] = coords[node_change  ]-(d_c1-dx*1.1);
								}else{
									coords[node_change  ] = coords[node_change  ]+(d_c1-dx*1.1);
								}
							}
						}else{
							for(int co = 0; co<4;co++){
								node_change = 3*elem->nodes[FaceNodesMap[4][co]].id;
								if(dy*0.5>d_c2){
									coords[node_change +1  ] = coords[node_change +1  ] - (d_c2-dy*1.1);
								}else{
									coords[node_change +1  ] = coords[node_change +1  ] + (d_c2-dy*1.1);
								}
							}
						}
						//move 7

					}else if(edge_list[11]&&edge_list[8]){
						int node_change = 3*elem->nodes[4].id;

						GtsPoint *p0 = gts_point_new(gts_point_class(),
								coords[node_change],
								coords[node_change + 1],
								coords[node_change + 2]);

						d_c1 = gts_point_distance(point[8], p0);
						d_c2 = gts_point_distance(point[11], p0);

						if(d_c1>d_c2){
							for(int co = 0; co<4;co++){
								node_change = 3*elem->nodes[FaceNodesMap[4][co]].id;
								if(dx*0.5>d_c1){
									coords[node_change  ] = coords[node_change  ]-(d_c1-dx*1.1);
								}else{
									coords[node_change  ] = coords[node_change  ]+(d_c1-dx*1.1);
								}
							}
						}else{
							for(int co = 0; co<4;co++){
								node_change = 3*elem->nodes[FaceNodesMap[4][co]].id;
								if(dy*0.5>d_c2){
									coords[node_change + 1] = coords[node_change +1  ] - (d_c2-dy*1.1);
								}else{
									coords[node_change + 1] = coords[node_change +1  ] + (d_c2-dy*1.1);
								}
							}
						}
						//move 4

					}else if( edge_list[0]&&edge_list[4] ){
						int node_change = 3*elem->nodes[0].id;

						GtsPoint *p0 = gts_point_new(gts_point_class(),
								coords[node_change],
								coords[node_change + 1],
								coords[node_change + 2]);

						d_c1 = gts_point_distance(point[0], p0);
						d_c2 = gts_point_distance(point[4], p0);

						if(d_c1>d_c2){
							for(int co = 0; co<4;co++){
								node_change = 3*elem->nodes[FaceNodesMap[2][co]].id;
								if(dx*0.5>d_c1){
									coords[node_change  ] = coords[node_change  ]-(d_c1-dx*1.1);
								}else{
									coords[node_change  ] = coords[node_change  ]+(d_c1-dx*1.1);
								}
							}
						}else{
							for(int co = 0; co<4;co++){
								node_change = 3*elem->nodes[FaceNodesMap[2][co]].id;
								if(dz*0.5>d_c2){
									coords[node_change + 2] = coords[node_change +2  ] - (d_c2-dz*1.1);
								}else{
									coords[node_change + 2] = coords[node_change +2  ] + (d_c2-dz*1.1);
								}
							}
						}
						//move 0

					}else if( edge_list[3]&&edge_list[4] ){
						int node_change = 3*elem->nodes[0].id;

						GtsPoint *p0 = gts_point_new(gts_point_class(),
								coords[node_change],
								coords[node_change + 1],
								coords[node_change + 2]);

						d_c1 = gts_point_distance(point[3], p0);
						d_c2 = gts_point_distance(point[4], p0);

						if(d_c1>d_c2){
							for(int co = 0; co<4;co++){
								node_change = 3*elem->nodes[FaceNodesMap[0][co]].id;
								if(dy*0.5>d_c1){
									coords[node_change + 1] = coords[node_change + 1]-(d_c1-dy*1.1);
								}else{
									coords[node_change + 1] = coords[node_change + 1]+(d_c1-dy*1.1);
								}
							}
						}else{
							for(int co = 0; co<4;co++){
								node_change = 3*elem->nodes[FaceNodesMap[0][co]].id;
								if(dz*0.5>d_c2){
									coords[node_change + 2] = coords[node_change +2  ] - (d_c2-dz*1.1);
								}else{
									coords[node_change + 2] = coords[node_change +2  ] + (d_c2-dz*1.1);
								}
							}
						}

						//move 0

					}else if(edge_list[0]&&edge_list[5]){
						int node_change = 3*elem->nodes[1].id;

						GtsPoint *p0 = gts_point_new(gts_point_class(),
								coords[node_change],
								coords[node_change + 1],
								coords[node_change + 2]);

						d_c1 = gts_point_distance(point[0], p0);
						d_c2 = gts_point_distance(point[5], p0);

						if(d_c1>d_c2){
							for(int co = 0; co<4;co++){
								node_change = 3*elem->nodes[FaceNodesMap[2][co]].id;
								if(dx*0.5>d_c1){
									coords[node_change  ] = coords[node_change  ]-(d_c1-dx*1.1);
								}else{
									coords[node_change  ] = coords[node_change  ]+(d_c1-dx*1.1);
								}
							}
						}else{
							for(int co = 0; co<4;co++){
								node_change = 3*elem->nodes[FaceNodesMap[2][co]].id;
								if(dz*0.5>d_c2){
									coords[node_change + 2] = coords[node_change +2  ] - (d_c2-dz*1.1);
								}else{
									coords[node_change + 2] = coords[node_change +2  ] + (d_c2-dz*1.1);
								}
							}
						}
						//move 1

					}else if(edge_list[1]&&edge_list[5]){
						int node_change = 3*elem->nodes[1].id;

						GtsPoint *p0 = gts_point_new(gts_point_class(),
								coords[node_change],
								coords[node_change + 1],
								coords[node_change + 2]);

						d_c1 = gts_point_distance(point[1], p0);
						d_c2 = gts_point_distance(point[5], p0);

						if(d_c1>d_c2){
							for(int co = 0; co<4;co++){
								node_change = 3*elem->nodes[FaceNodesMap[1][co]].id;
								if(dy*0.5>d_c1){
									coords[node_change + 1] = coords[node_change + 1]-(d_c1-dy*1.1);
								}else{
									coords[node_change + 1] = coords[node_change + 1]+(d_c1-dy*1.1);
								}
							}
						}else{
							for(int co = 0; co<4;co++){
								node_change = 3*elem->nodes[FaceNodesMap[1][co]].id;
								if(dz*0.5>d_c2){
									coords[node_change + 2] = coords[node_change +2  ] - (d_c2-dz*1.1);
								}else{
									coords[node_change + 2] = coords[node_change +2  ] + (d_c2-dz*1.1);
								}
							}
						}
						//move 1

					}else if(edge_list[1]&&edge_list[6]){
						int node_change = 3*elem->nodes[2].id;

						GtsPoint *p0 = gts_point_new(gts_point_class(),
								coords[node_change],
								coords[node_change + 1],
								coords[node_change + 2]);

						d_c1 = gts_point_distance(point[1], p0);
						d_c2 = gts_point_distance(point[6], p0);

						if(d_c1>d_c2){
							for(int co = 0; co<4;co++){
								node_change = 3*elem->nodes[FaceNodesMap[1][co]].id;
								if(dy*0.5>d_c1){
									coords[node_change + 1] = coords[node_change + 1]-(d_c1-dy*1.1);
								}else{
									coords[node_change + 1] = coords[node_change + 1]+(d_c1-dy*1.1);
								}
							}
						}else{
							for(int co = 0; co<4;co++){
								node_change = 3*elem->nodes[FaceNodesMap[1][co]].id;
								if(dz*0.5>d_c2){
									coords[node_change + 2] = coords[node_change +2  ] - (d_c2-dz*1.1);
								}else{
									coords[node_change + 2] = coords[node_change +2  ] + (d_c2-dz*1.1);
								}
							}
						}
						//move 2

					}else if(edge_list[2]&&edge_list[6]){
						int node_change = 3*elem->nodes[2].id;

						GtsPoint *p0 = gts_point_new(gts_point_class(),
								coords[node_change],
								coords[node_change + 1],
								coords[node_change + 2]);

						d_c1 = gts_point_distance(point[2], p0);
						d_c2 = gts_point_distance(point[6], p0);

						if(d_c1>d_c2){
							for(int co = 0; co<4;co++){
								node_change = 3*elem->nodes[FaceNodesMap[3][co]].id;
								if(dx*0.5>d_c1){
									coords[node_change  ] = coords[node_change  ]-(d_c1-dx*1.1);
								}else{
									coords[node_change  ] = coords[node_change  ]+(d_c1-dx*1.1);
								}
							}
						}else{
							for(int co = 0; co<4;co++){
								node_change = 3*elem->nodes[FaceNodesMap[3][co]].id;
								if(dz*0.5>d_c2){
									coords[node_change + 2] = coords[node_change +2  ] - (d_c2-dz*1.1);
								}else{
									coords[node_change + 2] = coords[node_change +2  ] + (d_c2-dz*1.1);
								}
							}
						}
						//move 2

					}else if(edge_list[2]&&edge_list[7]){
						int node_change = 3*elem->nodes[3].id;

						GtsPoint *p0 = gts_point_new(gts_point_class(),
								coords[node_change],
								coords[node_change + 1],
								coords[node_change + 2]);

						d_c1 = gts_point_distance(point[2], p0);
						d_c2 = gts_point_distance(point[7], p0);

						if(d_c1>d_c2){
							for(int co = 0; co<4;co++){
								node_change = 3*elem->nodes[FaceNodesMap[3][co]].id;
								if(dx*0.5>d_c1){
									coords[node_change  ] = coords[node_change  ]-(d_c1-dx*1.1);
								}else{
									coords[node_change  ] = coords[node_change  ]+(d_c1-dx*1.1);
								}
							}
						}else{
							for(int co = 0; co<4;co++){
								node_change = 3*elem->nodes[FaceNodesMap[3][co]].id;
								if(dz*0.5>d_c2){
									coords[node_change + 2] = coords[node_change +2  ] - (d_c2-dz*1.1);
								}else{
									coords[node_change + 2] = coords[node_change +2  ] + (d_c2-dz*1.1);
								}
							}
						}
						//move 3

					}else if(edge_list[3]&&edge_list[7]){
						int node_change = 3*elem->nodes[3].id;

						GtsPoint *p0 = gts_point_new(gts_point_class(),
								coords[node_change],
								coords[node_change + 1],
								coords[node_change + 2]);

						d_c1 = gts_point_distance(point[3], p0);
						d_c2 = gts_point_distance(point[7], p0);

						if(d_c1>d_c2){
							for(int co = 0; co<4;co++){
								node_change = 3*elem->nodes[FaceNodesMap[0][co]].id;
								if(dy*0.5>d_c1){
									coords[node_change + 1] = coords[node_change + 1]-(d_c1-dy*1.1);
								}else{
									coords[node_change + 1] = coords[node_change + 1]+(d_c1-dy*1.1);
								}
							}
						}else{
							for(int co = 0; co<4;co++){
								node_change = 3*elem->nodes[FaceNodesMap[0][co]].id;
								if(dz*0.5>d_c2){
									coords[node_change + 2] = coords[node_change +2  ] - (d_c2-dz*1.1);
								}else{
									coords[node_change + 2] = coords[node_change +2  ] + (d_c2-dz*1.1);
								}
							}
						}
						//move 3

					}else if(edge_list[8]&&edge_list[4]){
						int node_change = 3*elem->nodes[4].id;

						GtsPoint *p0 = gts_point_new(gts_point_class(),
								coords[node_change],
								coords[node_change + 1],
								coords[node_change + 2]);

						d_c1 = gts_point_distance(point[8], p0);
						d_c2 = gts_point_distance(point[4], p0);

						if(d_c1>d_c2){
							for(int co = 0; co<4;co++){
								node_change = 3*elem->nodes[FaceNodesMap[2][co]].id;
								if(dx*0.5>d_c1){
									coords[node_change  ] = coords[node_change  ]-(d_c1-dx*1.1);
								}else{
									coords[node_change  ] = coords[node_change  ]+(d_c1-dx*1.1);
								}
							}
						}else{
							for(int co = 0; co<4;co++){
								node_change = 3*elem->nodes[FaceNodesMap[2][co]].id;
								if(dz*0.5>d_c2){
									coords[node_change + 2] = coords[node_change +2  ] - (d_c2-dz*1.1);
								}else{
									coords[node_change + 2] = coords[node_change +2  ] + (d_c2-dz*1.1);
								}
							}
						}
						//move 4

					}else if(edge_list[11]&&edge_list[4]){
						int node_change = 3*elem->nodes[4].id;

						GtsPoint *p0 = gts_point_new(gts_point_class(),
								coords[node_change],
								coords[node_change + 1],
								coords[node_change + 2]);

						d_c1 = gts_point_distance(point[11], p0);
						d_c2 = gts_point_distance(point[4], p0);

						if(d_c1>d_c2){
							for(int co = 0; co<4;co++){
								node_change = 3*elem->nodes[FaceNodesMap[0][co]].id;
								if(dy*0.5>d_c1){
									coords[node_change + 1] = coords[node_change + 1]-(d_c1-dy*1.1);
								}else{
									coords[node_change + 1] = coords[node_change + 1]+(d_c1-dy*1.1);
								}
							}
						}else{
							for(int co = 0; co<4;co++){
								node_change = 3*elem->nodes[FaceNodesMap[0][co]].id;
								if(dz*0.5>d_c2){
									coords[node_change + 2] = coords[node_change +2  ] - (d_c2-dz*1.1);
								}else{
									coords[node_change + 2] = coords[node_change +2  ] + (d_c2-dz*1.1);
								}
							}
						}
						//move 4

					}else if(edge_list[8]&&edge_list[5]){
						int node_change = 3*elem->nodes[5].id;

						GtsPoint *p0 = gts_point_new(gts_point_class(),
								coords[node_change],
								coords[node_change + 1],
								coords[node_change + 2]);

						d_c1 = gts_point_distance(point[8], p0);
						d_c2 = gts_point_distance(point[5], p0);

						if(d_c1>d_c2){
							for(int co = 0; co<4;co++){
								node_change = 3*elem->nodes[FaceNodesMap[2][co]].id;
								if(dx*0.5>d_c1){
									coords[node_change  ] = coords[node_change  ]-(d_c1-dx*1.1);
								}else{
									coords[node_change  ] = coords[node_change  ]+(d_c1-dx*1.1);
								}
							}
						}else{
							for(int co = 0; co<4;co++){
								node_change = 3*elem->nodes[FaceNodesMap[2][co]].id;
								if(dz*0.5>d_c2){
									coords[node_change + 2] = coords[node_change +2  ] - (d_c2-dz*1.1);
								}else{
									coords[node_change + 2] = coords[node_change +2  ] + (d_c2-dz*1.1);
								}
							}
						}
						//move 5

					}else if(edge_list[9]&&edge_list[5]){
						int node_change = 3*elem->nodes[5].id;

						GtsPoint *p0 = gts_point_new(gts_point_class(),
								coords[node_change],
								coords[node_change + 1],
								coords[node_change + 2]);

						d_c1 = gts_point_distance(point[9], p0);
						d_c2 = gts_point_distance(point[5], p0);

						if(d_c1>d_c2){
							for(int co = 0; co<4;co++){
								node_change = 3*elem->nodes[FaceNodesMap[1][co]].id;
								if(dy*0.5>d_c1){
									coords[node_change + 1] = coords[node_change + 1]-(d_c1-dy*1.1);
								}else{
									coords[node_change + 1] = coords[node_change + 1]+(d_c1-dy*1.1);
								}
							}
						}else{
							for(int co = 0; co<4;co++){
								node_change = 3*elem->nodes[FaceNodesMap[1][co]].id;
								if(dz*0.5>d_c2){
									coords[node_change + 2] = coords[node_change +2  ] - (d_c2-dz*1.1);
								}else{
									coords[node_change + 2] = coords[node_change +2  ] + (d_c2-dz*1.1);
								}
							}
						}
						//move 5

					}else if(edge_list[9]&&edge_list[6]){
						int node_change = 3*elem->nodes[6].id;

						GtsPoint *p0 = gts_point_new(gts_point_class(),
								coords[node_change],
								coords[node_change + 1],
								coords[node_change + 2]);

						d_c1 = gts_point_distance(point[9], p0);
						d_c2 = gts_point_distance(point[6], p0);

						if(d_c1>d_c2){
							for(int co = 0; co<4;co++){
								node_change = 3*elem->nodes[FaceNodesMap[1][co]].id;
								if(dy*0.5>d_c1){
									coords[node_change + 1] = coords[node_change + 1]-(d_c1-dy*1.1);
								}else{
									coords[node_change + 1] = coords[node_change + 1]+(d_c1-dy*1.1);
								}
							}
						}else{
							for(int co = 0; co<4;co++){
								node_change = 3*elem->nodes[FaceNodesMap[1][co]].id;
								if(dz*0.5>d_c2){
									coords[node_change + 2] = coords[node_change +2  ] - (d_c2-dz*1.1);
								}else{
									coords[node_change + 2] = coords[node_change +2  ] + (d_c2-dz*1.1);
								}
							}
						}
						//move 6

					}else if(edge_list[10]&&edge_list[6]){
						int node_change = 3*elem->nodes[6].id;

						GtsPoint *p0 = gts_point_new(gts_point_class(),
								coords[node_change],
								coords[node_change + 1],
								coords[node_change + 2]);

						d_c1 = gts_point_distance(point[10], p0);
						d_c2 = gts_point_distance(point[6], p0);

						if(d_c1>d_c2){
							for(int co = 0; co<4;co++){
								node_change = 3*elem->nodes[FaceNodesMap[3][co]].id;
								if(dx*0.5>d_c1){
									coords[node_change  ] = coords[node_change  ]-(d_c1-dx*1.1);
								}else{
									coords[node_change  ] = coords[node_change  ]+(d_c1-dx*1.1);
								}
							}
						}else{
							for(int co = 0; co<4;co++){
								node_change = 3*elem->nodes[FaceNodesMap[3][co]].id;
								if(dz*0.5>d_c2){
									coords[node_change + 2] = coords[node_change +2  ] - (d_c2-dz*1.1);
								}else{
									coords[node_change + 2] = coords[node_change +2  ] + (d_c2-dz*1.1);
								}
							}
						}
						//move 6

					}else if(edge_list[10]&&edge_list[7]){
						int node_change = 3*elem->nodes[7].id;

						GtsPoint *p0 = gts_point_new(gts_point_class(),
								coords[node_change],
								coords[node_change + 1],
								coords[node_change + 2]);

						d_c1 = gts_point_distance(point[10], p0);
						d_c2 = gts_point_distance(point[7], p0);

						if(d_c1>d_c2){
							for(int co = 0; co<4;co++){
								node_change = 3*elem->nodes[FaceNodesMap[3][co]].id;
								if(dx*0.5>d_c1){
									coords[node_change  ] = coords[node_change  ]-(d_c1-dx*1.1);
								}else{
									coords[node_change  ] = coords[node_change  ]+(d_c1-dx*1.1);
								}
							}
						}else{
							for(int co = 0; co<4;co++){
								node_change = 3*elem->nodes[FaceNodesMap[3][co]].id;
								if(dz*0.5>d_c2){
									coords[node_change + 2] = coords[node_change +2  ] - (d_c2-dz*1.1);
								}else{
									coords[node_change + 2] = coords[node_change +2  ] + (d_c2-dz*1.1);
								}
							}
						}
						//move 7

					}else if(edge_list[11]&&edge_list[7]){
						int node_change = 3*elem->nodes[7].id;

						GtsPoint *p0 = gts_point_new(gts_point_class(),
								coords[node_change],
								coords[node_change + 1],
								coords[node_change + 2]);

						d_c1 = gts_point_distance(point[11], p0);
						d_c2 = gts_point_distance(point[7], p0);

						if(d_c1>d_c2){
							for(int co = 0; co<4;co++){
								node_change = 3*elem->nodes[FaceNodesMap[0][co]].id;
								if(dy*0.5>d_c1){
									coords[node_change + 1] = coords[node_change + 1]-(d_c1-dy*1.1);
								}else{
									coords[node_change + 1] = coords[node_change + 1]+(d_c1-dy*1.1);
								}
							}
						}else{
							for(int co = 0; co<4;co++){
								node_change = 3*elem->nodes[FaceNodesMap[0][co]].id;
								if(dz*0.5>d_c2){
									coords[node_change + 2] = coords[node_change +2  ] - (d_c2-dz*1.1);
								}else{
									coords[node_change + 2] = coords[node_change +2  ] + (d_c2-dz*1.1);
								}
							}
						}
						//move 7
					}


					CheckTemplate(mesh, coords, element_ids_local,false);
					elem = (octant_t*) sc_array_index(&mesh->elements, elements_ids[iel]);

					fprintf(fdbg,"Pad local:    ");
					for(int co = 0; co < element_ids_local.size(); co++){
						octant_t *h = (octant_t*) sc_array_index(&mesh->elements, element_ids_local[co]);
						var_aux_1[co]=h->pad;
						fprintf(fdbg,"%d ",h->pad);
					}
					fprintf(fdbg,"\n");

					if(elem->pad!=4 && elem->pad!=-10){
						flag2=false;
						for(int co = 0; co < element_ids_local.size(); co++){
							if(element_ids_local[co]!=elements_ids[iel]){
								if(var_aux[co]==var_aux_1[co]){
									flag2=false;
								}else{
									if(var_aux_1[co]==4){
										flag2=true;
										break;
									}else if(var_aux_1[co]==-10){
										flag2=true;
										break;
									}
								}
							}
						}
					}
				}

				if(!flag2){
					flag1=false;
					flag2=false;
					printf("To livreeeee rand tipo 4 movi minha fesse todinha :) El:%d, n_iter:%d\n",elements_ids[iel],n_iter);
					break;
				}

				if( flag2 ){
					for(int n_nodes = 0; n_nodes<8;n_nodes++){
						int node = 3*elem->nodes[n_nodes].id;
						coords[node]   = coords_orgi[3*n_nodes];
						coords[node+1] = coords_orgi[3*n_nodes+1];
						coords[node+2] = coords_orgi[3*n_nodes+2];
					}
				}


#endif
				n_iter++;
				v_1= 1;
			}

			el_4++;

		}else if(elem->pad==-10){

			v_1 = 0.75;
			n_iter = 0;
			bool flag1 = false;
			bool flag2 = true;

			fprintf(fdbg,"El: %d not handle\n",elements_ids[iel]);
			//printf("El:%d\n",elements_ids[iel]);
			element_ids_local.clear();
			SideEL(mesh, elem->x, elem->y, elem->z ,element_ids_local);

			fprintf(fdbg,"El: %d Case type 10\n",elements_ids[iel]);
			fprintf(fdbg,"El:");
			for(int co = 0; co<element_ids_local.size();co++){
				fprintf(fdbg,"%d ",element_ids_local[co]);
				if(element_ids_local[co]==elements_ids[iel]){ref=co;}
			}
			fprintf(fdbg,"\n");

			fprintf(fdbg,"REF:%d \n",ref);

			fprintf(fdbg,"Pad Original: ");
			for(int c = 0; c < element_ids_local.size(); c++){
				octant_t *h = (octant_t*) sc_array_index(&mesh->elements, element_ids_local[c]);
				var_aux[c]=h->pad;
				fprintf(fdbg,"%d ",var_aux[c]);
			}
			fprintf(fdbg,"\n");

			for(int n_nodes = 0; n_nodes<8;n_nodes++){
				octant_t *elem = (octant_t*) sc_array_index(&mesh->elements, elements_ids[iel]);
				int node = 3*elem->nodes[n_nodes].id;
				coords_orgi[3*n_nodes  ] = coords[node  ];
				coords_orgi[3*n_nodes+1] = coords[node+1];
				coords_orgi[3*n_nodes+2] = coords[node+2];
			}

			while(n_iter<iter_max && flag1){
				template2Rand(mesh, coords, elem, v_1);
				CheckTemplate(mesh, coords, element_ids_local,false);
				elem = (octant_t*) sc_array_index(&mesh->elements, elements_ids[iel]);

				fprintf(fdbg,"Pad local:    ");
				for(int co = 0; co < element_ids_local.size(); co++){
					octant_t *h = (octant_t*) sc_array_index(&mesh->elements, element_ids_local[co]);
					var_aux_1[co]=h->pad;
					fprintf(fdbg,"%d ",h->pad);
				}
				fprintf(fdbg,"\n");

				if(elem->pad!=4 && elem->pad!=-10){
					flag2=false;
					for(int co = 0; co < element_ids_local.size(); co++){
						if(element_ids_local[co]!=elements_ids[iel]){
							if(var_aux[co]==var_aux_1[co]){
								flag2=false;
							}else{
								if(var_aux_1[co]==4){
									flag2=true;
									break;
								}else if(var_aux_1[co]==-10){
									flag2=true;
									break;
								}
							}
						}
					}
				}

				if(!flag2){
					flag1=false;
					flag2=false;
					printf("To livreeeee rand tipo 10 manuuu El:%d, n_iter:%d\n",elements_ids[iel],n_iter);
					break;
				}

				if( flag2 ){
					for(int n_nodes = 0; n_nodes<8;n_nodes++){
						int node = 3*elem->nodes[n_nodes].id;
						coords[node]   = coords_orgi[3*n_nodes];
						coords[node+1] = coords_orgi[3*n_nodes+1];
						coords[node+2] = coords_orgi[3*n_nodes+2];
					}
				}
				n_iter++;
			}
			nao_sei++;

		}

	}

	fclose(fdbg);
#endif
}

/*
 1. Para cada octante,

 1.a vefificar se duas faces paralelas nao sao interceptadas pela
 superficie.
   => Aplica-se o template 1: divide elemento ao meio gerando dois novos elementos

 1.b: quatro faces inteceptadas e as restantes nao sao paralelas
   => Aplica-se o template 2:
           - Adicionar ponto na intersecao seguindo a projecao diagonal do ponto da esquina do elemento,
           -  Tres novos elementos são criados.
            Obs.: Possibilidade de formacao de triangulos.  Solucao: mover pontos.

 1.c: tres faces interceptadas: Veja template3 hexMesh/doc.
 */

void ApplyTemplate(hexa_tree_t* mesh, std::vector<double>& coords, std::vector<int>& elements_ids) {
	bool clamped = true;
	bool face_intecepted[6];
	int conn_p[8];
	int original_conn[8];
	FILE * fdbg;

	int Edge2GNode[12][2]={0};
	int Edge2GNode_s[12][2]={0};
	int Edge2GNode_v[4][2]={0};

	GtsSegment * segments[12]={0};
	GtsSegment * segments_s[12]={0};
	GtsSegment * segments_v[4]={0};

	GtsPoint * point[12]={NULL};
	GtsPoint * point_s[12]={NULL};
	GtsPoint * point_v[4]={NULL};
	bool edge_list[12]={false};
	bool edge_list_s[12]={false};
	bool edge_list_v[4]={false};
	int ed_cont = 0;

	fdbg = fopen("intercepted_faces.dbg", "w");

	sc_hash_array_t* hash_nodes = sc_hash_array_new(sizeof (node_t), edge_hash_fn, edge_equal_fn, &clamped);

	for (int iel = 0; iel < elements_ids.size(); ++iel) {

		octant_t *elem = (octant_t*) sc_array_index(&mesh->elements, elements_ids[iel]);

		for (int i = 0; i < 8; i++) original_conn[i] = elem->nodes[i].id;

		for (int edge = 0; edge < 12; ++edge) {
			point[edge] = NULL;
			int node1 = elem->nodes[EdgeVerticesMap[edge][0]].id;
			int node2 = elem->nodes[EdgeVerticesMap[edge][1]].id;

			Edge2GNode[edge][0] = node1 <= node2 ? node1 : node2;
			Edge2GNode[edge][1] = node1 >= node2 ? node1 : node2;

			GtsVertex *v1 = gts_vertex_new(gts_vertex_class(), coords[node1 * 3], coords[node1 * 3 + 1], coords[node1 * 3 + 2]);
			GtsVertex *v2 = gts_vertex_new(gts_vertex_class(), coords[node2 * 3], coords[node2 * 3 + 1], coords[node2 * 3 + 2]);

			segments[edge] = gts_segment_new(gts_segment_class(), v1, v2);
			GtsBBox *sb = gts_bbox_segment(gts_bbox_class(), segments[edge]);
			GSList* list = gts_bb_tree_overlap(mesh->gdata.bbt, sb);
			edge_list[edge] = false;
			if (list == NULL) continue;
			while (list) {
				GtsBBox *b = GTS_BBOX(list->data);
				point[edge] = SegmentTriangleIntersection(segments[edge], GTS_TRIANGLE(b->bounded));
				if (point[edge]) {
					edge_list[edge] = true;
					ed_cont++;
					break;
				}
				list = list->next;
			}
		}

		//check the diagonals in the surface and find the intersections
		for (int edge = 0; edge < 12; ++edge) {
			point_s[edge] = NULL;
			int node1 = elem->nodes[EdgeVerticesMap_surf_diagonal[edge][0]].id;
			int node2 = elem->nodes[EdgeVerticesMap_surf_diagonal[edge][1]].id;

			Edge2GNode_s[edge][0] = node1 <= node2 ? node1 : node2;
			Edge2GNode_s[edge][1] = node1 >= node2 ? node1 : node2;

			GtsVertex *v1 = gts_vertex_new(gts_vertex_class(), coords[node1 * 3], coords[node1 * 3 + 1], coords[node1 * 3 + 2]);
			GtsVertex *v2 = gts_vertex_new(gts_vertex_class(), coords[node2 * 3], coords[node2 * 3 + 1], coords[node2 * 3 + 2]);

			segments_s[edge] = gts_segment_new(gts_segment_class(), v1, v2);
			GtsBBox *sb = gts_bbox_segment(gts_bbox_class(), segments_s[edge]);
			GSList* list = gts_bb_tree_overlap(mesh->gdata.bbt, sb);
			edge_list_s[edge]=false;
			if (list == NULL) continue;
			while (list) {
				GtsBBox *b = GTS_BBOX(list->data);
				point_s[edge] = SegmentTriangleIntersection(segments_s[edge], GTS_TRIANGLE(b->bounded));
				if (point_s[edge]) {
					edge_list_s[edge]=true;
					break;
				}
				list = list->next;
			}
		}

		//check the diagonals in the volume and find the intersections
		for (int edge = 0; edge < 4; ++edge) {
			point_v[edge] = NULL;
			int node1 = elem->nodes[EdgeVerticesMap_vol_diagonal[edge][0]].id;
			int node2 = elem->nodes[EdgeVerticesMap_vol_diagonal[edge][1]].id;

			Edge2GNode_v[edge][0] = node1 <= node2 ? node1 : node2;
			Edge2GNode_v[edge][1] = node1 >= node2 ? node1 : node2;

			GtsVertex *v1 = gts_vertex_new(gts_vertex_class(), coords[node1 * 3], coords[node1 * 3 + 1], coords[node1 * 3 + 2]);
			GtsVertex *v2 = gts_vertex_new(gts_vertex_class(), coords[node2 * 3], coords[node2 * 3 + 1], coords[node2 * 3 + 2]);

			segments_v[edge] = gts_segment_new(gts_segment_class(), v1, v2);
			GtsBBox *sb = gts_bbox_segment(gts_bbox_class(), segments_v[edge]);
			GSList* list = gts_bb_tree_overlap(mesh->gdata.bbt, sb);
			edge_list_v[edge]=false;
			if (list == NULL) continue;
			while (list) {
				GtsBBox *b = GTS_BBOX(list->data);
				point_v[edge] = SegmentTriangleIntersection(segments_v[edge], GTS_TRIANGLE(b->bounded));
				if (point_v[edge]) {
					edge_list_v[edge]=true;
					break;
				}
				list = list->next;
			}
		}

		//check parallel faces
		for (int face = 0; face < 6; ++face) {
			face_intecepted[face] = false;
			for (int fe = 0; fe < 4; ++fe) {
				int edge = FaceEdgesMap[face][fe];
				if (point[edge] != NULL) {
					face_intecepted[face] = true;
					break;
				}
			}
			if (face_intecepted[face]) continue;
		}

		fprintf(fdbg, "Faces:  Elem. %d: %d %d %d %d %d %d\n", elements_ids[iel], face_intecepted[0],
				face_intecepted[1], face_intecepted[2],
				face_intecepted[3], face_intecepted[4], face_intecepted[5]);
		fprintf(fdbg, "Edges:  Elem. %d: %d %d %d %d %d %d %d %d %d %d %d %d %d\n", elements_ids[iel], edge_list[0],
				edge_list[1], edge_list[2], edge_list[3],
				edge_list[4], edge_list[5], edge_list[6],
				edge_list[7], edge_list[8], edge_list[9],
				edge_list[10], edge_list[11], ed_cont);

		int n_parallel_faces = (face_intecepted[0] && face_intecepted[1]) +
				(face_intecepted[2] && face_intecepted[3]) +
				(face_intecepted[4] && face_intecepted[5]);



		if(elem->pad == 22){
			//(edge_list[0]) && (!edge_list[1]) && (edge_list[2]) && (!edge_list[3]) &&
			//(!edge_list[4]) && (!edge_list[5]) && (!edge_list[6]) && (!edge_list[7]) &&
			//(edge_list[8]) && (!edge_list[9]) && (edge_list[10]) && (!edge_list[11])


		}else if(elem->pad == 23){
			//(!edge_list[0]) && (edge_list[1]) && (!edge_list[2]) && (edge_list[3]) &&
			//(!edge_list[4]) && (!edge_list[5]) && (!edge_list[6]) && (!edge_list[7]) &&
			//(!edge_list[8]) && (edge_list[9]) && (!edge_list[10]) && (edge_list[11])


		}else if(elem->pad == 24){
			//(!edge_list[0]) && (!edge_list[1]) && (!edge_list[2]) && (!edge_list[3]) &&
			//(edge_list[4]) && (edge_list[5]) && (edge_list[6]) && (edge_list[7]) &&
			//(!edge_list[8]) && (!edge_list[9]) && (!edge_list[10]) && (!edge_list[11])

			int cut_edge[4];

			cut_edge[0] = 4;
			cut_edge[1] = 5;
			cut_edge[2] = 6;
			cut_edge[3] = 7;

			GtsPoint *p0 = point[0];
			GtsPoint *p1 = point[0];
			GtsPoint *p2 = point[0];
			GtsPoint *p3 = point[0];
			GtsPoint *p4 = point[0];
			GtsPoint *p5 = point[0];
			GtsPoint *p6 = point[0];
			GtsPoint *p7 = point[0];

			if(p0==NULL) exit(1);
			if(p1==NULL) exit(1);
			if(p2==NULL) exit(1);
			if(p3==NULL) exit(1);
			if(p4==NULL) exit(1);
			if(p5==NULL) exit(1);
			if(p6==NULL) exit(1);
			if(p7==NULL) exit(1);
			/*
			conn_p[0] = AddPointOnEdge(Edge2GNode[cut_edge[0]], hash_nodes, mesh->local_n_nodes, p0, coords);
			conn_p[1] = AddPointOnEdge(Edge2GNode[cut_edge[0]], hash_nodes, mesh->local_n_nodes, p1, coords);
			conn_p[2] = AddPointOnEdge(Edge2GNode[cut_edge[1]], hash_nodes, mesh->local_n_nodes, p2, coords);
			conn_p[3] = AddPointOnEdge(Edge2GNode[cut_edge[1]], hash_nodes, mesh->local_n_nodes, p3, coords);
			conn_p[4] = AddPointOnEdge(Edge2GNode[cut_edge[2]], hash_nodes, mesh->local_n_nodes, p4, coords);
			conn_p[5] = AddPointOnEdge(Edge2GNode[cut_edge[2]], hash_nodes, mesh->local_n_nodes, p5, coords);
			conn_p[6] = AddPointOnEdge(Edge2GNode[cut_edge[3]], hash_nodes, mesh->local_n_nodes, p6, coords);
			conn_p[7] = AddPointOnEdge(Edge2GNode[cut_edge[3]], hash_nodes, mesh->local_n_nodes, p7, coords);
			 */
			octant_t *elem1 = (octant_t*) sc_array_index(&mesh->elements, elements_ids[iel]);

			elem1->nodes[0].id = original_conn[0];
			elem1->nodes[1].id = original_conn[1];
			elem1->nodes[2].id = original_conn[2];
			elem1->nodes[3].id = original_conn[3];

			elem1->nodes[4].id = conn_p[0];
			elem1->nodes[5].id = conn_p[1];
			elem1->nodes[6].id = conn_p[2];
			elem1->nodes[7].id = conn_p[3];

			elem1->level = elem->level;
			elem1->ref = 1;

			octant_t* elem2 = (octant_t*) sc_array_push(&mesh->elements);

			elem2->nodes[0].id = conn_p[0];
			elem2->nodes[1].id = conn_p[0];
			elem2->nodes[2].id = conn_p[0];
			elem2->nodes[3].id = conn_p[0];

			elem2->nodes[4].id = conn_p[0];
			elem2->nodes[5].id = conn_p[1];
			elem2->nodes[6].id = conn_p[2];
			elem2->nodes[7].id = conn_p[3];

			elem2->pad = elem->pad;
			elem2->level = elem->level;
			elem2->ref = 1;

			octant_t* elem3 = (octant_t*) sc_array_push(&mesh->elements);

			elem3->nodes[0].id = original_conn[0];
			elem3->nodes[1].id = original_conn[1];
			elem3->nodes[2].id = original_conn[2];
			elem3->nodes[3].id = original_conn[3];

			elem3->nodes[4].id = conn_p[0];
			elem3->nodes[5].id = conn_p[1];
			elem3->nodes[6].id = conn_p[2];
			elem3->nodes[7].id = conn_p[3];

			elem3->pad = elem->pad;
			elem3->level = elem->level;
			elem3->ref = 1;


		}


#if 0
		// verificação dos elementos
		if(elem->pad==1){
			if (n_parallel_faces == 2 && ed_cont == 4) {
				// Check template 1.
				int cut_edge[4];
				int connec_order_in[8];
				int connec_order_out1[8];
				int connec_order_out2[8];

				if (    (!edge_list[0]) && (edge_list[1]) && (!edge_list[2]) && (edge_list[3]) &&
						(!edge_list[4]) && (!edge_list[5]) && (!edge_list[6]) && (!edge_list[7]) &&
						(!edge_list[8]) && (edge_list[9]) && (!edge_list[10]) && (edge_list[11]) ) {

					cut_edge[0] = 3;
					cut_edge[1] = 1;
					cut_edge[2] = 11;
					cut_edge[3] = 9;

					connec_order_in[0] = 0;
					connec_order_in[1] = 1;
					connec_order_in[2] = 4;
					connec_order_in[3] = 5;
					connec_order_in[4] = 2;
					connec_order_in[5] = 3;
					connec_order_in[6] = 6;
					connec_order_in[7] = 7;

					connec_order_out1[0] = 0;
					connec_order_out1[1] = 1;
					connec_order_out1[2] = 4;
					connec_order_out1[3] = 5;
					connec_order_out1[4] = 3;
					connec_order_out1[5] = 2;
					connec_order_out1[6] = 7;
					connec_order_out1[7] = 6;

					connec_order_out2[0] = 0;
					connec_order_out2[1] = 1;
					connec_order_out2[2] = 4;
					connec_order_out2[3] = 5;
					connec_order_out2[4] = 2;
					connec_order_out2[5] = 3;
					connec_order_out2[6] = 6;
					connec_order_out2[7] = 7;

					elem->pad = 1;

				} else if (     (edge_list[0]) && (!edge_list[1]) && (edge_list[2]) && (!edge_list[3]) &&
						(!edge_list[4]) && (!edge_list[5]) && (!edge_list[6]) && (!edge_list[7]) &&
						(edge_list[8]) && (!edge_list[9]) && (edge_list[10]) && (!edge_list[11]) ) {

					cut_edge[0] = 0;
					cut_edge[1] = 2;
					cut_edge[2] = 8;
					cut_edge[3] = 10;

					connec_order_in[0] = 0;
					connec_order_in[1] = 3;
					connec_order_in[2] = 4;
					connec_order_in[3] = 7;
					connec_order_in[4] = 1;
					connec_order_in[5] = 2;
					connec_order_in[6] = 5;
					connec_order_in[7] = 6;

					connec_order_out1[0] = 0;
					connec_order_out1[1] = 3;
					connec_order_out1[2] = 4;
					connec_order_out1[3] = 7;
					connec_order_out1[4] = 1;
					connec_order_out1[5] = 2;
					connec_order_out1[6] = 5;
					connec_order_out1[7] = 6;

					connec_order_out2[0] = 0;
					connec_order_out2[1] = 3;
					connec_order_out2[2] = 4;
					connec_order_out2[3] = 7;
					connec_order_out2[4] = 1;
					connec_order_out2[5] = 2;
					connec_order_out2[6] = 5;
					connec_order_out2[7] = 6;

					elem->pad = 1;

				} else if (     (!edge_list[0]) && (!edge_list[1]) && (!edge_list[2]) && (!edge_list[3]) &&
						(edge_list[4]) && (edge_list[5]) && (edge_list[6]) && (edge_list[7]) &&
						(!edge_list[8]) && (!edge_list[9]) && (!edge_list[10]) && (!edge_list[11]) ) {

					cut_edge[0] = 4;
					cut_edge[1] = 5;
					cut_edge[2] = 6;
					cut_edge[3] = 7;

					connec_order_in[0] = 0;
					connec_order_in[1] = 1;
					connec_order_in[2] = 2;
					connec_order_in[3] = 3;
					connec_order_in[4] = 4;
					connec_order_in[5] = 5;
					connec_order_in[6] = 6;
					connec_order_in[7] = 7;

					connec_order_out1[0] = 0;
					connec_order_out1[1] = 1;
					connec_order_out1[2] = 2;
					connec_order_out1[3] = 3;
					connec_order_out1[4] = 4;
					connec_order_out1[5] = 5;
					connec_order_out1[6] = 6;
					connec_order_out1[7] = 7;

					connec_order_out2[0] = 0;
					connec_order_out2[1] = 1;
					connec_order_out2[2] = 2;
					connec_order_out2[3] = 3;
					connec_order_out2[4] = 4;
					connec_order_out2[5] = 5;
					connec_order_out2[6] = 6;
					connec_order_out2[7] = 7;

					elem->pad = 1;

				} else {
					elem->pad = -1;
					printf("Warning: Bad Element: %d, case 1\n",elements_ids[iel]);
					printf("Edges:  Elem. %d: %d %d %d %d %d %d %d %d %d %d %d %d %d\n", elements_ids[iel], edge_list[0],
							edge_list[1], edge_list[2], edge_list[3],
							edge_list[4], edge_list[5], edge_list[6],
							edge_list[7], edge_list[8], edge_list[9],
							edge_list[10], edge_list[11], ed_cont);

					continue;
				}

				GtsPoint *p0 = point[cut_edge[0]];
				GtsPoint *p1 = point[cut_edge[1]];
				GtsPoint *p2 = point[cut_edge[2]];
				GtsPoint *p3 = point[cut_edge[3]];

				//g_assert(p0 != NULL);
				//g_assert(p1 != NULL);
				//g_assert(p2 != NULL);
				//g_assert(p3 != NULL);

				if(p0==NULL) exit(1);
				if(p1==NULL) exit(1);
				if(p2==NULL) exit(1);
				if(p3==NULL) exit(1);


				conn_p[0] = AddPointOnEdge(Edge2GNode[cut_edge[0]], hash_nodes, mesh->local_n_nodes, p0, coords);
				conn_p[1] = AddPointOnEdge(Edge2GNode[cut_edge[1]], hash_nodes, mesh->local_n_nodes, p1, coords);
				conn_p[2] = AddPointOnEdge(Edge2GNode[cut_edge[2]], hash_nodes, mesh->local_n_nodes, p2, coords);
				conn_p[3] = AddPointOnEdge(Edge2GNode[cut_edge[3]], hash_nodes, mesh->local_n_nodes, p3, coords);

				octant_t *elem1 = (octant_t*) sc_array_index(&mesh->elements, elements_ids[iel]);

				elem1->nodes[connec_order_out1[0]].id = original_conn[connec_order_in[0]];
				elem1->nodes[connec_order_out1[1]].id = original_conn[connec_order_in[1]];
				elem1->nodes[connec_order_out1[2]].id = original_conn[connec_order_in[2]];
				elem1->nodes[connec_order_out1[3]].id = original_conn[connec_order_in[3]];

				elem1->nodes[connec_order_out1[4]].id = conn_p[0];
				elem1->nodes[connec_order_out1[5]].id = conn_p[1];
				elem1->nodes[connec_order_out1[6]].id = conn_p[2];
				elem1->nodes[connec_order_out1[7]].id = conn_p[3];

				octant_t* elem2 = (octant_t*) sc_array_push(&mesh->elements);

				elem2->nodes[connec_order_out2[0]].id = conn_p[0];
				elem2->nodes[connec_order_out2[1]].id = conn_p[1];
				elem2->nodes[connec_order_out2[2]].id = conn_p[2];
				elem2->nodes[connec_order_out2[3]].id = conn_p[3];

				elem2->nodes[connec_order_out2[4]].id = original_conn[connec_order_in[4]];
				elem2->nodes[connec_order_out2[5]].id = original_conn[connec_order_in[5]];
				elem2->nodes[connec_order_out2[6]].id = original_conn[connec_order_in[6]];
				elem2->nodes[connec_order_out2[7]].id = original_conn[connec_order_in[7]];

				elem2->pad = elem->pad;
				elem2->level = elem->level;


			}
		}else if(elem->pad==2){

			if (n_parallel_faces == 1 && ed_cont==4) {
				// Check template 2.

				int conn_t2[6];
				int cut_edge[4];
				int cut_edge_s[2];
				int connec_order_in1[8];
				int connec_order_out1[8];
				int connec_order_in2[8];
				int connec_order_out2[8];
				int connec_order_in3[8];
				int connec_order_out3[8];

				//Edge 0
				if ( (!edge_list[0]) && (edge_list[1]) && (!edge_list[2]) && (edge_list[3]) &&
						(edge_list[4]) && (edge_list[5]) && (!edge_list[6]) && (!edge_list[7]) &&
						(!edge_list[8]) && (!edge_list[9]) && (!edge_list[10]) && (!edge_list[11]) &&
						(edge_list_s[4]) && (edge_list_s[6]) ) {
					//printf("Entrou na 0!\n");

					cut_edge[0] = 3;
					cut_edge[1] = 1;
					cut_edge[2] = 4;
					cut_edge[3] = 5;

					cut_edge_s[0] = 4;
					cut_edge_s[1] = 6;

					connec_order_in1[0] = 2;
					connec_order_in1[1] = 3;
					connec_order_in1[2] = 6;
					connec_order_in1[3] = 7;
					connec_order_in1[4] = 0;
					connec_order_in1[5] = 1;
					connec_order_in1[6] = 4;
					connec_order_in1[7] = 5;

					connec_order_out1[0] = 2;
					connec_order_out1[1] = 3;
					connec_order_out1[2] = 6;
					connec_order_out1[3] = 7;
					connec_order_out1[4] = 0;
					connec_order_out1[5] = 1;
					connec_order_out1[6] = 4;
					connec_order_out1[7] = 5;

					connec_order_in2[0] = 0;
					connec_order_in2[1] = 1;
					connec_order_in2[2] = 2;
					connec_order_in2[3] = 3;
					connec_order_in2[4] = 4;
					connec_order_in2[5] = 5;
					connec_order_in2[6] = 6;
					connec_order_in2[7] = 7;

					connec_order_out2[0] = 2;
					connec_order_out2[1] = 3;
					connec_order_out2[2] = 5;
					connec_order_out2[3] = 4;
					connec_order_out2[4] = 4;
					connec_order_out2[5] = 5;
					connec_order_out2[6] = 6;
					connec_order_out2[7] = 7;

					connec_order_in3[0] = 0;
					connec_order_in3[1] = 1;
					connec_order_in3[2] = 2;
					connec_order_in3[3] = 3;
					connec_order_in3[4] = 4;
					connec_order_in3[5] = 5;
					connec_order_in3[6] = 6;
					connec_order_in3[7] = 7;

					connec_order_out3[0] = 0;
					connec_order_out3[1] = 1;
					connec_order_out3[2] = 1;
					connec_order_out3[3] = 0;
					connec_order_out3[4] = 2;
					connec_order_out3[5] = 3;
					connec_order_out3[6] = 5;
					connec_order_out3[7] = 4;

					elem->pad = 2;

				}//Edge 1
				else if ((edge_list[0]) && (!edge_list[1]) && (edge_list[2]) && (!edge_list[3]) &&
						(!edge_list[4]) && (edge_list[5]) && (edge_list[6]) && (!edge_list[7]) &&
						(!edge_list[8]) && (!edge_list[9]) && (!edge_list[10]) && (!edge_list[11])&&
						(edge_list_s[1]) && (edge_list_s[3]) ) {
					//printf("Entrou na 1!\n");

					cut_edge[0] = 0;
					cut_edge[1] = 2;
					cut_edge[2] = 5;
					cut_edge[3] = 6;

					cut_edge_s[0] = 1;
					cut_edge_s[1] = 3;

					connec_order_in1[0] = 0;
					connec_order_in1[1] = 3;
					connec_order_in1[2] = 4;
					connec_order_in1[3] = 7;
					connec_order_in1[4] = 1;
					connec_order_in1[5] = 2;
					connec_order_in1[6] = 5;
					connec_order_in1[7] = 6;

					connec_order_out1[0] = 0;
					connec_order_out1[1] = 3;
					connec_order_out1[2] = 4;
					connec_order_out1[3] = 7;
					connec_order_out1[4] = 0;
					connec_order_out1[5] = 1;
					connec_order_out1[6] = 4;
					connec_order_out1[7] = 5;

					connec_order_in2[0] = 0;
					connec_order_in2[1] = 1;
					connec_order_in2[2] = 2;
					connec_order_in2[3] = 3;
					connec_order_in2[4] = 4;
					connec_order_in2[5] = 5;
					connec_order_in2[6] = 6;
					connec_order_in2[7] = 7;

					connec_order_out2[0] = 4;
					connec_order_out2[1] = 2;
					connec_order_out2[2] = 3;
					connec_order_out2[3] = 5;
					connec_order_out2[4] = 4;
					connec_order_out2[5] = 5;
					connec_order_out2[6] = 6;
					connec_order_out2[7] = 7;

					connec_order_in3[0] = 1;
					connec_order_in3[1] = 2;
					connec_order_in3[2] = 0;
					connec_order_in3[3] = 3;
					connec_order_in3[4] = 4;
					connec_order_in3[5] = 5;
					connec_order_in3[6] = 6;
					connec_order_in3[7] = 7;

					connec_order_out3[0] = 1;
					connec_order_out3[1] = 2;
					connec_order_out3[2] = 0;
					connec_order_out3[3] = 1;
					connec_order_out3[4] = 4;
					connec_order_out3[5] = 2;
					connec_order_out3[6] = 3;
					connec_order_out3[7] = 5;

					elem->pad = 2;

				}//Edge 2
				else if ((!edge_list[0]) && (edge_list[1]) && (!edge_list[2]) && (edge_list[3]) &&
						(!edge_list[4]) && (!edge_list[5]) && (edge_list[6]) && (edge_list[7]) &&
						(!edge_list[8]) && (!edge_list[9]) && (!edge_list[10]) && (!edge_list[11])&&
						(edge_list_s[5]) && (edge_list_s[7]) ) {
					//printf("Entrou na 2!\n");

					cut_edge[0] = 3;
					cut_edge[1] = 1;
					cut_edge[2] = 7;
					cut_edge[3] = 6;

					cut_edge_s[0] = 5;
					cut_edge_s[1] = 7;

					connec_order_in1[0] = 0;
					connec_order_in1[1] = 1;
					connec_order_in1[2] = 4;
					connec_order_in1[3] = 5;
					connec_order_in1[4] = 2;
					connec_order_in1[5] = 3;
					connec_order_in1[6] = 6;
					connec_order_in1[7] = 7;

					connec_order_out1[0] = 0;
					connec_order_out1[1] = 1;
					connec_order_out1[2] = 4;
					connec_order_out1[3] = 5;
					connec_order_out1[4] = 1;
					connec_order_out1[5] = 0;
					connec_order_out1[6] = 5;
					connec_order_out1[7] = 4;

					connec_order_in2[0] = 0;
					connec_order_in2[1] = 1;
					connec_order_in2[2] = 2;
					connec_order_in2[3] = 3;
					connec_order_in2[4] = 4;
					connec_order_in2[5] = 5;
					connec_order_in2[6] = 6;
					connec_order_in2[7] = 7;

					connec_order_out2[0] = 4;
					connec_order_out2[1] = 5;
					connec_order_out2[2] = 3;
					connec_order_out2[3] = 2;
					connec_order_out2[4] = 4;
					connec_order_out2[5] = 5;
					connec_order_out2[6] = 6;
					connec_order_out2[7] = 7;

					connec_order_in3[0] = 2;
					connec_order_in3[1] = 3;
					connec_order_in3[2] = 0;
					connec_order_in3[3] = 1;
					connec_order_in3[4] = 4;
					connec_order_in3[5] = 5;
					connec_order_in3[6] = 6;
					connec_order_in3[7] = 7;

					connec_order_out3[0] = 2;
					connec_order_out3[1] = 3;
					connec_order_out3[2] = 0;
					connec_order_out3[3] = 1;
					connec_order_out3[4] = 4;
					connec_order_out3[5] = 5;
					connec_order_out3[6] = 3;
					connec_order_out3[7] = 2;

					elem->pad = 2;

				}//Edge 3
				else if ((edge_list[0]) && (!edge_list[1]) && (edge_list[2]) && (!edge_list[3]) &&
						(edge_list[4]) && (!edge_list[5]) && (!edge_list[6]) && (edge_list[7]) &&
						(!edge_list[8]) && (!edge_list[9]) && (!edge_list[10]) && (!edge_list[11])&&
						(edge_list_s[0]) && (edge_list_s[2])) {
					//printf("Entrou na 3!\n");

					cut_edge[0] = 0;
					cut_edge[1] = 2;
					cut_edge[2] = 4;
					cut_edge[3] = 7;

					cut_edge_s[0] = 0;
					cut_edge_s[1] = 2;

					connec_order_in1[0] = 1;
					connec_order_in1[1] = 2;
					connec_order_in1[2] = 5;
					connec_order_in1[3] = 6;
					connec_order_in1[4] = 0;
					connec_order_in1[5] = 3;
					connec_order_in1[6] = 4;
					connec_order_in1[7] = 7;

					connec_order_out1[0] = 1;
					connec_order_out1[1] = 2;
					connec_order_out1[2] = 5;
					connec_order_out1[3] = 6;
					connec_order_out1[4] = 0;
					connec_order_out1[5] = 1;
					connec_order_out1[6] = 4;
					connec_order_out1[7] = 5;

					connec_order_in2[0] = 0;
					connec_order_in2[1] = 1;
					connec_order_in2[2] = 2;
					connec_order_in2[3] = 3;
					connec_order_in2[4] = 4;
					connec_order_in2[5] = 5;
					connec_order_in2[6] = 6;
					connec_order_in2[7] = 7;

					connec_order_out2[0] = 2;
					connec_order_out2[1] = 4;
					connec_order_out2[2] = 5;
					connec_order_out2[3] = 3;
					connec_order_out2[4] = 4;
					connec_order_out2[5] = 5;
					connec_order_out2[6] = 6;
					connec_order_out2[7] = 7;

					connec_order_in3[0] = 0;
					connec_order_in3[1] = 3;
					connec_order_in3[2] = 1;
					connec_order_in3[3] = 2;
					connec_order_in3[4] = 4;
					connec_order_in3[5] = 5;
					connec_order_in3[6] = 6;
					connec_order_in3[7] = 7;

					connec_order_out3[0] = 0;
					connec_order_out3[1] = 3;
					connec_order_out3[2] = 0;
					connec_order_out3[3] = 1;
					connec_order_out3[4] = 2;
					connec_order_out3[5] = 4;
					connec_order_out3[6] = 5;
					connec_order_out3[7] = 3;

					elem->pad = 2;

				}//Edge 4
				else if ((edge_list[0]) && (!edge_list[1]) && (!edge_list[2]) && (edge_list[3]) &&
						(!edge_list[4]) && (!edge_list[5]) && (!edge_list[6]) && (!edge_list[7]) &&
						(edge_list[8]) && (!edge_list[9]) && (!edge_list[10]) && (edge_list[11])&&
						(edge_list_s[8]) && (edge_list_s[10])) {
					//printf("Entrou na 4!\n");

					cut_edge[0] = 0;
					cut_edge[1] = 3;
					cut_edge[2] = 11;
					cut_edge[3] = 8;

					cut_edge_s[0] = 10;
					cut_edge_s[1] = 8;

					connec_order_in1[0] = 1;
					connec_order_in1[1] = 2;
					connec_order_in1[2] = 5;
					connec_order_in1[3] = 6;
					connec_order_in1[4] = 0;
					connec_order_in1[5] = 3;
					connec_order_in1[6] = 4;
					connec_order_in1[7] = 7;

					connec_order_out1[0] = 1;
					connec_order_out1[1] = 2;
					connec_order_out1[2] = 5;
					connec_order_out1[3] = 6;
					connec_order_out1[4] = 0;
					connec_order_out1[5] = 4;
					connec_order_out1[6] = 3;
					connec_order_out1[7] = 5;

					connec_order_in2[0] = 0;
					connec_order_in2[1] = 1;
					connec_order_in2[2] = 4;
					connec_order_in2[3] = 5;
					connec_order_in2[4] = 2;
					connec_order_in2[5] = 3;
					connec_order_in2[6] = 6;
					connec_order_in2[7] = 7;

					connec_order_out2[0] = 1;
					connec_order_out2[1] = 4;
					connec_order_out2[2] = 2;
					connec_order_out2[3] = 5;
					connec_order_out2[4] = 2;
					connec_order_out2[5] = 3;
					connec_order_out2[6] = 6;
					connec_order_out2[7] = 7;

					connec_order_in3[0] = 0;
					connec_order_in3[1] = 4;
					connec_order_in3[2] = 1;
					connec_order_in3[3] = 2;
					connec_order_in3[4] = 3;
					connec_order_in3[5] = 5;
					connec_order_in3[6] = 6;
					connec_order_in3[7] = 7;

					connec_order_out3[0] = 0;
					connec_order_out3[1] = 4;
					connec_order_out3[2] = 0;
					connec_order_out3[3] = 4;
					connec_order_out3[4] = 1;
					connec_order_out3[5] = 3;
					connec_order_out3[6] = 5;
					connec_order_out3[7] = 2;

					elem->pad = 2;

				}//Edge 5
				else if ((edge_list[0]) && (edge_list[1]) && (!edge_list[2]) && (!edge_list[3]) &&
						(!edge_list[4]) && (!edge_list[5]) && (!edge_list[6]) && (!edge_list[7]) &&
						(edge_list[8]) && (edge_list[9]) && (!edge_list[10]) && (!edge_list[11])&&
						(edge_list_s[9]) && (edge_list_s[11])) {
					//printf("Entrou na 5!\n");

					cut_edge[0] = 1;
					cut_edge[1] = 0;
					cut_edge[2] = 8;
					cut_edge[3] = 9;

					cut_edge_s[0] = 11;
					cut_edge_s[1] = 9;

					connec_order_in1[0] = 0;
					connec_order_in1[1] = 3;
					connec_order_in1[2] = 4;
					connec_order_in1[3] = 7;
					connec_order_in1[4] = 1;
					connec_order_in1[5] = 2;
					connec_order_in1[6] = 5;
					connec_order_in1[7] = 6;

					connec_order_out1[0] = 0;
					connec_order_out1[1] = 3;
					connec_order_out1[2] = 4;
					connec_order_out1[3] = 7;
					connec_order_out1[4] = 1;
					connec_order_out1[5] = 4;
					connec_order_out1[6] = 2;
					connec_order_out1[7] = 5;

					connec_order_in2[0] = 0;
					connec_order_in2[1] = 1;
					connec_order_in2[2] = 4;
					connec_order_in2[3] = 5;
					connec_order_in2[4] = 2;
					connec_order_in2[5] = 3;
					connec_order_in2[6] = 6;
					connec_order_in2[7] = 7;

					connec_order_out2[0] = 4;
					connec_order_out2[1] = 0;
					connec_order_out2[2] = 5;
					connec_order_out2[3] = 3;
					connec_order_out2[4] = 2;
					connec_order_out2[5] = 3;
					connec_order_out2[6] = 6;
					connec_order_out2[7] = 7;

					connec_order_in3[0] = 1;
					connec_order_in3[1] = 5;
					connec_order_in3[2] = 0;
					connec_order_in3[3] = 2;
					connec_order_in3[4] = 3;
					connec_order_in3[5] = 4;
					connec_order_in3[6] = 6;
					connec_order_in3[7] = 7;

					connec_order_out3[0] = 1;
					connec_order_out3[1] = 5;
					connec_order_out3[2] = 1;
					connec_order_out3[3] = 0;
					connec_order_out3[4] = 4;
					connec_order_out3[5] = 2;
					connec_order_out3[6] = 3;
					connec_order_out3[7] = 5;

					elem->pad = 2;

				}//Edge 6
				else if ((!edge_list[0]) && (edge_list[1]) && (edge_list[2]) && (!edge_list[3]) &&
						(!edge_list[4]) && (!edge_list[5]) && (!edge_list[6]) && (!edge_list[7]) &&
						(!edge_list[8]) && (edge_list[9]) && (edge_list[10]) && (!edge_list[11])&&
						(edge_list_s[8]) && (edge_list_s[10])) {
					//printf("Entrou na 6!\n");

					cut_edge[0] = 2;
					cut_edge[1] = 10;
					cut_edge[2] = 1;
					cut_edge[3] = 9;

					cut_edge_s[0] = 10;
					cut_edge_s[1] = 8;

					connec_order_in1[0] = 0;
					connec_order_in1[1] = 1;
					connec_order_in1[2] = 4;
					connec_order_in1[3] = 5;
					connec_order_in1[4] = 2;
					connec_order_in1[5] = 3;
					connec_order_in1[6] = 6;
					connec_order_in1[7] = 7;

					connec_order_out1[0] = 0;
					connec_order_out1[1] = 1;
					connec_order_out1[2] = 4;
					connec_order_out1[3] = 5;
					connec_order_out1[4] = 2;
					connec_order_out1[5] = 4;
					connec_order_out1[6] = 3;
					connec_order_out1[7] = 5;

					connec_order_in2[0] = 1;
					connec_order_in2[1] = 2;
					connec_order_in2[2] = 5;
					connec_order_in2[3] = 6;
					connec_order_in2[4] = 0;
					connec_order_in2[5] = 3;
					connec_order_in2[6] = 4;
					connec_order_in2[7] = 7;

					connec_order_out2[0] = 4;
					connec_order_out2[1] = 0;
					connec_order_out2[2] = 5;
					connec_order_out2[3] = 1;
					connec_order_out2[4] = 0;
					connec_order_out2[5] = 3;
					connec_order_out2[6] = 4;
					connec_order_out2[7] = 7;

					connec_order_in3[0] = 2;
					connec_order_in3[1] = 6;
					connec_order_in3[2] = 0;
					connec_order_in3[3] = 1;
					connec_order_in3[4] = 3;
					connec_order_in3[5] = 4;
					connec_order_in3[6] = 5;
					connec_order_in3[7] = 7;

					connec_order_out3[0] = 2;
					connec_order_out3[1] = 6;
					connec_order_out3[2] = 4;
					connec_order_out3[3] = 2;
					connec_order_out3[4] = 0;
					connec_order_out3[5] = 5;
					connec_order_out3[6] = 3;
					connec_order_out3[7] = 1;

					elem->pad = 2;

				}//Edge 7
				else if ((!edge_list[0]) && (!edge_list[1]) && (edge_list[2]) && (edge_list[3]) &&
						(!edge_list[4]) && (!edge_list[5]) && (!edge_list[6]) && (!edge_list[7]) &&
						(!edge_list[8]) && (!edge_list[9]) && (edge_list[10]) && (edge_list[11])&&
						(edge_list_s[9]) && (edge_list_s[11])) {
					//printf("Entrou na 7!\n");

					cut_edge[0] = 3;
					cut_edge[1] = 11;
					cut_edge[2] = 2;
					cut_edge[3] = 10;

					cut_edge_s[0] = 11;
					cut_edge_s[1] = 9;

					connec_order_in1[0] = 0;
					connec_order_in1[1] = 1;
					connec_order_in1[2] = 4;
					connec_order_in1[3] = 5;
					connec_order_in1[4] = 2;
					connec_order_in1[5] = 3;
					connec_order_in1[6] = 6;
					connec_order_in1[7] = 7;

					connec_order_out1[0] = 0;
					connec_order_out1[1] = 1;
					connec_order_out1[2] = 4;
					connec_order_out1[3] = 5;
					connec_order_out1[4] = 4;
					connec_order_out1[5] = 0;
					connec_order_out1[6] = 5;
					connec_order_out1[7] = 1;

					connec_order_in2[0] = 0;
					connec_order_in2[1] = 3;
					connec_order_in2[2] = 4;
					connec_order_in2[3] = 7;
					connec_order_in2[4] = 1;
					connec_order_in2[5] = 2;
					connec_order_in2[6] = 5;
					connec_order_in2[7] = 6;

					connec_order_out2[0] = 4;
					connec_order_out2[1] = 2;
					connec_order_out2[2] = 5;
					connec_order_out2[3] = 3;
					connec_order_out2[4] = 1;
					connec_order_out2[5] = 2;
					connec_order_out2[6] = 5;
					connec_order_out2[7] = 6;

					connec_order_in3[0] = 3;
					connec_order_in3[1] = 7;
					connec_order_in3[2] = 0;
					connec_order_in3[3] = 1;
					connec_order_in3[4] = 2;
					connec_order_in3[5] = 4;
					connec_order_in3[6] = 5;
					connec_order_in3[7] = 6;

					connec_order_out3[0] = 3;
					connec_order_out3[1] = 7;
					connec_order_out3[2] = 0;
					connec_order_out3[3] = 4;
					connec_order_out3[4] = 2;
					connec_order_out3[5] = 1;
					connec_order_out3[6] = 5;
					connec_order_out3[7] = 3;

					elem->pad = 2;

				}//Edge 8
				else if ((!edge_list[0]) && (!edge_list[1]) && (!edge_list[2]) && (!edge_list[3]) &&
						(edge_list[4]) && (edge_list[5]) && (!edge_list[6]) && (!edge_list[7]) &&
						(!edge_list[8]) && (edge_list[9]) && (!edge_list[10]) && (edge_list[11])&&
						(edge_list_s[5]) && (edge_list_s[7])) {
					//printf("Entrou na 8!\n");

					cut_edge[0] = 11;
					cut_edge[1] = 9;
					cut_edge[2] = 4;
					cut_edge[3] = 5;

					cut_edge_s[0] = 5;
					cut_edge_s[1] = 7;

					connec_order_in1[0] = 2;
					connec_order_in1[1] = 3;
					connec_order_in1[2] = 6;
					connec_order_in1[3] = 7;
					connec_order_in1[4] = 0;
					connec_order_in1[5] = 1;
					connec_order_in1[6] = 4;
					connec_order_in1[7] = 5;

					connec_order_out1[0] = 2;
					connec_order_out1[1] = 3;
					connec_order_out1[2] = 6;
					connec_order_out1[3] = 7;
					connec_order_out1[4] = 4;
					connec_order_out1[5] = 5;
					connec_order_out1[6] = 0;
					connec_order_out1[7] = 1;

					connec_order_in2[0] = 4;
					connec_order_in2[1] = 5;
					connec_order_in2[2] = 6;
					connec_order_in2[3] = 7;
					connec_order_in2[4] = 0;
					connec_order_in2[5] = 1;
					connec_order_in2[6] = 2;
					connec_order_in2[7] = 3;

					connec_order_out2[0] = 2;
					connec_order_out2[1] = 3;
					connec_order_out2[2] = 5;
					connec_order_out2[3] = 4;
					connec_order_out2[4] = 0;
					connec_order_out2[5] = 1;
					connec_order_out2[6] = 2;
					connec_order_out2[7] = 3;

					connec_order_in3[0] = 4;
					connec_order_in3[1] = 5;
					connec_order_in3[2] = 0;
					connec_order_in3[3] = 1;
					connec_order_in3[4] = 2;
					connec_order_in3[5] = 3;
					connec_order_in3[6] = 6;
					connec_order_in3[7] = 7;

					connec_order_out3[0] = 4;
					connec_order_out3[1] = 5;
					connec_order_out3[2] = 2;
					connec_order_out3[3] = 3;
					connec_order_out3[4] = 5;
					connec_order_out3[5] = 4;
					connec_order_out3[6] = 1;
					connec_order_out3[7] = 0;

					elem->pad = 2;

				}//Edge 9
				else if ((!edge_list[0]) && (!edge_list[1]) && (!edge_list[2]) && (!edge_list[3]) &&
						(!edge_list[4]) && (edge_list[5]) && (edge_list[6]) && (!edge_list[7]) &&
						(edge_list[8]) && (!edge_list[9]) && (edge_list[10]) && (!edge_list[11])&&
						(edge_list_s[0]) && (edge_list_s[2])) {
					//printf("Entrou na 9!\n");

					cut_edge[0] = 8;
					cut_edge[1] = 10;
					cut_edge[2] = 5;
					cut_edge[3] = 6;

					cut_edge_s[0] = 0;
					cut_edge_s[1] = 2;

					connec_order_in1[0] = 0;
					connec_order_in1[1] = 3;
					connec_order_in1[2] = 4;
					connec_order_in1[3] = 7;
					connec_order_in1[4] = 1;
					connec_order_in1[5] = 2;
					connec_order_in1[6] = 5;
					connec_order_in1[7] = 6;

					connec_order_out1[0] = 0;
					connec_order_out1[1] = 3;
					connec_order_out1[2] = 4;
					connec_order_out1[3] = 7;
					connec_order_out1[4] = 4;
					connec_order_out1[5] = 5;
					connec_order_out1[6] = 0;
					connec_order_out1[7] = 1;

					connec_order_in2[0] = 4;
					connec_order_in2[1] = 5;
					connec_order_in2[2] = 6;
					connec_order_in2[3] = 7;
					connec_order_in2[4] = 0;
					connec_order_in2[5] = 1;
					connec_order_in2[6] = 2;
					connec_order_in2[7] = 3;

					connec_order_out2[0] = 4;
					connec_order_out2[1] = 2;
					connec_order_out2[2] = 3;
					connec_order_out2[3] = 5;
					connec_order_out2[4] = 0;
					connec_order_out2[5] = 1;
					connec_order_out2[6] = 2;
					connec_order_out2[7] = 3;

					connec_order_in3[0] = 5;
					connec_order_in3[1] = 6;
					connec_order_in3[2] = 0;
					connec_order_in3[3] = 1;
					connec_order_in3[4] = 2;
					connec_order_in3[5] = 3;
					connec_order_in3[6] = 4;
					connec_order_in3[7] = 7;

					connec_order_out3[0] = 5;
					connec_order_out3[1] = 6;
					connec_order_out3[2] = 4;
					connec_order_out3[3] = 2;
					connec_order_out3[4] = 3;
					connec_order_out3[5] = 5;
					connec_order_out3[6] = 0;
					connec_order_out3[7] = 1;

					elem->pad = 2;

				}//Edge 10
				else if ((!edge_list[0]) && (!edge_list[1]) && (!edge_list[2]) && (!edge_list[3]) &&
						(!edge_list[4]) && (!edge_list[5]) && (edge_list[6]) && (edge_list[7]) &&
						(!edge_list[8]) && (edge_list[9]) && (!edge_list[10]) && (edge_list[11])&&
						(edge_list_s[4]) && (edge_list_s[6])) {
					//printf("Entrou na 10!\n");

					cut_edge[0] = 9;
					cut_edge[1] = 11;
					cut_edge[2] = 6;
					cut_edge[3] = 7;

					cut_edge_s[0] = 6;
					cut_edge_s[1] = 4;

					connec_order_in1[0] = 0;
					connec_order_in1[1] = 1;
					connec_order_in1[2] = 4;
					connec_order_in1[3] = 5;
					connec_order_in1[4] = 2;
					connec_order_in1[5] = 3;
					connec_order_in1[6] = 6;
					connec_order_in1[7] = 7;

					connec_order_out1[0] = 0;
					connec_order_out1[1] = 1;
					connec_order_out1[2] = 4;
					connec_order_out1[3] = 5;
					connec_order_out1[4] = 4;
					connec_order_out1[5] = 5;
					connec_order_out1[6] = 0;
					connec_order_out1[7] = 1;

					connec_order_in2[0] = 4;
					connec_order_in2[1] = 5;
					connec_order_in2[2] = 6;
					connec_order_in2[3] = 7;
					connec_order_in2[4] = 0;
					connec_order_in2[5] = 1;
					connec_order_in2[6] = 2;
					connec_order_in2[7] = 3;

					connec_order_out2[0] = 5;
					connec_order_out2[1] = 4;
					connec_order_out2[2] = 2;
					connec_order_out2[3] = 3;
					connec_order_out2[4] = 0;
					connec_order_out2[5] = 1;
					connec_order_out2[6] = 2;
					connec_order_out2[7] = 3;

					connec_order_in3[0] = 6;
					connec_order_in3[1] = 7;
					connec_order_in3[2] = 0;
					connec_order_in3[3] = 1;
					connec_order_in3[4] = 2;
					connec_order_in3[5] = 3;
					connec_order_in3[6] = 4;
					connec_order_in3[7] = 5;

					connec_order_out3[0] = 6;
					connec_order_out3[1] = 7;
					connec_order_out3[2] = 5;
					connec_order_out3[3] = 4;
					connec_order_out3[4] = 2;
					connec_order_out3[5] = 3;
					connec_order_out3[6] = 1;
					connec_order_out3[7] = 0;

					elem->pad = 2;

				}//Edge 11
				else if ((!edge_list[0]) && (!edge_list[1]) && (!edge_list[2]) && (!edge_list[3]) &&
						(edge_list[4]) && (!edge_list[5]) && (!edge_list[6]) && (edge_list[7]) &&
						(edge_list[8]) && (!edge_list[9]) && (edge_list[10]) && (!edge_list[11])&&
						(edge_list_s[1]) && (edge_list_s[3])) {
					//printf("Entrou na 11!\n");

					cut_edge[0] = 8;
					cut_edge[1] = 10;
					cut_edge[2] = 4;
					cut_edge[3] = 7;

					cut_edge_s[0] = 1;
					cut_edge_s[1] = 3;

					connec_order_in1[0] = 1;
					connec_order_in1[1] = 2;
					connec_order_in1[2] = 5;
					connec_order_in1[3] = 6;
					connec_order_in1[4] = 0;
					connec_order_in1[5] = 3;
					connec_order_in1[6] = 4;
					connec_order_in1[7] = 7;

					connec_order_out1[0] = 1;
					connec_order_out1[1] = 2;
					connec_order_out1[2] = 5;
					connec_order_out1[3] = 6;
					connec_order_out1[4] = 4;
					connec_order_out1[5] = 5;
					connec_order_out1[6] = 0;
					connec_order_out1[7] = 1;

					connec_order_in2[0] = 4;
					connec_order_in2[1] = 5;
					connec_order_in2[2] = 6;
					connec_order_in2[3] = 7;
					connec_order_in2[4] = 0;
					connec_order_in2[5] = 1;
					connec_order_in2[6] = 2;
					connec_order_in2[7] = 3;

					connec_order_out2[0] = 2;
					connec_order_out2[1] = 4;
					connec_order_out2[2] = 5;
					connec_order_out2[3] = 3;
					connec_order_out2[4] = 0;
					connec_order_out2[5] = 1;
					connec_order_out2[6] = 2;
					connec_order_out2[7] = 3;

					connec_order_in3[0] = 4;
					connec_order_in3[1] = 7;
					connec_order_in3[2] = 0;
					connec_order_in3[3] = 1;
					connec_order_in3[4] = 2;
					connec_order_in3[5] = 3;
					connec_order_in3[6] = 5;
					connec_order_in3[7] = 6;

					connec_order_out3[0] = 4;
					connec_order_out3[1] = 7;
					connec_order_out3[2] = 2;
					connec_order_out3[3] = 4;
					connec_order_out3[4] = 5;
					connec_order_out3[5] = 3;
					connec_order_out3[6] = 0;
					connec_order_out3[7] = 1;

					elem->pad = 2;

				} else {
					elem->pad = -1;
					printf("Warning: Bad Element: %d, case 2\n",elements_ids[iel]);
					printf("Edges:  Elem. %d: %d %d %d %d %d %d %d %d %d %d %d %d %d\n", elements_ids[iel], edge_list[0],
							edge_list[1], edge_list[2], edge_list[3],
							edge_list[4], edge_list[5], edge_list[6],
							edge_list[7], edge_list[8], edge_list[9],
							edge_list[10], edge_list[11], ed_cont);
					continue;
				}

				GtsPoint *p0 = point[cut_edge[0]];
				GtsPoint *p1 = point[cut_edge[1]];
				GtsPoint *p2 = point[cut_edge[2]];
				GtsPoint *p3 = point[cut_edge[3]];

				GtsPoint *p4 = point_s[cut_edge_s[0]];
				GtsPoint *p5 = point_s[cut_edge_s[1]];

				//g_assert(p0 != NULL);
				//g_assert(p1 != NULL);
				//g_assert(p2 != NULL);
				//g_assert(p3 != NULL);
				//g_assert(p4 != NULL);
				//g_assert(p5 != NULL);

				if(p0==NULL) exit(1);
				if(p1==NULL) exit(2);
				if(p2==NULL) exit(3);
				if(p3==NULL) exit(4);
				if(p4==NULL) exit(5);
				if(p5==NULL) exit(6);

				conn_t2[0] = AddPointOnEdge(Edge2GNode[cut_edge[0]], hash_nodes, mesh->local_n_nodes, p0, coords);
				conn_t2[1] = AddPointOnEdge(Edge2GNode[cut_edge[1]], hash_nodes, mesh->local_n_nodes, p1, coords);
				conn_t2[2] = AddPointOnEdge(Edge2GNode[cut_edge[2]], hash_nodes, mesh->local_n_nodes, p2, coords);
				conn_t2[3] = AddPointOnEdge(Edge2GNode[cut_edge[3]], hash_nodes, mesh->local_n_nodes, p3, coords);

				// add 2 extra points in the surface
				conn_t2[4] = AddPointOnEdge(Edge2GNode_s[cut_edge_s[0]], hash_nodes, mesh->local_n_nodes, p4, coords);
				conn_t2[5] = AddPointOnEdge(Edge2GNode_s[cut_edge_s[1]], hash_nodes, mesh->local_n_nodes, p5, coords);

				octant_t *elem1 = (octant_t*) sc_array_index(&mesh->elements, elements_ids[iel]);

				elem1->nodes[connec_order_in1[0]].id = original_conn[connec_order_out1[0]];
				elem1->nodes[connec_order_in1[1]].id = original_conn[connec_order_out1[1]];
				elem1->nodes[connec_order_in1[2]].id = original_conn[connec_order_out1[2]];
				elem1->nodes[connec_order_in1[3]].id = original_conn[connec_order_out1[3]];

				elem1->nodes[connec_order_in1[4]].id = conn_t2[connec_order_out1[4]];
				elem1->nodes[connec_order_in1[5]].id = conn_t2[connec_order_out1[5]];
				elem1->nodes[connec_order_in1[6]].id = conn_t2[connec_order_out1[6]];
				elem1->nodes[connec_order_in1[7]].id = conn_t2[connec_order_out1[7]];

				octant_t* elem2 = (octant_t*) sc_array_push(&mesh->elements);

				elem2->nodes[connec_order_in2[0]].id = conn_t2[connec_order_out2[0]];
				elem2->nodes[connec_order_in2[1]].id = conn_t2[connec_order_out2[1]];
				elem2->nodes[connec_order_in2[2]].id = conn_t2[connec_order_out2[2]];
				elem2->nodes[connec_order_in2[3]].id = conn_t2[connec_order_out2[3]];

				elem2->nodes[connec_order_in2[4]].id = original_conn[connec_order_out2[4]];
				elem2->nodes[connec_order_in2[5]].id = original_conn[connec_order_out2[5]];
				elem2->nodes[connec_order_in2[6]].id = original_conn[connec_order_out2[6]];
				elem2->nodes[connec_order_in2[7]].id = original_conn[connec_order_out2[7]];

				elem2->pad = elem->pad;
				elem2->level = elem->level;

				octant_t* elem3 = (octant_t*) sc_array_push(&mesh->elements);

				elem3->nodes[connec_order_in3[0]].id = original_conn[connec_order_out3[0]];
				elem3->nodes[connec_order_in3[1]].id = original_conn[connec_order_out3[1]];

				elem3->nodes[connec_order_in3[2]].id = conn_t2[connec_order_out3[2]];
				elem3->nodes[connec_order_in3[3]].id = conn_t2[connec_order_out3[3]];
				elem3->nodes[connec_order_in3[4]].id = conn_t2[connec_order_out3[4]];
				elem3->nodes[connec_order_in3[5]].id = conn_t2[connec_order_out3[5]];
				elem3->nodes[connec_order_in3[6]].id = conn_t2[connec_order_out3[6]];
				elem3->nodes[connec_order_in3[7]].id = conn_t2[connec_order_out3[7]];

				elem3->pad = elem->pad;
				elem3->level = elem->level;
			}
		}else if(elem->pad==3){

			if (n_parallel_faces == 0 && ed_cont == 3) {
				// Check template 3.

				int conn_t2[7];
				int cut_edge[3];
				int cut_edge_s[3];
				int cut_edge_v[1];
				int connec_order_in1[8];
				int connec_order_out1[8];
				int connec_order_in2[8];
				int connec_order_out2[8];
				int connec_order_in3[8];
				int connec_order_out3[8];
				int connec_order_in4[8];
				int connec_order_out4[8];


				//Corner 0
				if ((edge_list[0]) && (!edge_list[1]) && (!edge_list[2]) && (edge_list[3]) &&
						(edge_list[4]) && (!edge_list[5]) && (!edge_list[6]) && (!edge_list[7]) &&
						(!edge_list[8]) && (!edge_list[9]) && (!edge_list[10]) && (!edge_list[11]) &&
						(edge_list_s[4]) && (edge_list_s[0]) && (edge_list_s[10]) &&
						(edge_list_v[0]) ) {
					//printf("Entrou no 0!\n");

					cut_edge[0] = 4;
					cut_edge[1] = 0;
					cut_edge[2] = 3;

					cut_edge_s[0] = 0;
					cut_edge_s[1] = 4;
					cut_edge_s[2] = 10;

					cut_edge_v[0] = 0;

					connec_order_in1[0] = 4;
					connec_order_in1[1] = 5;
					connec_order_in1[2] = 6;
					connec_order_in1[3] = 7;
					connec_order_in1[4] = 0;
					connec_order_in1[5] = 1;
					connec_order_in1[6] = 2;
					connec_order_in1[7] = 3;

					connec_order_out1[0] = 4;
					connec_order_out1[1] = 5;
					connec_order_out1[2] = 6;
					connec_order_out1[3] = 7;
					connec_order_out1[4] = 0;
					connec_order_out1[5] = 3;
					connec_order_out1[6] = 6;
					connec_order_out1[7] = 4;

					connec_order_in2[0] = 0;
					connec_order_in2[1] = 3;
					connec_order_in2[2] = 4;
					connec_order_in2[3] = 7;
					connec_order_in2[4] = 1;
					connec_order_in2[5] = 2;
					connec_order_in2[6] = 5;
					connec_order_in2[7] = 6;

					connec_order_out2[0] = 1;
					connec_order_out2[1] = 5;
					connec_order_out2[2] = 3;
					connec_order_out2[3] = 6;
					connec_order_out2[4] = 1;
					connec_order_out2[5] = 2;
					connec_order_out2[6] = 5;
					connec_order_out2[7] = 6;

					connec_order_in3[0] = 0;
					connec_order_in3[1] = 1;
					connec_order_in3[2] = 4;
					connec_order_in3[3] = 5;
					connec_order_in3[4] = 2;
					connec_order_in3[5] = 3;
					connec_order_in3[6] = 6;
					connec_order_in3[7] = 7;

					connec_order_out3[0] = 2;
					connec_order_out3[1] = 5;
					connec_order_out3[2] = 4;
					connec_order_out3[3] = 6;
					connec_order_out3[4] = 2;
					connec_order_out3[5] = 3;
					connec_order_out3[6] = 6;
					connec_order_out3[7] = 7;

					connec_order_in4[0] = 0;
					connec_order_in4[1] = 1;
					connec_order_in4[2] = 2;
					connec_order_in4[3] = 3;
					connec_order_in4[4] = 4;
					connec_order_in4[5] = 5;
					connec_order_in4[6] = 6;
					connec_order_in4[7] = 7;

					connec_order_out4[0] = 0;
					connec_order_out4[1] = 1;
					connec_order_out4[2] = 5;
					connec_order_out4[3] = 2;
					connec_order_out4[4] = 0;
					connec_order_out4[5] = 3;
					connec_order_out4[6] = 6;
					connec_order_out4[7] = 4;

					elem->pad = 3;

				}

				//Corner 1
				else if ((edge_list[0]) && (edge_list[1]) && (!edge_list[2]) && (!edge_list[3]) &&
						(!edge_list[4]) && (edge_list[5]) && (!edge_list[6]) && (!edge_list[7]) &&
						(!edge_list[8]) && (!edge_list[9]) && (!edge_list[10]) && (!edge_list[11])&&
						(edge_list_s[1]) && (edge_list_s[6]) && (edge_list_s[11])  &&
						(edge_list_v[1]) ) {
					//printf("Entrou no 1!\n");

					cut_edge[0] = 5;
					cut_edge[1] = 1;
					cut_edge[2] = 0;

					cut_edge_s[0] = 6;
					cut_edge_s[1] = 1;
					cut_edge_s[2] = 11;

					cut_edge_v[0] = 1;

					connec_order_in1[0] = 4;
					connec_order_in1[1] = 5;
					connec_order_in1[2] = 6;
					connec_order_in1[3] = 7;
					connec_order_in1[4] = 0;
					connec_order_in1[5] = 1;
					connec_order_in1[6] = 2;
					connec_order_in1[7] = 3;

					connec_order_out1[0] = 4;
					connec_order_out1[1] = 5;
					connec_order_out1[2] = 6;
					connec_order_out1[3] = 7;
					connec_order_out1[4] = 4;
					connec_order_out1[5] = 0;
					connec_order_out1[6] = 3;
					connec_order_out1[7] = 6;

					connec_order_in2[0] = 0;
					connec_order_in2[1] = 1;
					connec_order_in2[2] = 4;
					connec_order_in2[3] = 5;
					connec_order_in2[4] = 2;
					connec_order_in2[5] = 3;
					connec_order_in2[6] = 6;
					connec_order_in2[7] = 7;

					connec_order_out2[0] = 5;
					connec_order_out2[1] = 1;
					connec_order_out2[2] = 6;
					connec_order_out2[3] = 3;
					connec_order_out2[4] = 2;
					connec_order_out2[5] = 3;
					connec_order_out2[6] = 6;
					connec_order_out2[7] = 7;

					connec_order_in3[0] = 1;
					connec_order_in3[1] = 2;
					connec_order_in3[2] = 5;
					connec_order_in3[3] = 6;
					connec_order_in3[4] = 0;
					connec_order_in3[5] = 3;
					connec_order_in3[6] = 4;
					connec_order_in3[7] = 7;

					connec_order_out3[0] = 2;
					connec_order_out3[1] = 5;
					connec_order_out3[2] = 4;
					connec_order_out3[3] = 6;
					connec_order_out3[4] = 0;
					connec_order_out3[5] = 3;
					connec_order_out3[6] = 4;
					connec_order_out3[7] = 7;

					connec_order_in4[0] = 1;
					connec_order_in4[1] = 0;
					connec_order_in4[2] = 2;
					connec_order_in4[3] = 3;
					connec_order_in4[4] = 4;
					connec_order_in4[5] = 5;
					connec_order_in4[6] = 6;
					connec_order_in4[7] = 7;

					connec_order_out4[0] = 1;
					connec_order_out4[1] = 2;
					connec_order_out4[2] = 1;
					connec_order_out4[3] = 5;
					connec_order_out4[4] = 4;
					connec_order_out4[5] = 0;
					connec_order_out4[6] = 3;
					connec_order_out4[7] = 6;

					elem->pad = 3;

				}

				//Corner 2
				else if ((!edge_list[0]) && (edge_list[1]) && (edge_list[2]) && (!edge_list[3]) &&
						(!edge_list[4]) && (!edge_list[5]) && (edge_list[6]) && (!edge_list[7]) &&
						(!edge_list[8]) && (!edge_list[9]) && (!edge_list[10]) && (!edge_list[11]) &&
						(edge_list_s[3]) && (edge_list_s[7]) && (edge_list_s[10])  &&
						(edge_list_v[2]) ) {
					//printf("Entrou no 2!\n");

					cut_edge[0] = 6;
					cut_edge[1] = 2;
					cut_edge[2] = 1;

					cut_edge_s[0] = 3;
					cut_edge_s[1] = 7;
					cut_edge_s[2] = 10;

					cut_edge_v[0] = 2;

					connec_order_in1[0] = 4;
					connec_order_in1[1] = 5;
					connec_order_in1[2] = 6;
					connec_order_in1[3] = 7;
					connec_order_in1[4] = 1;
					connec_order_in1[5] = 2;
					connec_order_in1[6] = 3;
					connec_order_in1[7] = 4;

					connec_order_out1[0] = 4;
					connec_order_out1[1] = 5;
					connec_order_out1[2] = 6;
					connec_order_out1[3] = 7;
					connec_order_out1[4] = 4;
					connec_order_out1[5] = 0;
					connec_order_out1[6] = 3;
					connec_order_out1[7] = 6;

					connec_order_in2[0] = 1;
					connec_order_in2[1] = 2;
					connec_order_in2[2] = 5;
					connec_order_in2[3] = 6;
					connec_order_in2[4] = 0;
					connec_order_in2[5] = 3;
					connec_order_in2[6] = 4;
					connec_order_in2[7] = 7;

					connec_order_out2[0] = 5;
					connec_order_out2[1] = 1;
					connec_order_out2[2] = 6;
					connec_order_out2[3] = 3;
					connec_order_out2[4] = 0;
					connec_order_out2[5] = 3;
					connec_order_out2[6] = 4;
					connec_order_out2[7] = 7;

					connec_order_in3[0] = 2;
					connec_order_in3[1] = 3;
					connec_order_in3[2] = 6;
					connec_order_in3[3] = 7;
					connec_order_in3[4] = 0;
					connec_order_in3[5] = 1;
					connec_order_in3[6] = 4;
					connec_order_in3[7] = 5;

					connec_order_out3[0] = 2;
					connec_order_out3[1] = 5;
					connec_order_out3[2] = 4;
					connec_order_out3[3] = 6;
					connec_order_out3[4] = 0;
					connec_order_out3[5] = 1;
					connec_order_out3[6] = 4;
					connec_order_out3[7] = 5;

					connec_order_in4[0] = 2;
					connec_order_in4[1] = 0;
					connec_order_in4[2] = 1;
					connec_order_in4[3] = 3;
					connec_order_in4[4] = 4;
					connec_order_in4[5] = 5;
					connec_order_in4[6] = 6;
					connec_order_in4[7] = 7;

					connec_order_out4[0] = 2;
					connec_order_out4[1] = 5;
					connec_order_out4[2] = 2;
					connec_order_out4[3] = 1;
					connec_order_out4[4] = 6;
					connec_order_out4[5] = 4;
					connec_order_out4[6] = 0;
					connec_order_out4[7] = 3;

					elem->pad = 3;

				}

				//Corner 3
				else if ((!edge_list[0]) && (!edge_list[1]) && (edge_list[2]) && (edge_list[3]) &&
						(!edge_list[4]) && (!edge_list[5]) && (!edge_list[6]) && (edge_list[7]) &&
						(!edge_list[8]) && (!edge_list[9]) && (!edge_list[10]) && (!edge_list[11]) &&
						(edge_list_s[5]) && (edge_list_s[2]) && (edge_list_s[11])  &&
						(edge_list_v[3]) ) {
					//printf("Entrou no 3!\n");

					cut_edge[0] = 7;
					cut_edge[1] = 3;
					cut_edge[2] = 2;

					cut_edge_s[0] = 5;
					cut_edge_s[1] = 2;
					cut_edge_s[2] = 11;

					cut_edge_v[0] = 3;

					connec_order_in1[0] = 4;
					connec_order_in1[1] = 5;
					connec_order_in1[2] = 6;
					connec_order_in1[3] = 7;
					connec_order_in1[4] = 0;
					connec_order_in1[5] = 1;
					connec_order_in1[6] = 2;
					connec_order_in1[7] = 3;

					connec_order_out1[0] = 4;
					connec_order_out1[1] = 5;
					connec_order_out1[2] = 6;
					connec_order_out1[3] = 7;
					connec_order_out1[4] = 3;
					connec_order_out1[5] = 6;
					connec_order_out1[6] = 4;
					connec_order_out1[7] = 0;

					connec_order_in2[0] = 2;
					connec_order_in2[1] = 3;
					connec_order_in2[2] = 6;
					connec_order_in2[3] = 7;
					connec_order_in2[4] = 0;
					connec_order_in2[5] = 1;
					connec_order_in2[6] = 4;
					connec_order_in2[7] = 5;

					connec_order_out2[0] = 5;
					connec_order_out2[1] = 1;
					connec_order_out2[2] = 6;
					connec_order_out2[3] = 3;
					connec_order_out2[4] = 0;
					connec_order_out2[5] = 1;
					connec_order_out2[6] = 4;
					connec_order_out2[7] = 5;

					connec_order_in3[0] = 0;
					connec_order_in3[1] = 3;
					connec_order_in3[2] = 4;
					connec_order_in3[3] = 7;
					connec_order_in3[4] = 1;
					connec_order_in3[5] = 2;
					connec_order_in3[6] = 5;
					connec_order_in3[7] = 6;

					connec_order_out3[0] = 5;
					connec_order_out3[1] = 2;
					connec_order_out3[2] = 6;
					connec_order_out3[3] = 4;
					connec_order_out3[4] = 1;
					connec_order_out3[5] = 2;
					connec_order_out3[6] = 5;
					connec_order_out3[7] = 6;

					connec_order_in4[0] = 3;
					connec_order_in4[1] = 0;
					connec_order_in4[2] = 1;
					connec_order_in4[3] = 2;
					connec_order_in4[4] = 4;
					connec_order_in4[5] = 5;
					connec_order_in4[6] = 6;
					connec_order_in4[7] = 7;

					connec_order_out4[0] = 3;
					connec_order_out4[1] = 1;
					connec_order_out4[2] = 5;
					connec_order_out4[3] = 2;
					connec_order_out4[4] = 3;
					connec_order_out4[5] = 6;
					connec_order_out4[6] = 4;
					connec_order_out4[7] = 0;

					elem->pad = 3;

				}


				//Corner 4
				else if ((!edge_list[0]) && (!edge_list[1]) && (!edge_list[2]) && (!edge_list[3]) &&
						(edge_list[4]) && (!edge_list[5]) && (!edge_list[6]) && (!edge_list[7]) &&
						(edge_list[8]) && (!edge_list[9]) && (!edge_list[10]) && (edge_list[11]) &&
						(edge_list_s[1]) && (edge_list_s[5]) && (edge_list_s[8])  &&
						(edge_list_v[2]) ) {
					//printf("Entrou no 4!\n");

					cut_edge[0] = 4;
					cut_edge[1] = 8;
					cut_edge[2] = 11;

					cut_edge_s[0] = 1;
					cut_edge_s[1] = 5;
					cut_edge_s[2] = 8;

					cut_edge_v[0] = 2;

					connec_order_in1[0] = 0;
					connec_order_in1[1] = 1;
					connec_order_in1[2] = 2;
					connec_order_in1[3] = 3;
					connec_order_in1[4] = 4;
					connec_order_in1[5] = 5;
					connec_order_in1[6] = 6;
					connec_order_in1[7] = 7;

					connec_order_out1[0] = 0;
					connec_order_out1[1] = 1;
					connec_order_out1[2] = 2;
					connec_order_out1[3] = 3;
					connec_order_out1[4] = 0;
					connec_order_out1[5] = 3;
					connec_order_out1[6] = 6;
					connec_order_out1[7] = 4;

					connec_order_in2[0] = 0;
					connec_order_in2[1] = 3;
					connec_order_in2[2] = 4;
					connec_order_in2[3] = 7;
					connec_order_in2[4] = 1;
					connec_order_in2[5] = 2;
					connec_order_in2[6] = 5;
					connec_order_in2[7] = 6;

					connec_order_out2[0] = 3;
					connec_order_out2[1] = 6;
					connec_order_out2[2] = 1;
					connec_order_out2[3] = 5;
					connec_order_out2[4] = 1;
					connec_order_out2[5] = 2;
					connec_order_out2[6] = 5;
					connec_order_out2[7] = 6;

					connec_order_in3[0] = 0;
					connec_order_in3[1] = 1;
					connec_order_in3[2] = 4;
					connec_order_in3[3] = 5;
					connec_order_in3[4] = 2;
					connec_order_in3[5] = 3;
					connec_order_in3[6] = 6;
					connec_order_in3[7] = 7;

					connec_order_out3[0] = 4;
					connec_order_out3[1] = 6;
					connec_order_out3[2] = 2;
					connec_order_out3[3] = 5;
					connec_order_out3[4] = 2;
					connec_order_out3[5] = 3;
					connec_order_out3[6] = 6;
					connec_order_out3[7] = 7;

					connec_order_in4[0] = 4;
					connec_order_in4[1] = 0;
					connec_order_in4[2] = 1;
					connec_order_in4[3] = 2;
					connec_order_in4[4] = 3;
					connec_order_in4[5] = 5;
					connec_order_in4[6] = 6;
					connec_order_in4[7] = 7;

					connec_order_out4[0] = 4;
					connec_order_out4[1] = 0;
					connec_order_out4[2] = 3;
					connec_order_out4[3] = 6;
					connec_order_out4[4] = 4;
					connec_order_out4[5] = 1;
					connec_order_out4[6] = 5;
					connec_order_out4[7] = 2;

					elem->pad = 3;

				}

				//Corner 5
				else if ((!edge_list[0]) && (!edge_list[1]) && (!edge_list[2]) && (!edge_list[3]) &&
						(!edge_list[4]) && (edge_list[5]) && (!edge_list[6]) && (!edge_list[7]) &&
						(edge_list[8]) && (edge_list[9]) && (!edge_list[10]) && (!edge_list[11]) &&
						(edge_list_s[7]) && (edge_list_s[0]) && (edge_list_s[9])  &&
						(edge_list_v[3]) ) {
					//printf("Entrou no 5!\n");

					cut_edge[0] = 5;
					cut_edge[1] = 9;
					cut_edge[2] = 8;

					cut_edge_s[0] = 7;
					cut_edge_s[1] = 0;
					cut_edge_s[2] = 9;

					cut_edge_v[0] = 3;

					connec_order_in1[0] = 0;
					connec_order_in1[1] = 1;
					connec_order_in1[2] = 2;
					connec_order_in1[3] = 3;
					connec_order_in1[4] = 4;
					connec_order_in1[5] = 5;
					connec_order_in1[6] = 6;
					connec_order_in1[7] = 7;

					connec_order_out1[0] = 0;
					connec_order_out1[1] = 1;
					connec_order_out1[2] = 2;
					connec_order_out1[3] = 3;
					connec_order_out1[4] = 4;
					connec_order_out1[5] = 0;
					connec_order_out1[6] = 3;
					connec_order_out1[7] = 6;

					connec_order_in2[0] = 0;
					connec_order_in2[1] = 1;
					connec_order_in2[2] = 4;
					connec_order_in2[3] = 5;
					connec_order_in2[4] = 2;
					connec_order_in2[5] = 3;
					connec_order_in2[6] = 6;
					connec_order_in2[7] = 7;

					connec_order_out2[0] = 6;
					connec_order_out2[1] = 3;
					connec_order_out2[2] = 5;
					connec_order_out2[3] = 1;
					connec_order_out2[4] = 2;
					connec_order_out2[5] = 3;
					connec_order_out2[6] = 6;
					connec_order_out2[7] = 7;

					connec_order_in3[0] = 1;
					connec_order_in3[1] = 2;
					connec_order_in3[2] = 5;
					connec_order_in3[3] = 6;
					connec_order_in3[4] = 0;
					connec_order_in3[5] = 3;
					connec_order_in3[6] = 4;
					connec_order_in3[7] = 7;

					connec_order_out3[0] = 4;
					connec_order_out3[1] = 6;
					connec_order_out3[2] = 2;
					connec_order_out3[3] = 5;
					connec_order_out3[4] = 0;
					connec_order_out3[5] = 3;
					connec_order_out3[6] = 4;
					connec_order_out3[7] = 7;

					connec_order_in4[0] = 5;
					connec_order_in4[1] = 0;
					connec_order_in4[2] = 1;
					connec_order_in4[3] = 2;
					connec_order_in4[4] = 3;
					connec_order_in4[5] = 4;
					connec_order_in4[6] = 6;
					connec_order_in4[7] = 7;

					connec_order_out4[0] = 5;
					connec_order_out4[1] = 4;
					connec_order_out4[2] = 0;
					connec_order_out4[3] = 3;
					connec_order_out4[4] = 6;
					connec_order_out4[5] = 2;
					connec_order_out4[6] = 1;
					connec_order_out4[7] = 5;

					elem->pad = 3;

				}

				//Corner 6
				else if ((!edge_list[0]) && (!edge_list[1]) && (!edge_list[2]) && (!edge_list[3]) &&
						(!edge_list[4]) && (!edge_list[5]) && (edge_list[6]) && (!edge_list[7]) &&
						(!edge_list[8]) && (edge_list[9]) && (edge_list[10]) && (!edge_list[11]) &&
						(edge_list_s[2]) && (edge_list_s[6]) && (edge_list_s[8])  &&
						(edge_list_v[0]) ) {
					//printf("Entrou no 6!\n");


					cut_edge[0] = 6;
					cut_edge[1] = 10;
					cut_edge[2] = 9;

					cut_edge_s[0] = 2;
					cut_edge_s[1] = 6;
					cut_edge_s[2] = 8;

					cut_edge_v[0] = 0;

					connec_order_in1[0] = 0;
					connec_order_in1[1] = 1;
					connec_order_in1[2] = 2;
					connec_order_in1[3] = 3;
					connec_order_in1[4] = 4;
					connec_order_in1[5] = 5;
					connec_order_in1[6] = 6;
					connec_order_in1[7] = 7;

					connec_order_out1[0] = 0;
					connec_order_out1[1] = 1;
					connec_order_out1[2] = 2;
					connec_order_out1[3] = 3;
					connec_order_out1[4] = 6;
					connec_order_out1[5] = 4;
					connec_order_out1[6] = 0;
					connec_order_out1[7] = 3;

					connec_order_in2[0] = 1;
					connec_order_in2[1] = 2;
					connec_order_in2[2] = 5;
					connec_order_in2[3] = 6;
					connec_order_in2[4] = 0;
					connec_order_in2[5] = 3;
					connec_order_in2[6] = 4;
					connec_order_in2[7] = 7;

					connec_order_out2[0] = 6;
					connec_order_out2[1] = 3;
					connec_order_out2[2] = 5;
					connec_order_out2[3] = 1;
					connec_order_out2[4] = 0;
					connec_order_out2[5] = 3;
					connec_order_out2[6] = 4;
					connec_order_out2[7] = 7;

					connec_order_in3[0] = 2;
					connec_order_in3[1] = 3;
					connec_order_in3[2] = 6;
					connec_order_in3[3] = 7;
					connec_order_in3[4] = 0;
					connec_order_in3[5] = 1;
					connec_order_in3[6] = 4;
					connec_order_in3[7] = 5;

					connec_order_out3[0] = 4;
					connec_order_out3[1] = 6;
					connec_order_out3[2] = 2;
					connec_order_out3[3] = 5;
					connec_order_out3[4] = 0;
					connec_order_out3[5] = 1;
					connec_order_out3[6] = 4;
					connec_order_out3[7] = 5;

					connec_order_in4[0] = 6;
					connec_order_in4[1] = 0;
					connec_order_in4[2] = 1;
					connec_order_in4[3] = 2;
					connec_order_in4[4] = 3;
					connec_order_in4[5] = 4;
					connec_order_in4[6] = 5;
					connec_order_in4[7] = 7;

					connec_order_out4[0] = 6;
					connec_order_out4[1] = 6;
					connec_order_out4[2] = 4;
					connec_order_out4[3] = 0;
					connec_order_out4[4] = 3;
					connec_order_out4[5] = 5;
					connec_order_out4[6] = 2;
					connec_order_out4[7] = 1;

					elem->pad = 3;

				}

				//Corner 7
				else if ((!edge_list[0]) && (!edge_list[1]) && (!edge_list[2]) && (!edge_list[3]) &&
						(!edge_list[4]) && (!edge_list[5]) && (!edge_list[6]) && (edge_list[7]) &&
						(!edge_list[8]) && (!edge_list[9]) && (edge_list[10]) && (edge_list[11]) &&
						(edge_list_s[4]) && (edge_list_s[3]) && (edge_list_s[9])  &&
						(edge_list_v[1]) ) {
					//printf("Entrou no 7!\n");


					cut_edge[0] = 7;
					cut_edge[1] = 11;
					cut_edge[2] = 10;

					cut_edge_s[0] = 4;
					cut_edge_s[1] = 3;
					cut_edge_s[2] = 9;

					cut_edge_v[0] = 1;

					connec_order_in1[0] = 0;
					connec_order_in1[1] = 1;
					connec_order_in1[2] = 2;
					connec_order_in1[3] = 3;
					connec_order_in1[4] = 4;
					connec_order_in1[5] = 5;
					connec_order_in1[6] = 6;
					connec_order_in1[7] = 7;

					connec_order_out1[0] = 0;
					connec_order_out1[1] = 1;
					connec_order_out1[2] = 2;
					connec_order_out1[3] = 3;
					connec_order_out1[4] = 3;
					connec_order_out1[5] = 6;
					connec_order_out1[6] = 4;
					connec_order_out1[7] = 0;

					connec_order_in2[0] = 2;
					connec_order_in2[1] = 3;
					connec_order_in2[2] = 6;
					connec_order_in2[3] = 7;
					connec_order_in2[4] = 0;
					connec_order_in2[5] = 1;
					connec_order_in2[6] = 4;
					connec_order_in2[7] = 5;

					connec_order_out2[0] = 6;
					connec_order_out2[1] = 3;
					connec_order_out2[2] = 5;
					connec_order_out2[3] = 1;
					connec_order_out2[4] = 0;
					connec_order_out2[5] = 1;
					connec_order_out2[6] = 4;
					connec_order_out2[7] = 5;

					connec_order_in3[0] = 0;
					connec_order_in3[1] = 3;
					connec_order_in3[2] = 4;
					connec_order_in3[3] = 7;
					connec_order_in3[4] = 1;
					connec_order_in3[5] = 2;
					connec_order_in3[6] = 5;
					connec_order_in3[7] = 6;

					connec_order_out3[0] = 6;
					connec_order_out3[1] = 4;
					connec_order_out3[2] = 5;
					connec_order_out3[3] = 2;
					connec_order_out3[4] = 1;
					connec_order_out3[5] = 2;
					connec_order_out3[6] = 5;
					connec_order_out3[7] = 6;

					connec_order_in4[0] = 7;
					connec_order_in4[1] = 0;
					connec_order_in4[2] = 1;
					connec_order_in4[3] = 2;
					connec_order_in4[4] = 3;
					connec_order_in4[5] = 4;
					connec_order_in4[6] = 5;
					connec_order_in4[7] = 6;

					connec_order_out4[0] = 7;
					connec_order_out4[1] = 3;
					connec_order_out4[2] = 6;
					connec_order_out4[3] = 4;
					connec_order_out4[4] = 0;
					connec_order_out4[5] = 1;
					connec_order_out4[6] = 5;
					connec_order_out4[7] = 2;

					elem->pad = 3;

				} else {
					elem->pad = -1;
					printf("Warning: Bad Element: %d, case 3\n",elements_ids[iel]);
					printf("Edges:  Elem. %d: %d %d %d %d %d %d %d %d %d %d %d %d %d\n", elements_ids[iel], edge_list[0],
							edge_list[1], edge_list[2], edge_list[3],
							edge_list[4], edge_list[5], edge_list[6],
							edge_list[7], edge_list[8], edge_list[9],
							edge_list[10], edge_list[11], ed_cont);
					continue;
				}

				GtsPoint *p0 = point[cut_edge[0]];
				GtsPoint *p1 = point[cut_edge[1]];
				GtsPoint *p2 = point[cut_edge[2]];

				GtsPoint *p3 = point_s[cut_edge_s[0]];
				GtsPoint *p4 = point_s[cut_edge_s[1]];
				GtsPoint *p5 = point_s[cut_edge_s[2]];

				GtsPoint *p6 = point_v[cut_edge_v[0]];

				//g_assert(p0 != NULL);
				//g_assert(p1 != NULL);
				//g_assert(p2 != NULL);
				//g_assert(p3 != NULL);
				//g_assert(p4 != NULL);
				//g_assert(p5 != NULL);
				//g_assert(p6 != NULL);

				if(p0==NULL) exit(1);
				if(p1==NULL) exit(2);
				if(p2==NULL) exit(3);
				if(p3==NULL) exit(4);
				if(p4==NULL) exit(5);
				if(p5==NULL) exit(6);
				if(p6==NULL) exit(7);

				conn_t2[0] = AddPointOnEdge(Edge2GNode[cut_edge[0]], hash_nodes, mesh->local_n_nodes, p0, coords);
				conn_t2[1] = AddPointOnEdge(Edge2GNode[cut_edge[1]], hash_nodes, mesh->local_n_nodes, p1, coords);
				conn_t2[2] = AddPointOnEdge(Edge2GNode[cut_edge[2]], hash_nodes, mesh->local_n_nodes, p2, coords);

				// add 3 extra points in the surface
				conn_t2[3] = AddPointOnEdge(Edge2GNode_s[cut_edge_s[0]], hash_nodes, mesh->local_n_nodes, p3, coords);
				conn_t2[4] = AddPointOnEdge(Edge2GNode_s[cut_edge_s[1]], hash_nodes, mesh->local_n_nodes, p4, coords);
				conn_t2[5] = AddPointOnEdge(Edge2GNode_s[cut_edge_s[2]], hash_nodes, mesh->local_n_nodes, p5, coords);

				// add 1 extra points in the volume
				conn_t2[6] = AddPointOnEdge(Edge2GNode_s[cut_edge_v[0]], hash_nodes, mesh->local_n_nodes, p6, coords);

				octant_t *elem1 = (octant_t*) sc_array_index(&mesh->elements, elements_ids[iel]);

				elem1->nodes[connec_order_in1[0]].id = original_conn[connec_order_out1[0]];
				elem1->nodes[connec_order_in1[1]].id = original_conn[connec_order_out1[1]];
				elem1->nodes[connec_order_in1[2]].id = original_conn[connec_order_out1[2]];
				elem1->nodes[connec_order_in1[3]].id = original_conn[connec_order_out1[3]];

				elem1->nodes[connec_order_in1[4]].id = conn_t2[connec_order_out1[4]];
				elem1->nodes[connec_order_in1[5]].id = conn_t2[connec_order_out1[5]];
				elem1->nodes[connec_order_in1[6]].id = conn_t2[connec_order_out1[6]];
				elem1->nodes[connec_order_in1[7]].id = conn_t2[connec_order_out1[7]];

				octant_t* elem2 = (octant_t*) sc_array_push(&mesh->elements);

				elem2->nodes[connec_order_in2[0]].id = conn_t2[connec_order_out2[0]];
				elem2->nodes[connec_order_in2[1]].id = conn_t2[connec_order_out2[1]];
				elem2->nodes[connec_order_in2[2]].id = conn_t2[connec_order_out2[2]];
				elem2->nodes[connec_order_in2[3]].id = conn_t2[connec_order_out2[3]];

				elem2->nodes[connec_order_in2[4]].id = original_conn[connec_order_out2[4]];
				elem2->nodes[connec_order_in2[5]].id = original_conn[connec_order_out2[5]];
				elem2->nodes[connec_order_in2[6]].id = original_conn[connec_order_out2[6]];
				elem2->nodes[connec_order_in2[7]].id = original_conn[connec_order_out2[7]];

				elem2->pad = elem->pad;
				elem2->level = elem->level;

				octant_t* elem3 = (octant_t*) sc_array_push(&mesh->elements);

				elem3->nodes[connec_order_in3[0]].id = conn_t2[connec_order_out3[0]];
				elem3->nodes[connec_order_in3[1]].id = conn_t2[connec_order_out3[1]];
				elem3->nodes[connec_order_in3[2]].id = conn_t2[connec_order_out3[2]];
				elem3->nodes[connec_order_in3[3]].id = conn_t2[connec_order_out3[3]];

				elem3->nodes[connec_order_in3[4]].id = original_conn[connec_order_out3[4]];
				elem3->nodes[connec_order_in3[5]].id = original_conn[connec_order_out3[5]];
				elem3->nodes[connec_order_in3[6]].id = original_conn[connec_order_out3[6]];
				elem3->nodes[connec_order_in3[7]].id = original_conn[connec_order_out3[7]];

				elem3->pad = elem->pad;
				elem3->level = elem->level;

				octant_t* elem4 = (octant_t*) sc_array_push(&mesh->elements);

				elem4->nodes[connec_order_in4[0]].id = original_conn[connec_order_out4[0]];

				elem4->nodes[connec_order_in4[1]].id = conn_t2[connec_order_out4[1]];
				elem4->nodes[connec_order_in4[2]].id = conn_t2[connec_order_out4[2]];
				elem4->nodes[connec_order_in4[3]].id = conn_t2[connec_order_out4[3]];
				elem4->nodes[connec_order_in4[4]].id = conn_t2[connec_order_out4[4]];
				elem4->nodes[connec_order_in4[5]].id = conn_t2[connec_order_out4[5]];
				elem4->nodes[connec_order_in4[6]].id = conn_t2[connec_order_out4[6]];
				elem4->nodes[connec_order_in4[7]].id = conn_t2[connec_order_out4[7]];

				elem4->pad = elem->pad;
				elem4->level = elem->level;
			}

		}

#endif

		for (int edge = 0; edge < 4; edge++) {
			if (point_v[edge]) gts_object_destroy(GTS_OBJECT(point_v[edge]));
			point_v[edge] = NULL;
		}

		for (int edge = 0; edge < 12; edge++) {
			if (point_s[edge]) gts_object_destroy(GTS_OBJECT(point_s[edge]));
			point_s[edge] = NULL;
		}

		for (int edge = 0; edge < 12; edge++) {
			if (point[edge]) gts_object_destroy(GTS_OBJECT(point[edge]));
			point[edge] = NULL;
		}

	}

	mesh->local_n_elements = mesh->elements.elem_count;
	//mesh->local_n_nodes = coords.size()/3;

	//MPI_Allreduce(&mesh->local_n_elements, &mesh->total_n_elements, 1, MPI_LONG, MPI_SUM, MPI_COMM_WORLD);
	//MPI_Allreduce(&mesh->local_n_nodes, &mesh->total_n_nodes, 1, MPI_LONG, MPI_SUM, MPI_COMM_WORLD);

	//	if (mesh->mpi_rank == 0) {
	//		printf("Total number of elements: %lld\n", mesh->local_n_elements);
	//		printf("Total number of nodes: %lld\n", mesh->local_n_nodes);
	//	}

	fclose(fdbg);
	sc_hash_array_destroy(hash_nodes);
}
