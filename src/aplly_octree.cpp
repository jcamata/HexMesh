
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
	a += (uint32_t) q->coord[2];
	sc_hash_final(a, b, c);
	return (unsigned) c;
}

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

void CopyPropEl(hexa_tree_t* mesh, int id, octant_t *elem1){

	octant_t *elem = (octant_t*) sc_array_index(&mesh->elements, id);

	elem1->level = elem->level+1;
	elem1->tem = elem->tem;
	elem1->pad = elem->pad;
	elem1->n_mat = elem->n_mat;
	elem1->pml_id = elem->pml_id;

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

	octant_t *elem = (octant_t*) sc_array_index(&mesh->elements, elements_ids[iel]);

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
				conn_p[ii] = AddPoint(var, hash_nodes, mesh->local_n_nodes, point[ii] , coords);
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

}

void ApplyTemplate2(hexa_tree_t* mesh, std::vector<double>& coords, std::vector<int>& elements_ids, int iel, sc_hash_array_t* hash_nodes, double step){

	octant_t *elem = (octant_t*) sc_array_index(&mesh->elements, elements_ids[iel]);

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
				conn_p[ii] = AddPoint(var, hash_nodes, mesh->local_n_nodes, point[ii] , coords);
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
}

void ApplyTemplate3(hexa_tree_t* mesh, std::vector<double>& coords, std::vector<int>& elements_ids, int iel, sc_hash_array_t* hash_nodes, double step){

	octant_t *elem = (octant_t*) sc_array_index(&mesh->elements, elements_ids[iel]);

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
				conn_p[ii] = AddPoint(var, hash_nodes, mesh->local_n_nodes, point[ii] , coords);
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
}

void ApplyTemplate4(hexa_tree_t* mesh, std::vector<double>& coords, std::vector<int>& elements_ids, int iel, sc_hash_array_t* hash_nodes, double step){

	octant_t *elem = (octant_t*) sc_array_index(&mesh->elements, elements_ids[iel]);

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
				double var[3];
				var[0] = point[ii]->x;
				var[1] = point[ii]->y;
				var[2] = point[ii]->z;
				conn_p[ii] = AddPoint(var, hash_nodes, mesh->local_n_nodes, point[ii] , coords);
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
}

void ApplyTemplate5(hexa_tree_t* mesh, std::vector<double>& coords, std::vector<int>& elements_ids, int iel, sc_hash_array_t* hash_nodes, double step){

	octant_t *elem = (octant_t*) sc_array_index(&mesh->elements, elements_ids[iel]);

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
				conn_p[ii] = AddPoint(var, hash_nodes, mesh->local_n_nodes, point[ii] , coords);
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

}

void ApplyTemplate6(hexa_tree_t* mesh, std::vector<double>& coords, std::vector<int>& elements_ids, int iel, sc_hash_array_t* hash_nodes, double step){

	octant_t *elem = (octant_t*) sc_array_index(&mesh->elements, elements_ids[iel]);

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
				conn_p[ii] = AddPoint(var, hash_nodes, mesh->local_n_nodes, point[ii] , coords);
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
}

void ApplyTemplate7(hexa_tree_t* mesh, std::vector<double>& coords, std::vector<int>& elements_ids, int iel, sc_hash_array_t* hash_nodes, double step){

	octant_t *elem = (octant_t*) sc_array_index(&mesh->elements, elements_ids[iel]);

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
				conn_p[ii] = AddPoint(var, hash_nodes, mesh->local_n_nodes, point[ii] , coords);
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

}

void ApplyTemplate8(hexa_tree_t* mesh, std::vector<double>& coords, std::vector<int>& elements_ids, int iel, sc_hash_array_t* hash_nodes, double step){

	octant_t *elem = (octant_t*) sc_array_index(&mesh->elements, elements_ids[iel]);

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
				conn_p[ii] = AddPoint(var, hash_nodes, mesh->local_n_nodes, point[ii] , coords);
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
}

void ApplyTemplate9(hexa_tree_t* mesh, std::vector<double>& coords, std::vector<int>& elements_ids, int iel, sc_hash_array_t* hash_nodes, double step){

	octant_t *elem = (octant_t*) sc_array_index(&mesh->elements, elements_ids[iel]);

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
				conn_p[ii] = AddPoint(var, hash_nodes, mesh->local_n_nodes, point[ii] , coords);
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
}

void ApplyTemplate10(hexa_tree_t* mesh, std::vector<double>& coords, std::vector<int>& elements_ids, int iel, sc_hash_array_t* hash_nodes, double step){

	octant_t *elem = (octant_t*) sc_array_index(&mesh->elements, elements_ids[iel]);

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
				conn_p[ii] = AddPoint(var, hash_nodes, mesh->local_n_nodes, point[ii] , coords);
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
}

void ApplyTemplate11(hexa_tree_t* mesh, std::vector<double>& coords, std::vector<int>& elements_ids, int iel, sc_hash_array_t* hash_nodes, double step){

	octant_t *elem = (octant_t*) sc_array_index(&mesh->elements, elements_ids[iel]);

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
					conn_p[i*16+ii*4+iii] = AddPoint(var, hash_nodes, mesh->local_n_nodes, point[i*16+ii*4+iii] , coords);
				}

				//fprintf(mesh->fdbg,"id do no: %d\n",conn_p[i*16+ii*4+iii]);
				double xxx = coords[3*conn_p[i*16+ii*4+iii]];
				double yyy = coords[3*conn_p[i*16+ii*4+iii]+1];
				double zzz = coords[3*conn_p[i*16+ii*4+iii]+2];
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
}
void ApplyOctreeTemplate(hexa_tree_t* mesh, std::vector<double>& coords, std::vector<int>& elements_ids) {

	bool clamped = true;
	sc_hash_array_t* hash_nodes = sc_hash_array_new(sizeof(node_t), edge_hash_fn, edge_equal_fn, &clamped);

	for(int n = 0;n<mesh->nodes.elem_count;n++){
		size_t position;
		node_t *r;
		node_t key;
		octant_node_t* node = (octant_node_t*) sc_array_index (&mesh->nodes, n);

		key.coord[0] = coords[node->id];
		key.coord[1] = coords[node->id+1];
		key.coord[2] = coords[node->id+2];
		key.node_id = node->id;

		r = (node_t*) sc_hash_array_insert_unique(hash_nodes, &key, &position);
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
	MPI_Allreduce(&mesh->local_n_elements, &mesh->total_n_elements, 1, MPI_LONG, MPI_SUM, MPI_COMM_WORLD);
	MPI_Allreduce(&mesh->local_n_nodes, &mesh->total_n_nodes, 1, MPI_LONG, MPI_SUM, MPI_COMM_WORLD);
	mesh->nodes.elem_count =  mesh->local_n_nodes;

	if (mesh->mpi_rank == 0) {
		printf("Total number of elements: %lld\n", mesh->local_n_elements);
		printf("Total number of nodes: %lld\n", mesh->local_n_nodes);
	}

	//redefine the connectivity data for write the vtk file
	mesh->part_nodes = NULL;
	mesh->part_nodes = (int*) malloc (mesh->local_n_nodes*sizeof(int));
	for(int i =0; i < mesh->local_n_nodes; i++)
		mesh->part_nodes[i] = mesh->mpi_rank;

	for(int i=0; i < mesh->comm_map.RecvFrom.elem_count; i++){
		message_t* m = (message_t*) sc_array_index(&mesh->comm_map.RecvFrom, i);
		for(int j =0; j < m->idxs.elem_count; j++){
			int32_t *id = (int32_t*) sc_array_index(&m->idxs,j);
			mesh->part_nodes[*id] = m->rank;
		}
	}
}
