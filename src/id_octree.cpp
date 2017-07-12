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


#include "hexa.h"
#include "refinement.h"
#include "hilbert.h"


typedef struct {
	bitmask_t coord[2];
	unsigned int edge_id;
} edge_id_t;

typedef struct {
	unsigned int id;
} name_t;

unsigned edge_hash(const void *v, const void *u) {
	const edge_id_t *q = (const edge_id_t*) v;
	uint32_t a, b, c;

	a = (uint32_t) q->coord[0];
	b = (uint32_t) q->coord[1];
	c = (uint32_t) 0;
	sc_hash_mix(a, b, c);
	sc_hash_final(a, b, c);
	return (unsigned) c;
}

int edge_equal(const void *v, const void *u, const void *w) {
	const edge_id_t *e1 = (const edge_id_t*) v;
	const edge_id_t *e2 = (const edge_id_t*) u;

	return (unsigned) ((e1->coord[0] == e2->coord[0]) && (e1->coord[1] == e2->coord[1]));

}

unsigned id_hash(const void *v, const void *u) {
	const name_t *q = (const name_t*) v;
	uint32_t a, b, c;

	a = (uint32_t) q->id;
	b = (uint32_t) 0;
	c = (uint32_t) 0;
	sc_hash_mix(a, b, c);
	sc_hash_final(a, b, c);
	return (unsigned) c;
}

int id_equal(const void *v, const void *u, const void *w) {
	const name_t *e1 = (const name_t*) v;
	const name_t *e2 = (const name_t*) u;

	return (unsigned) ((e1->id == e2->id));

}

unsigned int edge_id(int* nodes, sc_hash_array_t* hash, int &npoints) {
	size_t position;
	edge_id_t *r;
	edge_id_t key;
	key.coord[0] = nodes[0];
	key.coord[1] = nodes[1];

	r = (edge_id_t*) sc_hash_array_insert_unique(hash, &key, &position);
	if (r != NULL) {
		r->coord[0] = key.coord[0];
		r->coord[1] = key.coord[1];
		r->coord[2] = key.coord[2];
		r->edge_id = npoints;
		npoints++;
		return r->edge_id;
	} else {
		r = (edge_id_t*) sc_array_index(&hash->a, position);
		return r->edge_id;
	}
}

void edge_add(int edge_id, sc_hash_array_t* hash_id ) {
	size_t position;
	name_t *r;
	name_t key;
        key.id = edge_id;

        r = (name_t*) sc_hash_array_insert_unique(hash_id, &key, &position);
        if(r != NULL)
        {
            //printf("dentro do add-edge");
            r->id = key.id;
        }
}

void IdentifyTemplate(hexa_tree_t* mesh, const std::vector<double>& coords, std::vector<int>& elements_ids){

	int el_0 = 0;
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

	for (int iel = 0; iel < elements_ids.size(); ++iel) {

		octant_t *elem = (octant_t*) sc_array_index(&mesh->elements, elements_ids[iel]);
		int ed_cont = 0;

		for (int edge = 0; edge < 12; ++edge) {
			if(elem->edge_ref[edge]){
				ed_cont++;
			}
		}

		//template 11
		elem->pad = 144;
		elem->tem = 11;

		// element verification
		// template 1
		if(ed_cont == 1){
			if ( (elem->edge_ref[0]) && (!elem->edge_ref[1]) && (!elem->edge_ref[2]) && (!elem->edge_ref[3]) &&
					(!elem->edge_ref[4]) && (!elem->edge_ref[5]) && (!elem->edge_ref[6]) && (!elem->edge_ref[7]) &&
					(!elem->edge_ref[8]) && (!elem->edge_ref[9]) && (!elem->edge_ref[10]) && (!elem->edge_ref[11])) {

				elem->pad = 10;
				elem->tem = 1;

			}
			else if ( (!elem->edge_ref[0]) && (elem->edge_ref[1]) && (!elem->edge_ref[2]) && (!elem->edge_ref[3]) &&
					(!elem->edge_ref[4]) && (!elem->edge_ref[5]) && (!elem->edge_ref[6]) && (!elem->edge_ref[7]) &&
					(!elem->edge_ref[8]) && (!elem->edge_ref[9]) && (!elem->edge_ref[10]) && (!elem->edge_ref[11])) {

				elem->pad = 11;
				elem->tem = 1;

			}
			else if ( (!elem->edge_ref[0]) && (!elem->edge_ref[1]) && (elem->edge_ref[2]) && (!elem->edge_ref[3]) &&
					(!elem->edge_ref[4]) && (!elem->edge_ref[5]) && (!elem->edge_ref[6]) && (!elem->edge_ref[7]) &&
					(!elem->edge_ref[8]) && (!elem->edge_ref[9]) && (!elem->edge_ref[10]) && (!elem->edge_ref[11])) {

				elem->pad = 12;
				elem->tem = 1;

			}
			else if ( (!elem->edge_ref[0]) && (!elem->edge_ref[1]) && (!elem->edge_ref[2]) && (elem->edge_ref[3]) &&
					(!elem->edge_ref[4]) && (!elem->edge_ref[5]) && (!elem->edge_ref[6]) && (!elem->edge_ref[7]) &&
					(!elem->edge_ref[8]) && (!elem->edge_ref[9]) && (!elem->edge_ref[10]) && (!elem->edge_ref[11])) {

				elem->pad = 13;
				elem->tem = 1;

			}
			else if ( (!elem->edge_ref[0]) && (!elem->edge_ref[1]) && (!elem->edge_ref[2]) && (!elem->edge_ref[3]) &&
					(elem->edge_ref[4]) && (!elem->edge_ref[5]) && (!elem->edge_ref[6]) && (!elem->edge_ref[7]) &&
					(!elem->edge_ref[8]) && (!elem->edge_ref[9]) && (!elem->edge_ref[10]) && (!elem->edge_ref[11])) {

				elem->pad = 14;
				elem->tem = 1;

			}
			else if ( (!elem->edge_ref[0]) && (!elem->edge_ref[1]) && (!elem->edge_ref[2]) && (!elem->edge_ref[3]) &&
					(!elem->edge_ref[4]) && (elem->edge_ref[5]) && (!elem->edge_ref[6]) && (!elem->edge_ref[7]) &&
					(!elem->edge_ref[8]) && (!elem->edge_ref[9]) && (!elem->edge_ref[10]) && (!elem->edge_ref[11])) {

				elem->pad = 15;
				elem->tem = 1;

			}
			else if ( (!elem->edge_ref[0]) && (!elem->edge_ref[1]) && (!elem->edge_ref[2]) && (!elem->edge_ref[3]) &&
					(!elem->edge_ref[4]) && (!elem->edge_ref[5]) && (elem->edge_ref[6]) && (!elem->edge_ref[7]) &&
					(!elem->edge_ref[8]) && (!elem->edge_ref[9]) && (!elem->edge_ref[10]) && (!elem->edge_ref[11])) {

				elem->pad = 16;
				elem->tem = 1;

			}
			else if ( (!elem->edge_ref[0]) && (!elem->edge_ref[1]) && (!elem->edge_ref[2]) && (!elem->edge_ref[3]) &&
					(!elem->edge_ref[4]) && (!elem->edge_ref[5]) && (!elem->edge_ref[6]) && (elem->edge_ref[7]) &&
					(!elem->edge_ref[8]) && (!elem->edge_ref[9]) && (!elem->edge_ref[10]) && (!elem->edge_ref[11])) {

				elem->pad = 17;
				elem->tem = 1;

			}
			else if ( (!elem->edge_ref[0]) && (!elem->edge_ref[1]) && (!elem->edge_ref[2]) && (!elem->edge_ref[3]) &&
					(!elem->edge_ref[4]) && (!elem->edge_ref[5]) && (!elem->edge_ref[6]) && (!elem->edge_ref[7]) &&
					(elem->edge_ref[8]) && (!elem->edge_ref[9]) && (!elem->edge_ref[10]) && (!elem->edge_ref[11])) {

				elem->pad = 18;
				elem->tem = 1;

			}
			else if ( (!elem->edge_ref[0]) && (!elem->edge_ref[1]) && (!elem->edge_ref[2]) && (!elem->edge_ref[3]) &&
					(!elem->edge_ref[4]) && (!elem->edge_ref[5]) && (!elem->edge_ref[6]) && (!elem->edge_ref[7]) &&
					(!elem->edge_ref[8]) && (elem->edge_ref[9]) && (!elem->edge_ref[10]) && (!elem->edge_ref[11])) {

				elem->pad = 19;
				elem->tem = 1;

			}
			else if ( (!elem->edge_ref[0]) && (!elem->edge_ref[1]) && (!elem->edge_ref[2]) && (!elem->edge_ref[3]) &&
					(!elem->edge_ref[4]) && (!elem->edge_ref[5]) && (!elem->edge_ref[6]) && (!elem->edge_ref[7]) &&
					(!elem->edge_ref[8]) && (!elem->edge_ref[9]) && (elem->edge_ref[10]) && (!elem->edge_ref[11])) {

				elem->pad = 20;
				elem->tem = 1;

			}
			else if ( (!elem->edge_ref[0]) && (!elem->edge_ref[1]) && (!elem->edge_ref[2]) && (!elem->edge_ref[3]) &&
					(!elem->edge_ref[4]) && (!elem->edge_ref[5]) && (!elem->edge_ref[6]) && (!elem->edge_ref[7]) &&
					(!elem->edge_ref[8]) && (!elem->edge_ref[9]) && (!elem->edge_ref[10]) && (elem->edge_ref[11])) {

				elem->pad = 21;
				elem->tem = 1;

			}
		}
		//template 2
		if(ed_cont == 4){
			if ( (elem->edge_ref[0]) && (!elem->edge_ref[1]) && (elem->edge_ref[2]) && (!elem->edge_ref[3]) &&
					(!elem->edge_ref[4]) && (!elem->edge_ref[5]) && (!elem->edge_ref[6]) && (!elem->edge_ref[7]) &&
					(elem->edge_ref[8]) && (!elem->edge_ref[9]) && (elem->edge_ref[10]) && (!elem->edge_ref[11])) {

				elem->pad = 22;
				elem->tem = 2;

			}
			else if ( (!elem->edge_ref[0]) && (elem->edge_ref[1]) && (!elem->edge_ref[2]) && (elem->edge_ref[3]) &&
					(!elem->edge_ref[4]) && (!elem->edge_ref[5]) && (!elem->edge_ref[6]) && (!elem->edge_ref[7]) &&
					(!elem->edge_ref[8]) && (elem->edge_ref[9]) && (!elem->edge_ref[10]) && (elem->edge_ref[11])) {

				elem->pad = 23;
				elem->tem = 2;

			}
			else if ( (!elem->edge_ref[0]) && (!elem->edge_ref[1]) && (!elem->edge_ref[2]) && (!elem->edge_ref[3]) &&
					(elem->edge_ref[4]) && (elem->edge_ref[5]) && (elem->edge_ref[6]) && (elem->edge_ref[7]) &&
					(!elem->edge_ref[8]) && (!elem->edge_ref[9]) && (!elem->edge_ref[10]) && (!elem->edge_ref[11])) {

				elem->pad = 24;
				elem->tem = 2;

			}
		}
		//template 3
		if(ed_cont == 2){
			if ( (elem->edge_ref[0]) && (!elem->edge_ref[1]) && (elem->edge_ref[2]) && (!elem->edge_ref[3]) &&
					(!elem->edge_ref[4]) && (!elem->edge_ref[5]) && (!elem->edge_ref[6]) && (!elem->edge_ref[7]) &&
					(!elem->edge_ref[8]) && (!elem->edge_ref[9]) && (!elem->edge_ref[10]) && (!elem->edge_ref[11])) {

				elem->pad = 25;
				elem->tem = 3;

			}
			else if ( (!elem->edge_ref[0]) && (elem->edge_ref[1]) && (!elem->edge_ref[2]) && (elem->edge_ref[3]) &&
					(!elem->edge_ref[4]) && (!elem->edge_ref[5]) && (!elem->edge_ref[6]) && (!elem->edge_ref[7]) &&
					(!elem->edge_ref[8]) && (!elem->edge_ref[9]) && (!elem->edge_ref[10]) && (!elem->edge_ref[11])) {

				elem->pad = 26;
				elem->tem = 3;

			}
			else if ( (!elem->edge_ref[0]) && (!elem->edge_ref[1]) && (!elem->edge_ref[2]) && (!elem->edge_ref[3]) &&
					(!elem->edge_ref[4]) && (!elem->edge_ref[5]) && (!elem->edge_ref[6]) && (!elem->edge_ref[7]) &&
					(elem->edge_ref[8]) && (!elem->edge_ref[9]) && (elem->edge_ref[10]) && (!elem->edge_ref[11])) {

				elem->pad = 27;
				elem->tem = 3;

			}
			else if ( (!elem->edge_ref[0]) && (!elem->edge_ref[1]) && (!elem->edge_ref[2]) && (!elem->edge_ref[3]) &&
					(!elem->edge_ref[4]) && (!elem->edge_ref[5]) && (!elem->edge_ref[6]) && (!elem->edge_ref[7]) &&
					(!elem->edge_ref[8]) && (elem->edge_ref[9]) && (!elem->edge_ref[10]) && (elem->edge_ref[11])) {

				elem->pad = 28;
				elem->tem = 3;

			}
			else if ( (!elem->edge_ref[0]) && (!elem->edge_ref[1]) && (!elem->edge_ref[2]) && (!elem->edge_ref[3]) &&
					(elem->edge_ref[4]) && (elem->edge_ref[5]) && (!elem->edge_ref[6]) && (!elem->edge_ref[7]) &&
					(!elem->edge_ref[8]) && (!elem->edge_ref[9]) && (!elem->edge_ref[10]) && (!elem->edge_ref[11])) {

				elem->pad = 29;
				elem->tem = 3;

			}
			else if ( (!elem->edge_ref[0]) && (!elem->edge_ref[1]) && (!elem->edge_ref[2]) && (!elem->edge_ref[3]) &&
					(!elem->edge_ref[4]) && (elem->edge_ref[5]) && (elem->edge_ref[6]) && (!elem->edge_ref[7]) &&
					(!elem->edge_ref[8]) && (!elem->edge_ref[9]) && (!elem->edge_ref[10]) && (!elem->edge_ref[11])) {

				elem->pad = 30;
				elem->tem = 3;

			}
			else if ( (!elem->edge_ref[0]) && (!elem->edge_ref[1]) && (!elem->edge_ref[2]) && (!elem->edge_ref[3]) &&
					(!elem->edge_ref[4]) && (!elem->edge_ref[5]) && (elem->edge_ref[6]) && (elem->edge_ref[7]) &&
					(!elem->edge_ref[8]) && (!elem->edge_ref[9]) && (!elem->edge_ref[10]) && (!elem->edge_ref[11])) {

				elem->pad = 31;
				elem->tem = 3;

			}
			else if ( (!elem->edge_ref[0]) && (!elem->edge_ref[1]) && (!elem->edge_ref[2]) && (!elem->edge_ref[3]) &&
					(elem->edge_ref[4]) && (!elem->edge_ref[5]) && (!elem->edge_ref[6]) && (elem->edge_ref[7]) &&
					(!elem->edge_ref[8]) && (!elem->edge_ref[9]) && (!elem->edge_ref[10]) && (!elem->edge_ref[11])) {

				elem->pad = 32;
				elem->tem = 3;

			}
			else if ( (elem->edge_ref[0]) && (!elem->edge_ref[1]) && (!elem->edge_ref[2]) && (!elem->edge_ref[3]) &&
					(!elem->edge_ref[4]) && (!elem->edge_ref[5]) && (!elem->edge_ref[6]) && (!elem->edge_ref[7]) &&
					(elem->edge_ref[8]) && (!elem->edge_ref[9]) && (!elem->edge_ref[10]) && (!elem->edge_ref[11])) {

				elem->pad = 33;
				elem->tem = 3;

			}
			else if ( (!elem->edge_ref[0]) && (elem->edge_ref[1]) && (!elem->edge_ref[2]) && (!elem->edge_ref[3]) &&
					(!elem->edge_ref[4]) && (!elem->edge_ref[5]) && (!elem->edge_ref[6]) && (!elem->edge_ref[7]) &&
					(!elem->edge_ref[8]) && (elem->edge_ref[9]) && (!elem->edge_ref[10]) && (!elem->edge_ref[11])) {

				elem->pad = 34;
				elem->tem = 3;

			}
			else if ( (!elem->edge_ref[0]) && (!elem->edge_ref[1]) && (elem->edge_ref[2]) && (!elem->edge_ref[3]) &&
					(!elem->edge_ref[4]) && (!elem->edge_ref[5]) && (!elem->edge_ref[6]) && (!elem->edge_ref[7]) &&
					(!elem->edge_ref[8]) && (!elem->edge_ref[9]) && (elem->edge_ref[10]) && (!elem->edge_ref[11])) {

				elem->pad = 35;
				elem->tem = 3;

			}
			else if ( (!elem->edge_ref[0]) && (!elem->edge_ref[1]) && (!elem->edge_ref[2]) && (elem->edge_ref[3]) &&
					(!elem->edge_ref[4]) && (!elem->edge_ref[5]) && (!elem->edge_ref[6]) && (!elem->edge_ref[7]) &&
					(!elem->edge_ref[8]) && (!elem->edge_ref[9]) && (!elem->edge_ref[10]) && (elem->edge_ref[11])) {

				elem->pad = 36;
				elem->tem = 3;

			}
		}
		//template 4
		if(ed_cont == 4){
			if ( (elem->edge_ref[0]) && (elem->edge_ref[1]) && (!elem->edge_ref[2]) && (!elem->edge_ref[3]) &&
					(!elem->edge_ref[4]) && (!elem->edge_ref[5]) && (!elem->edge_ref[6]) && (!elem->edge_ref[7]) &&
					(elem->edge_ref[8]) && (elem->edge_ref[9]) && (!elem->edge_ref[10]) && (!elem->edge_ref[11])) {

				elem->pad = 37;
				elem->tem = 4;

			}
			else if ( (!elem->edge_ref[0]) && (elem->edge_ref[1]) && (elem->edge_ref[2]) && (!elem->edge_ref[3]) &&
					(!elem->edge_ref[4]) && (!elem->edge_ref[5]) && (!elem->edge_ref[6]) && (!elem->edge_ref[7]) &&
					(!elem->edge_ref[8]) && (elem->edge_ref[9]) && (elem->edge_ref[10]) && (!elem->edge_ref[11])) {

				elem->pad = 38;
				elem->tem = 4;

			}
			else if ( (!elem->edge_ref[0]) && (!elem->edge_ref[1]) && (elem->edge_ref[2]) && (elem->edge_ref[3]) &&
					(!elem->edge_ref[4]) && (!elem->edge_ref[5]) && (!elem->edge_ref[6]) && (!elem->edge_ref[7]) &&
					(!elem->edge_ref[8]) && (!elem->edge_ref[9]) && (elem->edge_ref[10]) && (elem->edge_ref[11])) {

				elem->pad = 39;
				elem->tem = 4;

			}
			else if ( (elem->edge_ref[0]) && (!elem->edge_ref[1]) && (!elem->edge_ref[2]) && (elem->edge_ref[3]) &&
					(!elem->edge_ref[4]) && (!elem->edge_ref[5]) && (!elem->edge_ref[6]) && (!elem->edge_ref[7]) &&
					(elem->edge_ref[8]) && (!elem->edge_ref[9]) && (!elem->edge_ref[10]) && (elem->edge_ref[11])) {

				elem->pad = 40;
				elem->tem = 4;

			}
			else if ( (elem->edge_ref[0]) && (!elem->edge_ref[1]) && (elem->edge_ref[2]) && (!elem->edge_ref[3]) &&
					(elem->edge_ref[4]) && (!elem->edge_ref[5]) && (!elem->edge_ref[6]) && (elem->edge_ref[7]) &&
					(!elem->edge_ref[8]) && (!elem->edge_ref[9]) && (!elem->edge_ref[10]) && (!elem->edge_ref[11])) {

				elem->pad = 41;
				elem->tem = 4;

			}
			else if ( (elem->edge_ref[0]) && (!elem->edge_ref[1]) && (elem->edge_ref[2]) && (!elem->edge_ref[3]) &&
					(!elem->edge_ref[4]) && (elem->edge_ref[5]) && (elem->edge_ref[6]) && (!elem->edge_ref[7]) &&
					(!elem->edge_ref[8]) && (!elem->edge_ref[9]) && (!elem->edge_ref[10]) && (!elem->edge_ref[11])) {

				elem->pad = 42;
				elem->tem = 4;

			}
			else if ( (!elem->edge_ref[0]) && (!elem->edge_ref[1]) && (!elem->edge_ref[2]) && (!elem->edge_ref[3]) &&
					(!elem->edge_ref[4]) && (elem->edge_ref[5]) && (elem->edge_ref[6]) && (!elem->edge_ref[7]) &&
					(elem->edge_ref[8]) && (!elem->edge_ref[9]) && (elem->edge_ref[10]) && (!elem->edge_ref[11])) {

				elem->pad = 43;
				elem->tem = 4;

			}
			else if ( (!elem->edge_ref[0]) && (!elem->edge_ref[1]) && (!elem->edge_ref[2]) && (!elem->edge_ref[3]) &&
					(elem->edge_ref[4]) && (!elem->edge_ref[5]) && (!elem->edge_ref[6]) && (elem->edge_ref[7]) &&
					(elem->edge_ref[8]) && (!elem->edge_ref[9]) && (elem->edge_ref[10]) && (!elem->edge_ref[11])) {

				elem->pad = 44;
				elem->tem = 4;

			}
			else if ( (!elem->edge_ref[0]) && (elem->edge_ref[1]) && (!elem->edge_ref[2]) && (elem->edge_ref[3]) &&
					(elem->edge_ref[4]) && (elem->edge_ref[5]) && (!elem->edge_ref[6]) && (!elem->edge_ref[7]) &&
					(!elem->edge_ref[8]) && (!elem->edge_ref[9]) && (!elem->edge_ref[10]) && (!elem->edge_ref[11])) {

				elem->pad = 45;
				elem->tem = 4;

			}
			else if ( (!elem->edge_ref[0]) && (elem->edge_ref[1]) && (!elem->edge_ref[2]) && (elem->edge_ref[3]) &&
					(!elem->edge_ref[4]) && (!elem->edge_ref[5]) && (elem->edge_ref[6]) && (elem->edge_ref[7]) &&
					(!elem->edge_ref[8]) && (!elem->edge_ref[9]) && (!elem->edge_ref[10]) && (!elem->edge_ref[11])) {

				elem->pad = 46;
				elem->tem = 4;

			}
			else if ( (!elem->edge_ref[0]) && (!elem->edge_ref[1]) && (!elem->edge_ref[2]) && (!elem->edge_ref[3]) &&
					(!elem->edge_ref[4]) && (!elem->edge_ref[5]) && (elem->edge_ref[6]) && (elem->edge_ref[7]) &&
					(!elem->edge_ref[8]) && (elem->edge_ref[9]) && (!elem->edge_ref[10]) && (elem->edge_ref[11])) {

				elem->pad = 47;
				elem->tem = 4;

			}
			else if ( (!elem->edge_ref[0]) && (!elem->edge_ref[1]) && (!elem->edge_ref[2]) && (!elem->edge_ref[3]) &&
					(elem->edge_ref[4]) && (elem->edge_ref[5]) && (!elem->edge_ref[6]) && (!elem->edge_ref[7]) &&
					(!elem->edge_ref[8]) && (elem->edge_ref[9]) && (!elem->edge_ref[10]) && (elem->edge_ref[11])) {

				elem->pad = 48;
				elem->tem = 4;

			}
		}
		//template 5
		if(ed_cont == 4 || ed_cont == 3 || ed_cont == 2){
			//template 5 com 4
			if ( (elem->edge_ref[0]) && (elem->edge_ref[1]) && (elem->edge_ref[2]) && (elem->edge_ref[3]) &&
					(!elem->edge_ref[4]) && (!elem->edge_ref[5]) && (!elem->edge_ref[6]) && (!elem->edge_ref[7]) &&
					(!elem->edge_ref[8]) && (!elem->edge_ref[9]) && (!elem->edge_ref[10]) && (!elem->edge_ref[11])) {

				elem->pad = 49;
				elem->tem = 5;

			}
			else if ( (!elem->edge_ref[0]) && (!elem->edge_ref[1]) && (!elem->edge_ref[2]) && (!elem->edge_ref[3]) &&
					(!elem->edge_ref[4]) && (!elem->edge_ref[5]) && (!elem->edge_ref[6]) && (!elem->edge_ref[7]) &&
					(elem->edge_ref[8]) && (elem->edge_ref[9]) && (elem->edge_ref[10]) && (elem->edge_ref[11])) {

				elem->pad = 50;
				elem->tem = 5;

			}
			else if ( (elem->edge_ref[0]) && (!elem->edge_ref[1]) && (!elem->edge_ref[2]) && (!elem->edge_ref[3]) &&
					(elem->edge_ref[4]) && (elem->edge_ref[5]) && (!elem->edge_ref[6]) && (!elem->edge_ref[7]) &&
					(elem->edge_ref[8]) && (!elem->edge_ref[9]) && (!elem->edge_ref[10]) && (!elem->edge_ref[11])) {

				elem->pad = 51;
				elem->tem = 5;

			}
			else if ( (!elem->edge_ref[0]) && (elem->edge_ref[1]) && (!elem->edge_ref[2]) && (!elem->edge_ref[3]) &&
					(!elem->edge_ref[4]) && (elem->edge_ref[5]) && (elem->edge_ref[6]) && (!elem->edge_ref[7]) &&
					(!elem->edge_ref[8]) && (elem->edge_ref[9]) && (!elem->edge_ref[10]) && (!elem->edge_ref[11])) {

				elem->pad = 52;
				elem->tem = 5;

			}
			else if ( (!elem->edge_ref[0]) && (!elem->edge_ref[1]) && (elem->edge_ref[2]) && (!elem->edge_ref[3]) &&
					(!elem->edge_ref[4]) && (!elem->edge_ref[5]) && (elem->edge_ref[6]) && (elem->edge_ref[7]) &&
					(!elem->edge_ref[8]) && (!elem->edge_ref[9]) && (elem->edge_ref[10]) && (!elem->edge_ref[11])) {

				elem->pad = 53;
				elem->tem = 5;

			}
			else if ( (!elem->edge_ref[0]) && (!elem->edge_ref[1]) && (!elem->edge_ref[2]) && (elem->edge_ref[3]) &&
					(elem->edge_ref[4]) && (!elem->edge_ref[5]) && (!elem->edge_ref[6]) && (elem->edge_ref[7]) &&
					(!elem->edge_ref[8]) && (!elem->edge_ref[9]) && (!elem->edge_ref[10]) && (elem->edge_ref[11])) {

				elem->pad = 54;
				elem->tem = 5;

			}


			//template 5 com 3
			else if ( (elem->edge_ref[0]) && (elem->edge_ref[1]) && (!elem->edge_ref[2]) && (elem->edge_ref[3]) &&
					(!elem->edge_ref[4]) && (!elem->edge_ref[5]) && (!elem->edge_ref[6]) && (!elem->edge_ref[7]) &&
					(!elem->edge_ref[8]) && (!elem->edge_ref[9]) && (!elem->edge_ref[10]) && (!elem->edge_ref[11])) {

				elem->pad = 55;
				elem->tem = 5;

			}
			else if ( (elem->edge_ref[0]) && (elem->edge_ref[1]) && (elem->edge_ref[2]) && (!elem->edge_ref[3]) &&
					(!elem->edge_ref[4]) && (!elem->edge_ref[5]) && (!elem->edge_ref[6]) && (!elem->edge_ref[7]) &&
					(!elem->edge_ref[8]) && (!elem->edge_ref[9]) && (!elem->edge_ref[10]) && (!elem->edge_ref[11])) {

				elem->pad = 56;
				elem->tem = 5;

			}
			else if ( (!elem->edge_ref[0]) && (elem->edge_ref[1]) && (elem->edge_ref[2]) && (elem->edge_ref[3]) &&
					(!elem->edge_ref[4]) && (!elem->edge_ref[5]) && (!elem->edge_ref[6]) && (!elem->edge_ref[7]) &&
					(!elem->edge_ref[8]) && (!elem->edge_ref[9]) && (!elem->edge_ref[10]) && (!elem->edge_ref[11])) {

				elem->pad = 57;
				elem->tem = 5;

			}
			else if ( (elem->edge_ref[0]) && (!elem->edge_ref[1]) && (elem->edge_ref[2]) && (elem->edge_ref[3]) &&
					(!elem->edge_ref[4]) && (!elem->edge_ref[5]) && (!elem->edge_ref[6]) && (!elem->edge_ref[7]) &&
					(!elem->edge_ref[8]) && (!elem->edge_ref[9]) && (!elem->edge_ref[10]) && (!elem->edge_ref[11])) {

				elem->pad = 58;
				elem->tem = 5;

			}
			else if ( (!elem->edge_ref[0]) && (!elem->edge_ref[1]) && (!elem->edge_ref[2]) && (!elem->edge_ref[3]) &&
					(!elem->edge_ref[4]) && (!elem->edge_ref[5]) && (!elem->edge_ref[6]) && (!elem->edge_ref[7]) &&
					(elem->edge_ref[8]) && (elem->edge_ref[9]) && (elem->edge_ref[10]) && (!elem->edge_ref[11])) {

				elem->pad = 59;
				elem->tem = 5;

			}
			else if ( (!elem->edge_ref[0]) && (!elem->edge_ref[1]) && (!elem->edge_ref[2]) && (!elem->edge_ref[3]) &&
					(!elem->edge_ref[4]) && (!elem->edge_ref[5]) && (!elem->edge_ref[6]) && (!elem->edge_ref[7]) &&
					(elem->edge_ref[8]) && (elem->edge_ref[9]) && (!elem->edge_ref[10]) && (elem->edge_ref[11])) {

				elem->pad = 60;
				elem->tem = 5;

			}
			else if ( (!elem->edge_ref[0]) && (!elem->edge_ref[1]) && (!elem->edge_ref[2]) && (!elem->edge_ref[3]) &&
					(!elem->edge_ref[4]) && (!elem->edge_ref[5]) && (!elem->edge_ref[6]) && (!elem->edge_ref[7]) &&
					(!elem->edge_ref[8]) && (elem->edge_ref[9]) && (elem->edge_ref[10]) && (elem->edge_ref[11])) {

				elem->pad = 61;
				elem->tem = 5;

			}
			else if ( (!elem->edge_ref[0]) && (!elem->edge_ref[1]) && (!elem->edge_ref[2]) && (!elem->edge_ref[3]) &&
					(!elem->edge_ref[4]) && (!elem->edge_ref[5]) && (!elem->edge_ref[6]) && (!elem->edge_ref[7]) &&
					(elem->edge_ref[8]) && (!elem->edge_ref[9]) && (elem->edge_ref[10]) && (elem->edge_ref[11])) {

				elem->pad = 62;
				elem->tem = 5;

			}
			else if ( (elem->edge_ref[0]) && (!elem->edge_ref[1]) && (!elem->edge_ref[2]) && (!elem->edge_ref[3]) &&
					(elem->edge_ref[4]) && (elem->edge_ref[5]) && (!elem->edge_ref[6]) && (!elem->edge_ref[7]) &&
					(!elem->edge_ref[8]) && (!elem->edge_ref[9]) && (!elem->edge_ref[10]) && (!elem->edge_ref[11])) {

				elem->pad = 63;
				elem->tem = 5;

			}
			else if ( (elem->edge_ref[0]) && (!elem->edge_ref[1]) && (!elem->edge_ref[2]) && (!elem->edge_ref[3]) &&
					(!elem->edge_ref[4]) && (elem->edge_ref[5]) && (!elem->edge_ref[6]) && (!elem->edge_ref[7]) &&
					(elem->edge_ref[8]) && (!elem->edge_ref[9]) && (!elem->edge_ref[10]) && (!elem->edge_ref[11])) {

				elem->pad = 64;
				elem->tem = 5;

			}
			else if ( (!elem->edge_ref[0]) && (!elem->edge_ref[1]) && (!elem->edge_ref[2]) && (!elem->edge_ref[3]) &&
					(elem->edge_ref[4]) && (elem->edge_ref[5]) && (!elem->edge_ref[6]) && (!elem->edge_ref[7]) &&
					(elem->edge_ref[8]) && (!elem->edge_ref[9]) && (!elem->edge_ref[10]) && (!elem->edge_ref[11])) {

				elem->pad = 65;
				elem->tem = 5;

			}
			else if ( (elem->edge_ref[0]) && (!elem->edge_ref[1]) && (!elem->edge_ref[2]) && (!elem->edge_ref[3]) &&
					(elem->edge_ref[4]) && (!elem->edge_ref[5]) && (!elem->edge_ref[6]) && (!elem->edge_ref[7]) &&
					(elem->edge_ref[8]) && (!elem->edge_ref[9]) && (!elem->edge_ref[10]) && (!elem->edge_ref[11])) {

				elem->pad = 66;
				elem->tem = 5;

			}
			else if ( (!elem->edge_ref[0]) && (elem->edge_ref[1]) && (!elem->edge_ref[2]) && (!elem->edge_ref[3]) &&
					(!elem->edge_ref[4]) && (elem->edge_ref[5]) && (elem->edge_ref[6]) && (!elem->edge_ref[7]) &&
					(!elem->edge_ref[8]) && (!elem->edge_ref[9]) && (!elem->edge_ref[10]) && (!elem->edge_ref[11])) {

				elem->pad = 67;
				elem->tem = 5;

			}
			else if ( (!elem->edge_ref[0]) && (elem->edge_ref[1]) && (!elem->edge_ref[2]) && (!elem->edge_ref[3]) &&
					(!elem->edge_ref[4]) && (elem->edge_ref[5]) && (!elem->edge_ref[6]) && (!elem->edge_ref[7]) &&
					(!elem->edge_ref[8]) && (elem->edge_ref[9]) && (!elem->edge_ref[10]) && (!elem->edge_ref[11])) {

				elem->pad = 68;
				elem->tem = 5;

			}
			else if ( (!elem->edge_ref[0]) && (!elem->edge_ref[1]) && (!elem->edge_ref[2]) && (!elem->edge_ref[3]) &&
					(!elem->edge_ref[4]) && (elem->edge_ref[5]) && (elem->edge_ref[6]) && (!elem->edge_ref[7]) &&
					(!elem->edge_ref[8]) && (elem->edge_ref[9]) && (!elem->edge_ref[10]) && (!elem->edge_ref[11])) {

				elem->pad = 69;
				elem->tem = 5;

			}
			else if ( (!elem->edge_ref[0]) && (elem->edge_ref[1]) && (!elem->edge_ref[2]) && (!elem->edge_ref[3]) &&
					(!elem->edge_ref[4]) && (!elem->edge_ref[5]) && (elem->edge_ref[6]) && (!elem->edge_ref[7]) &&
					(!elem->edge_ref[8]) && (elem->edge_ref[9]) && (!elem->edge_ref[10]) && (!elem->edge_ref[11])) {

				elem->pad = 70;
				elem->tem = 5;

			}
			else if ( (!elem->edge_ref[0]) && (!elem->edge_ref[1]) && (elem->edge_ref[2]) && (!elem->edge_ref[3]) &&
					(!elem->edge_ref[4]) && (!elem->edge_ref[5]) && (elem->edge_ref[6]) && (elem->edge_ref[7]) &&
					(!elem->edge_ref[8]) && (!elem->edge_ref[9]) && (!elem->edge_ref[10]) && (!elem->edge_ref[11])) {

				elem->pad = 71;
				elem->tem = 5;

			}
			else if ( (!elem->edge_ref[0]) && (!elem->edge_ref[1]) && (elem->edge_ref[2]) && (!elem->edge_ref[3]) &&
					(!elem->edge_ref[4]) && (!elem->edge_ref[5]) && (elem->edge_ref[6]) && (!elem->edge_ref[7]) &&
					(!elem->edge_ref[8]) && (!elem->edge_ref[9]) && (elem->edge_ref[10]) && (!elem->edge_ref[11])) {

				elem->pad = 72;
				elem->tem = 5;

			}
			else if ( (!elem->edge_ref[0]) && (!elem->edge_ref[1]) && (!elem->edge_ref[2]) && (!elem->edge_ref[3]) &&
					(!elem->edge_ref[4]) && (!elem->edge_ref[5]) && (elem->edge_ref[6]) && (elem->edge_ref[7]) &&
					(!elem->edge_ref[8]) && (!elem->edge_ref[9]) && (elem->edge_ref[10]) && (!elem->edge_ref[11])) {

				elem->pad = 73;
				elem->tem = 5;

			}
			else if ( (!elem->edge_ref[0]) && (!elem->edge_ref[1]) && (elem->edge_ref[2]) && (!elem->edge_ref[3]) &&
					(!elem->edge_ref[4]) && (!elem->edge_ref[5]) && (!elem->edge_ref[6]) && (elem->edge_ref[7]) &&
					(!elem->edge_ref[8]) && (!elem->edge_ref[9]) && (elem->edge_ref[10]) && (!elem->edge_ref[11])) {

				elem->pad = 74;
				elem->tem = 5;

			}
			else if ( (!elem->edge_ref[0]) && (!elem->edge_ref[1]) && (!elem->edge_ref[2]) && (elem->edge_ref[3]) &&
					(elem->edge_ref[4]) && (!elem->edge_ref[5]) && (!elem->edge_ref[6]) && (elem->edge_ref[7]) &&
					(!elem->edge_ref[8]) && (!elem->edge_ref[9]) && (!elem->edge_ref[10]) && (!elem->edge_ref[11])) {

				elem->pad = 75;
				elem->tem = 5;

			}
			else if ( (!elem->edge_ref[0]) && (!elem->edge_ref[1]) && (!elem->edge_ref[2]) && (elem->edge_ref[3]) &&
					(elem->edge_ref[4]) && (!elem->edge_ref[5]) && (!elem->edge_ref[6]) && (!elem->edge_ref[7]) &&
					(!elem->edge_ref[8]) && (!elem->edge_ref[9]) && (!elem->edge_ref[10]) && (elem->edge_ref[11])) {

				elem->pad = 76;
				elem->tem = 5;

			}
			else if ( (!elem->edge_ref[0]) && (!elem->edge_ref[1]) && (!elem->edge_ref[2]) && (!elem->edge_ref[3]) &&
					(elem->edge_ref[4]) && (!elem->edge_ref[5]) && (!elem->edge_ref[6]) && (elem->edge_ref[7]) &&
					(!elem->edge_ref[8]) && (!elem->edge_ref[9]) && (!elem->edge_ref[10]) && (elem->edge_ref[11])) {

				elem->pad = 77;
				elem->tem = 5;

			}
			else if ( (!elem->edge_ref[0]) && (!elem->edge_ref[1]) && (!elem->edge_ref[2]) && (elem->edge_ref[3]) &&
					(!elem->edge_ref[4]) && (!elem->edge_ref[5]) && (!elem->edge_ref[6]) && (elem->edge_ref[7]) &&
					(!elem->edge_ref[8]) && (!elem->edge_ref[9]) && (!elem->edge_ref[10]) && (elem->edge_ref[11])) {

				elem->pad = 78;
				elem->tem = 5;

			}

			// templete 5 com 2
			else if ( (elem->edge_ref[0]) && (elem->edge_ref[1]) && (!elem->edge_ref[2]) && (!elem->edge_ref[3]) &&
					(!elem->edge_ref[4]) && (!elem->edge_ref[5]) && (!elem->edge_ref[6]) && (!elem->edge_ref[7]) &&
					(!elem->edge_ref[8]) && (!elem->edge_ref[9]) && (!elem->edge_ref[10]) && (!elem->edge_ref[11])) {

				elem->pad = 79;
				elem->tem = 5;

			}
			else if ( (!elem->edge_ref[0]) && (elem->edge_ref[1]) && (elem->edge_ref[2]) && (!elem->edge_ref[3]) &&
					(!elem->edge_ref[4]) && (!elem->edge_ref[5]) && (!elem->edge_ref[6]) && (!elem->edge_ref[7]) &&
					(!elem->edge_ref[8]) && (!elem->edge_ref[9]) && (!elem->edge_ref[10]) && (!elem->edge_ref[11])) {

				elem->pad = 80;
				elem->tem = 5;

			}
			else if ( (!elem->edge_ref[0]) && (!elem->edge_ref[1]) && (elem->edge_ref[2]) && (elem->edge_ref[3]) &&
					(!elem->edge_ref[4]) && (!elem->edge_ref[5]) && (!elem->edge_ref[6]) && (!elem->edge_ref[7]) &&
					(!elem->edge_ref[8]) && (!elem->edge_ref[9]) && (!elem->edge_ref[10]) && (!elem->edge_ref[11])) {

				elem->pad = 81;
				elem->tem = 5;

			}
			else if ( (elem->edge_ref[0]) && (!elem->edge_ref[1]) && (!elem->edge_ref[2]) && (elem->edge_ref[3]) &&
					(!elem->edge_ref[4]) && (!elem->edge_ref[5]) && (!elem->edge_ref[6]) && (!elem->edge_ref[7]) &&
					(!elem->edge_ref[8]) && (!elem->edge_ref[9]) && (!elem->edge_ref[10]) && (!elem->edge_ref[11])) {

				elem->pad = 82;
				elem->tem = 5;

			}
			else if ( (!elem->edge_ref[0]) && (!elem->edge_ref[1]) && (!elem->edge_ref[2]) && (!elem->edge_ref[3]) &&
					(!elem->edge_ref[4]) && (!elem->edge_ref[5]) && (!elem->edge_ref[6]) && (!elem->edge_ref[7]) &&
					(elem->edge_ref[8]) && (elem->edge_ref[9]) && (!elem->edge_ref[10]) && (!elem->edge_ref[11])) {

				elem->pad = 83;
				elem->tem = 5;

			}
			else if ( (!elem->edge_ref[0]) && (!elem->edge_ref[1]) && (!elem->edge_ref[2]) && (!elem->edge_ref[3]) &&
					(!elem->edge_ref[4]) && (!elem->edge_ref[5]) && (!elem->edge_ref[6]) && (!elem->edge_ref[7]) &&
					(!elem->edge_ref[8]) && (elem->edge_ref[9]) && (elem->edge_ref[10]) && (!elem->edge_ref[11])) {

				elem->pad = 84;
				elem->tem = 5;

			}
			else if ( (!elem->edge_ref[0]) && (!elem->edge_ref[1]) && (!elem->edge_ref[2]) && (!elem->edge_ref[3]) &&
					(!elem->edge_ref[4]) && (!elem->edge_ref[5]) && (!elem->edge_ref[6]) && (!elem->edge_ref[7]) &&
					(!elem->edge_ref[8]) && (!elem->edge_ref[9]) && (elem->edge_ref[10]) && (elem->edge_ref[11])) {

				elem->pad = 85;
				elem->tem = 5;

			}
			else if ( (!elem->edge_ref[0]) && (!elem->edge_ref[1]) && (!elem->edge_ref[2]) && (!elem->edge_ref[3]) &&
					(!elem->edge_ref[4]) && (!elem->edge_ref[5]) && (!elem->edge_ref[6]) && (!elem->edge_ref[7]) &&
					(elem->edge_ref[8]) && (!elem->edge_ref[9]) && (!elem->edge_ref[10]) && (elem->edge_ref[11])) {

				elem->pad = 86;
				elem->tem = 5;

			}
			else if ( (elem->edge_ref[0]) && (!elem->edge_ref[1]) && (!elem->edge_ref[2]) && (!elem->edge_ref[3]) &&
					(elem->edge_ref[4]) && (!elem->edge_ref[5]) && (!elem->edge_ref[6]) && (!elem->edge_ref[7]) &&
					(!elem->edge_ref[8]) && (!elem->edge_ref[9]) && (!elem->edge_ref[10]) && (!elem->edge_ref[11])) {

				elem->pad = 87;
				elem->tem = 5;

			}
			else if ( (elem->edge_ref[0]) && (!elem->edge_ref[1]) && (!elem->edge_ref[2]) && (!elem->edge_ref[3]) &&
					(!elem->edge_ref[4]) && (elem->edge_ref[5]) && (!elem->edge_ref[6]) && (!elem->edge_ref[7]) &&
					(!elem->edge_ref[8]) && (!elem->edge_ref[9]) && (!elem->edge_ref[10]) && (!elem->edge_ref[11])) {

				elem->pad = 88;
				elem->tem = 5;

			}
			else if ( (!elem->edge_ref[0]) && (!elem->edge_ref[1]) && (!elem->edge_ref[2]) && (!elem->edge_ref[3]) &&
					(!elem->edge_ref[4]) && (elem->edge_ref[5]) && (!elem->edge_ref[6]) && (!elem->edge_ref[7]) &&
					(elem->edge_ref[8]) && (!elem->edge_ref[9]) && (!elem->edge_ref[10]) && (!elem->edge_ref[11])) {

				elem->pad = 89;
				elem->tem = 5;

			}
			else if ( (!elem->edge_ref[0]) && (!elem->edge_ref[1]) && (!elem->edge_ref[2]) && (!elem->edge_ref[3]) &&
					(elem->edge_ref[4]) && (!elem->edge_ref[5]) && (!elem->edge_ref[6]) && (!elem->edge_ref[7]) &&
					(elem->edge_ref[8]) && (!elem->edge_ref[9]) && (!elem->edge_ref[10]) && (!elem->edge_ref[11])) {

				elem->pad = 90;
				elem->tem = 5;

			}
			else if ( (!elem->edge_ref[0]) && (!elem->edge_ref[1]) && (!elem->edge_ref[2]) && (!elem->edge_ref[3]) &&
					(!elem->edge_ref[4]) && (elem->edge_ref[5]) && (!elem->edge_ref[6]) && (!elem->edge_ref[7]) &&
					(!elem->edge_ref[8]) && (elem->edge_ref[9]) && (!elem->edge_ref[10]) && (!elem->edge_ref[11])) {

				elem->pad = 91;
				elem->tem = 5;

			}
			else if ( (!elem->edge_ref[0]) && (!elem->edge_ref[1]) && (!elem->edge_ref[2]) && (!elem->edge_ref[3]) &&
					(!elem->edge_ref[4]) && (!elem->edge_ref[5]) && (elem->edge_ref[6]) && (!elem->edge_ref[7]) &&
					(!elem->edge_ref[8]) && (elem->edge_ref[9]) && (!elem->edge_ref[10]) && (!elem->edge_ref[11])) {

				elem->pad = 92;
				elem->tem = 5;

			}
			else if ( (!elem->edge_ref[0]) && (elem->edge_ref[1]) && (!elem->edge_ref[2]) && (!elem->edge_ref[3]) &&
					(!elem->edge_ref[4]) && (!elem->edge_ref[5]) && (elem->edge_ref[6]) && (!elem->edge_ref[7]) &&
					(!elem->edge_ref[8]) && (!elem->edge_ref[9]) && (!elem->edge_ref[10]) && (!elem->edge_ref[11])) {

				elem->pad = 93;
				elem->tem = 5;

			}
			else if ( (!elem->edge_ref[0]) && (elem->edge_ref[1]) && (!elem->edge_ref[2]) && (!elem->edge_ref[3]) &&
					(!elem->edge_ref[4]) && (elem->edge_ref[5]) && (!elem->edge_ref[6]) && (!elem->edge_ref[7]) &&
					(!elem->edge_ref[8]) && (!elem->edge_ref[9]) && (!elem->edge_ref[10]) && (!elem->edge_ref[11])) {

				elem->pad = 94;
				elem->tem = 5;

			}
			else if ( (!elem->edge_ref[0]) && (!elem->edge_ref[1]) && (elem->edge_ref[2]) && (!elem->edge_ref[3]) &&
					(!elem->edge_ref[4]) && (!elem->edge_ref[5]) && (elem->edge_ref[6]) && (!elem->edge_ref[7]) &&
					(!elem->edge_ref[8]) && (!elem->edge_ref[9]) && (!elem->edge_ref[10]) && (!elem->edge_ref[11])) {

				elem->pad = 95;
				elem->tem = 5;

			}
			else if ( (!elem->edge_ref[0]) && (!elem->edge_ref[1]) && (!elem->edge_ref[2]) && (!elem->edge_ref[3]) &&
					(!elem->edge_ref[4]) && (!elem->edge_ref[5]) && (elem->edge_ref[6]) && (!elem->edge_ref[7]) &&
					(!elem->edge_ref[8]) && (!elem->edge_ref[9]) && (elem->edge_ref[10]) && (!elem->edge_ref[11])) {

				elem->pad = 96;
				elem->tem = 5;

			}
			else if ( (!elem->edge_ref[0]) && (!elem->edge_ref[1]) && (!elem->edge_ref[2]) && (!elem->edge_ref[3]) &&
					(!elem->edge_ref[4]) && (!elem->edge_ref[5]) && (!elem->edge_ref[6]) && (elem->edge_ref[7]) &&
					(!elem->edge_ref[8]) && (!elem->edge_ref[9]) && (elem->edge_ref[10]) && (!elem->edge_ref[11])) {

				elem->pad = 97;
				elem->tem = 5;

			}
			else if ( (!elem->edge_ref[0]) && (!elem->edge_ref[1]) && (elem->edge_ref[2]) && (!elem->edge_ref[3]) &&
					(!elem->edge_ref[4]) && (!elem->edge_ref[5]) && (!elem->edge_ref[6]) && (elem->edge_ref[7]) &&
					(!elem->edge_ref[8]) && (!elem->edge_ref[9]) && (!elem->edge_ref[10]) && (!elem->edge_ref[11])) {

				elem->pad = 98;
				elem->tem = 5;

			}
			else if ( (!elem->edge_ref[0]) && (!elem->edge_ref[1]) && (!elem->edge_ref[2]) && (elem->edge_ref[3]) &&
					(elem->edge_ref[4]) && (!elem->edge_ref[5]) && (!elem->edge_ref[6]) && (!elem->edge_ref[7]) &&
					(!elem->edge_ref[8]) && (!elem->edge_ref[9]) && (!elem->edge_ref[10]) && (!elem->edge_ref[11])) {

				elem->pad = 99;
				elem->tem = 5;

			}
			else if ( (!elem->edge_ref[0]) && (!elem->edge_ref[1]) && (!elem->edge_ref[2]) && (!elem->edge_ref[3]) &&
					(elem->edge_ref[4]) && (!elem->edge_ref[5]) && (elem->edge_ref[6]) && (!elem->edge_ref[7]) &&
					(!elem->edge_ref[8]) && (!elem->edge_ref[9]) && (!elem->edge_ref[10]) && (!elem->edge_ref[11])) {

				elem->pad = 100;
				elem->tem = 5;

			}
			else if ( (!elem->edge_ref[0]) && (!elem->edge_ref[1]) && (!elem->edge_ref[2]) && (!elem->edge_ref[3]) &&
					(!elem->edge_ref[4]) && (!elem->edge_ref[5]) && (elem->edge_ref[6]) && (elem->edge_ref[7]) &&
					(!elem->edge_ref[8]) && (!elem->edge_ref[9]) && (!elem->edge_ref[10]) && (!elem->edge_ref[11])) {

				elem->pad = 101;
				elem->tem = 5;

			}
			else if ( (!elem->edge_ref[0]) && (!elem->edge_ref[1]) && (!elem->edge_ref[2]) && (elem->edge_ref[3]) &&
					(!elem->edge_ref[4]) && (!elem->edge_ref[5]) && (!elem->edge_ref[6]) && (elem->edge_ref[7]) &&
					(!elem->edge_ref[8]) && (!elem->edge_ref[9]) && (!elem->edge_ref[10]) && (!elem->edge_ref[11])) {

				elem->pad = 102;
				elem->tem = 5;

			}
		}
		//template 6
		if(ed_cont == 6){
			if ( (elem->edge_ref[0]) && (elem->edge_ref[1]) && (!elem->edge_ref[2]) && (elem->edge_ref[3]) &&
					(!elem->edge_ref[4]) && (!elem->edge_ref[5]) && (!elem->edge_ref[6]) && (!elem->edge_ref[7]) &&
					(elem->edge_ref[8]) && (elem->edge_ref[9]) && (!elem->edge_ref[10]) && (elem->edge_ref[11])) {

				elem->pad = 103;
				elem->tem = 6;

			}
			else if ( (elem->edge_ref[0]) && (elem->edge_ref[1]) && (elem->edge_ref[2]) && (!elem->edge_ref[3]) &&
					(!elem->edge_ref[4]) && (!elem->edge_ref[5]) && (!elem->edge_ref[6]) && (!elem->edge_ref[7]) &&
					(elem->edge_ref[8]) && (elem->edge_ref[9]) && (elem->edge_ref[10]) && (!elem->edge_ref[11])) {

				elem->pad = 104;
				elem->tem = 6;

			}
			else if ( (!elem->edge_ref[0]) && (elem->edge_ref[1]) && (elem->edge_ref[2]) && (elem->edge_ref[3]) &&
					(!elem->edge_ref[4]) && (!elem->edge_ref[5]) && (!elem->edge_ref[6]) && (!elem->edge_ref[7]) &&
					(!elem->edge_ref[8]) && (elem->edge_ref[9]) && (elem->edge_ref[10]) && (elem->edge_ref[11])) {

				elem->pad = 105;
				elem->tem = 6;

			}
			else if ( (elem->edge_ref[0]) && (!elem->edge_ref[1]) && (elem->edge_ref[2]) && (elem->edge_ref[3]) &&
					(!elem->edge_ref[4]) && (!elem->edge_ref[5]) && (!elem->edge_ref[6]) && (!elem->edge_ref[7]) &&
					(elem->edge_ref[8]) && (!elem->edge_ref[9]) && (elem->edge_ref[10]) && (elem->edge_ref[11])) {

				elem->pad = 106;
				elem->tem = 6;

			}
			else if ( (elem->edge_ref[0]) && (!elem->edge_ref[1]) && (elem->edge_ref[2]) && (!elem->edge_ref[3]) &&
					(!elem->edge_ref[4]) && (elem->edge_ref[5]) && (elem->edge_ref[6]) && (!elem->edge_ref[7]) &&
					(elem->edge_ref[8]) && (!elem->edge_ref[9]) && (elem->edge_ref[10]) && (!elem->edge_ref[11])) {

				elem->pad = 107;
				elem->tem = 6;

			}
			else if ( (!elem->edge_ref[0]) && (!elem->edge_ref[1]) && (!elem->edge_ref[2]) && (!elem->edge_ref[3]) &&
					(elem->edge_ref[4]) && (elem->edge_ref[5]) && (elem->edge_ref[6]) && (elem->edge_ref[7]) &&
					(elem->edge_ref[8]) && (!elem->edge_ref[9]) && (elem->edge_ref[10]) && (!elem->edge_ref[11])) {

				elem->pad = 108;
				elem->tem = 6;

			}
			else if ( (elem->edge_ref[0]) && (!elem->edge_ref[1]) && (elem->edge_ref[2]) && (!elem->edge_ref[3]) &&
					(elem->edge_ref[4]) && (!elem->edge_ref[5]) && (!elem->edge_ref[6]) && (elem->edge_ref[7]) &&
					(elem->edge_ref[8]) && (!elem->edge_ref[9]) && (elem->edge_ref[10]) && (!elem->edge_ref[11])) {

				elem->pad = 109;
				elem->tem = 6;

			}
			else if ( (elem->edge_ref[0]) && (!elem->edge_ref[1]) && (elem->edge_ref[2]) && (!elem->edge_ref[3]) &&
					(elem->edge_ref[4]) && (elem->edge_ref[5]) && (elem->edge_ref[6]) && (elem->edge_ref[7]) &&
					(!elem->edge_ref[8]) && (!elem->edge_ref[9]) && (!elem->edge_ref[10]) && (!elem->edge_ref[11])) {

				elem->pad = 110;
				elem->tem = 6;

			}
			else if ( (!elem->edge_ref[0]) && (elem->edge_ref[1]) && (!elem->edge_ref[2]) && (elem->edge_ref[3]) &&
					(elem->edge_ref[4]) && (elem->edge_ref[5]) && (elem->edge_ref[6]) && (elem->edge_ref[7]) &&
					(!elem->edge_ref[8]) && (!elem->edge_ref[9]) && (!elem->edge_ref[10]) && (!elem->edge_ref[11])) {

				elem->pad = 111;
				elem->tem = 6;

			}
			else if ( (!elem->edge_ref[0]) && (elem->edge_ref[1]) && (!elem->edge_ref[2]) && (elem->edge_ref[3]) &&
					(elem->edge_ref[4]) && (elem->edge_ref[5]) && (!elem->edge_ref[6]) && (!elem->edge_ref[7]) &&
					(!elem->edge_ref[8]) && (elem->edge_ref[9]) && (!elem->edge_ref[10]) && (elem->edge_ref[11])) {

				elem->pad = 112;
				elem->tem = 6;

			}
			else if ( (!elem->edge_ref[0]) && (!elem->edge_ref[1]) && (!elem->edge_ref[2]) && (!elem->edge_ref[3]) &&
					(elem->edge_ref[4]) && (elem->edge_ref[5]) && (elem->edge_ref[6]) && (elem->edge_ref[7]) &&
					(!elem->edge_ref[8]) && (elem->edge_ref[9]) && (!elem->edge_ref[10]) && (elem->edge_ref[11])) {

				elem->pad = 113;
				elem->tem = 6;

			}
			else if ( (!elem->edge_ref[0]) && (elem->edge_ref[1]) && (!elem->edge_ref[2]) && (elem->edge_ref[3]) &&
					(!elem->edge_ref[4]) && (!elem->edge_ref[5]) && (elem->edge_ref[6]) && (elem->edge_ref[7]) &&
					(!elem->edge_ref[8]) && (elem->edge_ref[9]) && (!elem->edge_ref[10]) && (elem->edge_ref[11])) {

				elem->pad = 114;
				elem->tem = 6;

			}
		}
		//template 7
		if(ed_cont == 3){
			if ( (elem->edge_ref[0]) && (!elem->edge_ref[1]) && (!elem->edge_ref[2]) && (elem->edge_ref[3]) &&
					(elem->edge_ref[4]) && (!elem->edge_ref[5]) && (!elem->edge_ref[6]) && (!elem->edge_ref[7]) &&
					(!elem->edge_ref[8]) && (!elem->edge_ref[9]) && (!elem->edge_ref[10]) && (!elem->edge_ref[11])) {

				elem->pad = 115;
				elem->tem = 7;

			}
			else if ( (elem->edge_ref[0]) && (elem->edge_ref[1]) && (!elem->edge_ref[2]) && (!elem->edge_ref[3]) &&
					(!elem->edge_ref[4]) && (elem->edge_ref[5]) && (!elem->edge_ref[6]) && (!elem->edge_ref[7]) &&
					(!elem->edge_ref[8]) && (!elem->edge_ref[9]) && (!elem->edge_ref[10]) && (!elem->edge_ref[11])) {

				elem->pad = 116;
				elem->tem = 7;

			}
			else if ( (!elem->edge_ref[0]) && (elem->edge_ref[1]) && (elem->edge_ref[2]) && (!elem->edge_ref[3]) &&
					(!elem->edge_ref[4]) && (!elem->edge_ref[5]) && (elem->edge_ref[6]) && (!elem->edge_ref[7]) &&
					(!elem->edge_ref[8]) && (!elem->edge_ref[9]) && (!elem->edge_ref[10]) && (!elem->edge_ref[11])) {

				elem->pad = 117;
				elem->tem = 7;

			}
			else if ( (!elem->edge_ref[0]) && (!elem->edge_ref[1]) && (elem->edge_ref[2]) && (elem->edge_ref[3]) &&
					(!elem->edge_ref[4]) && (!elem->edge_ref[5]) && (!elem->edge_ref[6]) && (elem->edge_ref[7]) &&
					(!elem->edge_ref[8]) && (!elem->edge_ref[9]) && (!elem->edge_ref[10]) && (!elem->edge_ref[11])) {

				elem->pad = 118;
				elem->tem = 7;

			}
			else if ( (!elem->edge_ref[0]) && (!elem->edge_ref[1]) && (!elem->edge_ref[2]) && (!elem->edge_ref[3]) &&
					(elem->edge_ref[4]) && (!elem->edge_ref[5]) && (!elem->edge_ref[6]) && (!elem->edge_ref[7]) &&
					(elem->edge_ref[8]) && (!elem->edge_ref[9]) && (!elem->edge_ref[10]) && (elem->edge_ref[11])) {

				elem->pad = 119;
				elem->tem = 7;

			}
			else if ( (!elem->edge_ref[0]) && (!elem->edge_ref[1]) && (!elem->edge_ref[2]) && (!elem->edge_ref[3]) &&
					(!elem->edge_ref[4]) && (elem->edge_ref[5]) && (!elem->edge_ref[6]) && (!elem->edge_ref[7]) &&
					(elem->edge_ref[8]) && (elem->edge_ref[9]) && (!elem->edge_ref[10]) && (!elem->edge_ref[11])) {

				elem->pad = 120;
				elem->tem = 7;

			}
			else if ( (!elem->edge_ref[0]) && (!elem->edge_ref[1]) && (!elem->edge_ref[2]) && (!elem->edge_ref[3]) &&
					(!elem->edge_ref[4]) && (!elem->edge_ref[5]) && (elem->edge_ref[6]) && (!elem->edge_ref[7]) &&
					(!elem->edge_ref[8]) && (elem->edge_ref[9]) && (elem->edge_ref[10]) && (!elem->edge_ref[11])) {

				elem->pad = 121;
				elem->tem = 7;

			}
			else if ( (!elem->edge_ref[0]) && (!elem->edge_ref[1]) && (!elem->edge_ref[2]) && (!elem->edge_ref[3]) &&
					(!elem->edge_ref[4]) && (!elem->edge_ref[5]) && (!elem->edge_ref[6]) && (elem->edge_ref[7]) &&
					(!elem->edge_ref[8]) && (!elem->edge_ref[9]) && (elem->edge_ref[10]) && (elem->edge_ref[11])) {

				elem->pad = 122;
				elem->tem = 7;

			}
		}
		//template 8
		if(ed_cont == 8){
			if ( (elem->edge_ref[0]) && (elem->edge_ref[1]) && (elem->edge_ref[2]) && (elem->edge_ref[3]) &&
					(!elem->edge_ref[4]) && (!elem->edge_ref[5]) && (!elem->edge_ref[6]) && (!elem->edge_ref[7]) &&
					(elem->edge_ref[8]) && (elem->edge_ref[9]) && (elem->edge_ref[10]) && (elem->edge_ref[11])) {

				elem->pad = 123;
				elem->tem = 8;

			}
			else if ( (elem->edge_ref[0]) && (!elem->edge_ref[1]) && (elem->edge_ref[2]) && (!elem->edge_ref[3]) &&
					(elem->edge_ref[4]) && (elem->edge_ref[5]) && (elem->edge_ref[6]) && (elem->edge_ref[7]) &&
					(elem->edge_ref[8]) && (!elem->edge_ref[9]) && (elem->edge_ref[10]) && (!elem->edge_ref[11])) {

				elem->pad = 124;
				elem->tem = 8;

			}
			else if ( (!elem->edge_ref[0]) && (elem->edge_ref[1]) && (!elem->edge_ref[2]) && (elem->edge_ref[3]) &&
					(elem->edge_ref[4]) && (elem->edge_ref[5]) && (elem->edge_ref[6]) && (elem->edge_ref[7]) &&
					(!elem->edge_ref[8]) && (elem->edge_ref[9]) && (!elem->edge_ref[10]) && (elem->edge_ref[11])) {

				elem->pad = 125;
				elem->tem = 8;

			}
		}
		//template 9
		if(ed_cont == 8){
			if ( (elem->edge_ref[0]) && (elem->edge_ref[1]) && (elem->edge_ref[2]) && (elem->edge_ref[3]) &&
					(elem->edge_ref[4]) && (elem->edge_ref[5]) && (elem->edge_ref[6]) && (elem->edge_ref[7]) &&
					(!elem->edge_ref[8]) && (!elem->edge_ref[9]) && (!elem->edge_ref[10]) && (!elem->edge_ref[11])) {

				elem->pad = 126;
				elem->tem = 9;

			}
			else if ( (!elem->edge_ref[0]) && (!elem->edge_ref[1]) && (!elem->edge_ref[2]) && (!elem->edge_ref[3]) &&
					(elem->edge_ref[4]) && (elem->edge_ref[5]) && (elem->edge_ref[6]) && (elem->edge_ref[7]) &&
					(elem->edge_ref[8]) && (elem->edge_ref[9]) && (elem->edge_ref[10]) && (elem->edge_ref[11])) {

				elem->pad = 127;
				elem->tem = 9;

			}
			else if ( (elem->edge_ref[0]) && (elem->edge_ref[1]) && (!elem->edge_ref[2]) && (elem->edge_ref[3]) &&
					(elem->edge_ref[4]) && (elem->edge_ref[5]) && (!elem->edge_ref[6]) && (!elem->edge_ref[7]) &&
					(elem->edge_ref[8]) && (elem->edge_ref[9]) && (!elem->edge_ref[10]) && (elem->edge_ref[11])) {

				elem->pad = 128;
				elem->tem = 9;

			}
			else if ( (!elem->edge_ref[0]) && (elem->edge_ref[1]) && (elem->edge_ref[2]) && (elem->edge_ref[3]) &&
					(!elem->edge_ref[4]) && (!elem->edge_ref[5]) && (elem->edge_ref[6]) && (elem->edge_ref[7]) &&
					(!elem->edge_ref[8]) && (elem->edge_ref[9]) && (elem->edge_ref[10]) && (elem->edge_ref[11])) {

				elem->pad = 129;
				elem->tem = 9;

			}
			else if ( (elem->edge_ref[0]) && (elem->edge_ref[1]) && (elem->edge_ref[2]) && (!elem->edge_ref[3]) &&
					(!elem->edge_ref[4]) && (elem->edge_ref[5]) && (elem->edge_ref[6]) && (!elem->edge_ref[7]) &&
					(elem->edge_ref[8]) && (elem->edge_ref[9]) && (elem->edge_ref[10]) && (!elem->edge_ref[11])) {

				elem->pad = 130;
				elem->tem = 9;

			}
			else if ( (elem->edge_ref[0]) && (!elem->edge_ref[1]) && (elem->edge_ref[2]) && (elem->edge_ref[3]) &&
					(elem->edge_ref[4]) && (!elem->edge_ref[5]) && (!elem->edge_ref[6]) && (elem->edge_ref[7]) &&
					(elem->edge_ref[8]) && (!elem->edge_ref[9]) && (elem->edge_ref[10]) && (elem->edge_ref[11])) {

				elem->pad = 131;
				elem->tem = 9;

			}
		}
		//template 10
		if(ed_cont == 5){
			if ( (elem->edge_ref[0]) && (elem->edge_ref[1]) && (!elem->edge_ref[2]) && (elem->edge_ref[3]) &&
					(elem->edge_ref[4]) && (elem->edge_ref[5]) && (!elem->edge_ref[6]) && (!elem->edge_ref[7]) &&
					(!elem->edge_ref[8]) && (!elem->edge_ref[9]) && (!elem->edge_ref[10]) && (!elem->edge_ref[11])) {

				elem->pad = 132;
				elem->tem = 10;

			}
			else if ( (elem->edge_ref[0]) && (elem->edge_ref[1]) && (elem->edge_ref[2]) && (!elem->edge_ref[3]) &&
					(!elem->edge_ref[4]) && (elem->edge_ref[5]) && (elem->edge_ref[6]) && (!elem->edge_ref[7]) &&
					(!elem->edge_ref[8]) && (!elem->edge_ref[9]) && (!elem->edge_ref[10]) && (!elem->edge_ref[11])) {

				elem->pad = 133;
				elem->tem = 10;

			}
			else if ( (!elem->edge_ref[0]) && (elem->edge_ref[1]) && (elem->edge_ref[2]) && (elem->edge_ref[3]) &&
					(!elem->edge_ref[4]) && (!elem->edge_ref[5]) && (elem->edge_ref[6]) && (elem->edge_ref[7]) &&
					(!elem->edge_ref[8]) && (!elem->edge_ref[9]) && (!elem->edge_ref[10]) && (!elem->edge_ref[11])) {

				elem->pad = 134;
				elem->tem = 10;

			}
			else if ( (elem->edge_ref[0]) && (!elem->edge_ref[1]) && (elem->edge_ref[2]) && (elem->edge_ref[3]) &&
					(elem->edge_ref[4]) && (!elem->edge_ref[5]) && (!elem->edge_ref[6]) && (elem->edge_ref[7]) &&
					(!elem->edge_ref[8]) && (!elem->edge_ref[9]) && (!elem->edge_ref[10]) && (!elem->edge_ref[11])) {

				elem->pad = 135;
				elem->tem = 10;

			}
			else if ( (elem->edge_ref[0]) && (!elem->edge_ref[1]) && (!elem->edge_ref[2]) && (elem->edge_ref[3]) &&
					(elem->edge_ref[4]) && (!elem->edge_ref[5]) && (!elem->edge_ref[6]) && (!elem->edge_ref[7]) &&
					(elem->edge_ref[8]) && (!elem->edge_ref[9]) && (!elem->edge_ref[10]) && (elem->edge_ref[11])) {

				elem->pad = 136;
				elem->tem = 10;

			}
			else if ( (elem->edge_ref[0]) && (elem->edge_ref[1]) && (!elem->edge_ref[2]) && (!elem->edge_ref[3]) &&
					(!elem->edge_ref[4]) && (elem->edge_ref[5]) && (!elem->edge_ref[6]) && (!elem->edge_ref[7]) &&
					(elem->edge_ref[8]) && (elem->edge_ref[9]) && (!elem->edge_ref[10]) && (!elem->edge_ref[11])) {

				elem->pad = 137;
				elem->tem = 10;

			}
			else if ( (!elem->edge_ref[0]) && (elem->edge_ref[1]) && (elem->edge_ref[2]) && (!elem->edge_ref[3]) &&
					(!elem->edge_ref[4]) && (!elem->edge_ref[5]) && (elem->edge_ref[6]) && (!elem->edge_ref[7]) &&
					(!elem->edge_ref[8]) && (elem->edge_ref[9]) && (elem->edge_ref[10]) && (!elem->edge_ref[11])) {

				elem->pad = 138;
				elem->tem = 10;

			}
			else if ( (!elem->edge_ref[0]) && (!elem->edge_ref[1]) && (elem->edge_ref[2]) && (elem->edge_ref[3]) &&
					(!elem->edge_ref[4]) && (!elem->edge_ref[5]) && (!elem->edge_ref[6]) && (elem->edge_ref[7]) &&
					(!elem->edge_ref[8]) && (!elem->edge_ref[9]) && (elem->edge_ref[10]) && (elem->edge_ref[11])) {

				elem->pad = 139;
				elem->tem = 10;

			}
			else if ( (!elem->edge_ref[0]) && (!elem->edge_ref[1]) && (!elem->edge_ref[2]) && (!elem->edge_ref[3]) &&
					(elem->edge_ref[4]) && (elem->edge_ref[5]) && (!elem->edge_ref[6]) && (!elem->edge_ref[7]) &&
					(elem->edge_ref[8]) && (elem->edge_ref[9]) && (!elem->edge_ref[10]) && (elem->edge_ref[11])) {

				elem->pad = 140;
				elem->tem = 10;

			}
			else if ( (!elem->edge_ref[0]) && (!elem->edge_ref[1]) && (!elem->edge_ref[2]) && (!elem->edge_ref[3]) &&
					(!elem->edge_ref[4]) && (elem->edge_ref[5]) && (elem->edge_ref[6]) && (!elem->edge_ref[7]) &&
					(elem->edge_ref[8]) && (elem->edge_ref[9]) && (elem->edge_ref[10]) && (!elem->edge_ref[11])) {

				elem->pad = 141;
				elem->tem = 10;

			}
			else if ( (!elem->edge_ref[0]) && (!elem->edge_ref[1]) && (!elem->edge_ref[2]) && (!elem->edge_ref[3]) &&
					(!elem->edge_ref[4]) && (!elem->edge_ref[5]) && (elem->edge_ref[6]) && (elem->edge_ref[7]) &&
					(!elem->edge_ref[8]) && (elem->edge_ref[9]) && (elem->edge_ref[10]) && (elem->edge_ref[11])) {

				elem->pad = 142;
				elem->tem = 10;

			}
			else if ( (!elem->edge_ref[0]) && (!elem->edge_ref[1]) && (!elem->edge_ref[2]) && (!elem->edge_ref[3]) &&
					(elem->edge_ref[4]) && (!elem->edge_ref[5]) && (!elem->edge_ref[6]) && (elem->edge_ref[7]) &&
					(elem->edge_ref[8]) && (!elem->edge_ref[9]) && (elem->edge_ref[10]) && (elem->edge_ref[11])) {

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

	if(true){
		printf("case_1: %d\n",el_1);
		printf("case_2: %d\n",el_2);
		printf("case_3: %d\n",el_3);
		printf("case_4: %d\n",el_4);
		printf("case_5: %d\n",el_5);
		printf("case_6: %d\n",el_6);
		printf("case_7: %d\n",el_7);
		printf("case_8: %d\n",el_8);
		printf("case_9: %d\n",el_9);
		printf("case_10: %d\n",el_10);
		printf("case_11: %d\n",el_11);
		printf("sum: %d\n",su);
		printf("case_0: %d\n",el_0);
	}

	//for debug, work just in serial!!!
#if 0
	FILE * fdbg;
	char filename[80];
	sprintf(filename, "Edge_id_%04d.txt", mesh->mpi_rank);
	fdbg = fopen(filename, "w");

	for(int iel= 0; iel < mesh->total_n_elements; iel++ ){
		octant_t *elem = (octant_t*) sc_array_index(&mesh->elements, iel);
		fprintf(fdbg,"El: %d\n",iel);
		for (int edge = 0; edge < 12; ++edge) {
			fprintf(fdbg,"%d ",elem->edge_id[edge]);
		}
		fprintf(fdbg,"\n");
		for (int edge = 0; edge < 12; ++edge) {
			fprintf(fdbg,"%d ",elem->edge_ref[edge]);
		}
		fprintf(fdbg,"\n");
	}
	fclose(fdbg);
#endif

}

void Edge_identification(hexa_tree_t* mesh, std::vector<int>& elements_ids, sc_hash_array_t* hash) {

	for (int iel = 0; iel < elements_ids.size(); ++iel) {

		octant_t *elem = (octant_t*) sc_array_index(&mesh->elements, elements_ids[iel]);

		//template 11
		if(elem->pad==144){
                    for (int edge = 0; edge < 12; ++edge) {
                        edge_add(elem->edge_id[edge], hash );
                            //edges_ids.push_back(elem->edge_id[edge]);
                    }
		}
		//template 10
		else if(elem->pad==143){
                    edge_add(elem->edge_id[4], hash );
                    edge_add(elem->edge_id[7], hash );
                    edge_add(elem->edge_id[8], hash );
                    edge_add(elem->edge_id[10], hash );
                    edge_add(elem->edge_id[11], hash );
			//edges_ids.push_back(elem->edge_id[4]);
			//edges_ids.push_back(elem->edge_id[7]);
			//edges_ids.push_back(elem->edge_id[8]);
			//edges_ids.push_back(elem->edge_id[10]);
			//edges_ids.push_back(elem->edge_id[11]);
		}else if(elem->pad==142){
                    edge_add(elem->edge_id[6], hash );
                    edge_add(elem->edge_id[7], hash );
                    edge_add(elem->edge_id[9], hash );
                    edge_add(elem->edge_id[10], hash );
                    edge_add(elem->edge_id[11], hash );
			//edges_ids.push_back(elem->edge_id[6]);
			//edges_ids.push_back(elem->edge_id[7]);
			//edges_ids.push_back(elem->edge_id[9]);
			//edges_ids.push_back(elem->edge_id[10]);
			//edges_ids.push_back(elem->edge_id[11]);
		}else if(elem->pad==141){
                    edge_add(elem->edge_id[5], hash );
                    edge_add(elem->edge_id[6], hash );
                    edge_add(elem->edge_id[8], hash );
                    edge_add(elem->edge_id[9], hash );
                    edge_add(elem->edge_id[10], hash );
			//edges_ids.push_back(elem->edge_id[5]);
			//edges_ids.push_back(elem->edge_id[6]);
			//edges_ids.push_back(elem->edge_id[8]);
			//edges_ids.push_back(elem->edge_id[9]);
			//edges_ids.push_back(elem->edge_id[10]);
		}else if(elem->pad==140){
                    edge_add(elem->edge_id[4], hash );
                    edge_add(elem->edge_id[5], hash );
                    edge_add(elem->edge_id[8], hash );
                    edge_add(elem->edge_id[9], hash );
                    edge_add(elem->edge_id[11], hash );
			//edges_ids.push_back(elem->edge_id[4]);
			//edges_ids.push_back(elem->edge_id[5]);
			//edges_ids.push_back(elem->edge_id[8]);
			//edges_ids.push_back(elem->edge_id[9]);
			//edges_ids.push_back(elem->edge_id[11]);
		}else if(elem->pad==139){
                    edge_add(elem->edge_id[2], hash );
                    edge_add(elem->edge_id[3], hash );
                    edge_add(elem->edge_id[7], hash );
                    edge_add(elem->edge_id[10], hash );
                    edge_add(elem->edge_id[11], hash );
			//edges_ids.push_back(elem->edge_id[2]);
			//edges_ids.push_back(elem->edge_id[3]);
			//edges_ids.push_back(elem->edge_id[7]);
			//edges_ids.push_back(elem->edge_id[10]);
			//edges_ids.push_back(elem->edge_id[11]);
		}else if(elem->pad==138){                    
                    edge_add(elem->edge_id[1], hash );
                    edge_add(elem->edge_id[2], hash );
                    edge_add(elem->edge_id[6], hash );
                    edge_add(elem->edge_id[9], hash );
                    edge_add(elem->edge_id[10], hash );
			//edges_ids.push_back(elem->edge_id[1]);
			//edges_ids.push_back(elem->edge_id[2]);
			//edges_ids.push_back(elem->edge_id[6]);
			//edges_ids.push_back(elem->edge_id[9]);
			//edges_ids.push_back(elem->edge_id[10]);
		}else if(elem->pad==137){
                    edge_add(elem->edge_id[0], hash );
                    edge_add(elem->edge_id[1], hash );
                    edge_add(elem->edge_id[5], hash );
                    edge_add(elem->edge_id[8], hash );
                    edge_add(elem->edge_id[9], hash );
			//edges_ids.push_back(elem->edge_id[0]);
			//edges_ids.push_back(elem->edge_id[1]);
			//edges_ids.push_back(elem->edge_id[5]);
			//edges_ids.push_back(elem->edge_id[8]);
			//edges_ids.push_back(elem->edge_id[9]);
		}else if(elem->pad==136){
                    edge_add(elem->edge_id[0], hash );
                    edge_add(elem->edge_id[3], hash );
                    edge_add(elem->edge_id[4], hash );
                    edge_add(elem->edge_id[8], hash );
                    edge_add(elem->edge_id[11], hash );
			//edges_ids.push_back(elem->edge_id[0]);
			//edges_ids.push_back(elem->edge_id[3]);
			//edges_ids.push_back(elem->edge_id[4]);
			//edges_ids.push_back(elem->edge_id[8]);
			//edges_ids.push_back(elem->edge_id[11]);
		}else if(elem->pad==135){
                    edge_add(elem->edge_id[0], hash );
                    edge_add(elem->edge_id[2], hash );
                    edge_add(elem->edge_id[3], hash );
                    edge_add(elem->edge_id[4], hash );
                    edge_add(elem->edge_id[7], hash );
			//edges_ids.push_back(elem->edge_id[0]);
			//edges_ids.push_back(elem->edge_id[2]);
			//edges_ids.push_back(elem->edge_id[3]);
			//edges_ids.push_back(elem->edge_id[4]);
			//edges_ids.push_back(elem->edge_id[7]);
		}else if(elem->pad==134){
                    edge_add(elem->edge_id[1], hash );
                    edge_add(elem->edge_id[2], hash );
                    edge_add(elem->edge_id[3], hash );
                    edge_add(elem->edge_id[6], hash );
                    edge_add(elem->edge_id[7], hash );
			//edges_ids.push_back(elem->edge_id[1]);
			//edges_ids.push_back(elem->edge_id[2]);
			//edges_ids.push_back(elem->edge_id[3]);
			//edges_ids.push_back(elem->edge_id[6]);
			//edges_ids.push_back(elem->edge_id[7]);
		}else if(elem->pad==133){
                    edge_add(elem->edge_id[0], hash );
                    edge_add(elem->edge_id[1], hash );
                    edge_add(elem->edge_id[2], hash );
                    edge_add(elem->edge_id[5], hash );
                    edge_add(elem->edge_id[6], hash );
			//edges_ids.push_back(elem->edge_id[0]);
			//edges_ids.push_back(elem->edge_id[1]);
			//edges_ids.push_back(elem->edge_id[2]);
			//edges_ids.push_back(elem->edge_id[5]);
			//edges_ids.push_back(elem->edge_id[6]);
		}else if(elem->pad==132){
                    edge_add(elem->edge_id[0], hash );
                    edge_add(elem->edge_id[1], hash );
                    edge_add(elem->edge_id[3], hash );
                    edge_add(elem->edge_id[4], hash );
                    edge_add(elem->edge_id[5], hash );
			//edges_ids.push_back(elem->edge_id[0]);
			//edges_ids.push_back(elem->edge_id[1]);
			//edges_ids.push_back(elem->edge_id[3]);
			//edges_ids.push_back(elem->edge_id[4]);
			//edges_ids.push_back(elem->edge_id[5]);
		}
		//template 9
		else if(elem->pad==131){
                    edge_add(elem->edge_id[0], hash );
                    edge_add(elem->edge_id[2], hash );
                    edge_add(elem->edge_id[3], hash );
                    edge_add(elem->edge_id[4], hash );
                    edge_add(elem->edge_id[7], hash );
                    edge_add(elem->edge_id[8], hash );
                    edge_add(elem->edge_id[10], hash );
                    edge_add(elem->edge_id[11], hash );
			//edges_ids.push_back(elem->edge_id[0]);
			//edges_ids.push_back(elem->edge_id[2]);
			//edges_ids.push_back(elem->edge_id[3]);
			//edges_ids.push_back(elem->edge_id[4]);
			//edges_ids.push_back(elem->edge_id[7]);
			//edges_ids.push_back(elem->edge_id[8]);
			//edges_ids.push_back(elem->edge_id[10]);
			//edges_ids.push_back(elem->edge_id[11]);
		}else if(elem->pad==130){
                    edge_add(elem->edge_id[0], hash );
                    edge_add(elem->edge_id[1], hash );
                    edge_add(elem->edge_id[3], hash );
                    edge_add(elem->edge_id[5], hash );
                    edge_add(elem->edge_id[6], hash );
                    edge_add(elem->edge_id[8], hash );
                    edge_add(elem->edge_id[9], hash );
                    edge_add(elem->edge_id[11], hash );
			//edges_ids.push_back(elem->edge_id[0]);
			//edges_ids.push_back(elem->edge_id[1]);
			//edges_ids.push_back(elem->edge_id[2]);
			//edges_ids.push_back(elem->edge_id[5]);
			//edges_ids.push_back(elem->edge_id[6]);
			//edges_ids.push_back(elem->edge_id[8]);
			//edges_ids.push_back(elem->edge_id[9]);
			//edges_ids.push_back(elem->edge_id[10]);
		}else if(elem->pad==129){
                    edge_add(elem->edge_id[1], hash );
                    edge_add(elem->edge_id[2], hash );
                    edge_add(elem->edge_id[3], hash );
                    edge_add(elem->edge_id[6], hash );
                    edge_add(elem->edge_id[7], hash );
                    edge_add(elem->edge_id[9], hash );
                    edge_add(elem->edge_id[10], hash );
                    edge_add(elem->edge_id[11], hash );
			//edges_ids.push_back(elem->edge_id[1]);
			//edges_ids.push_back(elem->edge_id[2]);
			//edges_ids.push_back(elem->edge_id[3]);
			//edges_ids.push_back(elem->edge_id[6]);
			//edges_ids.push_back(elem->edge_id[7]);
			//edges_ids.push_back(elem->edge_id[9]);
			//edges_ids.push_back(elem->edge_id[10]);
			//edges_ids.push_back(elem->edge_id[11]);
		}else if(elem->pad==128){
                    edge_add(elem->edge_id[0], hash );
                    edge_add(elem->edge_id[1], hash );
                    edge_add(elem->edge_id[3], hash );
                    edge_add(elem->edge_id[5], hash );
                    edge_add(elem->edge_id[8], hash );
                    edge_add(elem->edge_id[9], hash );
                    edge_add(elem->edge_id[9], hash );
                    edge_add(elem->edge_id[11], hash );
			//edges_ids.push_back(elem->edge_id[0]);
			//edges_ids.push_back(elem->edge_id[1]);
			//edges_ids.push_back(elem->edge_id[3]);
			//edges_ids.push_back(elem->edge_id[4]);
			//edges_ids.push_back(elem->edge_id[5]);
			//edges_ids.push_back(elem->edge_id[8]);
			//edges_ids.push_back(elem->edge_id[9]);
			//edges_ids.push_back(elem->edge_id[11]);
		}else if(elem->pad==127){
                    edge_add(elem->edge_id[4], hash );
                    edge_add(elem->edge_id[5], hash );
                    edge_add(elem->edge_id[6], hash );
                    edge_add(elem->edge_id[7], hash );
                    edge_add(elem->edge_id[8], hash );
                    edge_add(elem->edge_id[9], hash );
                    edge_add(elem->edge_id[10], hash );
                    edge_add(elem->edge_id[11], hash );
			//edges_ids.push_back(elem->edge_id[4]);
			//edges_ids.push_back(elem->edge_id[5]);
			//edges_ids.push_back(elem->edge_id[6]);
			//edges_ids.push_back(elem->edge_id[7]);
			//edges_ids.push_back(elem->edge_id[8]);
			//edges_ids.push_back(elem->edge_id[9]);
			//edges_ids.push_back(elem->edge_id[10]);
			//edges_ids.push_back(elem->edge_id[11]);
		}else if(elem->pad==126){
                    edge_add(elem->edge_id[0], hash );
                    edge_add(elem->edge_id[1], hash );
                    edge_add(elem->edge_id[2], hash );
                    edge_add(elem->edge_id[3], hash );
                    edge_add(elem->edge_id[4], hash );
                    edge_add(elem->edge_id[5], hash );
                    edge_add(elem->edge_id[6], hash );
                    edge_add(elem->edge_id[7], hash );
			//edges_ids.push_back(elem->edge_id[0]);
			//edges_ids.push_back(elem->edge_id[1]);
			//edges_ids.push_back(elem->edge_id[2]);
			//edges_ids.push_back(elem->edge_id[3]);
			//edges_ids.push_back(elem->edge_id[4]);
			//edges_ids.push_back(elem->edge_id[5]);
			//edges_ids.push_back(elem->edge_id[6]);
			//edges_ids.push_back(elem->edge_id[7]);
		}
		//template 8
		else if(elem->pad==125){
                    edge_add(elem->edge_id[1], hash );
                    edge_add(elem->edge_id[3], hash );
                    edge_add(elem->edge_id[4], hash );
                    edge_add(elem->edge_id[5], hash );
                    edge_add(elem->edge_id[6], hash );
                    edge_add(elem->edge_id[7], hash );
                    edge_add(elem->edge_id[9], hash );
                    edge_add(elem->edge_id[11], hash );
			//edges_ids.push_back(elem->edge_id[1]);
			//edges_ids.push_back(elem->edge_id[3]);
			//edges_ids.push_back(elem->edge_id[4]);
			//edges_ids.push_back(elem->edge_id[5]);
			//edges_ids.push_back(elem->edge_id[6]);
			//edges_ids.push_back(elem->edge_id[7]);
			//edges_ids.push_back(elem->edge_id[9]);
			//edges_ids.push_back(elem->edge_id[11]);
		}else if(elem->pad==124){
                    edge_add(elem->edge_id[0], hash );
                    edge_add(elem->edge_id[2], hash );
                    edge_add(elem->edge_id[4], hash );
                    edge_add(elem->edge_id[5], hash );
                    edge_add(elem->edge_id[6], hash );
                    edge_add(elem->edge_id[7], hash );
                    edge_add(elem->edge_id[8], hash );
                    edge_add(elem->edge_id[10], hash );
			//edges_ids.push_back(elem->edge_id[0]);
			//edges_ids.push_back(elem->edge_id[2]);
			//edges_ids.push_back(elem->edge_id[4]);
			//edges_ids.push_back(elem->edge_id[5]);
			//edges_ids.push_back(elem->edge_id[6]);
			//edges_ids.push_back(elem->edge_id[7]);
			//edges_ids.push_back(elem->edge_id[8]);
			//edges_ids.push_back(elem->edge_id[10]);
		}else if(elem->pad==123){
                    edge_add(elem->edge_id[0], hash );
                    edge_add(elem->edge_id[1], hash );
                    edge_add(elem->edge_id[2], hash );
                    edge_add(elem->edge_id[3], hash );
                    edge_add(elem->edge_id[8], hash );
                    edge_add(elem->edge_id[9], hash );
                    edge_add(elem->edge_id[10], hash );
                    edge_add(elem->edge_id[11], hash );
			//edges_ids.push_back(elem->edge_id[0]);
			//edges_ids.push_back(elem->edge_id[1]);
			//edges_ids.push_back(elem->edge_id[2]);
			//edges_ids.push_back(elem->edge_id[3]);
			//edges_ids.push_back(elem->edge_id[8]);
			//edges_ids.push_back(elem->edge_id[9]);
			//edges_ids.push_back(elem->edge_id[10]);
			//edges_ids.push_back(elem->edge_id[11]);
		}
		//template 7
		else if(elem->pad==122){
                    edge_add(elem->edge_id[7], hash );
                    edge_add(elem->edge_id[10], hash );
                    edge_add(elem->edge_id[11], hash );
			//edges_ids.push_back(elem->edge_id[7]);
			//edges_ids.push_back(elem->edge_id[10]);
			//edges_ids.push_back(elem->edge_id[11]);
		}else if(elem->pad==121){
                    edge_add(elem->edge_id[6], hash );
                    edge_add(elem->edge_id[9], hash );
                    edge_add(elem->edge_id[10], hash );
			//edges_ids.push_back(elem->edge_id[6]);
			//edges_ids.push_back(elem->edge_id[9]);
			//edges_ids.push_back(elem->edge_id[10]);
		}else if(elem->pad==120){
                    edge_add(elem->edge_id[5], hash );
                    edge_add(elem->edge_id[8], hash );
                    edge_add(elem->edge_id[9], hash );
			//edges_ids.push_back(elem->edge_id[5]);
			//edges_ids.push_back(elem->edge_id[8]);
			//edges_ids.push_back(elem->edge_id[9]);
		}else if(elem->pad==119){
                    edge_add(elem->edge_id[4], hash );
                    edge_add(elem->edge_id[8], hash );
                    edge_add(elem->edge_id[11], hash );
			//edges_ids.push_back(elem->edge_id[4]);
			//edges_ids.push_back(elem->edge_id[8]);
			//edges_ids.push_back(elem->edge_id[11]);
		}else if(elem->pad==118){
                    edge_add(elem->edge_id[2], hash );
                    edge_add(elem->edge_id[3], hash );
                    edge_add(elem->edge_id[7], hash );
			//edges_ids.push_back(elem->edge_id[2]);
			//edges_ids.push_back(elem->edge_id[3]);
			//edges_ids.push_back(elem->edge_id[7]);
		}else if(elem->pad==117){
                    edge_add(elem->edge_id[1], hash );
                    edge_add(elem->edge_id[2], hash );
                    edge_add(elem->edge_id[6], hash );
			//edges_ids.push_back(elem->edge_id[1]);
			//edges_ids.push_back(elem->edge_id[2]);
			//edges_ids.push_back(elem->edge_id[6]);
		}else if(elem->pad==116){
                    edge_add(elem->edge_id[0], hash );
                    edge_add(elem->edge_id[1], hash );
                    edge_add(elem->edge_id[5], hash );
			//edges_ids.push_back(elem->edge_id[0]);
			//edges_ids.push_back(elem->edge_id[1]);
			//edges_ids.push_back(elem->edge_id[5]);
		}else if(elem->pad==115){
                    edge_add(elem->edge_id[0], hash );
                    edge_add(elem->edge_id[3], hash );
                    edge_add(elem->edge_id[4], hash );
			//edges_ids.push_back(elem->edge_id[0]);
			//edges_ids.push_back(elem->edge_id[3]);
			//edges_ids.push_back(elem->edge_id[4]);
		}
		//template 6
		else if(elem->pad==114){
                    edge_add(elem->edge_id[1], hash );
                    edge_add(elem->edge_id[3], hash );
                    edge_add(elem->edge_id[6], hash );                    
                    edge_add(elem->edge_id[7], hash );
                    edge_add(elem->edge_id[9], hash );
                    edge_add(elem->edge_id[11], hash );
			//edges_ids.push_back(elem->edge_id[1]);
			//edges_ids.push_back(elem->edge_id[3]);
			//edges_ids.push_back(elem->edge_id[6]);
			//edges_ids.push_back(elem->edge_id[7]);
			//edges_ids.push_back(elem->edge_id[9]);
			//edges_ids.push_back(elem->edge_id[11]);
		}else if(elem->pad==113){
                    edge_add(elem->edge_id[4], hash );
                    edge_add(elem->edge_id[5], hash );
                    edge_add(elem->edge_id[6], hash );                    
                    edge_add(elem->edge_id[7], hash );
                    edge_add(elem->edge_id[9], hash );
                    edge_add(elem->edge_id[11], hash );
			//edges_ids.push_back(elem->edge_id[4]);
			//edges_ids.push_back(elem->edge_id[5]);
			//edges_ids.push_back(elem->edge_id[6]);
			//edges_ids.push_back(elem->edge_id[7]);
			//edges_ids.push_back(elem->edge_id[9]);
			//edges_ids.push_back(elem->edge_id[11]);
		}else if(elem->pad==112){
                    edge_add(elem->edge_id[1], hash );
                    edge_add(elem->edge_id[3], hash );
                    edge_add(elem->edge_id[4], hash );                    
                    edge_add(elem->edge_id[5], hash );
                    edge_add(elem->edge_id[8], hash );
                    edge_add(elem->edge_id[11], hash );
			//edges_ids.push_back(elem->edge_id[1]);
			//edges_ids.push_back(elem->edge_id[3]);
			//edges_ids.push_back(elem->edge_id[4]);
			//edges_ids.push_back(elem->edge_id[5]);
			//edges_ids.push_back(elem->edge_id[8]);
			//edges_ids.push_back(elem->edge_id[11]);
		}else if(elem->pad==111){
                    edge_add(elem->edge_id[1], hash );
                    edge_add(elem->edge_id[3], hash );
                    edge_add(elem->edge_id[4], hash );                    
                    edge_add(elem->edge_id[5], hash );
                    edge_add(elem->edge_id[6], hash );
                    edge_add(elem->edge_id[7], hash );
			//edges_ids.push_back(elem->edge_id[1]);
			//edges_ids.push_back(elem->edge_id[3]);
			//edges_ids.push_back(elem->edge_id[4]);
			//edges_ids.push_back(elem->edge_id[5]);
			//edges_ids.push_back(elem->edge_id[6]);
			//edges_ids.push_back(elem->edge_id[7]);
		}else if(elem->pad==110){
                    edge_add(elem->edge_id[0], hash );
                    edge_add(elem->edge_id[2], hash );
                    edge_add(elem->edge_id[4], hash );                    
                    edge_add(elem->edge_id[5], hash );
                    edge_add(elem->edge_id[6], hash );
                    edge_add(elem->edge_id[7], hash );
			//edges_ids.push_back(elem->edge_id[0]);
			//edges_ids.push_back(elem->edge_id[2]);
			//edges_ids.push_back(elem->edge_id[4]);
			//edges_ids.push_back(elem->edge_id[5]);
			//edges_ids.push_back(elem->edge_id[6]);
			//edges_ids.push_back(elem->edge_id[7]);
		}else if(elem->pad==109){
                    edge_add(elem->edge_id[0], hash );
                    edge_add(elem->edge_id[2], hash );
                    edge_add(elem->edge_id[4], hash );                    
                    edge_add(elem->edge_id[7], hash );
                    edge_add(elem->edge_id[8], hash );
                    edge_add(elem->edge_id[10], hash );
			//edges_ids.push_back(elem->edge_id[0]);
			//edges_ids.push_back(elem->edge_id[2]);
			//edges_ids.push_back(elem->edge_id[4]);
			//edges_ids.push_back(elem->edge_id[7]);
			//edges_ids.push_back(elem->edge_id[8]);
			//edges_ids.push_back(elem->edge_id[10]);
		}else if(elem->pad==108){
                    edge_add(elem->edge_id[4], hash );
                    edge_add(elem->edge_id[5], hash );
                    edge_add(elem->edge_id[6], hash );                    
                    edge_add(elem->edge_id[7], hash );
                    edge_add(elem->edge_id[8], hash );
                    edge_add(elem->edge_id[10], hash );
			//edges_ids.push_back(elem->edge_id[4]);
			//edges_ids.push_back(elem->edge_id[5]);
			//edges_ids.push_back(elem->edge_id[6]);
			//edges_ids.push_back(elem->edge_id[7]);
			//edges_ids.push_back(elem->edge_id[8]);
			//edges_ids.push_back(elem->edge_id[10]);
		}else if(elem->pad==107){
                    edge_add(elem->edge_id[0], hash );
                    edge_add(elem->edge_id[2], hash );
                    edge_add(elem->edge_id[5], hash );                    
                    edge_add(elem->edge_id[6], hash );
                    edge_add(elem->edge_id[8], hash );
                    edge_add(elem->edge_id[10], hash );
			//edges_ids.push_back(elem->edge_id[0]);
			//edges_ids.push_back(elem->edge_id[2]);
			//edges_ids.push_back(elem->edge_id[5]);
			//edges_ids.push_back(elem->edge_id[6]);
			//edges_ids.push_back(elem->edge_id[8]);
			//edges_ids.push_back(elem->edge_id[10]);
		}else if(elem->pad==106){
                    edge_add(elem->edge_id[0], hash );
                    edge_add(elem->edge_id[2], hash );
                    edge_add(elem->edge_id[3], hash );                    
                    edge_add(elem->edge_id[8], hash );
                    edge_add(elem->edge_id[10], hash );
                    edge_add(elem->edge_id[11], hash );
			//edges_ids.push_back(elem->edge_id[0]);
			//edges_ids.push_back(elem->edge_id[2]);
			//edges_ids.push_back(elem->edge_id[3]);
			//edges_ids.push_back(elem->edge_id[8]);
			//edges_ids.push_back(elem->edge_id[10]);
			//edges_ids.push_back(elem->edge_id[11]);
		}else if(elem->pad==105){
                    edge_add(elem->edge_id[1], hash );
                    edge_add(elem->edge_id[2], hash );
                    edge_add(elem->edge_id[3], hash );                    
                    edge_add(elem->edge_id[9], hash );
                    edge_add(elem->edge_id[10], hash );
                    edge_add(elem->edge_id[11], hash );
			//edges_ids.push_back(elem->edge_id[1]);
			//edges_ids.push_back(elem->edge_id[2]);
			//edges_ids.push_back(elem->edge_id[3]);
			//edges_ids.push_back(elem->edge_id[9]);
			//edges_ids.push_back(elem->edge_id[10]);
			//edges_ids.push_back(elem->edge_id[11]);
		}else if(elem->pad==104){
                    edge_add(elem->edge_id[0], hash );
                    edge_add(elem->edge_id[1], hash );
                    edge_add(elem->edge_id[2], hash );                    
                    edge_add(elem->edge_id[8], hash );
                    edge_add(elem->edge_id[9], hash );
                    edge_add(elem->edge_id[10], hash );
			//edges_ids.push_back(elem->edge_id[0]);
			//edges_ids.push_back(elem->edge_id[1]);
			//edges_ids.push_back(elem->edge_id[2]);
			//edges_ids.push_back(elem->edge_id[8]);
			//edges_ids.push_back(elem->edge_id[9]);
			//edges_ids.push_back(elem->edge_id[10]);
		}else if(elem->pad==103){
                    edge_add(elem->edge_id[0], hash );
                    edge_add(elem->edge_id[1], hash );
                    edge_add(elem->edge_id[3], hash );                    
                    edge_add(elem->edge_id[8], hash );
                    edge_add(elem->edge_id[9], hash );
                    edge_add(elem->edge_id[11], hash );
			//edges_ids.push_back(elem->edge_id[0]);
			//edges_ids.push_back(elem->edge_id[1]);
			//edges_ids.push_back(elem->edge_id[3]);
			//edges_ids.push_back(elem->edge_id[8]);
			//edges_ids.push_back(elem->edge_id[9]);
			//edges_ids.push_back(elem->edge_id[11]);
		}
		//template 5
		else if(elem->pad==102){
                    edge_add(elem->edge_id[3], hash );
                    edge_add(elem->edge_id[7], hash );
			//edges_ids.push_back(elem->edge_id[3]);
			//edges_ids.push_back(elem->edge_id[7]);
		}else if(elem->pad==101){
                    edge_add(elem->edge_id[6], hash );
                    edge_add(elem->edge_id[7], hash );
			//edges_ids.push_back(elem->edge_id[6]);
			//edges_ids.push_back(elem->edge_id[7]);
		}else if(elem->pad==100){
                    edge_add(elem->edge_id[4], hash );
                    edge_add(elem->edge_id[6], hash );
			//edges_ids.push_back(elem->edge_id[4]);
			//edges_ids.push_back(elem->edge_id[6]);
		}else if(elem->pad==99){
                    edge_add(elem->edge_id[3], hash );
                    edge_add(elem->edge_id[4], hash );
			//edges_ids.push_back(elem->edge_id[3]);
			//edges_ids.push_back(elem->edge_id[4]);
		}else if(elem->pad==98){
                    edge_add(elem->edge_id[2], hash );
                    edge_add(elem->edge_id[7], hash );
			//edges_ids.push_back(elem->edge_id[2]);
			//edges_ids.push_back(elem->edge_id[7]);
		}else if(elem->pad==97){
                    edge_add(elem->edge_id[7], hash );
                    edge_add(elem->edge_id[10], hash );
			//edges_ids.push_back(elem->edge_id[7]);
			//edges_ids.push_back(elem->edge_id[10]);
		}else if(elem->pad==96){
                    edge_add(elem->edge_id[6], hash );
                    edge_add(elem->edge_id[10], hash );
			//edges_ids.push_back(elem->edge_id[6]);
			//edges_ids.push_back(elem->edge_id[10]);
		}else if(elem->pad==95){
                    edge_add(elem->edge_id[2], hash );
                    edge_add(elem->edge_id[6], hash );
			//edges_ids.push_back(elem->edge_id[2]);
			//edges_ids.push_back(elem->edge_id[6]);
		}else if(elem->pad==94){
                    edge_add(elem->edge_id[1], hash );
                    edge_add(elem->edge_id[5], hash );
			//edges_ids.push_back(elem->edge_id[1]);
			//edges_ids.push_back(elem->edge_id[5]);
		}else if(elem->pad==93){
                    edge_add(elem->edge_id[1], hash );
                    edge_add(elem->edge_id[6], hash );
			//edges_ids.push_back(elem->edge_id[1]);
			//edges_ids.push_back(elem->edge_id[6]);
		}else if(elem->pad==92){
                    edge_add(elem->edge_id[6], hash );
                    edge_add(elem->edge_id[9], hash );
			//edges_ids.push_back(elem->edge_id[6]);
			//edges_ids.push_back(elem->edge_id[9]);
		}else if(elem->pad==91){
                    edge_add(elem->edge_id[5], hash );
                    edge_add(elem->edge_id[9], hash );
			//edges_ids.push_back(elem->edge_id[5]);
			//edges_ids.push_back(elem->edge_id[9]);
		}else if(elem->pad==90){
                    edge_add(elem->edge_id[0], hash );
                    edge_add(elem->edge_id[4], hash );
			//edges_ids.push_back(elem->edge_id[0]);
			//edges_ids.push_back(elem->edge_id[4]);
		}else if(elem->pad==89){
                    edge_add(elem->edge_id[0], hash );
                    edge_add(elem->edge_id[5], hash );
			//edges_ids.push_back(elem->edge_id[0]);
			//edges_ids.push_back(elem->edge_id[5]);
		}else if(elem->pad==88){
                    edge_add(elem->edge_id[5], hash );
                    edge_add(elem->edge_id[8], hash );
			//edges_ids.push_back(elem->edge_id[5]);
			//edges_ids.push_back(elem->edge_id[8]);
		}else if(elem->pad==87){
                    edge_add(elem->edge_id[8], hash );
                    edge_add(elem->edge_id[4], hash );
			//edges_ids.push_back(elem->edge_id[8]);
			//edges_ids.push_back(elem->edge_id[4]);
		}else if(elem->pad==86){
                    edge_add(elem->edge_id[8], hash );
                    edge_add(elem->edge_id[11], hash );
			//edges_ids.push_back(elem->edge_id[8]);
			//edges_ids.push_back(elem->edge_id[11]);
		}else if(elem->pad==85){
                    edge_add(elem->edge_id[10], hash );
                    edge_add(elem->edge_id[11], hash );
			//edges_ids.push_back(elem->edge_id[10]);
			//edges_ids.push_back(elem->edge_id[11]);
		}else if(elem->pad==84){
                    edge_add(elem->edge_id[9], hash );
                    edge_add(elem->edge_id[10], hash );
			//edges_ids.push_back(elem->edge_id[9]);
			//edges_ids.push_back(elem->edge_id[10]);
		}else if(elem->pad==83){
                    edge_add(elem->edge_id[8], hash );
                    edge_add(elem->edge_id[9], hash );
			//edges_ids.push_back(elem->edge_id[8]);
			//edges_ids.push_back(elem->edge_id[9]);
		}else if(elem->pad==82){
                    edge_add(elem->edge_id[0], hash );
                    edge_add(elem->edge_id[3], hash );
			//edges_ids.push_back(elem->edge_id[0]);
			//edges_ids.push_back(elem->edge_id[3]);
		}else if(elem->pad==81){
                    edge_add(elem->edge_id[2], hash );
                    edge_add(elem->edge_id[3], hash );
			//edges_ids.push_back(elem->edge_id[2]);
			//edges_ids.push_back(elem->edge_id[3]);
		}else if(elem->pad==80){
                    edge_add(elem->edge_id[1], hash );
                    edge_add(elem->edge_id[2], hash );
			//edges_ids.push_back(elem->edge_id[1]);
			//edges_ids.push_back(elem->edge_id[2]);
		}else if(elem->pad==79){
                    edge_add(elem->edge_id[0], hash );
                    edge_add(elem->edge_id[1], hash );
			//edges_ids.push_back(elem->edge_id[0]);
			//edges_ids.push_back(elem->edge_id[1]);
		}else if(elem->pad==78){
                    edge_add(elem->edge_id[3], hash );
                    edge_add(elem->edge_id[7], hash );
                    edge_add(elem->edge_id[11], hash );
			//edges_ids.push_back(elem->edge_id[3]);
			//edges_ids.push_back(elem->edge_id[7]);
			//edges_ids.push_back(elem->edge_id[11]);
		}else if(elem->pad==77){
                    edge_add(elem->edge_id[4], hash );
                    edge_add(elem->edge_id[7], hash );
                    edge_add(elem->edge_id[11], hash );
			//edges_ids.push_back(elem->edge_id[4]);
			//edges_ids.push_back(elem->edge_id[7]);
			//edges_ids.push_back(elem->edge_id[11]);
		}else if(elem->pad==76){
                    edge_add(elem->edge_id[3], hash );
                    edge_add(elem->edge_id[4], hash );
                    edge_add(elem->edge_id[11], hash );
			//edges_ids.push_back(elem->edge_id[3]);
			//edges_ids.push_back(elem->edge_id[4]);
			//edges_ids.push_back(elem->edge_id[11]);
		}else if(elem->pad==75){
                    edge_add(elem->edge_id[3], hash );
                    edge_add(elem->edge_id[4], hash );
                    edge_add(elem->edge_id[7], hash );
			//edges_ids.push_back(elem->edge_id[3]);
			//edges_ids.push_back(elem->edge_id[4]);
			//edges_ids.push_back(elem->edge_id[7]);
		}else if(elem->pad==74){
                    edge_add(elem->edge_id[2], hash );
                    edge_add(elem->edge_id[7], hash );
                    edge_add(elem->edge_id[10], hash );
			//edges_ids.push_back(elem->edge_id[2]);
			//edges_ids.push_back(elem->edge_id[7]);
			//edges_ids.push_back(elem->edge_id[10]);
		}else if(elem->pad==73){
                    edge_add(elem->edge_id[6], hash );
                    edge_add(elem->edge_id[7], hash );
                    edge_add(elem->edge_id[10], hash );
			//edges_ids.push_back(elem->edge_id[6]);
			//edges_ids.push_back(elem->edge_id[7]);
			//edges_ids.push_back(elem->edge_id[10]);
		}else if(elem->pad==72){
                    edge_add(elem->edge_id[2], hash );
                    edge_add(elem->edge_id[6], hash );
                    edge_add(elem->edge_id[10], hash );
			//edges_ids.push_back(elem->edge_id[2]);
			//edges_ids.push_back(elem->edge_id[6]);
			//edges_ids.push_back(elem->edge_id[10]);
		}else if(elem->pad==71){
                    edge_add(elem->edge_id[2], hash );
                    edge_add(elem->edge_id[6], hash );
                    edge_add(elem->edge_id[7], hash );
			//edges_ids.push_back(elem->edge_id[2]);
			//edges_ids.push_back(elem->edge_id[6]);
			//edges_ids.push_back(elem->edge_id[7]);
		}else if(elem->pad==70){
                    edge_add(elem->edge_id[1], hash );
                    edge_add(elem->edge_id[6], hash );
                    edge_add(elem->edge_id[9], hash );
			//edges_ids.push_back(elem->edge_id[1]);
			//edges_ids.push_back(elem->edge_id[6]);
			//edges_ids.push_back(elem->edge_id[9]);
		}else if(elem->pad==69){
                    edge_add(elem->edge_id[5], hash );
                    edge_add(elem->edge_id[6], hash );
                    edge_add(elem->edge_id[9], hash );
			//edges_ids.push_back(elem->edge_id[5]);
			//edges_ids.push_back(elem->edge_id[6]);
			//edges_ids.push_back(elem->edge_id[9]);
		}else if(elem->pad==68){
                    edge_add(elem->edge_id[1], hash );
                    edge_add(elem->edge_id[5], hash );
                    edge_add(elem->edge_id[9], hash );
			//edges_ids.push_back(elem->edge_id[1]);
			//edges_ids.push_back(elem->edge_id[5]);
			//edges_ids.push_back(elem->edge_id[9]);
		}else if(elem->pad==67){
                    edge_add(elem->edge_id[1], hash );
                    edge_add(elem->edge_id[5], hash );
                    edge_add(elem->edge_id[6], hash );
			//edges_ids.push_back(elem->edge_id[1]);
			//edges_ids.push_back(elem->edge_id[5]);
			//edges_ids.push_back(elem->edge_id[6]);
		}else if(elem->pad==66){
                    edge_add(elem->edge_id[0], hash );
                    edge_add(elem->edge_id[4], hash );
                    edge_add(elem->edge_id[8], hash );
			//edges_ids.push_back(elem->edge_id[0]);
			//edges_ids.push_back(elem->edge_id[4]);
			//edges_ids.push_back(elem->edge_id[8]);
		}else if(elem->pad==65){
                    edge_add(elem->edge_id[4], hash );
                    edge_add(elem->edge_id[5], hash );
                    edge_add(elem->edge_id[8], hash );
			//edges_ids.push_back(elem->edge_id[4]);
			//edges_ids.push_back(elem->edge_id[5]);
			//edges_ids.push_back(elem->edge_id[8]);
		}else if(elem->pad==64){
                    edge_add(elem->edge_id[0], hash );
                    edge_add(elem->edge_id[5], hash );
                    edge_add(elem->edge_id[8], hash );
			//edges_ids.push_back(elem->edge_id[0]);
			//edges_ids.push_back(elem->edge_id[5]);
			//edges_ids.push_back(elem->edge_id[8]);
		}else if(elem->pad==63){
                    edge_add(elem->edge_id[0], hash );
                    edge_add(elem->edge_id[4], hash );
                    edge_add(elem->edge_id[5], hash );
			//edges_ids.push_back(elem->edge_id[0]);
			//edges_ids.push_back(elem->edge_id[4]);
			//edges_ids.push_back(elem->edge_id[5]);
		}else if(elem->pad==62){
                    edge_add(elem->edge_id[8], hash );
                    edge_add(elem->edge_id[10], hash );
                    edge_add(elem->edge_id[11], hash );
			//edges_ids.push_back(elem->edge_id[8]);
			//edges_ids.push_back(elem->edge_id[10]);
			//edges_ids.push_back(elem->edge_id[11]);
		}else if(elem->pad==61){
                    edge_add(elem->edge_id[9], hash );
                    edge_add(elem->edge_id[10], hash );
                    edge_add(elem->edge_id[11], hash );
			//edges_ids.push_back(elem->edge_id[9]);
			//edges_ids.push_back(elem->edge_id[10]);
			//edges_ids.push_back(elem->edge_id[11]);
		}else if(elem->pad==60){
                    edge_add(elem->edge_id[8], hash );
                    edge_add(elem->edge_id[9], hash );
                    edge_add(elem->edge_id[11], hash );
			//edges_ids.push_back(elem->edge_id[8]);
			//edges_ids.push_back(elem->edge_id[9]);
			//edges_ids.push_back(elem->edge_id[11]);
		}else if(elem->pad==59){
                    edge_add(elem->edge_id[8], hash );
                    edge_add(elem->edge_id[9], hash );
                    edge_add(elem->edge_id[10], hash );
			//edges_ids.push_back(elem->edge_id[8]);
			//edges_ids.push_back(elem->edge_id[9]);
			//edges_ids.push_back(elem->edge_id[10]);
		}else if(elem->pad==58){
                    edge_add(elem->edge_id[0], hash );
                    edge_add(elem->edge_id[2], hash );
                    edge_add(elem->edge_id[11], hash );
			//edges_ids.push_back(elem->edge_id[0]);
			//edges_ids.push_back(elem->edge_id[2]);
			//edges_ids.push_back(elem->edge_id[3]);
		}else if(elem->pad==57){
                    edge_add(elem->edge_id[1], hash );
                    edge_add(elem->edge_id[2], hash );
                    edge_add(elem->edge_id[3], hash );
			//edges_ids.push_back(elem->edge_id[1]);
			//edges_ids.push_back(elem->edge_id[2]);
			//edges_ids.push_back(elem->edge_id[3]);
		}else if(elem->pad==56){
                    edge_add(elem->edge_id[0], hash );
                    edge_add(elem->edge_id[1], hash );
                    edge_add(elem->edge_id[2], hash );
			//edges_ids.push_back(elem->edge_id[0]);
			//edges_ids.push_back(elem->edge_id[1]);
			//edges_ids.push_back(elem->edge_id[2]);
		}else if(elem->pad==55){
                    edge_add(elem->edge_id[0], hash );
                    edge_add(elem->edge_id[1], hash );
                    edge_add(elem->edge_id[3], hash );
			//edges_ids.push_back(elem->edge_id[0]);
			//edges_ids.push_back(elem->edge_id[1]);
			//edges_ids.push_back(elem->edge_id[3]);
		}else if(elem->pad==54){
                    edge_add(elem->edge_id[3], hash );
                    edge_add(elem->edge_id[4], hash );
                    edge_add(elem->edge_id[7], hash );
                    edge_add(elem->edge_id[11], hash );
			//edges_ids.push_back(elem->edge_id[3]);
			//edges_ids.push_back(elem->edge_id[4]);
			//edges_ids.push_back(elem->edge_id[7]);
			//edges_ids.push_back(elem->edge_id[11]);
		}else if(elem->pad==53){
                    edge_add(elem->edge_id[2], hash );
                    edge_add(elem->edge_id[6], hash );
                    edge_add(elem->edge_id[7], hash );
                    edge_add(elem->edge_id[10], hash );
                        //edges_ids.push_back(elem->edge_id[2]);
			//edges_ids.push_back(elem->edge_id[6]);
			//edges_ids.push_back(elem->edge_id[7]);
			//edges_ids.push_back(elem->edge_id[10]);
		}else if(elem->pad==52){
                    edge_add(elem->edge_id[1], hash );
                    edge_add(elem->edge_id[5], hash );
                    edge_add(elem->edge_id[6], hash );
                    edge_add(elem->edge_id[9], hash );
			//edges_ids.push_back(elem->edge_id[1]);
			//edges_ids.push_back(elem->edge_id[5]);
			//edges_ids.push_back(elem->edge_id[6]);
			//edges_ids.push_back(elem->edge_id[9]);
		}else if(elem->pad==51){
                    edge_add(elem->edge_id[0], hash );
                    edge_add(elem->edge_id[4], hash );
                    edge_add(elem->edge_id[5], hash );
                    edge_add(elem->edge_id[8], hash );
			//edges_ids.push_back(elem->edge_id[0]);
			//edges_ids.push_back(elem->edge_id[4]);
			//edges_ids.push_back(elem->edge_id[5]);
			//edges_ids.push_back(elem->edge_id[8]);
		}else if(elem->pad==50){
                    edge_add(elem->edge_id[8], hash );
                    edge_add(elem->edge_id[9], hash );
                    edge_add(elem->edge_id[10], hash );
                    edge_add(elem->edge_id[11], hash );
			//edges_ids.push_back(elem->edge_id[8]);
			//edges_ids.push_back(elem->edge_id[9]);
			//edges_ids.push_back(elem->edge_id[10]);
			//edges_ids.push_back(elem->edge_id[11]);
		}else if(elem->pad==49){
                    edge_add(elem->edge_id[0], hash );
                    edge_add(elem->edge_id[1], hash );
                    edge_add(elem->edge_id[2], hash );
                    edge_add(elem->edge_id[3], hash );
			//edges_ids.push_back(elem->edge_id[0]);
			//edges_ids.push_back(elem->edge_id[1]);
			//edges_ids.push_back(elem->edge_id[2]);
			//edges_ids.push_back(elem->edge_id[3]);
		}
		//template 4
		else if(elem->pad==48){
                    edge_add(elem->edge_id[4], hash );
                    edge_add(elem->edge_id[5], hash );
                    edge_add(elem->edge_id[9], hash );
                    edge_add(elem->edge_id[11], hash );
			//edges_ids.push_back(elem->edge_id[4]);
			//edges_ids.push_back(elem->edge_id[5]);
			//edges_ids.push_back(elem->edge_id[9]);
			//edges_ids.push_back(elem->edge_id[11]);
		}else if(elem->pad==47){
                    edge_add(elem->edge_id[6], hash );
                    edge_add(elem->edge_id[7], hash );
                    edge_add(elem->edge_id[9], hash );
                    edge_add(elem->edge_id[11], hash );
			//edges_ids.push_back(elem->edge_id[6]);
			//edges_ids.push_back(elem->edge_id[7]);
			//edges_ids.push_back(elem->edge_id[9]);
			//edges_ids.push_back(elem->edge_id[11]);
		}else if(elem->pad==46){
                    edge_add(elem->edge_id[1], hash );
                    edge_add(elem->edge_id[3], hash );
                    edge_add(elem->edge_id[6], hash );
                    edge_add(elem->edge_id[7], hash );
			//edges_ids.push_back(elem->edge_id[1]);
			//edges_ids.push_back(elem->edge_id[3]);
			//edges_ids.push_back(elem->edge_id[6]);
			//edges_ids.push_back(elem->edge_id[7]);
		}else if(elem->pad==45){
                    edge_add(elem->edge_id[1], hash );
                    edge_add(elem->edge_id[3], hash );
                    edge_add(elem->edge_id[4], hash );
                    edge_add(elem->edge_id[5], hash );
			//edges_ids.push_back(elem->edge_id[1]);
			//edges_ids.push_back(elem->edge_id[3]);
			//edges_ids.push_back(elem->edge_id[4]);
			//edges_ids.push_back(elem->edge_id[5]);
		}else if(elem->pad==44){
                    edge_add(elem->edge_id[4], hash );
                    edge_add(elem->edge_id[7], hash );
                    edge_add(elem->edge_id[8], hash );
                    edge_add(elem->edge_id[10], hash );
			//edges_ids.push_back(elem->edge_id[4]);
			//edges_ids.push_back(elem->edge_id[7]);
			//edges_ids.push_back(elem->edge_id[8]);
			//edges_ids.push_back(elem->edge_id[10]);
		}else if(elem->pad==43){
                    edge_add(elem->edge_id[5], hash );
                    edge_add(elem->edge_id[6], hash );
                    edge_add(elem->edge_id[8], hash );
                    edge_add(elem->edge_id[10], hash );
			//edges_ids.push_back(elem->edge_id[5]);
			//edges_ids.push_back(elem->edge_id[6]);
			//edges_ids.push_back(elem->edge_id[8]);
			//edges_ids.push_back(elem->edge_id[10]);
		}else if(elem->pad==42){
                    edge_add(elem->edge_id[0], hash );
                    edge_add(elem->edge_id[2], hash );
                    edge_add(elem->edge_id[5], hash );
                    edge_add(elem->edge_id[6], hash );
			//edges_ids.push_back(elem->edge_id[0]);
			//edges_ids.push_back(elem->edge_id[2]);
			//edges_ids.push_back(elem->edge_id[5]);
			//edges_ids.push_back(elem->edge_id[6]);
		}else if(elem->pad==41){
                    edge_add(elem->edge_id[0], hash );
                    edge_add(elem->edge_id[2], hash );
                    edge_add(elem->edge_id[4], hash );
                    edge_add(elem->edge_id[7], hash );
			//edges_ids.push_back(elem->edge_id[0]);
			//edges_ids.push_back(elem->edge_id[2]);
			//edges_ids.push_back(elem->edge_id[4]);
			//edges_ids.push_back(elem->edge_id[7]);
		}else if(elem->pad==40){
                    edge_add(elem->edge_id[0], hash );
                    edge_add(elem->edge_id[3], hash );
                    edge_add(elem->edge_id[8], hash );
                    edge_add(elem->edge_id[11], hash );
			//edges_ids.push_back(elem->edge_id[0]);
			//edges_ids.push_back(elem->edge_id[3]);
			//edges_ids.push_back(elem->edge_id[8]);
			//edges_ids.push_back(elem->edge_id[11]);
		}else if(elem->pad==39){
                    edge_add(elem->edge_id[2], hash );
                    edge_add(elem->edge_id[3], hash );
                    edge_add(elem->edge_id[10], hash );
                    edge_add(elem->edge_id[11], hash );
			//edges_ids.push_back(elem->edge_id[2]);
			//edges_ids.push_back(elem->edge_id[3]);
			//edges_ids.push_back(elem->edge_id[10]);
			//edges_ids.push_back(elem->edge_id[11]);
		}else if(elem->pad==38){
                    edge_add(elem->edge_id[1], hash );
                    edge_add(elem->edge_id[2], hash );
                    edge_add(elem->edge_id[9], hash );
                    edge_add(elem->edge_id[10], hash );
			//edges_ids.push_back(elem->edge_id[1]);
			//edges_ids.push_back(elem->edge_id[2]);
			//edges_ids.push_back(elem->edge_id[9]);
			//edges_ids.push_back(elem->edge_id[10]);
		}else if(elem->pad==37){
                    edge_add(elem->edge_id[0], hash );
                    edge_add(elem->edge_id[1], hash );
                    edge_add(elem->edge_id[8], hash );
                    edge_add(elem->edge_id[9], hash );
			//edges_ids.push_back(elem->edge_id[0]);
			//edges_ids.push_back(elem->edge_id[1]);
			//edges_ids.push_back(elem->edge_id[8]);
			//edges_ids.push_back(elem->edge_id[9]);
		}
		//template 3
		else if(elem->pad==36){
                    edge_add(elem->edge_id[3], hash );
                    edge_add(elem->edge_id[11], hash );
			//edges_ids.push_back(elem->edge_id[3]);
			//edges_ids.push_back(elem->edge_id[11]);
		}else if(elem->pad==35){
                    edge_add(elem->edge_id[1], hash );
                    edge_add(elem->edge_id[10], hash );
			//edges_ids.push_back(elem->edge_id[1]);
			//edges_ids.push_back(elem->edge_id[10]);
		}else if(elem->pad==34){
                    edge_add(elem->edge_id[1], hash );
                    edge_add(elem->edge_id[9], hash );
			//edges_ids.push_back(elem->edge_id[1]);
			//edges_ids.push_back(elem->edge_id[9]);
		}else if(elem->pad==33){
                    edge_add(elem->edge_id[0], hash );
                    edge_add(elem->edge_id[8], hash );
			//edges_ids.push_back(elem->edge_id[0]);
			//edges_ids.push_back(elem->edge_id[8]);
		}else if(elem->pad==32){
                    edge_add(elem->edge_id[4], hash );
                    edge_add(elem->edge_id[7], hash );
			//edges_ids.push_back(elem->edge_id[4]);
			//edges_ids.push_back(elem->edge_id[7]);
		}else if(elem->pad==31){
                    edge_add(elem->edge_id[6], hash );
                    edge_add(elem->edge_id[7], hash );
			//edges_ids.push_back(elem->edge_id[6]);
			//edges_ids.push_back(elem->edge_id[7]);
		}else if(elem->pad==30){
                    edge_add(elem->edge_id[5], hash );
                    edge_add(elem->edge_id[6], hash );
			//edges_ids.push_back(elem->edge_id[5]);
			//edges_ids.push_back(elem->edge_id[6]);
		}else if(elem->pad==29){
                    edge_add(elem->edge_id[4], hash );
                    edge_add(elem->edge_id[5], hash );
			//edges_ids.push_back(elem->edge_id[4]);
			//edges_ids.push_back(elem->edge_id[5]);
		}else if(elem->pad==28){
                    edge_add(elem->edge_id[9], hash );
                    edge_add(elem->edge_id[11], hash );
			//edges_ids.push_back(elem->edge_id[9]);
			//edges_ids.push_back(elem->edge_id[11]);
		}else if(elem->pad==27){
                    edge_add(elem->edge_id[8], hash );
                    edge_add(elem->edge_id[10], hash );
			//edges_ids.push_back(elem->edge_id[8]);
			//edges_ids.push_back(elem->edge_id[10]);
		}else if(elem->pad==26){
                    edge_add(elem->edge_id[1], hash );
                    edge_add(elem->edge_id[3], hash );
			//edges_ids.push_back(elem->edge_id[1]);
			//edges_ids.push_back(elem->edge_id[3]);
		}else if(elem->pad==25){
                    edge_add(elem->edge_id[0], hash );
                    edge_add(elem->edge_id[2], hash );
			//edges_ids.push_back(elem->edge_id[0]);
			//edges_ids.push_back(elem->edge_id[2]);
		}
		//template 2
		else if(elem->pad==24){
                    edge_add(elem->edge_id[4], hash );
                    edge_add(elem->edge_id[5], hash );
                    edge_add(elem->edge_id[6], hash );
                    edge_add(elem->edge_id[7], hash );
			//edges_ids.push_back(elem->edge_id[4]);
			//edges_ids.push_back(elem->edge_id[5]);
			//edges_ids.push_back(elem->edge_id[6]);
			//edges_ids.push_back(elem->edge_id[7]);
		}else if(elem->pad==23){
                    edge_add(elem->edge_id[1], hash );
                    edge_add(elem->edge_id[3], hash );
                    edge_add(elem->edge_id[9], hash );
                    edge_add(elem->edge_id[11], hash );
			//edges_ids.push_back(elem->edge_id[1]);
			//edges_ids.push_back(elem->edge_id[3]);
			//edges_ids.push_back(elem->edge_id[9]);
			//edges_ids.push_back(elem->edge_id[11]);
		}else if(elem->pad==22){
                    edge_add(elem->edge_id[0], hash );
                    edge_add(elem->edge_id[2], hash );
                    edge_add(elem->edge_id[8], hash );
                    edge_add(elem->edge_id[10], hash );
			//edges_ids.push_back(elem->edge_id[0]);
			//edges_ids.push_back(elem->edge_id[2]);
			//edges_ids.push_back(elem->edge_id[8]);
			//edges_ids.push_back(elem->edge_id[10]);
		}
		//template 1
		else if(elem->pad==21){
                    edge_add(elem->edge_id[11], hash );
			//edges_ids.push_back(elem->edge_id[11]);
		}else if(elem->pad==20){
                    edge_add(elem->edge_id[10], hash );
			//edges_ids.push_back(elem->edge_id[10]);
		}else if(elem->pad==19){
                    edge_add(elem->edge_id[9], hash );
			//edges_ids.push_back(elem->edge_id[9]);
		}else if(elem->pad==18){
                    edge_add(elem->edge_id[8], hash );
			//edges_ids.push_back(elem->edge_id[8]);
		}else if(elem->pad==17){
                    edge_add(elem->edge_id[7], hash );
			//edges_ids.push_back(elem->edge_id[7]);
		}else if(elem->pad==16){
                    edge_add(elem->edge_id[6], hash );
			//edges_ids.push_back(elem->edge_id[6]);
		}else if(elem->pad==15){
                    edge_add(elem->edge_id[5], hash );
			//edges_ids.push_back(elem->edge_id[5]);
		}else if(elem->pad==14){
                    edge_add(elem->edge_id[4], hash );
			//edges_ids.push_back(elem->edge_id[4]);
		}else if(elem->pad==13){
                    edge_add(elem->edge_id[3], hash );
			//edges_ids.push_back(elem->edge_id[3]);
		}else if(elem->pad==12){
                    edge_add(elem->edge_id[2], hash );
			//edges_ids.push_back(elem->edge_id[2]);
		}else if(elem->pad==11){
                    edge_add(elem->edge_id[1], hash );
			//edges_ids.push_back(elem->edge_id[1]);
		}else if(elem->pad==10){
                    edge_add(elem->edge_id[0], hash );
			//edges_ids.push_back(elem->edge_id[0]);
		}
	}

	//TODO check the time of the two implementations
	//printf(" Edges intercepted: %d\n", edges_ids.size());

	//remove duplicates edges
	//std::set<int> s;
	//unsigned size = edges_ids.size();
	//for( unsigned i = 0; i < size; ++i ) s.insert( edges_ids[i] );
	//edges_ids.assign( s.begin(), s.end() );
	//printf("Cleaned Edges intercepted: %d\n", edges_ids.size());

	//std::sort( edges_ids.begin(), edges_ids.end() );
	//edges_ids.erase( std::unique( edges_ids.begin(), edges_ids.end() ), edges_ids.end() );
	//printf("Cleaned Edges intercepted: %d\n", edges_ids.size());

}

void Edge_propagation(hexa_tree_t* mesh, std::vector<int>& elements_ids, sc_hash_array_t* hash_id) {

    size_t* position;
    
	for(int iel= 0; iel < mesh->total_n_elements; iel++ ){
		octant_t *elem = (octant_t*) sc_array_index(&mesh->elements, iel);

		for (int edge = 0; edge < 12; ++edge) {
                    bool out = false;
                    out =  sc_hash_array_lookup(hash_id, &elem->edge_id[edge], position);	
                    
                    if(out){
                        elem->edge_ref[edge]=true;
                        elem->pad = -1;
                        elements_ids.push_back(iel);
                    }
                    /*
			if(binary_search(edges_ids.begin(), edges_ids.end(), elem->edge_id[edge])){
				elem->edge_ref[edge]=true;
				elem->pad = -1;
				elements_ids.push_back(iel);
			}
                    */
		}
	}
	//TODO segfault in MPI
	//cleaning the element vector
	std::sort( elements_ids.begin(), elements_ids.end() );
	elements_ids.erase( std::unique( elements_ids.begin(), elements_ids.end() ), elements_ids.end() );
}

void CheckOctreeTemplate(hexa_tree_t* mesh, const std::vector<double>& coords, std::vector<int>& elements_ids, bool flag) {

	bool clamped = true;

	std::vector<unsigned int> edges_ids;

	sc_array_t *elements = &mesh->elements;
	int npoints = 0;

	FILE * fdbg;
	char filename[80];
	sprintf(filename, "Edge_deb_%04d.txt", mesh->mpi_rank);
	fdbg = fopen(filename, "w");

	sc_hash_array_t* hash_edge = sc_hash_array_new(sizeof (edge_id_t), edge_hash, edge_equal, &clamped);
        
        sc_hash_array_t* hash_id = sc_hash_array_new(sizeof (name_t), id_hash, id_equal, &clamped);

	for (int iel = 0; iel < elements->elem_count; ++iel) {

		octant_t *elem = (octant_t*) sc_array_index(&mesh->elements, iel);

		fprintf(fdbg,"Element: %d\n",iel);

		//initialization for the edge id
		for (int edge = 0; edge < 12; ++edge) {
			int Edge2GNode[2];

			int node1 = elem->nodes[EdgeVerticesMap[edge][0]].id;
			int node2 = elem->nodes[EdgeVerticesMap[edge][1]].id;

			Edge2GNode[0] = node1 <= node2 ? node1 : node2;
			Edge2GNode[1] = node1 >= node2 ? node1 : node2;

			elem->edge_id[edge] = edge_id(Edge2GNode, hash_edge, npoints);
			elem->edge_ref[edge] = false;
			fprintf(fdbg,"Edge: %d, global: %d %d, id: %d\n",edge, Edge2GNode[0],Edge2GNode[1], elem->edge_id[edge]);

		}

	}
	fclose(fdbg);

	//initialization for the node color
	for (int iel = 0; iel < elements->elem_count; ++iel) {
		octant_t *elem = (octant_t*) sc_array_index(&mesh->elements, iel);
		for(int j = 0; j < 8; j++) {
			octant_node_t* node = &elem->nodes[j];
			node->color=-1;
		}
	}


	for (int iel = 0; iel < elements_ids.size(); ++iel) {

		octant_t *elem = (octant_t*) sc_array_index(&mesh->elements, elements_ids[iel]);

		GtsSegment * segments[12]={0};
		GtsPoint * point[12]={NULL};
		int ed_cont = 0;

		for (int edge = 0; edge < 12; ++edge) {
			point[edge] = NULL;

			int node1 = elem->nodes[EdgeVerticesMap[edge][0]].id;
			int node2 = elem->nodes[EdgeVerticesMap[edge][1]].id;

			GtsVertex *v1 = gts_vertex_new(gts_vertex_class(), coords[node1 * 3], coords[node1 * 3 + 1], coords[node1 * 3 + 2]);
			GtsVertex *v2 = gts_vertex_new(gts_vertex_class(), coords[node2 * 3], coords[node2 * 3 + 1], coords[node2 * 3 + 2]);

			segments[edge] = gts_segment_new(gts_segment_class(), v1, v2);
			GtsBBox *sb = gts_bbox_segment(gts_bbox_class(), segments[edge]);
			GSList* list = gts_bb_tree_overlap(mesh->gdata.bbt, sb);
			if (list == NULL) continue;
			while (list) {
				GtsBBox *b = GTS_BBOX(list->data);
				point[edge] = SegmentTriangleIntersection(segments[edge], GTS_TRIANGLE(b->bounded));
				if (point[edge]) {
					elem->edge_ref[edge] = true;
					ed_cont++;
					break;
				}
				list = list->next;
			}
		}

		//Bounding box intercepted
		if(elem->pad == -1 && ed_cont == 0){
			elements_ids.erase(elements_ids.begin() + iel);
			elem->pad = 0;
			for (int edge = 0; edge < 12; ++edge) {
				elem->edge_ref[edge] = false;
			}
		}

		//clean points
		for (int edge = 0; edge < 12; edge++) {
			if (point[edge]) gts_object_destroy(GTS_OBJECT(point[edge]));
			//if (segments[edge]) gts_object_destroy(GTS_OBJECT(segments[edge]));
			point[edge] = NULL;
		}
	}

        //TODO change the for to while, check convergence criteria
	for (int i = 0; i < 100; i++){
		printf("numero %d\n",i);
		printf(" Elements ref: %d\n", elements_ids.size());
		IdentifyTemplate(mesh, coords, elements_ids);
		Edge_identification( mesh, elements_ids, hash_id);
                printf("Tamanho da hash: %d\n",hash_id->a.elem_count);
		Edge_propagation( mesh, elements_ids, hash_id);
	}
	printf(" Elements ref: %d\n", elements_ids.size());
	IdentifyTemplate(mesh, coords, elements_ids);

}
