#include <gts.h>
#include <glib.h>
#include <vector>
#include <iostream>
using namespace std;
#include <set>
#include <algorithm>
#include <assert.h>
#include <sc.h>
#include <sc_io.h>
#include <sc_containers.h>


#include "hexa.h"
#include "refinement.h"

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

void IdentifyTemplate(hexa_tree_t* mesh, std::vector<int>& elements_ids){

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
    if(false){
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
        MPI_Irecv(&recvbuf[offset], 2*m->idxs.elem_count, MPI_INT, m->rank,0,MPI_COMM_WORLD, &requests[c]);
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
            //TODO fix this bug
            bool out =  sc_hash_array_lookup(hash_edge_ref, &key, &position); 
            //bool out = 0;
            //for(int j=0; j< hash_edge_ref->a.elem_count; j++){
            //    octant_edge_t* edge = (octant_edge_t*) sc_array_index(&hash_edge_ref->a,j);        
            //    if(edge->id==*mm){
            //        out = 1;
            //    }
            //}
            sendbuf[(j+offset)] = (long long) *mm;
            sendbuf[2*(j+offset)+1] = (long long)   out; 
        }
        MPI_Isend(&sendbuf[offset], 2*m->idxs.elem_count, MPI_INT, m->rank,0,MPI_COMM_WORLD, &requests[c]);
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
                if(recvbuf[2*(j+offset)+1]){
                        edge_add(recvbuf[2*(j+offset)], hash_edge_ref );
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
            //TODO fix this bug
            bool out =  sc_hash_array_lookup(hash_edge_ref, &key, &position); 
            //bool out = 0;
            //for(int j=0; j< hash_edge_ref->a.elem_count; j++){
            //   octant_edge_t* edge = (octant_edge_t*) sc_array_index(&hash_edge_ref->a,j);        
            //    if(edge->id==*mm){
            //        out = 1;
            //    }
            //}
            sendbuf[(j+offset)] = (long long) *mm;
            sendbuf[2*(j+offset)+1] = (long long)   out; 
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
    for(int i = 0; i < mesh->comm_map_edge.SendTo.elem_count; ++i) {
        message_t *m = (message_t*) sc_array_index(&mesh->comm_map_edge.SendTo, i);
        for(int j = 0; j < m->idxs.elem_count; ++j){
                if(recvbuf[2*(j+offset)+1]){
                        edge_add(recvbuf[2*(j+offset)], hash_edge_ref );
            }
        }
        offset += 2*m->idxs.elem_count;
    }
    
    
    //
    for (int iel = 0; iel < elements_ids.size(); ++iel) {
        size_t position;
        octant_t *elem = (octant_t*) sc_array_index(&mesh->elements, elements_ids[iel]);
        for(int ied=0; ied<12; ied++){
            bool out = sc_hash_array_lookup(hash_edge_ref, &elem->edge[ied].id, &position);
            elem->edge[ied].ref = out;
        }
    }
    
#ifdef HEXA_DEBUG_
    if(1){
        fprintf(mesh->fdbg ,"after comm\n");

        fprintf(mesh->fdbg,"\n");
        for(int iel= 0; iel < mesh->elements.elem_count; iel++ ){
            octant_t *elem = (octant_t*) sc_array_index(&mesh->elements, iel);
            fprintf(mesh->fdbg,"El: %d template:%d pad:%d\n",iel,elem->tem,elem->pad);
            for (int edge = 0; edge < 12; ++edge) {
                fprintf(mesh->fdbg,"%lld ",elem->edge[edge].id);
            }
            fprintf(mesh->fdbg,"\n");
            for (int edge = 0; edge < 12; ++edge) {
                fprintf(mesh->fdbg,"%d ",elem->edge[edge].ref);
            }
            fprintf(mesh->fdbg,"\n");
        }
    }
#endif
    
    free(&recvbuf[0]);
    free(&sendbuf[0]);
    free(requests);
    free(statuses);
    
}

void Edge_propagation(hexa_tree_t* mesh, std::vector<int>& elements_ids, sc_hash_array_t* hash_edge_ref) {

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

void CheckOctreeTemplate(hexa_tree_t* mesh, const std::vector<double>& coords, std::vector<int>& elements_ids, bool flag) {

    bool clamped = true;
    int npoints = 0;

    //sc_array_t *elements = &mesh->elements;        
    sc_hash_array_t* hash_edge_ref = sc_hash_array_new(sizeof (octant_edge_t), id_hash, id_equal, &clamped);

    for (int iel = 0; iel < elements_ids.size(); ++iel) {

        octant_t *elem = (octant_t*) sc_array_index(&mesh->elements, elements_ids[iel]);

        GtsSegment * segments[12]={0};
        GtsPoint * point[12]={NULL};
        int ed_cont = 0;

        for (int edge = 0; edge < 12; ++edge) {
            elem->edge[edge].ref = false;
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
                        elem->edge[edge].ref = true;
                        edge_add(elem->edge[edge].id, hash_edge_ref );
                        elem->pad = -1;
                        ed_cont++;
                        break;
                }
                list = list->next;
            }
        }

/*
        //TODO check this! Bounding box intercepted
        if(elem->pad == -1 && ed_cont == 0){
                elements_ids.erase(elements_ids.begin() + iel);
                elem->pad = 0;
                for (int edge = 0; edge < 12; ++edge) {
                        elem->edge[edge].ref = false;
                }
        }
*/
        //clean points
        for (int edge = 0; edge < 12; edge++) {
                if (point[edge]) gts_object_destroy(GTS_OBJECT(point[edge]));
                //if (segments[edge]) gts_object_destroy(GTS_OBJECT(segments[edge]));
                point[edge] = NULL;
        }
    }
/*
#ifdef HEXA_DEBUG_
    fprintf(mesh->fdbg,"\nInitial state\n");
    for(int iel= 0; iel < mesh->elements.elem_count; iel++ ){
        octant_t *elem = (octant_t*) sc_array_index(&mesh->elements, iel);
        fprintf(mesh->fdbg,"El: %d template:%d pad:%d\n",iel,elem->tem,elem->pad);
        for (int edge = 0; edge < 12; ++edge) {
            fprintf(mesh->fdbg,"%lld ",elem->edge[edge].id);
        }
        fprintf(mesh->fdbg,"\n");
        for (int edge = 0; edge < 12; ++edge) {
            fprintf(mesh->fdbg,"%d ",elem->edge[edge].ref);
        }
        fprintf(mesh->fdbg,"\n");
    }
#endif
    */
    //work in serial
    for (int i = 0; i < 100; i++){
        printf(" Elements ref: %d\n", elements_ids.size());
        IdentifyTemplate(mesh, elements_ids);
        Edge_identification( mesh, elements_ids, hash_edge_ref);
        Edge_propagation (mesh, elements_ids, hash_edge_ref);
    }
    printf(" Elements ref: %d\n", elements_ids.size());
    IdentifyTemplate(mesh, elements_ids);
    
    
    /*
    for (int i = 0; i < 10; i++){
            //printf("numero %d, in proc:%d\n",i,mesh->mpi_rank);
            //printf(" Elements ref: %d\n", elements_ids.size());
            Edge_comunication(mesh, elements_ids, hash_edge_ref);
            Edge_propagation (mesh, elements_ids, hash_edge_ref);
    }
    printf(" Elements ref: %d\n", elements_ids.size());
    //IdentifyTemplate(mesh, elements_ids);
    //Edge_identification( mesh, elements_ids, hash_edge_ref);

#ifdef HEXA_DEBUG_
    fprintf(mesh->fdbg ,"Final state\n");
    for(int iel= 0; iel < mesh->elements.elem_count; iel++ ){
        octant_t *elem = (octant_t*) sc_array_index(&mesh->elements, iel);
        fprintf(mesh->fdbg,"El: %d template:%d pad:%d\n",iel,elem->tem,elem->pad);
        for (int edge = 0; edge < 12; ++edge) {
            fprintf(mesh->fdbg,"%d ",elem->edge[edge].id);
        }
        fprintf(mesh->fdbg,"\n");
        for (int edge = 0; edge < 12; ++edge) {
            fprintf(mesh->fdbg,"%d ",elem->edge[edge].ref);
        }
        fprintf(mesh->fdbg,"\n");
    }
#endif
   
   //sc_hash_array_rip(hash_edge_ref,&mesh->edges_ref);
*/
}
