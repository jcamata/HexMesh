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

unsigned pillow_hash_fn (const void *v, const void *u)
{
	const pillow_t *q = (const pillow_t*) v;
	uint32_t            a, b, c;

	a = (uint32_t) q->x;
	b = (uint32_t) q->y;
	c = (uint32_t) q->z;
	sc_hash_mix (a, b, c);
	sc_hash_final (a, b, c);
	return (unsigned) c;
}

int pillow_equal_fn (const void *v1, const void *v2, const void *u)
{
	const pillow_t *q1 = (const pillow_t*) v1;
	const pillow_t *q2 = (const pillow_t*) v2;
	return (q1->x == q2->x && q1->y == q2->y && q1->z == q2->z);
}
/*
typedef struct aux_teste
{
	octant_node_t orig;
	octant_node_t pill;
}aux_teste_t;

unsigned aux_teste_fn (const void *v, const void *u)
{
	const aux_teste_t *q = (const aux_teste_t*) v;
	uint32_t            a, b, c;

	a = (uint32_t) q->orig.x;
	b = (uint32_t) q->orig.y;
	c = (uint32_t) q->orig.z;
	sc_hash_mix (a, b, c);
	sc_hash_final (a, b, c);
	return (unsigned) c;
}

int aux_teste_equal_fn (const void *v1, const void *v2, const void *u)
{
	const aux_teste_t *q1 = (const aux_teste_t*) v1;
	const aux_teste_t *q2 = (const aux_teste_t*) v2;
	return (q1->orig.x == q2->orig.x && q1->orig.y == q2->orig.y && q1->orig.z == q2->orig.z);
}
 */
int AddPoint(hexa_tree_t* mesh, sc_hash_array_t* hash_nodes, GtsPoint *p, std::vector<double> &coords, int x, int y, int z)
{
	size_t position;
	octant_node_t *r;
	octant_node_t key;

	key.x = x;
	key.y = y;
	key.z = z;

	r = (octant_node_t*) sc_hash_array_insert_unique(hash_nodes, &key, &position);

	if (r != NULL)
	{
		r->x = x;
		r->y = y;
		r->z = z;
		r->id = mesh->nodes.elem_count;
		octant_node_t* n = (octant_node_t*) sc_array_push(&mesh->nodes);
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
		r = (octant_node_t*) sc_array_index(&hash_nodes->a, position);
		return r->id;
	}
}

GtsPoint* LinearMapHex(const double* cord_in_ref, const double* cord_in_x, const double* cord_in_y, const double* cord_in_z)
{

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

void RedoNodeMapping(hexa_tree_t* mesh)
{

	int factor = 12;
	//just the multiplication
	// it allow us add int points in the mesh
	// keeping a structured mesh
	for(int iel = 0; iel <mesh->elements.elem_count; iel++)
	{
		octant_t * elem = (octant_t*) sc_array_index (&mesh->elements, iel);
		for(int ino = 0; ino < 8; ino++)
		{
			elem->nodes[ino].x = factor*elem->nodes[ino].x;
			elem->nodes[ino].y = factor*elem->nodes[ino].y;
			elem->nodes[ino].z = factor*elem->nodes[ino].z;
		}
		elem->x = 4*elem->x;
		elem->y = 4*elem->y;
		elem->z = 4*elem->z;
	}

	for(int ino = 0; ino <mesh->nodes.elem_count; ino++)
	{
		octant_node_t * node = (octant_node_t*) sc_array_index (&mesh->nodes, ino);
		node->x = factor*node->x;
		node->y = factor*node->y;
		node->z = factor*node->z;
	}
	mesh->x_start = factor*mesh->x_start;
	mesh->y_start = factor*mesh->y_start;
	mesh->x_end = factor*mesh->x_end;
	mesh->y_end = factor*mesh->y_end;

	mesh->ncellx = 4*mesh->ncellx;
	mesh->ncelly = 4*mesh->ncelly;
	mesh->max_z = 4*mesh->max_z;
}

void CopyPropEl(hexa_tree_t* mesh, int id, octant_t *elem1)
{

	octant_t *elem = (octant_t*) sc_array_index(&mesh->elements, id);

	elem1->level = -1;
	elem1->tem = elem->tem;
	elem1->pad = elem->pad;
	elem1->n_mat = elem->n_mat;
	elem1->pml_id = elem->pml_id;
	elem1->father = id;
	elem1->boundary = elem->boundary;

	for(int ino = 0; ino < 8; ino++)
	{
		elem1->nodes[ino].color = elem->nodes[ino].color;
		elem1->nodes[ino].fixed = elem->nodes[ino].fixed;
		elem1->nodes[ino].x = elem->nodes[ino].x;
		elem1->nodes[ino].y = elem->nodes[ino].y;
		elem1->nodes[ino].z = elem->nodes[ino].z;
	}

	for(int iedge = 0; iedge < 12; iedge++)
	{
		elem1->edge[iedge].coord[0] = elem->edge[iedge].coord[0];
		elem1->edge[iedge].coord[1] = elem->edge[iedge].coord[1];
		elem1->edge[iedge].id = elem->edge[iedge].id;
		elem1->edge[iedge].ref = elem->edge[iedge].ref;
	}

	for(int isurf = 0; isurf < 6; isurf++)
	{
		elem1->surf[isurf].ext = elem->surf[isurf].ext;
	}

	elem1->x=elem->x;
	elem1->y=elem->y;
	elem1->z=elem->z;
}

void Pillowing(hexa_tree_t* mesh, std::vector<double>& coords, std::vector<int>& nodes_b_mat)
{

	bool deb = false;
	bool clamped = true;

	//noeud hash
	sc_hash_array_t*   hash_nodes  = (sc_hash_array_t *)sc_hash_array_new(sizeof(octant_node_t), node_hash_fn , node_equal_fn, &clamped);
	for(int iel = 0; iel < mesh->elements.elem_count; iel++)
	{
		octant_t* elem = (octant_t*) sc_array_index(&mesh->elements, iel);
		for (int ino = 0; ino < 8; ino++)
		{
			size_t position;
			octant_node_t *r;
			octant_node_t key;
			key.x = elem->nodes[ino].x;
			key.y = elem->nodes[ino].y;
			key.z = elem->nodes[ino].z;
			r = (octant_node_t*) sc_hash_array_insert_unique(hash_nodes, &key, &position);
			if(r!=NULL){
				r->x = elem->nodes[ino].x;
				r->y = elem->nodes[ino].y;
				r->z = elem->nodes[ino].z;
				r->id = elem->nodes[ino].id;
			}
		}
	}

	//faire noeud hash nodes_b_mat
	sc_hash_array_t*   hash_b_mat  = (sc_hash_array_t *)sc_hash_array_new(sizeof(octant_node_t), node_hash_fn, node_equal_fn, &clamped);
	for(int ino = 0; ino < nodes_b_mat.size(); ino++)
	{
		size_t position;
		octant_node_t *r;
		octant_node_t key;
		octant_node_t* node = (octant_node_t*) sc_array_index (&mesh->nodes, nodes_b_mat[ino]);
		key.x = node->x;
		key.y = node->y;
		key.z = node->z;

		r = (octant_node_t*) sc_hash_array_insert_unique(hash_b_mat, &key, &position);
		if(r!=NULL){
			r->x = key.x;
			r->y = key.y;
			r->z = key.z;
			r->id = node->id;
		}
	}

	//mesh->oct.elem_count
	for(int ioc = 0; ioc < mesh->oct.elem_count; ioc++){
		octree_t* oct = (octree_t*) sc_array_index(&mesh->oct, ioc);

		//générer le tableau pour le pillow
		sc_hash_array_t* pillow = (sc_hash_array_t*) sc_hash_array_new(sizeof(pillow_t),pillow_hash_fn,pillow_equal_fn,&clamped);
		for(int iel = 0; iel < 8; iel++){
			octant_t* elem = (octant_t*) sc_array_index(&mesh->elements,oct->id[iel]);

			for(int ino = 0; ino < 8; ino++){
				//check si le noeud est sur la hahs
				//hash_b_mat
				size_t position;
				octant_node_t key;
				key.x = elem->nodes[ino].x;
				key.y = elem->nodes[ino].y;
				key.z = elem->nodes[ino].z;
				key.id = elem->nodes[ino].id;

				bool lnode = sc_hash_array_lookup(hash_b_mat, &key, &position);

				if(lnode)
				{
					if(deb) printf("J'ai trouvée le noeud %d (nombre:%d) élémént %d (nombre:%d) du octree %d\n",ino,key.id,iel,elem->id,ioc);

					//on commance a ajouter les noeuds dans la hahs de pillow
					pillow_t keyP;
					size_t positionP;
					keyP.x = elem->nodes[ino].x;
					keyP.y = elem->nodes[ino].y;
					keyP.z = elem->nodes[ino].z;
					pillow_t* p = (pillow_t*) sc_hash_array_insert_unique(pillow, &keyP, &positionP);
					if(p!=NULL)
					{
						p->id = elem->nodes[ino].id;
						p->a = -1;
						p->b = -1;
						p->x = keyP.x;
						p->y = keyP.y;
						p->z = keyP.z;
						p->elem[0] = elem->id;
						p->list_elem = 1;
						p->pa = false;
						p->pb = false;
						// on doit gerer les surfaces que sont fixes...
						p->list_face[0] = 0;
						for(int isurf = 0; isurf < 3; isurf++)
						{
							bool surf = true;
							if(deb) printf("Essaye de la surface %d dans le élément %d nombre du noeud %d",VertexSurfMap[ino][isurf], iel,elem->nodes[ino].id);
							for(int ive = 0; ive < 4; ive++)
							{
								size_t positionSurfVert;
								octant_node_t SurfVert;
								int refnode = FaceNodesMap[VertexSurfMap[ino][isurf]][ive];
								SurfVert.x = elem->nodes[refnode].x;
								SurfVert.y = elem->nodes[refnode].y;
								SurfVert.z = elem->nodes[refnode].z;
								bool lvert = sc_hash_array_lookup(hash_b_mat, &SurfVert, &positionSurfVert);
								if(!lvert) surf = false;
							}
							if(deb) printf("surface nombre:%d\n",surf);
							if(surf) p->face[0][p->list_face[0]] = VertexSurfMap[ino][isurf];
							if(surf) p->list_face[0]++;
						}
					}
					else
					{
						pillow_t* p1 = (pillow_t*) sc_array_index(&pillow->a,positionP);
						p1->elem[p1->list_elem] = elem->id;
						p1->list_face[p1->list_elem] = 0;

						for(int isurf = 0; isurf < 3; isurf++)
						{
							bool surf = true;
							if(deb) printf("Essaye de la surface %d dans le élément %d nombre du noeud %d",VertexSurfMap[ino][isurf], iel,elem->nodes[ino].id);
							for(int ive = 0; ive < 4; ive++)
							{
								size_t positionSurfVert;
								octant_node_t SurfVert;
								int refnode = FaceNodesMap[VertexSurfMap[ino][isurf]][ive];
								SurfVert.x = elem->nodes[refnode].x;
								SurfVert.y = elem->nodes[refnode].y;
								SurfVert.z = elem->nodes[refnode].z;
								bool lvert = sc_hash_array_lookup(hash_b_mat, &SurfVert, &positionSurfVert);
								if(!lvert) surf = false;
							}
							if(deb) printf("surface nombre:%d\n",surf);
							if(surf) p1->face[p1->list_elem][p1->list_face[p1->list_elem]] = VertexSurfMap[ino][isurf];
							if(surf) p1->list_face[p1->list_elem]++;
						}
						p1->list_elem++;
					}
				}
			}
		}

		//rétrécir octree/pillow
		sc_array_t toto;
		sc_array_init(&toto, sizeof(octant_t));
		for(int iel = 0; iel < 8; iel++)
		{
			octant_t* elemOrig = (octant_t*) sc_array_index(&mesh->elements, oct->id[iel]);
			octant_t * elem = (octant_t*) sc_array_push(&toto);

			hexa_element_copy(elemOrig,elem);

			double ref_in_x[8],ref_in_y[8],ref_in_z[8];
			//on garde la connectivite original...
			int connec[8];
			for (int ino = 0; ino < 8; ino++)
			{
				connec[ino] = elem->nodes[ino].id;
				ref_in_x[ino] = coords[3*elem->nodes[ino].id+0];
				ref_in_y[ino] = coords[3*elem->nodes[ino].id+1];
				ref_in_z[ino] = coords[3*elem->nodes[ino].id+2];
			}

			for(int isurf = 0; isurf < 6; isurf++){
				bool surf = true;
				int new_nodes[4] = {-1};

				//check if the surface should be pillowed
				for(int ive = 0; ive < 4; ive++)
				{
					size_t positionSurfVert;
					octant_node_t SurfVert;
					SurfVert.x = elem->nodes[FaceNodesMap[isurf][ive]].x;
					SurfVert.y = elem->nodes[FaceNodesMap[isurf][ive]].y;
					SurfVert.z = elem->nodes[FaceNodesMap[isurf][ive]].z;
					//printf("%d %d %d\n",SurfVert.x,SurfVert.y,SurfVert.z);
					bool lvert = sc_hash_array_lookup(hash_b_mat, &SurfVert, &positionSurfVert);
					if(!lvert) surf = false;
				}

				//pillow the surface
				if(surf){
					//printf("El:%d iel:%d mat:%d face %d\n",elem->id,iel,elem->n_mat,isurf);
					//criar os nos
					for(int ive = 0; ive < 4; ive++)
					{
						//inicializa
						int x = elem->nodes[FaceNodesMap[isurf][ive]].x;
						int y = elem->nodes[FaceNodesMap[isurf][ive]].y;
						int z = elem->nodes[FaceNodesMap[isurf][ive]].z;
						//printf("Node original: %d %d %d\n",x,y,z);
						//printf("%d %d %d\n",x,y,z);
						//achar o pillow dele
						size_t position;
						pillow_t vert;
						vert.x = x;
						vert.y = y;
						vert.z = z;
						bool lvert = sc_hash_array_lookup(pillow, &vert, &position);
						pillow_t* p = (pillow_t*) sc_array_index(&pillow->a,position);

						bool flag[6];
						for(int i = 0; i<6;i++){
							flag[i] = true;
						}
						for(int ielp = 0; ielp < p->list_elem; ielp++){
							octant_t* el = (octant_t*) sc_array_index(&mesh->elements,p->elem[ielp]);
							if(el->n_mat == elem->n_mat && p->list_face[ielp] != 0){
								//printf("eu o elemento %d sou irmãozinho do elemento %d pq tenho o memso material\n",el->id,elem->id);
								//printf("Tenho %d faces para alterar\n",p->list_face[ielp]);
								for(int iisu = 0; iisu < p->list_face[ielp]; iisu++){
									//printf("face: %d\n",p->face[ielp][iisu]);
									//printf("%d %d %d %d %d %d\n",flag[0],flag[1],flag[2],flag[3],flag[4],flag[5]);
									if(p->face[ielp][iisu]==0 && flag[0]){
										//printf("entrei no 0\n");
										x += 6;
										flag[0] = false;
									}else if(p->face[ielp][iisu]==1 && flag[1]){
										//printf("entrei no 1\n");
										x -= 6;
										flag[1] = false;
									}else if(p->face[ielp][iisu]==2 && flag[2]){
										//printf("entrei no 2\n");
										y += 6;
										flag[2] = false;
									}else if(p->face[ielp][iisu]==3 && flag[3]){
										//printf("entrei no 3\n");
										y -= 6;
										flag[3] = false;
									}else if(p->face[ielp][iisu]==4 && flag[4]){
										//printf("entrei no 4\n");
										z += 6;
										flag[4] = false;
									}else if(p->face[ielp][iisu]==5 && flag[5]){
										//printf("entrei no 5\n");
										z -= 6;
										flag[5] = false;
									}
									//printf("%d %d %d %d %d %d\n",flag[0],flag[1],flag[2],flag[3],flag[4],flag[5]);
								}
							}
						}
						//printf("I found the node: %d for this node I have %d faces\n",p->id,p->list_face[iel]);

						//printf("Node Alterado: %d %d %d\n",x,y,z);
						//printf("%d %d %d\n",x,y,z);

						//ajouter sur la hash de noeuds
						octant_node_t key;
						key.x = x;
						key.y = y;
						key.z = z;
						octant_node_t* ra = (octant_node_t*) sc_hash_array_insert_unique(hash_nodes, &key, &position);
						if(ra!=NULL){
							ra->x = x;
							ra->y = y;
							ra->z = z;
							ra->id = hash_nodes->a.elem_count-1;
							new_nodes[ive] = ra->id;
							ra->fixed = 0;
							ra->color = 0;
							if(elem->n_mat == 0) p->a = ra->id;
							if(elem->n_mat == 1) p->b = ra->id;
							if(elem->n_mat == 0) p->pa = true;
							if(elem->n_mat == 1) p->pb = true;
							if(ra->id == -1)
							{
								printf("ra->id == -1 id:%d",ra->id);
							}
							double cord_in_ref[3];
							cord_in_ref[0] = (x - elem->nodes[0].x)/6 -1;
							cord_in_ref[1] = (y - elem->nodes[0].y)/6 -1;
							cord_in_ref[2] = (z - elem->nodes[0].z)/6 -1;

							GtsPoint* ref_point = LinearMapHex(cord_in_ref, ref_in_x, ref_in_y, ref_in_z);
							coords.push_back(ref_point->x);
							coords.push_back(ref_point->y);
							coords.push_back(ref_point->z);
						}else{
							octant_node_t* ra = (octant_node_t*) sc_array_index(&hash_nodes->a,position);
							new_nodes[ive] = ra->id;
							if(elem->n_mat == 0) p->a = ra->id;
							if(elem->n_mat == 1) p->b = ra->id;
							if(ra->id == -1)
							{
								printf("no else ra->id == -1 id:%d",ra->id);
							}
						}
					}

					//faire les éléments
					if(true)
					{
						//printf("entrei para fazer elementos\n");
						//créer l'élément
						octant_t* pelem = (octant_t*) sc_array_push(&mesh->elements);
						pelem->id = mesh->elements.elem_count-1;
						pelem->n_mat = elem->n_mat;
						pelem->pad = elem->pad;
						pelem->x = elem->x;
						pelem->y = elem->y;
						pelem->z = elem->z;

						//printf("  Face:%d\n",aux[isurf]);
						for(int ino = 0; ino < 4; ino++)
						{
							pelem->nodes[FaceNodesMap[isurf][ino]].id = connec[FaceNodesMap[isurf][ino]];
							pelem->nodes[FaceNodesMap[isurf][ino]].fixed = elem->nodes[FaceNodesMap[isurf][ino]].fixed;
							pelem->nodes[FaceNodesMap_inv[isurf][ino]].id= new_nodes[ino];

							pelem->nodes[FaceNodesMap[isurf][ino]].fixed= 0;
							octant_node_t* node = (octant_node_t*) sc_array_index(&hash_nodes->a,connec[FaceNodesMap[isurf][ino]]);
							pelem->nodes[FaceNodesMap[isurf][ino]].x = node->x;
							pelem->nodes[FaceNodesMap[isurf][ino]].y = node->y;
							pelem->nodes[FaceNodesMap[isurf][ino]].z = node->z;
							node = (octant_node_t*) sc_array_index(&hash_nodes->a,new_nodes[ino]);
							pelem->nodes[FaceNodesMap_inv[isurf][ino]].x = node->x;
							pelem->nodes[FaceNodesMap_inv[isurf][ino]].y = node->y;
							pelem->nodes[FaceNodesMap_inv[isurf][ino]].z = node->z;

							elemOrig->nodes[FaceNodesMap[isurf][ino]].id = new_nodes[ino];

							if(new_nodes[ino]==-1){
								printf("New node == -1");
							}
						}
					}
				}
			}
		}

		//element update
		//TODO bug here when I have an element tha there are no surface only a vertex or an edge
		for(int ive = 0; ive < pillow->a.elem_count; ive++){
			pillow_t* p = (pillow_t*) sc_array_index(&pillow->a,ive);
			for(int iel = 0; iel < p->list_elem; iel++){
				octant_t* el = (octant_t*) sc_array_index(&mesh->elements,p->elem[iel]);
				if(p->list_face[iel] == 0){
					//printf("para o pillow %d, ive:%d, sou o elemento %d\n",p->id,ive,el->id);
					for(int ino = 0; ino < 8; ino++){
						if(el->nodes[ino].id == p->id){
							int ori = el->nodes[ino].id;
							//printf("Sou o ino:%d e sou o no:%d\n",ino,p->id);
							//printf("o id do no pode ser %d ou %d\n",p->a,p->b);
							if(el->n_mat == 0){
								octant_node_t* node = (octant_node_t*) sc_array_index(&hash_nodes->a,p->a);
								el->nodes[ino].id = p->a;
								el->nodes[ino].x = node->x;
								el->nodes[ino].y = node->y;
								el->nodes[ino].z = node->z;
							}
							if(el->n_mat == 1){
								octant_node_t* node = (octant_node_t*) sc_array_index(&hash_nodes->a,p->b);
								//printf("nesse caso sou o %d\n",node->id);
								el->nodes[ino].id = p->b;
								el->nodes[ino].x = node->x;
								el->nodes[ino].y = node->y;
								el->nodes[ino].z = node->z;
							}
							if(el->nodes[ino].id == -1){
								printf("Me acharam antes %d, p->a: %d p->b: %d\n",el->id,p->a,p->b);
								printf("ino: %d pid:%d elnodeid: %d or %d\n",ino,p->id,el->nodes[ino].id,ori);
								el->nodes[ino].id = ori;

								printf("Sou o octree numero %d\n",ioc);
								for(int ied = 0; ied < 12; ied++){
									printf("%d ",oct->edge[ied]);
								}
								printf("\n");
							}
						}
					}
				}
			}
		}


		//TODO check here the bug...
		//sc_array_destroy(&el_copy);
		for(int iel = 0; iel < 0; iel++)
		{
			octant_t* elem = (octant_t*) sc_array_index(&mesh->elements,oct->id[iel]);
			for(int ino = 0; ino < 8; ino ++)
			{
				octant_node_t* node = (octant_node_t*) sc_array_index(&mesh->nodes,elem->nodes[ino].id);
				elem->nodes[ino].x = node->x;
				elem->nodes[ino].y = node->y;
				elem->nodes[ino].z = node->z;
			}
		}
	}
	//faire la atualisation de mon tableau de noeuds
	sc_array_reset(&mesh->nodes);
	sc_hash_array_rip(hash_nodes,&mesh->nodes);
}

void SurfaceIdentification(hexa_tree_t* mesh, std::vector<double>& coords)
{

	//assign a color for the nodes...

	// free node
	// color= 0 && fixed = 0;

	//exterior surface fixed nodes
	//fixed = 1 && color = 1

	//exterior global nodes
	//color && fixed = -1;
	// x- = 2 x+ = 3
	// y- = 4 y+ = 5
	// z- = 6 z+ = 7

	//exterior local nodes
	//color && fixed = -1;
	// x- = -2 x+ = -3
	// y- = -4 y+ = -5
	// z- = -6 z+ = -7

	//exterior global edges
	//color && fixed = -1;
	// z-y- = 0+11
	// z-x+ = 1+11
	// z-y+ = 2+11
	// z-x- = 3+11

	// x-y- = 4+11
	// x+y- = 5+11
	// x+y+ = 6+11
	// x-y+ = 7+11

	// z+y- = 8+11
	// z+x+ = 9+11
	// z+y+ = 10+11
	// z+x- = 11+11

	//exterior global vertex
	//color && fixed = -1;
	//x-y-z- = 30
	//x+y-z- = 31
	//x+y+z- = 32
	//x-y+z- = 33

	//x-y-z+ = 34
	//x+y-z+ = 35
	//x+y+z+ = 36
	//x-y+z+ = 37

	//exterior local edges
	//color && fixed = -1;
	// z-y- = -0-11
	// z-x+ = -1-11
	// z-y+ = -2-11
	// z-x- = -3-11

	// x-y- = -4-11
	// x+y- = -5-11
	// x+y+ = -6-11
	// x-y+ = -7-11

	// z+y- = -8-11
	// z+x+ = -9-11
	// z+y+ = -10-11
	// z+x- = -11-11

	//exterior local vertex
	//color && fixed = -1;
	//x-y-z- = -30
	//x+y-z- = -31
	//x+y+z- = -32
	//x-y+z- = -33

	//x-y-z+ = -34
	//x+y-z+ = -35
	//x+y+z+ = -36
	//x-y+z+ = -37

	bool deb = false;
	bool clamped = true;
	//vertex hash
	sc_hash_array_t*	  vertex_hash  = (sc_hash_array_t *)sc_hash_array_new(sizeof (octant_vertex_t), vertex_hash_id, vertex_equal_id, &clamped);


	//fazendo vertex hash & assign the color to the local nodes
	for(int iel = 0; iel < mesh->elements.elem_count ; iel++)
	{
		size_t  position;
		octant_t *elem = (octant_t*) sc_array_index(&mesh->elements, iel);
		for (int ino = 0; ino < 8; ino++)
		{
			//build the hash
			octant_vertex_t key;
			key.id = elem->nodes[ino].id;
			octant_vertex_t* vert = (octant_vertex_t*) sc_hash_array_insert_unique(vertex_hash, &key, &position);
			if(vert != NULL)
			{
				vert->id = elem->nodes[ino].id;
				vert->list_elem = 1;
				vert->elem[vert->list_elem-1] = elem->id;
			}
			else
			{
				vert = (octant_vertex_t*) sc_array_index(&vertex_hash->a, position);
				vert->elem[vert->list_elem] = elem->id;
				vert->list_elem++;
			}


			//setting free the free nodes
			if(elem->nodes[ino].fixed == 0) elem->nodes[ino].color = 0;

			//assign the color for the local nodes...

			//surface and edges
			if(elem->nodes[ino].x == mesh->x_start)
			{
				elem->nodes[ino].fixed = -1;
				elem->nodes[ino].color = -2;
			}

			if(elem->nodes[ino].x == mesh->x_end)
			{
				elem->nodes[ino].fixed = -1;
				elem->nodes[ino].color = -3;
			}

			if(elem->nodes[ino].y == mesh->y_start)
			{
				elem->nodes[ino].fixed = -1;
				elem->nodes[ino].color = -4;
			}

			if(elem->nodes[ino].y == mesh->y_end)
			{
				elem->nodes[ino].fixed = -1;
				elem->nodes[ino].color = -5;
			}

			if(elem->nodes[ino].z == 0)
			{
				elem->nodes[ino].fixed = -1;
				elem->nodes[ino].color = -6;
			}

			if(elem->nodes[ino].z == 3*mesh->max_z)
			{
				elem->nodes[ino].fixed = -1;
				elem->nodes[ino].color = -7;
			}

			// z-y- = -0-11
			if(elem->nodes[ino].z == 0 && elem->nodes[ino].y == mesh->y_start)
			{
				elem->nodes[ino].fixed = -1;
				elem->nodes[ino].color = -0-11;
			}
			// z-x+ = -1-11
			if(elem->nodes[ino].z == 0 && elem->nodes[ino].x == mesh->x_end)
			{
				elem->nodes[ino].fixed = -1;
				elem->nodes[ino].color = -1-11;
			}
			// z-y+ = -2-11
			if(elem->nodes[ino].z == 0 && elem->nodes[ino].y == mesh->y_end)
			{
				elem->nodes[ino].fixed = -1;
				elem->nodes[ino].color = -2-11;
			}
			// z-x- = -3-11
			if(elem->nodes[ino].z == 0 && elem->nodes[ino].x == mesh->x_start)
			{
				elem->nodes[ino].fixed = -1;
				elem->nodes[ino].color = -3-11;
			}


			// x-y- = -4-11
			if(elem->nodes[ino].x == mesh->x_start && elem->nodes[ino].y == mesh->y_start)
			{
				elem->nodes[ino].fixed = -1;
				elem->nodes[ino].color = -4-11;
			}
			// x+y- = -5-11
			if(elem->nodes[ino].x == mesh->x_end && elem->nodes[ino].y == mesh->y_start)
			{
				elem->nodes[ino].fixed = -1;
				elem->nodes[ino].color = -5-11;
			}
			// x+y+ = -6-11
			if(elem->nodes[ino].x == mesh->x_end && elem->nodes[ino].y == mesh->y_end)
			{
				elem->nodes[ino].fixed = -1;
				elem->nodes[ino].color = -6-11;
			}
			// x-y+ = -7-11
			if(elem->nodes[ino].x == mesh->x_start && elem->nodes[ino].y == mesh->y_end)
			{
				elem->nodes[ino].fixed = -1;
				elem->nodes[ino].color = -7-11;
			}

			// z+y- = -8-11
			if(elem->nodes[ino].z == 3*mesh->max_z && elem->nodes[ino].y == mesh->y_start)
			{
				elem->nodes[ino].fixed = -1;
				elem->nodes[ino].color = -8-11;
			}
			// z+x+ = -9-11
			if(elem->nodes[ino].z == 3*mesh->max_z && elem->nodes[ino].x == mesh->x_end)
			{
				elem->nodes[ino].fixed = -1;
				elem->nodes[ino].color = -9-11;
			}
			// z+y+ = -10-11
			if(elem->nodes[ino].z == 3*mesh->max_z && elem->nodes[ino].y == mesh->y_end)
			{
				elem->nodes[ino].fixed = -1;
				elem->nodes[ino].color = -10-11;
			}
			// z+x- = -11-11
			if(elem->nodes[ino].z == 3*mesh->max_z && elem->nodes[ino].x == mesh->x_start)
			{
				elem->nodes[ino].fixed = -1;
				elem->nodes[ino].color = -11-11;
			}

			// x-y-z- = -30
			if(elem->nodes[ino].x == mesh->x_start && elem->nodes[ino].y == mesh->y_start && elem->nodes[ino].z == 0 )
			{
				elem->nodes[ino].fixed = -1;
				elem->nodes[ino].color = -30;
			}
			// x+y-z- = -31
			if(elem->nodes[ino].x == mesh->x_end && elem->nodes[ino].y == mesh->y_start && elem->nodes[ino].z == 0 )
			{
				elem->nodes[ino].fixed = -1;
				elem->nodes[ino].color = -31;
			}
			// x+y+z- = -32
			if(elem->nodes[ino].x == mesh->x_end && elem->nodes[ino].y == mesh->y_end && elem->nodes[ino].z == 0 )
			{
				elem->nodes[ino].fixed = -1;
				elem->nodes[ino].color = -32;
			}
			// x-y+z- = -33
			if(elem->nodes[ino].x == mesh->x_start && elem->nodes[ino].y == mesh->y_end && elem->nodes[ino].z == 0 )
			{
				elem->nodes[ino].fixed = -1;
				elem->nodes[ino].color = -33;
			}
			// x-y-z+ = -34
			if(elem->nodes[ino].x == mesh->x_start && elem->nodes[ino].y == mesh->y_start && elem->nodes[ino].z == 3*mesh->max_z )
			{
				elem->nodes[ino].fixed = -1;
				elem->nodes[ino].color = -34;
			}
			// x+y-z+ = -35
			if(elem->nodes[ino].x == mesh->x_end && elem->nodes[ino].y == mesh->y_start && elem->nodes[ino].z == 3*mesh->max_z )
			{
				elem->nodes[ino].fixed = -1;
				elem->nodes[ino].color = -35;
			}
			// x+y+z+ = -36
			if(elem->nodes[ino].x == mesh->x_end && elem->nodes[ino].y == mesh->y_end && elem->nodes[ino].z == 3*mesh->max_z )
			{
				elem->nodes[ino].fixed = -1;
				elem->nodes[ino].color = -36;
			}
			// x-y+z+ = -37
			if(elem->nodes[ino].x == mesh->x_start && elem->nodes[ino].y == mesh->y_end && elem->nodes[ino].z == 3*mesh->max_z )
			{
				elem->nodes[ino].fixed = -1;
				elem->nodes[ino].color = -37;
			}
		}
	}

	sc_array_init(&mesh->outsurf, sizeof(octant_t));

	sc_array_t toto;
	sc_array_init(&toto, sizeof(octant_t));

	//id global exterior surface
	for(int iel = 0; iel < mesh->elements.elem_count; iel++)
	{
		octant_t* elem = (octant_t*) sc_array_index(&mesh->elements, iel);
		//octant_t * elem = (octant_t*) sc_array_push(&toto);
		//hexa_element_copy(elemOrig,elem);

		int isurf;
		isurf = 0;
		elem->surf[isurf].ext = false;
		if(elem->nodes[FaceNodesMap[isurf][0]].x == 0 && elem->nodes[FaceNodesMap[isurf][1]].x == 0 &&
				elem->nodes[FaceNodesMap[isurf][2]].x == 0 && elem->nodes[FaceNodesMap[isurf][3]].x == 0)
		{
			elem->surf[isurf].ext = true;
			if(deb) for(int ino = 0; ino < 4; ino++) mesh->part_nodes[elem->nodes[FaceNodesMap[isurf][ino]].id] = isurf+20;
		}

		isurf = 1;
		elem->surf[isurf].ext = false;
		if(elem->nodes[FaceNodesMap[isurf][0]].x == 3*mesh->ncellx && elem->nodes[FaceNodesMap[isurf][1]].x == 3*mesh->ncellx &&
				elem->nodes[FaceNodesMap[isurf][2]].x == 3*mesh->ncellx && elem->nodes[FaceNodesMap[isurf][3]].x == 3*mesh->ncellx)
		{
			elem->surf[isurf].ext = true;
			if(deb) for(int ino = 0; ino < 4; ino++) mesh->part_nodes[elem->nodes[FaceNodesMap[isurf][ino]].id] = isurf+20;
		}

		isurf = 2;
		elem->surf[isurf].ext = false;
		if(elem->nodes[FaceNodesMap[isurf][0]].y == 0 && elem->nodes[FaceNodesMap[isurf][1]].y == 0 &&
				elem->nodes[FaceNodesMap[isurf][2]].y == 0 && elem->nodes[FaceNodesMap[isurf][3]].y == 0)
		{
			elem->surf[isurf].ext = true;
			if(deb) for(int ino = 0; ino < 4; ino++) mesh->part_nodes[elem->nodes[FaceNodesMap[isurf][ino]].id] = isurf+20;
		}

		isurf = 3;
		elem->surf[isurf].ext = false;
		if(elem->nodes[FaceNodesMap[isurf][0]].y == 3*mesh->ncelly && elem->nodes[FaceNodesMap[isurf][1]].y == 3*mesh->ncelly &&
				elem->nodes[FaceNodesMap[isurf][2]].y == 3*mesh->ncelly && elem->nodes[FaceNodesMap[isurf][3]].y == 3*mesh->ncelly)
		{
			elem->surf[isurf].ext = true;
			if(deb) for(int ino = 0; ino < 4; ino++) mesh->part_nodes[elem->nodes[FaceNodesMap[isurf][ino]].id] = isurf+20;
		}

		isurf = 5;
		elem->surf[isurf].ext = false;
		if(elem->nodes[FaceNodesMap[isurf][0]].z == 3*mesh->max_z && elem->nodes[FaceNodesMap[isurf][1]].z == 3*mesh->max_z &&
				elem->nodes[FaceNodesMap[isurf][2]].z == 3*mesh->max_z && elem->nodes[FaceNodesMap[isurf][3]].z == 3*mesh->max_z)
		{
			elem->surf[isurf].ext = true;
			if(deb) for(int ino = 0; ino < 4; ino++) mesh->part_nodes[elem->nodes[FaceNodesMap[isurf][ino]].id] = isurf+20;
		}

		isurf = 4;
		elem->surf[isurf].ext = false;
		if(elem->nodes[FaceNodesMap[isurf][0]].z == 0 && elem->nodes[FaceNodesMap[isurf][1]].z == 0 &&
				elem->nodes[FaceNodesMap[isurf][2]].z == 0 && elem->nodes[FaceNodesMap[isurf][3]].z == 0)
		{
			elem->surf[isurf].ext = true;
			if(deb) for(int ino = 0; ino < 4; ino++) mesh->part_nodes[elem->nodes[FaceNodesMap[isurf][ino]].id] = isurf+20;
		}

		bool aux = false;
		for(int isurf = 0; isurf < 6; isurf++) if(elem->surf[isurf].ext) aux = true;
		if(aux)
		{
			octant_t* elem1 = (octant_t*) sc_array_push(&mesh->outsurf);

			elem1->level = -1;
			elem1->id = elem->id;
			elem1->tem = elem->tem;
			elem1->pad = elem->pad;
			elem1->n_mat = elem->n_mat;
			elem1->pml_id = elem->pml_id;
			elem1->father = elem->id;
			elem1->boundary = elem->boundary;
			elem1->x=elem->x;
			elem1->y=elem->y;
			elem1->z=elem->z;

			for(int ino = 0; ino < 8; ino++)
			{
				elem1->nodes[ino].color = elem->nodes[ino].color;
				elem1->nodes[ino].fixed = 0;

				elem1->nodes[ino].fixed = elem->nodes[ino].fixed;

				elem1->nodes[ino].id = elem->nodes[ino].id;
				elem1->nodes[ino].x = elem->nodes[ino].x;
				elem1->nodes[ino].y = elem->nodes[ino].y;
				elem1->nodes[ino].z = elem->nodes[ino].z;
			}

			for(int iedge = 0; iedge < 12; iedge++)
			{
				elem1->edge[iedge].coord[0] = elem->edge[iedge].coord[0];
				elem1->edge[iedge].coord[1] = elem->edge[iedge].coord[1];
				elem1->edge[iedge].id = elem->edge[iedge].id;
				elem1->edge[iedge].ref = false;
			}

			for(int isurf = 0; isurf < 6; isurf++) elem1->surf[isurf].ext = elem->surf[isurf].ext;
		}

		sc_array_reset(&toto);
	}

	//id global exterior edges
	for(int iel = 0; iel < mesh->outsurf.elem_count; iel++)
	{
		octant_t* elem = (octant_t*) sc_array_index (&mesh->outsurf, iel);

		//edge 0
		int iedge;
		iedge = 0;
		if(elem->nodes[EdgeVerticesMap[iedge][0]].y == 0 && elem->nodes[EdgeVerticesMap[iedge][1]].y == 0)
		{
			if(elem->nodes[EdgeVerticesMap[iedge][0]].z == 0 && elem->nodes[EdgeVerticesMap[iedge][1]].z == 0)
			{
				elem->edge[iedge].ref = true;
				if(deb) mesh->part_nodes[elem->nodes[EdgeVerticesMap[iedge][0]].id] = iedge;
				if(deb) mesh->part_nodes[elem->nodes[EdgeVerticesMap[iedge][1]].id] = iedge;
			}
		}

		iedge = 1;
		if(elem->nodes[EdgeVerticesMap[iedge][0]].x == 3*mesh->ncellx && elem->nodes[EdgeVerticesMap[iedge][1]].x == 3*mesh->ncellx)
		{
			if(elem->nodes[EdgeVerticesMap[iedge][0]].z == 0 && elem->nodes[EdgeVerticesMap[iedge][1]].z == 0)
			{
				elem->edge[iedge].ref = true;
				if(deb) mesh->part_nodes[elem->nodes[EdgeVerticesMap[iedge][0]].id] = iedge;
				if(deb) mesh->part_nodes[elem->nodes[EdgeVerticesMap[iedge][1]].id] = iedge;
			}
		}

		iedge = 2;
		if(elem->nodes[EdgeVerticesMap[iedge][0]].y == 3*mesh->ncelly && elem->nodes[EdgeVerticesMap[iedge][1]].y == 3*mesh->ncelly)
		{
			if(elem->nodes[EdgeVerticesMap[iedge][0]].z == 0 && elem->nodes[EdgeVerticesMap[iedge][1]].z == 0)
			{
				elem->edge[iedge].ref = true;
				if(deb) mesh->part_nodes[elem->nodes[EdgeVerticesMap[iedge][0]].id] = iedge;
				if(deb) mesh->part_nodes[elem->nodes[EdgeVerticesMap[iedge][1]].id] = iedge;
			}
		}

		iedge = 3;
		if(elem->nodes[EdgeVerticesMap[iedge][0]].x == 0 && elem->nodes[EdgeVerticesMap[iedge][1]].x == 0)
		{
			if(elem->nodes[EdgeVerticesMap[iedge][0]].z == 0 && elem->nodes[EdgeVerticesMap[iedge][1]].z == 0)
			{
				elem->edge[iedge].ref = true;
				if(deb) mesh->part_nodes[elem->nodes[EdgeVerticesMap[iedge][0]].id] = iedge;
				if(deb) mesh->part_nodes[elem->nodes[EdgeVerticesMap[iedge][1]].id] = iedge;
			}
		}

		iedge = 4;
		if(elem->nodes[EdgeVerticesMap[iedge][0]].x == 0 && elem->nodes[EdgeVerticesMap[iedge][1]].x == 0)
		{
			if(elem->nodes[EdgeVerticesMap[iedge][0]].y == 0 && elem->nodes[EdgeVerticesMap[iedge][1]].y == 0)
			{
				elem->edge[iedge].ref = true;
				if(deb) mesh->part_nodes[elem->nodes[EdgeVerticesMap[iedge][0]].id] = iedge;
				if(deb) mesh->part_nodes[elem->nodes[EdgeVerticesMap[iedge][1]].id] = iedge;
			}
		}

		iedge = 5;
		if(elem->nodes[EdgeVerticesMap[iedge][0]].x == 3*mesh->ncellx && elem->nodes[EdgeVerticesMap[iedge][1]].x == 3*mesh->ncellx)
		{
			if(elem->nodes[EdgeVerticesMap[iedge][0]].y == 0 && elem->nodes[EdgeVerticesMap[iedge][1]].y == 0)
			{
				elem->edge[iedge].ref = true;
				if(deb) mesh->part_nodes[elem->nodes[EdgeVerticesMap[iedge][0]].id] = iedge;
				if(deb) mesh->part_nodes[elem->nodes[EdgeVerticesMap[iedge][1]].id] = iedge;
			}
		}

		iedge = 6;
		if(elem->nodes[EdgeVerticesMap[iedge][0]].y == 3*mesh->ncelly && elem->nodes[EdgeVerticesMap[iedge][1]].y == 3*mesh->ncelly)
		{
			if(elem->nodes[EdgeVerticesMap[iedge][0]].x == 3*mesh->ncellx && elem->nodes[EdgeVerticesMap[iedge][1]].x == 3*mesh->ncellx)
			{
				elem->edge[iedge].ref = true;
				if(deb) mesh->part_nodes[elem->nodes[EdgeVerticesMap[iedge][0]].id] = iedge;
				if(deb) mesh->part_nodes[elem->nodes[EdgeVerticesMap[iedge][1]].id] = iedge;
			}
		}

		iedge = 7;
		if(elem->nodes[EdgeVerticesMap[iedge][0]].x == 0 && elem->nodes[EdgeVerticesMap[iedge][1]].x == 0)
		{
			if(elem->nodes[EdgeVerticesMap[iedge][0]].y == 3*mesh->ncelly && elem->nodes[EdgeVerticesMap[iedge][1]].y == 3*mesh->ncelly)
			{
				elem->edge[iedge].ref = true;
				if(deb) mesh->part_nodes[elem->nodes[EdgeVerticesMap[iedge][0]].id] = iedge;
				if(deb) mesh->part_nodes[elem->nodes[EdgeVerticesMap[iedge][1]].id] = iedge;
			}
		}

		iedge = 8;
		if(elem->nodes[EdgeVerticesMap[iedge][0]].y == 0 && elem->nodes[EdgeVerticesMap[iedge][1]].y == 0)
		{
			if(elem->nodes[EdgeVerticesMap[iedge][0]].z == 3*mesh->max_z && elem->nodes[EdgeVerticesMap[iedge][1]].z == 3*mesh->max_z)
			{
				elem->edge[iedge].ref = true;
				if(deb) mesh->part_nodes[elem->nodes[EdgeVerticesMap[iedge][0]].id] = iedge;
				if(deb) mesh->part_nodes[elem->nodes[EdgeVerticesMap[iedge][1]].id] = iedge;
			}
		}

		iedge = 9;
		if(elem->nodes[EdgeVerticesMap[iedge][0]].x == 3*mesh->ncellx && elem->nodes[EdgeVerticesMap[iedge][1]].x == 3*mesh->ncellx)
		{
			if(elem->nodes[EdgeVerticesMap[iedge][0]].z == 3*mesh->max_z && elem->nodes[EdgeVerticesMap[iedge][1]].z == 3*mesh->max_z)
			{
				elem->edge[iedge].ref = true;
				if(deb) mesh->part_nodes[elem->nodes[EdgeVerticesMap[iedge][0]].id] = iedge;
				if(deb) mesh->part_nodes[elem->nodes[EdgeVerticesMap[iedge][1]].id] = iedge;
			}
		}

		iedge = 10;
		if(elem->nodes[EdgeVerticesMap[iedge][0]].y == 3*mesh->ncelly && elem->nodes[EdgeVerticesMap[iedge][1]].y == 3*mesh->ncelly)
		{
			if(elem->nodes[EdgeVerticesMap[iedge][0]].z == 3*mesh->max_z && elem->nodes[EdgeVerticesMap[iedge][1]].z == 3*mesh->max_z)
			{
				elem->edge[iedge].ref = true;
				if(deb) mesh->part_nodes[elem->nodes[EdgeVerticesMap[iedge][0]].id] = iedge;
				if(deb) mesh->part_nodes[elem->nodes[EdgeVerticesMap[iedge][1]].id] = iedge;
			}
		}

		iedge = 11;
		if(elem->nodes[EdgeVerticesMap[iedge][0]].x == 0 && elem->nodes[EdgeVerticesMap[iedge][1]].x == 0)
		{
			if(elem->nodes[EdgeVerticesMap[iedge][0]].z == 3*mesh->max_z && elem->nodes[EdgeVerticesMap[iedge][1]].z == 3*mesh->max_z)
			{
				elem->edge[iedge].ref = true;
				if(deb) mesh->part_nodes[elem->nodes[EdgeVerticesMap[iedge][0]].id] = iedge;
				if(deb) mesh->part_nodes[elem->nodes[EdgeVerticesMap[iedge][1]].id] = iedge;
			}
		}
	}

	if(deb)
	{
		//debug
		for(int iel = 0; iel < mesh->outsurf.elem_count; iel++)
		{
			octant_t* elem = (octant_t*) sc_array_index (&mesh->outsurf, iel);
			for(int ino = 0; ino < 8; ino++)
			{
				if(elem->nodes[ino].fixed == 1) mesh->part_nodes[elem->nodes[ino].id] = 1;
			}
		}

		for(int iel = 0; iel < mesh->outsurf.elem_count; iel++)
		{
			octant_t* elem = (octant_t*) sc_array_index (&mesh->outsurf, iel);

			printf("Sou o elemento: %d\n",elem->id);
			printf("Surface: \n",elem->id);
			for(int isurf = 0; isurf < 6; isurf++) printf("%s ", elem->surf[isurf].ext ? "true" : "false");

			printf("\nEdge: \n",elem->id);
			for(int isurf = 0; isurf < 12; isurf++) printf("%d ",elem->edge[isurf].ref);

			printf("\nNode: \n",elem->id);
			for(int isurf = 0; isurf < 8; isurf++) printf("%d ",elem->nodes[isurf].fixed);
			printf("\n",elem->id);

		}
	}
}

void PillowingInterface(hexa_tree_t* mesh, std::vector<double>& coords, std::vector<int>& nodes_b_mat)
{
	int elem_old = mesh->elements.elem_count;
	int nodes_old = mesh->nodes.elem_count;
	bool clamped = true;
	fprintf(mesh->profile,"Time inside PillowingInterface\n");

	//redo the mapping in the nodes
	auto start = std::chrono::steady_clock::now( );
	printf("     Redo Node Mapping\n");
	RedoNodeMapping(mesh);
	auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>( std::chrono::steady_clock::now( ) - start );
	fprintf(mesh->profile,"    Redo the mapping in the nodes %lld millisecond(s).\n",elapsed.count());
	//std::cout << "Redo the mapping in the nodes "<< elapsed.count() <<" millisecond(s)."<< std::endl;

	//Make the pillow
	start = std::chrono::steady_clock::now( );
	printf("     Pillow Layer\n");
	Pillowing(mesh, coords,nodes_b_mat);
	fprintf(mesh->profile,"    Time in PillowLayer %lld millisecond(s).\n",elapsed.count());
	//std::cout << "Time SurfaceIdentification "<< elapsed.count() <<" millisecond(s)."<< std::endl;

	//Identify the global and local boundaries
	start = std::chrono::steady_clock::now( );
	printf("     Surface Identification\n");
	SurfaceIdentification(mesh, coords);
	fprintf(mesh->profile,"    Time in SurfaceIdentification %lld millisecond(s).\n",elapsed.count());
	//std::cout << "Time SurfaceIdentification "<< elapsed.count() <<" millisecond(s)."<< std::endl;

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
		mesh->part_nodes[ino] = mesh->mpi_rank;
	}
	for(int ino = 0; ino < nodes_b_mat.size();ino++)
	{
		mesh->part_nodes[nodes_b_mat[ino]] = 1;
	}
}
