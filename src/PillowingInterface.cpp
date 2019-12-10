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

typedef struct pillow
{
	uint64_t id;
	uint64_t a;
	uint64_t b;
	int mata;
	int matb;
	//	uint64_t elem_a;
	//	uint64_t elem_b;
	int32_t x,y,z;
	int8_t  list_elem;
	int8_t  list_face[8];
	int32_t  elem[8]; // we must check this values
	int8_t   face[8][3]; // we must check this values
} pillow_t;

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


void Pillowing(hexa_tree_t* mesh, std::vector<double>& coords, std::vector<int>& nodes_b_mat){

	bool deb = false;
	bool clamped = true;
	bool createPillow = true;

	//criando a hash de nos para evitar nos duplicados no pillowing
	sc_hash_array_t*   hash_nodes  = (sc_hash_array_t *)sc_hash_array_new(sizeof(octant_node_t), node_hash_fn , node_equal_fn, &clamped);
	//hash nodes
	for(int ino = 0; ino < mesh->nodes.elem_count; ino++)
	{
		size_t position;
		octant_node_t *r;
		octant_node_t key;
		octant_node_t* node = (octant_node_t*) sc_array_index (&mesh->nodes, ino);
		key.x = node->x;
		key.y = node->y;
		key.z = node->z;
		r = (octant_node_t*) sc_hash_array_insert_unique(hash_nodes, &key, &position);
		if(r!=NULL){
			r->x = node->x;
			r->y = node->y;
			r->z = node->z;
			r->id = node->id;
		}
	}

	//assert(hash_nodes->a.elem_count == mesh->nodes.elem_count);
	//hash nodes nodes_b_mat
	sc_hash_array_t*   hash_b_mat  = (sc_hash_array_t *)sc_hash_array_new(sizeof(octant_node_t), node_hash_fn, node_equal_fn, &clamped);

	//fazendo hash nos nodes_b_mat
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

	//fazendo o pillow
	for(int ioc = 0; ioc < mesh->oct.elem_count; ioc++)
	{
		octree_t* oct = (octree_t*) sc_array_index(&mesh->oct,ioc);
		sc_hash_array_t*   pillow  = (sc_hash_array_t *)sc_hash_array_new(sizeof(pillow_t), pillow_hash_fn, pillow_equal_fn, &clamped);

		//en train de faire la hash du pillowing
		for(int iel = 0; iel < 8; iel++)
		{
			octant_t* elem = (octant_t*) sc_array_index(&mesh->elements, oct->id[iel]);
			for (int ino = 0; ino < 8; ino++)
			{
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

		if(false)
		{
			printf("Nombre de éléments dans Pillow:%d\n", pillow->a.elem_count);
			for(int ino = 0; ino < pillow->a.elem_count; ino ++)
			{
				pillow_t* p = (pillow_t*) sc_array_index(&pillow->a, ino);
				printf("    Noeud %d x:%d y:%d z:%d\n",p->id,p->x,p->y,p->z);
				printf("        les éléments que partagent le noeud sont: ");
				for(int i = 0; i < p->list_elem; i++)
				{
					printf("%d ",p->elem[i]);
				}
				printf("\n");
				printf("        les surfaces sur les éléments que partagent le noeud sont: \n");
				for(int i = 0; i < p->list_elem; i++)
				{
					printf("        pour le élément:%d Surface:",p->elem[i]);
					for(int j = 0; j < p->list_face[i]; j++) printf(" %d",p->face[i][j]);
					printf("\n");
				}
			}
		}

		bool move_node = false;
		//faire la géneration des noeuds
		for(int ino = 0; ino < pillow->a.elem_count; ino ++)
		{
			pillow_t* p = (pillow_t*) sc_array_index(&pillow->a, ino);
			if(deb)printf("    Noeud %d x:%d y:%d z:%d\n",p->id,p->x,p->y,p->z);

			bool saux[3] = {false,false,false};
			for(int iel = 0; iel < p->list_elem; iel++)
			{
				for(int isurf = 0; isurf < p->list_face[iel]; isurf++)
				{
					if(p->face[iel][isurf] == 0 || p->face[iel][isurf] == 1)
					{
						saux[0] = true;
					}
					if(p->face[iel][isurf] == 2 || p->face[iel][isurf] == 3)
					{
						saux[1] = true;
					}
					if(p->face[iel][isurf] == 4 || p->face[iel][isurf] == 5)
					{
						saux[2] = true;
					}
				}
			}
			std::vector<int>auxsurf;
			for(int i = 0; i < 3; i++) if(saux[i])auxsurf.push_back(2*i);

			if(deb) printf("%d %d %d %d %d\n",p->id,auxsurf.size(),saux[0],saux[1],saux[2]);
			if(auxsurf.size()!=0)
			{
				//le noeud base:
				int32_t x = p->x;
				int32_t y = p->y;
				int32_t z = p->z;
				//noeud a
				int32_t xa = p->x;
				int32_t ya = p->y;
				int32_t za = p->z;
				//noeud b
				int32_t xb = p->x;
				int32_t yb = p->y;
				int32_t zb = p->z;

				if(auxsurf.size() == 1)
				{
					if(auxsurf[0] == 0)
					{
						xa += 6;
						xb -= 6;
					}
					if(auxsurf[0] == 2)
					{
						ya += 6;
						yb -= 6;
					}
					if(auxsurf[0] == 4)
					{
						za += 6;
						zb -= 6;
					}
				}
				if(auxsurf.size() == 2)
				{
					//on a deux possibilité à chaque surface
					std::vector<int> aux;
					for(int iel = 0; iel < p->list_elem; iel ++)
					{
						if(p->list_face[iel] == 2)
						{
							aux.push_back(p->face[iel][0]);
							aux.push_back(p->face[iel][1]);
							break;
						}
					}
					std::sort(aux.begin(), aux.end());

					//surface x et y
					if(auxsurf[0] == 0 && auxsurf[1] == 2)
					{
						//diagonal noeud 1-3
						if((aux[0] == 1 && aux[1] == 2) || (aux[0] == 0 && aux[1] == 3))
						{
							xa += 6;
							ya -= 6;
							xb -= 6;
							yb += 6;
						}
						//diagonal noeud 0-2
						if((aux[0] == 1 && aux[1] == 3) || (aux[0] == 0 && aux[1] == 2))
						{
							xa += 6;
							ya += 6;
							xb -= 6;
							yb -= 6;
						}
					}

					//surface x et z
					if(auxsurf[0] == 0 && auxsurf[1] == 4)
					{
						//diagonal noeud 0-5
						if((aux[0] == 1 && aux[1] == 5) || (aux[0] == 0 && aux[1] == 4))
						{
							xa += 6;
							za += 6;
							xb -= 6;
							zb -= 6;
						}
						//diagonal noeud 1-4
						if((aux[0] == 1 && aux[1] == 4) || (aux[0] == 0 && aux[1] == 5))
						{
							xa += 6;
							za -= 6;
							xb -= 6;
							zb += 6;
						}
					}

					//surface y et z
					if(auxsurf[0] == 2 && auxsurf[1] == 4)
					{
						//diagonal noeud 0-7
						if((aux[0] == 2 && aux[1] == 4) || (aux[0] == 3 && aux[1] == 5))
						{
							za -= 6;
							ya -= 6;
							zb += 6;
							yb += 6;
						}
						//diagonal noeud 4-3
						if((aux[0] == 2 && aux[1] == 5) || (aux[0] == 3 && aux[1] == 4))
						{
							za += 6;
							ya -= 6;
							zb -= 6;
							yb += 6;
						}
					}

				}
				if(auxsurf.size() == 3)
				{
					//on a quatre possibilité
					std::vector<int> aux;
					for(int iel = 0; iel < p->list_elem; iel ++)
					{
						if(p->list_face[iel] == 3)
						{
							aux.push_back(p->face[iel][0]);
							aux.push_back(p->face[iel][1]);
							aux.push_back(p->face[iel][2]);
							break;
						}
					}
					//on doit regarder si on a que 2 surface...
					if(aux.size() == 0)
					{
						if(deb)printf("Je suis ici, je suis le noeud %d\n",p->id);
						int elemId = 0;
						for(int iel = 0; iel < p->list_elem; iel ++)
						{
							if(p->list_face[iel] == 2){
								if(p->face[iel][1] == 4 || p->face[iel][1] == 5) continue;
								elemId = p->elem[iel];
								aux.push_back(p->face[iel][0]);
								aux.push_back(p->face[iel][1]);
								if(deb)printf("%d %d\n",aux[0],aux[1]);
								//break;
							}
						}
						for(int iel = 0; iel < 8; iel++)
						{
							if(oct->id[iel] == elemId)
							{
								if(iel < 4)
								{
									aux.push_back(5);
								}
								else
								{
									aux.push_back(4);
								}
							}
						}
						move_node = true;
					}
					std::sort(aux.begin(), aux.end());

					if(deb) printf("%d %d %d\n",aux[0],aux[1],aux[2]);

					//diagonal 0-6
					if((aux[0] == 1 && aux[1] == 3 && aux[2] == 5) ||
							(aux[0] == 0 && aux[1] == 2 && aux[2] == 4))
					{
						if(deb)printf("Diagonal 0-6\n");
						xa += 6;
						ya += 6;
						za += 6;
						xb -= 6;
						yb -= 6;
						zb -= 6;
					}
					//diagonal 1-7
					if((aux[0] == 0 && aux[1] == 3 && aux[2] == 5) ||
							(aux[0] == 1 && aux[1] == 2 && aux[2] == 4))
					{
						if(deb)printf("Diagonal 1-7\n");
						xa -= 6;
						ya += 6;
						za += 6;
						xb += 6;
						yb -= 6;
						zb -= 6;
					}
					//diagonal 2-4
					if((aux[0] == 0 && aux[1] == 2 && aux[2] == 5) ||
							(aux[0] == 1 && aux[1] == 3 && aux[2] == 4))
					{
						if(deb)printf("Diagonal 2-4\n");
						xa -= 6;
						ya -= 6;
						za += 6;
						xb += 6;
						yb += 6;
						zb -= 6;
					}
					//diagonal 3-5
					if((aux[0] == 1 && aux[1] == 2 && aux[2] == 5) ||
							(aux[0] == 0 && aux[1] == 3 && aux[2] == 4))
					{
						if(deb)printf("Diagonal 3-5\n");
						xa += 6;
						ya -= 6;
						za += 6;
						xb -= 6;
						yb += 6;
						zb -= 6;
					}

				}

				if(deb)printf("Node id %d Surface %d\n",p->id,auxsurf.size());
				if(deb)printf("xa %d ya %d za %d\n",xa,ya,za);
				if(deb)printf("xb %d yb %d zb %d\n",xb,yb,zb);

				//ajouter sur la hash de noeuds
				size_t position;
				octant_node_t key;
				key.x = xa;
				key.y = ya;
				key.z = za;
				octant_node_t* ra = (octant_node_t*) sc_hash_array_insert_unique(hash_nodes, &key, &position);
				if(ra!=NULL){
					ra->x = xa;
					ra->y = ya;
					ra->z = za;
					ra->id = hash_nodes->a.elem_count-1;
					p->a = ra->id;
					ra->fixed = 0;
					ra->color = 0;
					coords.push_back(0);
					coords.push_back(0);
					coords.push_back(0);
				}else{
					octant_node_t* ra = (octant_node_t*) sc_array_index(&hash_nodes->a,position);
					p->a = ra->id;
				}

				key.x = xb;
				key.y = yb;
				key.z = zb;
				octant_node_t* rb = (octant_node_t*) sc_hash_array_insert_unique(hash_nodes, &key, &position);
				if(rb!=NULL){
					rb->x = xb;
					rb->y = yb;
					rb->z = zb;
					rb->id = hash_nodes->a.elem_count-1;
					p->b = rb->id;
					rb->fixed = 0;
					rb->color = 0;
					coords.push_back(0);
					coords.push_back(0);
					coords.push_back(0);
				}else{
					octant_node_t* rb = (octant_node_t*) sc_array_index(&hash_nodes->a,position);
					p->b = rb->id;
				}
			}

			//trouver les coords
			octant_node_t * a = (octant_node_t*) sc_array_index(&hash_nodes->a,p->a);
			octant_node_t * b = (octant_node_t*) sc_array_index(&hash_nodes->a,p->b);
			bool ina[8];
			bool inb[8];
			for(int iel = 0; iel < 8; iel++)
			{
				octant_t* elem = (octant_t*) sc_array_index(&mesh->elements, oct->id[iel]);
				ina[iel] = false;
				inb[iel] = false;
				if(elem->nodes[0].x <= a->x && elem->nodes[6].x >= a->x &&
						elem->nodes[0].y <= a->y && elem->nodes[6].y >= a->y &&
						elem->nodes[0].z <= a->z && elem->nodes[6].z >= a->z) ina[iel] = true;
				if(elem->nodes[0].x <= b->x && elem->nodes[6].x >= b->x &&
						elem->nodes[0].y <= b->y && elem->nodes[6].y >= b->y &&
						elem->nodes[0].z <= b->z && elem->nodes[6].z >= b->z) inb[iel] = true;

				if(ina[iel] || inb[iel])
				{
					if(deb)printf("%d %d\n",elem->id,p->id);
					if(deb)printf("%d %d %d %d %s\n",a->id,a->x,a->y,a->z, ina[iel] ? "true" : "false");
					if(deb)printf("%d %d %d %d %s\n",b->id,b->x,b->y,b->z, inb[iel] ? "true" : "false");
					if(deb)printf("%d %d %d %d %d %d\n",elem->nodes[0].x,elem->nodes[6].x,elem->nodes[0].y,elem->nodes[6].y,
							elem->nodes[0].z,elem->nodes[6].z);

					double cord_in_ref[3];
					int nnode;
					double ref_in_x[8],ref_in_y[8],ref_in_z[8];
					int color = 0;
					int colorc = 0;
					for (int ino = 0; ino < 8; ino++)
					{
						ref_in_x[ino]=coords[3*elem->nodes[ino].id+0];
						ref_in_y[ino]=coords[3*elem->nodes[ino].id+1];
						ref_in_z[ino]=coords[3*elem->nodes[ino].id+2];
						if(elem->nodes[ino].color != 0)
						{
							color += elem->nodes[ino].color;
							colorc++;
						}
					}

					if(ina[iel])
					{
						cord_in_ref[0] = (a->x - elem->nodes[0].x)/6 -1;
						cord_in_ref[1] = (a->y - elem->nodes[0].y)/6 -1;
						cord_in_ref[2] = (a->z - elem->nodes[0].z)/6 -1;
						nnode = p->a;
						p->mata = elem->n_mat;
						a->color = int(color)/int(colorc);
						if(deb)printf("Id do no %d, aqui achei no a:%d no elemento %d que tem o material %d\n",
								p->id,p->a,elem->id,elem->n_mat);
					}
					if(inb[iel])
					{
						cord_in_ref[0] = (b->x - elem->nodes[0].x)/6 -1;
						cord_in_ref[1] = (b->y - elem->nodes[0].y)/6 -1;
						cord_in_ref[2] = (b->z - elem->nodes[0].z)/6 -1;
						nnode = p->b;
						p->matb = elem->n_mat;
						b->color = int(color)/int(colorc);
						if(deb)printf("Id do no %d, aqui achei no b:%d no elemento %d que tem o material %d\n",
								p->id,p->b,elem->id,elem->n_mat);
					}
					GtsPoint* ref_point = LinearMapHex(cord_in_ref, ref_in_x, ref_in_y, ref_in_z);

					coords[3*nnode + 0] = ref_point->x;
					coords[3*nnode + 1] = ref_point->y;
					coords[3*nnode + 2] = ref_point->z;
				}
			}
		}

		//faire la criation des éléments pour le pillowing
		sc_array_t toto;
		sc_array_init(&toto, sizeof(octant_t));
		for(int iel = 0; iel < 8; iel++)
		{
			octant_t* elemOrig = (octant_t*) sc_array_index(&mesh->elements, oct->id[iel]);
			octant_t * elem = (octant_t*) sc_array_push(&toto);

			hexa_element_copy(elemOrig,elem);

			//on doit voir si le élément demande de extrusion...
			bool surf[6];
			std::vector<int> aux;
			int surf_count = 0;
			for(int isurf = 0; isurf < 6; isurf++)
			{
				surf[isurf] = true;
				for(int ino = 0; ino < 4; ino++)
				{
					size_t position;
					octant_node_t key;
					key.x = elem->nodes[FaceNodesMap[isurf][ino]].x;
					key.y = elem->nodes[FaceNodesMap[isurf][ino]].y;
					key.z = elem->nodes[FaceNodesMap[isurf][ino]].z;

					bool lnode = sc_hash_array_lookup(hash_b_mat, &key, &position);
					if(!lnode) surf[isurf] = false;
				}
				if(surf[isurf]) surf_count++;
				if(surf[isurf]) aux.push_back(isurf);
			}

			if(deb) printf("Para o elemento %d eu vou extrudar %d elementos\n",elem->id,aux.size());

			//on garde la connectivite original...
			int connec[8];
			for(int ino = 0; ino <8; ino++){
				connec[ino] = elem->nodes[ino].id;
			}

			//extrusion se fait à chaque surface du élément
			for(int isurf = 0; isurf < aux.size(); isurf++)
			{
				int new_nodes[4];
				for(int ino = 0; ino < 4; ino++)
				{
					size_t position;
					pillow_t key;
					key.x = elem->nodes[FaceNodesMap[aux[isurf]][ino]].x;
					key.y = elem->nodes[FaceNodesMap[aux[isurf]][ino]].y;
					key.z = elem->nodes[FaceNodesMap[aux[isurf]][ino]].z;

					bool lnode = sc_hash_array_lookup(pillow, &key, &position);
					pillow_t* p = (pillow_t*) sc_array_index(&pillow->a,position);
					if(deb) printf("Id do no %d, mata:%d matb:%d elem mat %d\n",p->id,p->mata,p->matb,elem->n_mat);
					if(elem->n_mat == p->mata)
					{
						new_nodes[ino] = p->a;
					}
					else
					{
						new_nodes[ino] = p->b;
					}
				}

				if(createPillow)
				{
					//créer l'élément
					octant_t* pelem;
					pelem = (octant_t*) sc_array_push(&mesh->elements);
					pelem->id = mesh->elements.elem_count-1;
					pelem->n_mat = elem->n_mat;
					pelem->pad = elem->pad;
					pelem->x = elem->x;
					pelem->y = elem->y;
					pelem->z = elem->z;

					//printf("  Face:%d\n",aux[isurf]);
					for(int ino = 0; ino < 4; ino++)
					{
						pelem->nodes[FaceNodesMap[aux[isurf]][ino]].id = connec[FaceNodesMap[aux[isurf]][ino]];
						pelem->nodes[FaceNodesMap[aux[isurf]][ino]].fixed = elem->nodes[FaceNodesMap[aux[isurf]][ino]].fixed;
						pelem->nodes[FaceNodesMap_inv[aux[isurf]][ino]].id= new_nodes[ino];

						elemOrig->nodes[FaceNodesMap[aux[isurf]][ino]].id = new_nodes[ino];

						pelem->nodes[FaceNodesMap[aux[isurf]][ino]].fixed= 0;
						octant_node_t* node = (octant_node_t*) sc_array_index(&hash_nodes->a,connec[FaceNodesMap[aux[isurf]][ino]]);
						pelem->nodes[FaceNodesMap[aux[isurf]][ino]].x = node->x;
						pelem->nodes[FaceNodesMap[aux[isurf]][ino]].y = node->y;
						pelem->nodes[FaceNodesMap[aux[isurf]][ino]].z = node->z;
						node = (octant_node_t*) sc_array_index(&hash_nodes->a,new_nodes[ino]);
						pelem->nodes[FaceNodesMap_inv[aux[isurf]][ino]].x = node->x;
						pelem->nodes[FaceNodesMap_inv[aux[isurf]][ino]].y = node->y;
						pelem->nodes[FaceNodesMap_inv[aux[isurf]][ino]].z = node->z;
					}
				}
			}

			if(aux.size() == 0 || move_node)
			{
				//chercher le noeuds pour trouver la bonne connectivite dans les éléments
				for (int ino = 0; ino < 8; ino++)
				{
					int ina = 0;
					int inb = 0;
					pillow_t keyP;
					size_t positionP;
					keyP.x = elem->nodes[ino].x;
					keyP.y = elem->nodes[ino].y;
					keyP.z = elem->nodes[ino].z;
					bool lnode = sc_hash_array_lookup(pillow, &keyP, &positionP);
					if(lnode)
					{
						pillow_t * p = (pillow_t*) sc_array_index(&pillow->a,positionP);

						octant_node_t * a = (octant_node_t*) sc_array_index(&hash_nodes->a,p->a);
						octant_node_t * b = (octant_node_t*) sc_array_index(&hash_nodes->a,p->b);
						if(deb) printf("%d %d %d %d %d %d\n",elem->id,p->id,p->a, a->id,p->b,b->id);

						if(createPillow)
						{
							if(elem->n_mat == p->mata)
							{
								elemOrig->nodes[ino].id = p->a;
								elemOrig->nodes[ino].x = a->x;
								elemOrig->nodes[ino].y = a->y;
								elemOrig->nodes[ino].z = a->z;
							}
							else
							{
								elemOrig->nodes[ino].id = p->b;
								elemOrig->nodes[ino].x = b->x;
								elemOrig->nodes[ino].y = b->y;
								elemOrig->nodes[ino].z = b->z;
							}
						}
					}
				}
			}
			sc_array_reset(&toto);
		}
		sc_hash_array_destroy(pillow);
	}

	//faire la atualisation de mon tableau de noeuds
	sc_array_reset(&mesh->nodes);
	sc_hash_array_rip(hash_nodes,&mesh->nodes);

	//update the x y z coords in mesh->elements
	for(int oc = 0; oc < mesh->oct.elem_count; oc++)
	{
		octree_t* oct = (octree_t*) sc_array_index(&mesh->oct,oc);
		for(int iel = 0; iel < 8; iel++)
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
}

void Pillowing_oct(hexa_tree_t* mesh, std::vector<double>& coords, std::vector<int>& nodes_b_mat){

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


	//faire le pillow mesh->oct.elem_count
	for(int ioc = 0; ioc < mesh->oct.elem_count; ioc++)
	{
		octree_t* oct = (octree_t*) sc_array_index(&mesh->oct,ioc);
		sc_hash_array_t*   pillow  = (sc_hash_array_t *)sc_hash_array_new(sizeof(pillow_t), pillow_hash_fn, pillow_equal_fn, &clamped);

		int elName = -1;
		//en train de faire la hash du pillowing
		for(int iel = 0; iel < 8; iel++)
		{
			octant_t* elem = (octant_t*) sc_array_index(&mesh->elements, oct->id[iel]);
			for (int ino = 0; ino < 8; ino++)
			{
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

			//trouver l'element avec 3 surfaces
			bool surf[6];
			std::vector<int> aux;
			for(int isurf = 0; isurf < 6; isurf++)
			{
				surf[isurf] = true;
				for(int ino = 0; ino < 4; ino++)
				{
					size_t position;
					octant_node_t key;
					key.x = elem->nodes[FaceNodesMap[isurf][ino]].x;
					key.y = elem->nodes[FaceNodesMap[isurf][ino]].y;
					key.z = elem->nodes[FaceNodesMap[isurf][ino]].z;

					bool lnode = sc_hash_array_lookup(hash_b_mat, &key, &position);
					if(!lnode) surf[isurf] = false;
				}
				if(surf[isurf]) aux.push_back(isurf);
			}

			if(aux.size() == 3) elName = iel;
		}

		sc_hash_array_t*   pauxhash  = (sc_hash_array_t *)sc_hash_array_new(sizeof(aux_teste_t),aux_teste_fn ,aux_teste_equal_fn , &clamped);

		//en train de bouger que les noeuds que on besoin pas faira la extrution
		octant_t* elemaux[8];
		for(int iel = 0; iel < 8; iel++)
		{
			elemaux[iel] = (octant_t*) sc_array_index(&mesh->elements, oct->id[iel]);
		}
		for(int iel = 0; iel < 8; iel++)
		{
			double cord_in_ref[3];

			//verifica quais dos 8 nos estao bloqueados
			int ncount = 0;
			int nnode[8];
			for(int ino = 0; ino < 8; ino++)
			{
				size_t position;
				octant_node_t key;
				key.x = elemaux[iel]->nodes[ino].x;
				key.y = elemaux[iel]->nodes[ino].y;
				key.z = elemaux[iel]->nodes[ino].z;
				key.id = elemaux[iel]->nodes[ino].id;

				bool lnode = sc_hash_array_lookup(hash_b_mat, &key, &position);

				if(lnode)
				{
					nnode[ncount] = ino;
					ncount++;
				}
			}

			//agora verifica se o no deve ser movido
			if(ncount == 1)
			{
				int x,y,z;
				int aaux;
				if(elName == 0)
				{
					x = elemaux[iel]->nodes[6].x + 6;
					y = elemaux[iel]->nodes[6].y + 6;
					z = elemaux[iel]->nodes[6].z + 6;
					aaux = 6;
				}
				if(elName == 1)
				{
					x = elemaux[iel]->nodes[7].x + 6;
					y = elemaux[iel]->nodes[7].y + 6;
					z = elemaux[iel]->nodes[7].z + 6;
					aaux = 7;
				}
				if(elName == 2)
				{
					x = elemaux[iel]->nodes[4].x + 6;
					y = elemaux[iel]->nodes[4].y + 6;
					z = elemaux[iel]->nodes[4].z + 6;
					aaux = 4;
				}
				if(elName == 3)
				{
					x = elemaux[iel]->nodes[5].x + 6;
					y = elemaux[iel]->nodes[5].y + 6;
					z = elemaux[iel]->nodes[5].z + 6;
					aaux = 5;
				}

				if(elName == 4)
				{
					x = elemaux[iel]->nodes[2].x + 6;
					y = elemaux[iel]->nodes[2].y + 6;
					z = elemaux[iel]->nodes[2].z + 6;
					aaux = 2;
				}
				if(elName == 5)
				{
					x = elemaux[iel]->nodes[3].x + 6;
					y = elemaux[iel]->nodes[3].y + 6;
					z = elemaux[iel]->nodes[3].z + 6;
					aaux = 3;
				}
				if(elName == 6)
				{
					x = elemaux[iel]->nodes[0].x + 6;
					y = elemaux[iel]->nodes[0].y + 6;
					z = elemaux[iel]->nodes[0].z + 6;
					aaux = 0;
				}
				if(elName == 7)
				{
					x = elemaux[iel]->nodes[1].x + 6;
					y = elemaux[iel]->nodes[1].y + 6;
					z = elemaux[iel]->nodes[1].z + 6;
					aaux = 1;
				}

				double ref_in_x[8],ref_in_y[8],ref_in_z[8];
				for (int ino = 0; ino < 8; ino++)
				{
					ref_in_x[ino] = coords[3*elemaux[aaux]->nodes[ino].id+0];
					ref_in_y[ino] = coords[3*elemaux[aaux]->nodes[ino].id+1];
					ref_in_z[ino] = coords[3*elemaux[aaux]->nodes[ino].id+2];
				}

				cord_in_ref[0] = 0;
				cord_in_ref[1] = 0;
				cord_in_ref[2] = 0;

				//ajouter sur la hash de noeuds
				size_t position;
				octant_node_t key;
				key.x = x;
				key.y = y;
				key.z = z;
				key.id = elemaux[iel]->nodes[nnode[0]].id;

				GtsPoint* ref_point = LinearMapHex(cord_in_ref, ref_in_x, ref_in_y, ref_in_z);
				int nodeID;
				octant_node_t* r = (octant_node_t*) sc_hash_array_insert_unique(hash_nodes, &key, &position);
				if(r!=NULL){
					r->x = x;
					r->y = y;
					r->z = z;
					r->id = hash_nodes->a.elem_count-1;
					nodeID = r->id;
					coords.push_back(ref_point->x);
					coords.push_back(ref_point->y);
					coords.push_back(ref_point->z);
				}else{
					octant_node_t* r = (octant_node_t*) sc_array_index(&hash_nodes->a,position);
					nodeID = r->id;
				}
				//colocando o novo no na hash auxiliar
				aux_teste_t keyaux;
				keyaux.orig.x = elemaux[iel]->nodes[nnode[0]].x;
				keyaux.orig.y = elemaux[iel]->nodes[nnode[0]].y;
				keyaux.orig.z = elemaux[iel]->nodes[nnode[0]].z;
				keyaux.orig.id = elemaux[iel]->nodes[nnode[0]].id;

				elemaux[iel]->nodes[nnode[0]].id = nodeID;
				elemaux[iel]->nodes[nnode[0]].x = x;
				elemaux[iel]->nodes[nnode[0]].y = y;
				elemaux[iel]->nodes[nnode[0]].z = z;

				aux_teste_t* paux = (aux_teste_t*) sc_hash_array_insert_unique(pauxhash, &keyaux, &position);
				if(paux!=NULL){
					paux->orig.x = keyaux.orig.x;
					paux->orig.y = keyaux.orig.y;
					paux->orig.z = keyaux.orig.z;
					paux->orig.id = keyaux.orig.id;

					paux->pill.x = x;
					paux->pill.y = y;
					paux->pill.z = z;
					paux->pill.id = nodeID;
				}
			}
		}

		//en train de bouger que les lignes que on besoin pas faira la extrution
		for(int iel = 0; iel < 8; iel++)
		{
			octant_t* elem = (octant_t*) sc_array_index(&mesh->elements, oct->id[iel]);
			double cord_in_ref[3];
			double cord_in_refa[3];
			double cord_in_refb[3];
			double ref_in_x[8],ref_in_y[8],ref_in_z[8];
			for (int ino = 0; ino < 8; ino++)
			{
				ref_in_x[ino] = coords[3*elem->nodes[ino].id+0];
				ref_in_y[ino] = coords[3*elem->nodes[ino].id+1];
				ref_in_z[ino] = coords[3*elem->nodes[ino].id+2];
			}

			//verifica quais dos 8 nos estao bloqueados
			int ncount = 0;
			int nnode[8];
			for(int ino = 0; ino < 8; ino++)
			{
				size_t position;
				octant_node_t key;
				key.x = elem->nodes[ino].x;
				key.y = elem->nodes[ino].y;
				key.z = elem->nodes[ino].z;
				key.id = elem->nodes[ino].id;

				bool lnode = sc_hash_array_lookup(hash_b_mat, &key, &position);

				if(lnode)
				{
					nnode[ncount] = ino;
					ncount++;
				}
			}

			//agora verifica se uma aresta deve ser movida
			if(ncount == 2)
			{
				int xa,xb;
				int ya,yb;
				int za,zb;

				//sempre vai ser interior
				//preciso achar o eixo...
				if((nnode[0] == 0 && nnode[1] == 1) ||
						(nnode[0] == 4 && nnode[1] == 5))
				{
					xa = elem->nodes[0].x;
					xb = elem->nodes[1].x;

					ya = elem->nodes[0].y + 6;
					yb = elem->nodes[0].y + 6;

					za = elem->nodes[0].z + 6;
					zb = elem->nodes[0].z + 6;

					cord_in_refa[0] = -1;
					cord_in_refa[1] = 0;
					cord_in_refa[2] = 0;

					cord_in_refb[0] = 1;
					cord_in_refb[1] = 0;
					cord_in_refb[2] = 0;
				}

				if((nnode[0] == 2 && nnode[1] == 3) ||
						(nnode[0] == 6 && nnode[1] == 7))
				{
					xa = elem->nodes[2].x;
					xb = elem->nodes[3].x;

					ya = elem->nodes[0].y + 6;
					yb = elem->nodes[0].y + 6;

					za = elem->nodes[0].z + 6;
					zb = elem->nodes[0].z + 6;

					cord_in_refa[0] = 1;
					cord_in_refa[1] = 0;
					cord_in_refa[2] = 0;

					cord_in_refb[0] = -1;
					cord_in_refb[1] = 0;
					cord_in_refb[2] = 0;
				}

				if((nnode[0] == 1 && nnode[1] == 2) ||
						(nnode[0] == 0 && nnode[1] == 3))
				{
					xa = elem->nodes[0].x + 6;
					xb = elem->nodes[0].x + 6;

					ya = elem->nodes[0].y;
					yb = elem->nodes[3].y;

					za = elem->nodes[0].z + 6;
					zb = elem->nodes[0].z + 6;

					cord_in_refa[0] = 0;
					cord_in_refa[1] = -1;
					cord_in_refa[2] = 0;

					cord_in_refb[0] = 0;
					cord_in_refb[1] = 1;
					cord_in_refb[2] = 0;
				}

				if((nnode[0] == 5 && nnode[1] == 6) ||
						(nnode[0] == 4 && nnode[1] == 7))
				{
					xa = elem->nodes[0].x + 6;
					xb = elem->nodes[0].x + 6;

					ya = elem->nodes[0].y;
					yb = elem->nodes[3].y;

					za = elem->nodes[0].z + 6;
					zb = elem->nodes[0].z + 6;

					cord_in_refa[0] = 0;
					cord_in_refa[1] = -1;
					cord_in_refa[2] = 0;

					cord_in_refb[0] = 0;
					cord_in_refb[1] = 1;
					cord_in_refb[2] = 0;
				}

				if((nnode[0] == 0 && nnode[1] == 4) ||
						(nnode[0] == 1 && nnode[1] == 5) ||
						(nnode[0] == 2 && nnode[1] == 6) ||
						(nnode[0] == 3 && nnode[1] == 7))
				{
					xa = elem->nodes[0].x + 6 ;
					xb = elem->nodes[0].x + 6;

					ya = elem->nodes[0].y + 6;
					yb = elem->nodes[0].y + 6;

					za = elem->nodes[0].z;
					zb = elem->nodes[4].z;

					cord_in_refa[0] = 0;
					cord_in_refa[1] = 0;
					cord_in_refa[2] = -1;

					cord_in_refb[0] = 0;
					cord_in_refb[1] = 0;
					cord_in_refb[2] = 1;
				}

				aux_teste_t keyt;
				size_t position;
				keyt.orig.x = elem->nodes[nnode[0]].x;
				keyt.orig.y = elem->nodes[nnode[0]].y;
				keyt.orig.z = elem->nodes[nnode[0]].z;

				bool lnode = sc_hash_array_lookup(pauxhash, &keyt, &position);

				if(lnode)
				{
					//printf("Entrei no a\n");
					aux_teste_t* p = (aux_teste_t*) sc_array_index(&pauxhash->a,position);
					xa = p->pill.x;
					ya = p->pill.y;
					za = p->pill.z;
				}

				//ajouter sur la hash de noeuds
				octant_node_t key;
				//primeiro dos dois nos
				key.x = xa;
				key.y = ya;
				key.z = za;
				key.id = elem->nodes[nnode[0]].id;

				GtsPoint* ref_point = LinearMapHex(cord_in_refa, ref_in_x, ref_in_y, ref_in_z);
				int nodeID;
				octant_node_t* r = (octant_node_t*) sc_hash_array_insert_unique(hash_nodes, &key, &position);
				if(r!=NULL){
					r->x = xa;
					r->y = ya;
					r->z = za;
					r->id = hash_nodes->a.elem_count-1;
					nodeID = r->id;
					coords.push_back(ref_point->x);
					coords.push_back(ref_point->y);
					coords.push_back(ref_point->z);
				}else{
					octant_node_t* r = (octant_node_t*) sc_array_index(&hash_nodes->a,position);
					nodeID = r->id;
				}
				//colocando o novo no na hash auxiliar
				aux_teste_t keyaux;
				keyaux.orig.x = elem->nodes[nnode[0]].x;
				keyaux.orig.y = elem->nodes[nnode[0]].y;
				keyaux.orig.z = elem->nodes[nnode[0]].z;
				keyaux.orig.id = elem->nodes[nnode[0]].id;

				elem->nodes[nnode[0]].id = nodeID;
				elem->nodes[nnode[0]].x = xa;
				elem->nodes[nnode[0]].y = ya;
				elem->nodes[nnode[0]].z = za;

				aux_teste_t* paux = (aux_teste_t*) sc_hash_array_insert_unique(pauxhash, &keyaux, &position);
				if(paux!=NULL){
					paux->orig.x = keyaux.orig.x;
					paux->orig.y = keyaux.orig.y;
					paux->orig.z = keyaux.orig.z;
					paux->orig.id = keyaux.orig.id;

					paux->pill.x = xa;
					paux->pill.y = ya;
					paux->pill.z = za;
					paux->pill.id = nodeID;
				}

				//agora o segundo
				keyt.orig.x = elem->nodes[nnode[1]].x;
				keyt.orig.y = elem->nodes[nnode[1]].y;
				keyt.orig.z = elem->nodes[nnode[1]].z;

				lnode = sc_hash_array_lookup(pauxhash, &keyt, &position);

				if(lnode)
				{
					//printf("Entrei no b\n");
					aux_teste_t* p = (aux_teste_t*) sc_array_index(&pauxhash->a,position);
					xb = p->pill.x;
					yb = p->pill.y;
					zb = p->pill.z;
				}

				key.x = xb;
				key.y = yb;
				key.z = zb;
				key.id = elem->nodes[nnode[1]].id;

				ref_point = LinearMapHex(cord_in_refb, ref_in_x, ref_in_y, ref_in_z);
				r = (octant_node_t*) sc_hash_array_insert_unique(hash_nodes, &key, &position);
				if(r!=NULL){
					r->x = xb;
					r->y = yb;
					r->z = zb;
					r->id = hash_nodes->a.elem_count-1;
					nodeID = r->id;
					coords.push_back(ref_point->x);
					coords.push_back(ref_point->y);
					coords.push_back(ref_point->z);
				}else{
					octant_node_t* r = (octant_node_t*) sc_array_index(&hash_nodes->a,position);
					nodeID = r->id;
				}
				//colocando o novo no na hash auxiliar
				keyaux.orig.x = elem->nodes[nnode[1]].x;
				keyaux.orig.y = elem->nodes[nnode[1]].y;
				keyaux.orig.z = elem->nodes[nnode[1]].z;
				keyaux.orig.id = elem->nodes[nnode[1]].id;

				elem->nodes[nnode[1]].id = nodeID;
				elem->nodes[nnode[1]].x = xb;
				elem->nodes[nnode[1]].y = yb;
				elem->nodes[nnode[1]].z = zb;

				paux = (aux_teste_t*) sc_hash_array_insert_unique(pauxhash, &keyaux, &position);
				if(paux!=NULL){
					paux->orig.x = keyaux.orig.x;
					paux->orig.y = keyaux.orig.y;
					paux->orig.z = keyaux.orig.z;
					paux->orig.id = keyaux.orig.id;

					paux->pill.x = xb;
					paux->pill.y = yb;
					paux->pill.z = zb;
					paux->pill.id = nodeID;
				}
			}

		}
		//printf("Tamanho da bagaca %d\n",pauxhash->a.elem_count);
		//for(int ino = 0; ino < pauxhash->a.elem_count; ino++)
		//{
		//		aux_teste_t* aux = (aux_teste_t*) sc_array_index(&pauxhash->a, ino);
		//	printf("%d %d %d %d\n",aux->orig.id,aux->orig.x,aux->orig.y,aux->orig.z);
		//}
		sc_array_t toto;
		sc_array_init(&toto, sizeof(octant_t));
		//en train de faire le pillowing
		for(int iel = 0; iel < 8; iel++)
		{
			octant_t* elem = (octant_t*) sc_array_index(&mesh->elements, oct->id[iel]);
			octant_t * elemcp = (octant_t*) sc_array_push(&toto);

			hexa_element_copy(elem,elemcp);
			double cord_in_ref[3];
			double cord_in_refa[3];
			double cord_in_refb[3];
			double ref_in_x[8],ref_in_y[8],ref_in_z[8];
			for (int ino = 0; ino < 8; ino++)
			{
				ref_in_x[ino] = coords[3*elem->nodes[ino].id+0];
				ref_in_y[ino] = coords[3*elem->nodes[ino].id+1];
				ref_in_z[ino] = coords[3*elem->nodes[ino].id+2];
			}

			//verifica quais dos 8 nos estao bloqueados
			int ncount = 0;
			int nnode[8];
			for(int ino = 0; ino < 8; ino++)
			{
				size_t position;
				octant_node_t key;
				key.x = elem->nodes[ino].x;
				key.y = elem->nodes[ino].y;
				key.z = elem->nodes[ino].z;
				key.id = elem->nodes[ino].id;

				bool lnode = sc_hash_array_lookup(hash_b_mat, &key, &position);

				if(lnode)
				{
					nnode[ncount] = ino;
					ncount++;
				}
			}

			//on garde la connectivite original...
			int connec[8];
			for(int ino = 0; ino < 8; ino++){
				connec[ino] = elemcp->nodes[ino].id;
			}

			//agora verifica as faces que devem ser extrudada
			if(ncount > 3)
			{
				//trouver la surface(s)
				std::vector<int> aux;
				bool surf[6];
				for(int isurf = 0; isurf < 6; isurf++)
				{
					surf[isurf] = true;
					for(int ino = 0; ino < 4; ino++)
					{
						size_t position;
						octant_node_t key;
						key.x = elem->nodes[FaceNodesMap[isurf][ino]].x;
						key.y = elem->nodes[FaceNodesMap[isurf][ino]].y;
						key.z = elem->nodes[FaceNodesMap[isurf][ino]].z;

						bool lnode = sc_hash_array_lookup(hash_b_mat, &key, &position);
						if(!lnode) surf[isurf] = false;
					}
					if(surf[isurf]) aux.push_back(isurf);
				}

				//achar os 4 vertices das faces
				for(int isurf = 0; isurf < aux.size(); isurf++)
				{
					std::vector<int> caux;
					//verificando o tipo de tratamento que sera dado na sequencia....
					//se caux tem o mesmo tamanho que aux tudo certo, so gerar os elementos
					//se caux é maior que aux a geracao do no precisa
					//ser diferente...
					for(int ino = 0; ino < 4; ino++)
					{
						size_t position;
						pillow_t key;
						key.x = elem->nodes[FaceNodesMap[aux[isurf]][ino]].x;
						key.y = elem->nodes[FaceNodesMap[aux[isurf]][ino]].y;
						key.z = elem->nodes[FaceNodesMap[aux[isurf]][ino]].z;
						key.id = elem->nodes[FaceNodesMap[aux[isurf]][ino]].id;

						bool lnode = sc_hash_array_lookup(pillow, &key, &position);
						pillow_t* p = (pillow_t*) sc_array_index(&pillow->a,position);

						for(int i = 0; i < p->list_elem; i++)
						{
							for(int j = 0; j < p->list_face[i]; j ++)
							{
								caux.push_back(p->face[i][j]);
								if(aux[isurf] == 0 || aux[isurf] == 1 && false)
								{
									if(p->face[i][j] != 0 && p->face[i][j] != 1)
									{
										//printf("Elemento %d, face %d no noh %d, achei o seguinte %d no elemento %d\n",elem->id,aux[isurf],FaceNodesMap[aux[isurf]][ino],p->face[i][j],p->elem[i]);
										caux.push_back(p->face[i][j]);
									}
								}
								if(aux[isurf] == 2 || aux[isurf] == 3 && false)
								{
									if(p->face[i][j] != 2 && p->face[i][j] != 3)
									{
										//printf("Elemento %d, face %d no noh %d, achei o seguinte %d no elemento %d\n",elem->id,aux[isurf],FaceNodesMap[aux[isurf]][ino],p->face[i][j],p->elem[i]);
										caux.push_back(p->face[i][j]);
									}
								}
								if(aux[isurf] == 4 || aux[isurf] == 5 && false)
								{
									if(p->face[i][j] != 4 && p->face[i][j] != 5)
									{
										//printf("Elemento %d, face %d no noh %d, achei o seguinte %d no elemento %d\n",elem->id,aux[isurf],FaceNodesMap[aux[isurf]][ino],p->face[i][j],p->elem[i]);
										caux.push_back(p->face[i][j]);
									}
								}
							}
						}
					}

					//processando caux para saber quantos unicos...
					std::sort(caux.begin(),caux.end());
					caux.erase(std::unique(caux.begin(),caux.end()),caux.end());
					int scount = 0;
					bool su[3] = {false,false,false};
					for(int ic = 0; ic < caux.size(); ic++)
					{
						if(caux[ic] == 0 || caux[ic] == 1) su[0] = true;
						if(caux[ic] == 2 || caux[ic] == 3) su[1] = true;
						if(caux[ic] == 4 || caux[ic] == 5) su[2] = true;
					}
					if(su[0]) scount++;
					if(su[1]) scount++;
					if(su[2]) scount++;
					bool exel = true;
					if(scount != aux.size()) exel = false;
					if(scount == aux.size()) exel = true;

					int new_id[4];
					for(int ino = 0; ino < 4; ino++)
					{
						int x,y,z;
						if(aux.size() == 1)
						{
							x = elem->nodes[FaceNodesMap[aux[isurf]][ino]].x;
							y = elem->nodes[FaceNodesMap[aux[isurf]][ino]].y;
							z = elem->nodes[FaceNodesMap[aux[isurf]][ino]].z;

							if(aux[0] == 0)
							{
								if(ino == 0)
								{
									x += 6;

									cord_in_ref[0] = 0;
									cord_in_ref[1] = -1;
									cord_in_ref[2] = -1;
								}
								if(ino == 1)
								{
									x += 6;

									cord_in_ref[0] = 0;
									cord_in_ref[1] = -1;
									cord_in_ref[2] = 1;
								}
								if(ino == 2)
								{
									x += 6;

									cord_in_ref[0] = 0;
									cord_in_ref[1] = 1;
									cord_in_ref[2] = 1;
								}
								if(ino == 3)
								{
									x += 6;

									cord_in_ref[0] = 0;
									cord_in_ref[1] = 1;
									cord_in_ref[2] = -1;
								}
							}
							if(aux[0] == 1)
							{
								if(ino == 0)
								{
									x -= 6;

									cord_in_ref[0] = 0;
									cord_in_ref[1] = -1;
									cord_in_ref[2] = -1;
								}
								if(ino == 1)
								{
									x -= 6;

									cord_in_ref[0] = 0;
									cord_in_ref[1] = -1;
									cord_in_ref[2] = 1;
								}
								if(ino == 2)
								{
									x -= 6;

									cord_in_ref[0] = 0;
									cord_in_ref[1] = 1;
									cord_in_ref[2] = 1;
								}
								if(ino == 3)
								{
									x -= 6;

									cord_in_ref[0] = 0;
									cord_in_ref[1] = 1;
									cord_in_ref[2] = -1;
								}
							}
							if(aux[0] == 2)
							{
								if(ino == 0)
								{
									y += 6;

									cord_in_ref[0] = 1;
									cord_in_ref[1] = 0;
									cord_in_ref[2] = -1;
								}
								if(ino == 1)
								{
									y += 6;

									cord_in_ref[0] = 1;
									cord_in_ref[1] = 0;
									cord_in_ref[2] = 1;
								}
								if(ino == 2)
								{
									y += 6;

									cord_in_ref[0] = -1;
									cord_in_ref[1] = 0;
									cord_in_ref[2] = 1;
								}
								if(ino == 3)
								{
									y += 6;

									cord_in_ref[0] = -1;
									cord_in_ref[1] = 0;
									cord_in_ref[2] = -1;
								}
							}
							if(aux[0] == 3)
							{
								if(ino == 0)
								{
									y -= 6;

									cord_in_ref[0] = 1;
									cord_in_ref[1] = 0;
									cord_in_ref[2] = -1;
								}
								if(ino == 1)
								{
									y -= 6;

									cord_in_ref[0] = 1;
									cord_in_ref[1] = 0;
									cord_in_ref[2] = 1;
								}
								if(ino == 2)
								{
									y -= 6;

									cord_in_ref[0] = -1;
									cord_in_ref[1] = 0;
									cord_in_ref[2] = 1;
								}
								if(ino == 3)
								{
									y -= 6;

									cord_in_ref[0] = -1;
									cord_in_ref[1] = 0;
									cord_in_ref[2] = -1;
								}
							}
							if(aux[0] == 4)
							{
								if(ino == 0)
								{
									z += 6;

									cord_in_ref[0] = 1;
									cord_in_ref[1] = -1;
									cord_in_ref[2] = 0;
								}
								if(ino == 1)
								{
									z += 6;

									cord_in_ref[0] = -1;
									cord_in_ref[1] = -1;
									cord_in_ref[2] = 0;
								}
								if(ino == 2)
								{
									z += 6;

									cord_in_ref[0] = -1;
									cord_in_ref[1] = 1;
									cord_in_ref[2] = 0;
								}
								if(ino == 3)
								{
									z += 6;

									cord_in_ref[0] = 1;
									cord_in_ref[1] = 1;
									cord_in_ref[2] = 0;
								}
							}
							if(aux[0] == 5)
							{
								if(ino == 0)
								{
									z -= 6;

									cord_in_ref[0] = 1;
									cord_in_ref[1] = -1;
									cord_in_ref[2] = 0;
								}
								if(ino == 1)
								{
									z -= 6;

									cord_in_ref[0] = -1;
									cord_in_ref[1] = -1;
									cord_in_ref[2] = 0;
								}
								if(ino == 2)
								{
									z -= 6;

									cord_in_ref[0] = -1;
									cord_in_ref[1] = 1;
									cord_in_ref[2] = 0;
								}
								if(ino == 3)
								{
									z -= 6;

									cord_in_ref[0] = 1;
									cord_in_ref[1] = 1;
									cord_in_ref[2] = 0;
								}
							}
						}
						else if(aux.size() == 2)
						{
							x = elem->nodes[FaceNodesMap[aux[isurf]][ino]].x;
							y = elem->nodes[FaceNodesMap[aux[isurf]][ino]].y;
							z = elem->nodes[FaceNodesMap[aux[isurf]][ino]].z;
							//edge 0
							if(aux[0] == 2 && aux[1] == 4)
							{
								if(aux[isurf] == 2)
								{
									if(ino == 0)
									{
										y += 6;
										z += 6;
										cord_in_ref[0] = 1;
										cord_in_ref[1] = 0;
										cord_in_ref[2] = 0;
									}
									if(ino == 1)
									{
										y += 6;
										cord_in_ref[0] = 1;
										cord_in_ref[1] = 0;
										cord_in_ref[2] = 1;
									}
									if(ino == 2)
									{
										y += 6;
										cord_in_ref[0] = -1;
										cord_in_ref[1] = 0;
										cord_in_ref[2] = 1;
									}
									if(ino == 3)
									{
										y += 6;
										z += 6;
										cord_in_ref[0] = -1;
										cord_in_ref[1] = 0;
										cord_in_ref[2] = 0;
									}
								}
								if(aux[isurf] == 4)
								{//1 0 3 2
									if(ino == 0)
									{
										y += 6;
										z += 6;
										cord_in_ref[0] = 1;
										cord_in_ref[1] = 0;
										cord_in_ref[2] = 0;
									}
									if(ino == 1)
									{
										y += 6;
										z += 6;
										cord_in_ref[0] = -1;
										cord_in_ref[1] = 0;
										cord_in_ref[2] = 0;
									}
									if(ino == 2)
									{
										z += 6;
										cord_in_ref[0] = -1;
										cord_in_ref[1] = 1;
										cord_in_ref[2] = 0;
									}
									if(ino == 3)
									{
										z += 6;
										cord_in_ref[0] = 1;
										cord_in_ref[1] = 1;
										cord_in_ref[2] = 0;
									}
								}
							}
							//edge 1
							if(aux[0] == 1 && aux[1] == 4)
							{
								if(aux[isurf] == 1)
								{
									if(ino == 0)
									{
										x -= 6;
										z += 6;
										cord_in_ref[0] = 0;
										cord_in_ref[1] = -1;
										cord_in_ref[2] = 0;
									}
									if(ino == 1)
									{
										x -= 6;
										cord_in_ref[0] = 0;
										cord_in_ref[1] = -1;
										cord_in_ref[2] = 1;
									}
									if(ino == 2)
									{
										x -= 6;
										cord_in_ref[0] = 0;
										cord_in_ref[1] = 1;
										cord_in_ref[2] = 1;
									}
									if(ino == 3)
									{
										x -= 6;
										z += 6;
										cord_in_ref[0] = 0;
										cord_in_ref[1] = 1;
										cord_in_ref[2] = 0;
									}
								}
								if(aux[isurf] == 4)
								{
									if(ino == 0)
									{
										x -= 6;
										z += 6;
										cord_in_ref[0] = 0;
										cord_in_ref[1] = -1;
										cord_in_ref[2] = 0;
									}
									if(ino == 1)
									{
										z += 6;
										cord_in_ref[0] = -1;
										cord_in_ref[1] = -1;
										cord_in_ref[2] = 0;
									}
									if(ino == 2)
									{
										z += 6;
										cord_in_ref[0] = -1;
										cord_in_ref[1] = 1;
										cord_in_ref[2] = 0;
									}
									if(ino == 3)
									{
										x -= 6;
										z += 6;
										cord_in_ref[0] = 0;
										cord_in_ref[1] = 1;
										cord_in_ref[2] = 0;
									}
								}
							}
							//edge 2
							if(aux[0] == 3 && aux[1] == 4)
							{
								if(aux[isurf] == 3)
								{
									if(ino == 0)
									{
										y -= 6;
										z += 6;
										cord_in_ref[0] = 1;
										cord_in_ref[1] = 0;
										cord_in_ref[2] = 0;
									}
									if(ino == 1)
									{
										y -= 6;
										cord_in_ref[0] = 1;
										cord_in_ref[1] = 0;
										cord_in_ref[2] = 1;
									}
									if(ino == 2)
									{
										y -= 6;
										cord_in_ref[0] = -1;
										cord_in_ref[1] = 0;
										cord_in_ref[2] = 1;
									}
									if(ino == 3)
									{
										y -= 6;
										z += 6;
										cord_in_ref[0] = -1;
										cord_in_ref[1] = 0;
										cord_in_ref[2] = 0;
									}
								}
								if(aux[isurf] == 4)
								{
									if(ino == 0)
									{
										z += 6;
										cord_in_ref[0] = 1;
										cord_in_ref[1] = -1;
										cord_in_ref[2] = 0;
									}
									if(ino == 1)
									{
										z += 6;
										cord_in_ref[0] = -1;
										cord_in_ref[1] = -1;
										cord_in_ref[2] = 0;
									}
									if(ino == 2)
									{
										y -= 6;
										z += 6;
										cord_in_ref[0] = 1;
										cord_in_ref[1] = 0;
										cord_in_ref[2] = 0;
									}
									if(ino == 3)
									{
										y -= 6;
										z += 6;
										cord_in_ref[0] = -1;
										cord_in_ref[1] = 0;
										cord_in_ref[2] = 0;
									}
								}
							}
							//edge 3
							if(aux[0] == 0 && aux[1] == 4)
							{
								if(aux[isurf] == 0)
								{//0,4,7,3
									if(ino == 0)
									{
										x += 6;
										z += 6;
										cord_in_ref[0] = 0;
										cord_in_ref[1] = -1;
										cord_in_ref[2] = 0;
									}
									if(ino == 1)
									{
										x += 6;
										cord_in_ref[0] = 0;
										cord_in_ref[1] = -1;
										cord_in_ref[2] = 1;
									}
									if(ino == 2)
									{
										x += 6;
										cord_in_ref[0] = 0;
										cord_in_ref[1] = 1;
										cord_in_ref[2] = 1;
									}
									if(ino == 3)
									{
										x += 6;
										z += 6;
										cord_in_ref[0] = 0;
										cord_in_ref[1] = 1;
										cord_in_ref[2] = 0;
									}
								}
								if(aux[isurf] == 4)
								{//1,0,3,2
									if(ino == 0)
									{
										z += 6;
										cord_in_ref[0] = 1;
										cord_in_ref[1] = -1;
										cord_in_ref[2] = 0;
									}

									if(ino == 1)
									{
										x += 6;
										z += 6;
										cord_in_ref[0] = 0;
										cord_in_ref[1] = -1;
										cord_in_ref[2] = 0;
									}
									if(ino == 2)
									{
										x += 6;
										z += 6;
										cord_in_ref[0] = 0;
										cord_in_ref[1] = 1;
										cord_in_ref[2] = 0;
									}

									if(ino == 3)
									{
										z += 6;
										cord_in_ref[0] = 1;
										cord_in_ref[1] = 1;
										cord_in_ref[2] = 0;
									}
								}
							}

							//edge 4
							if(aux[0] == 0 && aux[1] == 2)
							{
								if(aux[isurf] == 0)
								{
									if(ino == 0)
									{
										x += 6;
										y += 6;
										cord_in_ref[0] = 0;
										cord_in_ref[1] = 0;
										cord_in_ref[2] = -1;
									}
									if(ino == 1)
									{
										x += 6;
										y += 6;
										cord_in_ref[0] = 0;
										cord_in_ref[1] = 0;
										cord_in_ref[2] = 1;
									}
									if(ino == 2)
									{
										x += 6;
										cord_in_ref[0] = 0;
										cord_in_ref[1] = 1;
										cord_in_ref[2] = 1;
									}
									if(ino == 3)
									{
										x += 6;
										cord_in_ref[0] = 0;
										cord_in_ref[1] = 1;
										cord_in_ref[2] = -1;
									}
								}
								if(aux[isurf] == 2)
								{
									if(ino == 0)
									{
										y += 6;
										cord_in_ref[0] = 1;
										cord_in_ref[1] = 0;
										cord_in_ref[2] = -1;
									}
									if(ino == 1)
									{
										y += 6;
										cord_in_ref[0] = 0;
										cord_in_ref[1] = 1;
										cord_in_ref[2] = 1;
									}
									if(ino == 2)
									{
										x += 6;
										y += 6;
										cord_in_ref[0] = 0;
										cord_in_ref[1] = 0;
										cord_in_ref[2] = 1;
									}
									if(ino == 3)
									{
										x += 6;
										y += 6;
										cord_in_ref[0] = 0;
										cord_in_ref[1] = 0;
										cord_in_ref[2] = -1;
									}
								}
							}
							//edge 5
							if(aux[0] == 1 && aux[1] == 2)
							{
								if(aux[isurf] == 1)
								{//1 5 6 2
									if(ino == 0)
									{
										x -= 6;
										y += 6;
										cord_in_ref[0] = 0;
										cord_in_ref[1] = 0;
										cord_in_ref[2] = -1;
									}
									if(ino == 1)
									{
										x -= 6;
										y += 6;
										cord_in_ref[0] = 0;
										cord_in_ref[1] = 0;
										cord_in_ref[2] = 1;
									}
									if(ino == 2)
									{
										x -= 6;
										cord_in_ref[0] = 0;
										cord_in_ref[1] = 1;
										cord_in_ref[2] = 1;
									}
									if(ino == 3)
									{
										x -= 6;
										cord_in_ref[0] = 0;
										cord_in_ref[1] = 1;
										cord_in_ref[2] = -1;
									}
								}
								if(aux[isurf] == 2)
								{
									if(ino == 0)
									{
										y += 6;
										x -= 6;
										cord_in_ref[0] = 0;
										cord_in_ref[1] = 0;
										cord_in_ref[2] = -1;
									}
									if(ino == 1)
									{
										y += 6;
										x -= 6;
										cord_in_ref[0] = 0;
										cord_in_ref[1] = 0;
										cord_in_ref[2] = 1;
									}
									if(ino == 2)
									{
										y += 6;
										cord_in_ref[0] = -1;
										cord_in_ref[1] = 0;
										cord_in_ref[2] = 1;
									}
									if(ino == 3)
									{
										y += 6;
										cord_in_ref[0] = -1;
										cord_in_ref[1] = 0;
										cord_in_ref[2] = -1;
									}
								}
							}
							//edge 6
							if(aux[0] == 1 && aux[1] == 3)
							{
								if(aux[isurf] == 1)
								{//1 5 6 2
									if(ino == 0)
									{
										x -= 6;
										cord_in_ref[0] = 0;
										cord_in_ref[1] = -1;
										cord_in_ref[2] = -1;
									}
									if(ino == 1)
									{
										x -= 6;
										cord_in_ref[0] = 0;
										cord_in_ref[1] = -1;
										cord_in_ref[2] = 1;
									}
									if(ino == 2)
									{
										x -= 6;
										y -= 6;
										cord_in_ref[0] = 0;
										cord_in_ref[1] = 0;
										cord_in_ref[2] = 1;
									}
									if(ino == 3)
									{
										x -= 6;
										y -= 6;
										cord_in_ref[0] = 0;
										cord_in_ref[1] = 0;
										cord_in_ref[2] = -1;
									}
								}
								if(aux[isurf] == 3)
								{//2 6 7 3
									if(ino == 0)
									{
										x -= 6;
										y -= 6;
										cord_in_ref[0] = 0;
										cord_in_ref[1] = 0;
										cord_in_ref[2] = -1;
									}
									if(ino == 1)
									{
										x -= 6;
										y -= 6;
										cord_in_ref[0] = 0;
										cord_in_ref[1] = 0;
										cord_in_ref[2] = 1;
									}
									if(ino == 2)
									{
										y -= 6;
										cord_in_ref[0] = -1;
										cord_in_ref[1] = 0;
										cord_in_ref[2] = 1;
									}
									if(ino == 3)
									{
										y -= 6;
										cord_in_ref[0] = -1;
										cord_in_ref[1] = 0;
										cord_in_ref[2] = -1;
									}
								}
							}
							//edge 7
							if(aux[0] == 0 && aux[1] == 3)
							{
								if(aux[isurf] == 0)
								{//0 4 7 3
									if(ino == 0)
									{
										x += 6;
										cord_in_ref[0] = 0;
										cord_in_ref[1] = -1;
										cord_in_ref[2] = -1;
									}
									if(ino == 1)
									{
										x += 6;
										cord_in_ref[0] = 0;
										cord_in_ref[1] = -1;
										cord_in_ref[2] = 1;
									}
									if(ino == 2)
									{
										x += 6;
										y += 6;
										cord_in_ref[0] = 0;
										cord_in_ref[1] = 0;
										cord_in_ref[2] = 1;
									}
									if(ino == 3)
									{
										x += 6;
										y += 6;
										cord_in_ref[0] = 0;
										cord_in_ref[1] = 0;
										cord_in_ref[2] = -1;
									}
								}
								if(aux[isurf] == 3)
								{//2 6 7 3
									if(ino == 0)
									{
										y -= 6;
										cord_in_ref[0] = 0;
										cord_in_ref[1] = 1;
										cord_in_ref[2] = -1;
									}
									if(ino == 1)
									{
										y -= 6;
										cord_in_ref[0] = 0;
										cord_in_ref[1] = 0;
										cord_in_ref[2] = -1;
									}
									if(ino == 2)
									{
										y -= 6;
										x += 6;
										cord_in_ref[0] = 0;
										cord_in_ref[1] = 0;
										cord_in_ref[2] = 1;
									}
									if(ino == 3)
									{
										y -= 6;
										x += 6;
										cord_in_ref[0] = 0;
										cord_in_ref[1] = 0;
										cord_in_ref[2] = -1;
									}
								}
							}

							//edge 8
							if(aux[0] == 2 && aux[1] == 5)
							{
								if(aux[isurf] == 2)
								{//1 5 4 0
									if(ino == 0)
									{
										y += 6;
										cord_in_ref[0] = 1;
										cord_in_ref[1] = 0;
										cord_in_ref[2] = 1;
									}
									if(ino == 1)
									{
										y += 6;
										z -= 6;
										cord_in_ref[0] = 1;
										cord_in_ref[1] = 0;
										cord_in_ref[2] = 0;
									}
									if(ino == 2)
									{
										y += 6;
										z -= 6;
										cord_in_ref[0] = -1;
										cord_in_ref[1] = 0;
										cord_in_ref[2] = 0;
									}
									if(ino == 3)
									{
										y += 6;
										cord_in_ref[0] = -1;
										cord_in_ref[1] = 0;
										cord_in_ref[2] = 1;
									}
								}
								if(aux[isurf] == 5)
								{//5 4 7 6
									if(ino == 0)
									{
										y += 6;
										z -= 6;
										cord_in_ref[0] = 1;
										cord_in_ref[1] = 0;
										cord_in_ref[2] = 0;
									}
									if(ino == 1)
									{
										y += 6;
										z -= 6;
										cord_in_ref[0] = -1;
										cord_in_ref[1] = 0;
										cord_in_ref[2] = 0;
									}
									if(ino == 2)
									{
										z -= 6;
										cord_in_ref[0] = -1;
										cord_in_ref[1] = 1;
										cord_in_ref[2] = 0;
									}
									if(ino == 3)
									{
										z -= 6;
										cord_in_ref[0] = 1;
										cord_in_ref[1] = 1;
										cord_in_ref[2] = 0;
									}
								}
							}
							//edge 9
							if(aux[0] == 1 && aux[1] == 5)
							{
								if(aux[isurf] == 1)
								{//1,5,6,2
									if(ino == 0)
									{
										x -= 6;
										cord_in_ref[0] = 0;
										cord_in_ref[1] = -1;
										cord_in_ref[2] = -1;
									}
									if(ino == 1)
									{
										x -= 6;
										z -= 6;
										cord_in_ref[0] = 0;
										cord_in_ref[1] = -1;
										cord_in_ref[2] = 0;
									}
									if(ino == 2)
									{
										x -= 6;
										z -= 6;
										cord_in_ref[0] = 0;
										cord_in_ref[1] = 1;
										cord_in_ref[2] = 0;
									}
									if(ino == 3)
									{
										x -= 6;
										cord_in_ref[0] = 0;
										cord_in_ref[1] = 1;
										cord_in_ref[2] = 1;
									}
								}
								if(aux[isurf] == 5)
								{
									if(ino == 0)
									{
										x -= 6;
										z -= 6;
										cord_in_ref[0] = 0;
										cord_in_ref[1] = -1;
										cord_in_ref[2] = 0;
									}
									if(ino == 1)
									{
										z -= 6;
										cord_in_ref[0] = -1;
										cord_in_ref[1] = -1;
										cord_in_ref[2] = 0;
									}
									if(ino == 2)
									{
										z -= 6;
										cord_in_ref[0] = -1;
										cord_in_ref[1] = 1;
										cord_in_ref[2] = 0;
									}
									if(ino == 3)
									{
										x -= 6;
										z -= 6;
										cord_in_ref[0] = 0;
										cord_in_ref[1] = 1;
										cord_in_ref[2] = 0;
									}
								}
							}
							//edge 10
							if(aux[0] == 3 && aux[1] == 5)
							{
								if(aux[isurf] == 3)
								{
									if(ino == 0)
									{
										y -= 6;
										cord_in_ref[0] = 1;
										cord_in_ref[1] = 0;
										cord_in_ref[2] = -1;
									}
									if(ino == 1)
									{
										y -= 6;
										z -= 6;
										cord_in_ref[0] = 1;
										cord_in_ref[1] = 0;
										cord_in_ref[2] = 0;
									}
									if(ino == 2)
									{
										y -= 6;
										z -= 6;
										cord_in_ref[0] = -1;
										cord_in_ref[1] = 0;
										cord_in_ref[2] = 0;
									}
									if(ino == 3)
									{
										y -= 6;
										cord_in_ref[0] = -1;
										cord_in_ref[1] = 0;
										cord_in_ref[2] = -1;
									}
								}
								if(aux[isurf] == 5)
								{
									if(ino == 1)
									{
										z -= 6;
										cord_in_ref[0] = -1;
										cord_in_ref[1] = 1;
										cord_in_ref[2] = 0;
									}
									if(ino == 0)
									{
										z -= 6;
										cord_in_ref[0] = 1;
										cord_in_ref[1] = 1;
										cord_in_ref[2] = 0;
									}
									if(ino == 2)
									{
										y -= 6;
										z -= 6;
										cord_in_ref[0] = -1;
										cord_in_ref[1] = 0;
										cord_in_ref[2] = 0;
									}
									if(ino == 3)
									{
										y -= 6;
										z -= 6;
										cord_in_ref[0] = 1;
										cord_in_ref[1] = 0;
										cord_in_ref[2] = 0;
									}
								}
							}
							//edge 11
							if(aux[0] == 0 && aux[1] == 5)
							{
								if(aux[isurf] == 0)
								{//0,4,7,3
									if(ino == 1)
									{
										x += 6;
										z -= 6;
										cord_in_ref[0] = 0;
										cord_in_ref[1] = -1;
										cord_in_ref[2] = 0;
									}
									if(ino == 0)
									{
										x += 6;
										cord_in_ref[0] = 0;
										cord_in_ref[1] = -1;
										cord_in_ref[2] = 1;
									}
									if(ino == 3)
									{
										x += 6;
										cord_in_ref[0] = 0;
										cord_in_ref[1] = 1;
										cord_in_ref[2] = 1;
									}
									if(ino == 2)
									{
										x += 6;
										z -= 6;
										cord_in_ref[0] = 0;
										cord_in_ref[1] = 1;
										cord_in_ref[2] = 0;
									}
								}
								if(aux[isurf] == 5)
								{//1,0,3,2
									if(ino == 0)
									{
										z -= 6;
										cord_in_ref[0] = 1;
										cord_in_ref[1] = -1;
										cord_in_ref[2] = 0;
									}

									if(ino == 1)
									{
										x += 6;
										z -= 6;
										cord_in_ref[0] = 0;
										cord_in_ref[1] = -1;
										cord_in_ref[2] = 0;
									}
									if(ino == 2)
									{
										x += 6;
										z -= 6;
										cord_in_ref[0] = 0;
										cord_in_ref[1] = 1;
										cord_in_ref[2] = 0;
									}

									if(ino == 3)
									{
										z -= 6;
										cord_in_ref[0] = 1;
										cord_in_ref[1] = 1;
										cord_in_ref[2] = 0;
									}
								}
							}


						}
						else if(aux.size() == 3)
						{
							x = elem->nodes[FaceNodesMap[aux[isurf]][ino]].x;
							y = elem->nodes[FaceNodesMap[aux[isurf]][ino]].y;
							z = elem->nodes[FaceNodesMap[aux[isurf]][ino]].z;
							//corner 0 teste 8
							if(aux[0] == 1 && aux[1] == 3 && aux[2] == 5)
							{
								if(aux[isurf] == 1)
								{//1 5 6 2
									if(ino == 0)
									{
										x -= 6;
										cord_in_ref[0] = 0;
										cord_in_ref[1] = -1;
										cord_in_ref[2] = -1;
									}
									if(ino == 1)
									{
										x -= 6;
										z -= 6;
										cord_in_ref[0] = 0;
										cord_in_ref[1] = -1;
										cord_in_ref[2] = 0;
									}
									if(ino == 2)
									{
										x -= 6;
										y -= 6;
										z -= 6;
										cord_in_ref[0] = 0;
										cord_in_ref[1] = 0;
										cord_in_ref[2] = 0;
									}
									if(ino == 3)
									{
										x -= 6;
										y -= 6;
										cord_in_ref[0] = 0;
										cord_in_ref[1] = 0;
										cord_in_ref[2] = -1;
									}
								}
								if(aux[isurf] == 3)
								{//2 6 7 3
									if(ino == 0)
									{
										x -= 6;
										y -= 6;
										cord_in_ref[0] = 0;
										cord_in_ref[1] = 0;
										cord_in_ref[2] = -1;
									}
									if(ino == 1)
									{
										x -= 6;
										y -= 6;
										z -= 6;
										cord_in_ref[0] = 0;
										cord_in_ref[1] = 0;
										cord_in_ref[2] = 0;
									}
									if(ino == 2)
									{
										y -= 6;
										z -= 6;
										cord_in_ref[0] = -1;
										cord_in_ref[1] = 0;
										cord_in_ref[2] = 0;
									}
									if(ino == 3)
									{
										y -= 6;
										cord_in_ref[0] = -1;
										cord_in_ref[1] = 0;
										cord_in_ref[2] = -1;
									}
								}
								if(aux[isurf] == 5)
								{//5 4 7 6
									if(ino == 0)
									{
										x -= 6;
										z -= 6;
										cord_in_ref[0] = 0;
										cord_in_ref[1] = -1;
										cord_in_ref[2] = 0;
									}
									if(ino == 1)
									{
										z -= 6;
										cord_in_ref[0] = -1;
										cord_in_ref[1] = -1;
										cord_in_ref[2] = 0;
									}
									if(ino == 2)
									{
										y -= 6;
										z -= 6;
										cord_in_ref[0] = -1;
										cord_in_ref[1] = 0;
										cord_in_ref[2] = 0;
									}
									if(ino == 3)
									{
										x -= 6;
										y -= 6;
										z -= 6;
										cord_in_ref[0] = 0;
										cord_in_ref[1] = 0;
										cord_in_ref[2] = 0;
									}
								}
							}
							//corner 1
							if(aux[0] == 0 && aux[1] == 3 && aux[2] == 5 && false)
							{
								if(aux[isurf] == 0)
								{//1 5 6 2
									if(ino == 0)
									{
										x += 6;
										cord_in_ref[0] = 0;
										cord_in_ref[1] = -1;
										cord_in_ref[2] = -1;
									}
									if(ino == 1)
									{
										x += 6;
										z -= 6;
										cord_in_ref[0] = 0;
										cord_in_ref[1] = -1;
										cord_in_ref[2] = 0;
									}
									if(ino == 2)
									{
										x += 6;
										y -= 6;
										z -= 6;
										cord_in_ref[0] = 0;
										cord_in_ref[1] = 0;
										cord_in_ref[2] = 0;
									}
									if(ino == 3)
									{
										x += 6;
										y -= 6;
										cord_in_ref[0] = 0;
										cord_in_ref[1] = 0;
										cord_in_ref[2] = -1;
									}
								}
								if(aux[isurf] == 3)
								{//2 6 7 3
									if(ino == 0)
									{
										x += 6;
										y -= 6;
										cord_in_ref[0] = 0;
										cord_in_ref[1] = 0;
										cord_in_ref[2] = -1;
									}
									if(ino == 1)
									{
										x += 6;
										y -= 6;
										z -= 6;
										cord_in_ref[0] = 0;
										cord_in_ref[1] = 0;
										cord_in_ref[2] = 0;
									}
									if(ino == 2)
									{
										y -= 6;
										z -= 6;
										cord_in_ref[0] = 1;
										cord_in_ref[1] = 0;
										cord_in_ref[2] = 0;
									}
									if(ino == 3)
									{
										y -= 6;
										cord_in_ref[0] = 1;
										cord_in_ref[1] = 0;
										cord_in_ref[2] = -1;
									}
								}
								if(aux[isurf] == 5)
								{//5 4 7 6
									if(ino == 0)
									{
										x += 6;
										z -= 6;
										cord_in_ref[0] = 0;
										cord_in_ref[1] = -1;
										cord_in_ref[2] = 0;
									}
									if(ino == 1)
									{
										z -= 6;
										cord_in_ref[0] = 1;
										cord_in_ref[1] = -1;
										cord_in_ref[2] = 0;
									}
									if(ino == 2)
									{
										y -= 6;
										z -= 6;
										cord_in_ref[0] = 1;
										cord_in_ref[1] = 0;
										cord_in_ref[2] = 0;
									}
									if(ino == 3)
									{
										x += 6;
										y -= 6;
										z -= 6;
										cord_in_ref[0] = 0;
										cord_in_ref[1] = 0;
										cord_in_ref[2] = 0;
									}
								}
							}
							//corner 2 teste10
							if(aux[0] == 0 && aux[1] == 2 && aux[2] == 5)
							{
								if(aux[isurf] == 0)
								{//0,4,7,3
									if(ino == 0)
									{
										x += 6;
										y += 6;
										cord_in_ref[0] = 0;
										cord_in_ref[1] = 0;
										cord_in_ref[2] = -1;
									}
									if(ino == 1)
									{
										x += 6;
										y += 6;
										z -= 6;
										cord_in_ref[0] = 0;
										cord_in_ref[1] = 0;
										cord_in_ref[2] = 0;
									}
									if(ino == 2)
									{
										x += 6;
										z -= 6;
										cord_in_ref[0] = 0;
										cord_in_ref[1] = 1;
										cord_in_ref[2] = 0;
									}
									if(ino == 3)
									{
										x += 6;
										cord_in_ref[0] = 0;
										cord_in_ref[1] = 1;
										cord_in_ref[2] = -1;
									}
								}
								if(aux[isurf] == 2)
								{//2 6 7 3
									if(ino == 0)
									{
										y += 6;
										cord_in_ref[0] = 1;
										cord_in_ref[1] = 0;
										cord_in_ref[2] = -1;
									}
									if(ino == 1)
									{
										y += 6;
										z -= 6;
										cord_in_ref[0] = 1;
										cord_in_ref[1] = 0;
										cord_in_ref[2] = 0;
									}
									if(ino == 2)
									{
										x += 6;
										y += 6;
										z -= 6;
										cord_in_ref[0] = 0;
										cord_in_ref[1] = 0;
										cord_in_ref[2] = 0;
									}
									if(ino == 3)
									{
										x += 6;
										y += 6;
										cord_in_ref[0] = 0;
										cord_in_ref[1] = 0;
										cord_in_ref[2] = -1;
									}
								}
								if(aux[isurf] == 5)
								{//5 4 7 6
									if(ino == 0)
									{
										z -= 6;
										y += 6;
										cord_in_ref[0] = 1;
										cord_in_ref[1] = 0;
										cord_in_ref[2] = 0;
									}
									if(ino == 1)
									{
										x += 6;
										z -= 6;
										y += 6;
										cord_in_ref[0] = 0;
										cord_in_ref[1] = 0;
										cord_in_ref[2] = 0;
									}
									if(ino == 2)
									{
										x += 6;
										z -= 6;
										cord_in_ref[0] = 0;
										cord_in_ref[1] = 1;
										cord_in_ref[2] = 0;
									}
									if(ino == 3)
									{
										z -= 6;
										cord_in_ref[0] = 1;
										cord_in_ref[1] = 1;
										cord_in_ref[2] = 0;
									}
								}
							}
							//corner 3 teste9
							if(aux[0] == 1 && aux[1] == 2 && aux[2] == 5)
							{
								if(aux[isurf] == 1)
								{//1 5 6 2
									if(ino == 0)
									{
										x -= 6;
										y += 6;
										cord_in_ref[0] = 0;
										cord_in_ref[1] = 0;
										cord_in_ref[2] = -1;
									}
									if(ino == 1)
									{
										x -= 6;
										z -= 6;
										y += 6;
										cord_in_ref[0] = 0;
										cord_in_ref[1] = 0;
										cord_in_ref[2] = 0;
									}
									if(ino == 2)
									{
										x -= 6;
										z -= 6;
										cord_in_ref[0] = 0;
										cord_in_ref[1] = 1;
										cord_in_ref[2] = 0;
									}
									if(ino == 3)
									{
										x -= 6;
										cord_in_ref[0] = 0;
										cord_in_ref[1] = 1;
										cord_in_ref[2] = -1;
									}
								}
								if(aux[isurf] == 2)
								{//2 6 7 3
									if(ino == 0)
									{
										x -= 6;
										y += 6;
										cord_in_ref[0] = 0;
										cord_in_ref[1] = 0;
										cord_in_ref[2] = -1;
									}
									if(ino == 1)
									{
										x -= 6;
										y += 6;
										z -= 6;
										cord_in_ref[0] = 0;
										cord_in_ref[1] = 0;
										cord_in_ref[2] = 0;
									}
									if(ino == 2)
									{
										y += 6;
										z -= 6;
										cord_in_ref[0] = -1;
										cord_in_ref[1] = 0;
										cord_in_ref[2] = 0;
									}
									if(ino == 3)
									{
										y += 6;
										cord_in_ref[0] = -1;
										cord_in_ref[1] = 0;
										cord_in_ref[2] = -1;
									}
								}
								if(aux[isurf] == 5)
								{//5 4 7 6
									if(ino == 0)
									{
										x -= 6;
										z -= 6;
										y += 6;
										cord_in_ref[0] = 0;
										cord_in_ref[1] = 0;
										cord_in_ref[2] = 0;
									}
									if(ino == 1)
									{
										z -= 6;
										y += 6;
										cord_in_ref[0] = -1;
										cord_in_ref[1] = 0;
										cord_in_ref[2] = 0;
									}
									if(ino == 2)
									{
										z -= 6;
										cord_in_ref[0] = -1;
										cord_in_ref[1] = 1;
										cord_in_ref[2] = 0;
									}
									if(ino == 3)
									{
										x -= 6;
										z -= 6;
										cord_in_ref[0] = 0;
										cord_in_ref[1] = 1;
										cord_in_ref[2] = 0;
									}
								}
							}

							//corner 4
							if(aux[0] == 1 && aux[1] == 3 && aux[2] == 4 && false)
							{
								if(aux[isurf] == 1)
								{//1 5 6 2
									if(ino == 0)
									{
										x -= 6;
										cord_in_ref[0] = 0;
										cord_in_ref[1] = -1;
										cord_in_ref[2] = 1;
									}
									if(ino == 1)
									{
										x -= 6;
										z += 6;
										cord_in_ref[0] = 0;
										cord_in_ref[1] = -1;
										cord_in_ref[2] = 0;
									}
									if(ino == 2)
									{
										x -= 6;
										y -= 6;
										z += 6;
										cord_in_ref[0] = 0;
										cord_in_ref[1] = 0;
										cord_in_ref[2] = 0;
									}
									if(ino == 3)
									{
										x -= 6;
										y -= 6;
										cord_in_ref[0] = 0;
										cord_in_ref[1] = 0;
										cord_in_ref[2] = 1;
									}
								}
								if(aux[isurf] == 3)
								{//2 6 7 3
									if(ino == 0)
									{
										x -= 6;
										y -= 6;
										cord_in_ref[0] = 0;
										cord_in_ref[1] = 0;
										cord_in_ref[2] = 1;
									}
									if(ino == 1)
									{
										x -= 6;
										y -= 6;
										z += 6;
										cord_in_ref[0] = 0;
										cord_in_ref[1] = 0;
										cord_in_ref[2] = 0;
									}
									if(ino == 2)
									{
										y -= 6;
										z += 6;
										cord_in_ref[0] = -1;
										cord_in_ref[1] = 0;
										cord_in_ref[2] = 0;
									}
									if(ino == 3)
									{
										y -= 6;
										cord_in_ref[0] = -1;
										cord_in_ref[1] = 0;
										cord_in_ref[2] = 1;
									}
								}
								if(aux[isurf] == 4)
								{//5 4 7 6
									if(ino == 0)
									{
										x -= 6;
										z += 6;
										cord_in_ref[0] = 0;
										cord_in_ref[1] = -1;
										cord_in_ref[2] = 0;
									}
									if(ino == 1)
									{
										z += 6;
										cord_in_ref[0] = -1;
										cord_in_ref[1] = -1;
										cord_in_ref[2] = 0;
									}
									if(ino == 2)
									{
										y -= 6;
										z += 6;
										cord_in_ref[0] = -1;
										cord_in_ref[1] = 0;
										cord_in_ref[2] = 0;
									}
									if(ino == 3)
									{
										x -= 6;
										y -= 6;
										z += 6;
										cord_in_ref[0] = 0;
										cord_in_ref[1] = 0;
										cord_in_ref[2] = 0;
									}
								}
							}
							//corner 5
							if(aux[0] == 0 && aux[1] == 3 && aux[2] == 4 && false)
							{
								if(aux[isurf] == 0)
								{//1 5 6 2
									if(ino == 0)
									{
										x += 6;
										cord_in_ref[0] = 0;
										cord_in_ref[1] = -1;
										cord_in_ref[2] = 1;
									}
									if(ino == 1)
									{
										x += 6;
										z += 6;
										cord_in_ref[0] = 0;
										cord_in_ref[1] = -1;
										cord_in_ref[2] = 0;
									}
									if(ino == 2)
									{
										x += 6;
										y -= 6;
										z += 6;
										cord_in_ref[0] = 0;
										cord_in_ref[1] = 0;
										cord_in_ref[2] = 0;
									}
									if(ino == 3)
									{
										x += 6;
										y -= 6;
										cord_in_ref[0] = 0;
										cord_in_ref[1] = 0;
										cord_in_ref[2] = 1;
									}
								}
								if(aux[isurf] == 3)
								{//2 6 7 3
									if(ino == 0)
									{
										x += 6;
										y -= 6;
										cord_in_ref[0] = 0;
										cord_in_ref[1] = 0;
										cord_in_ref[2] = 1;
									}
									if(ino == 1)
									{
										x += 6;
										y -= 6;
										z += 6;
										cord_in_ref[0] = 0;
										cord_in_ref[1] = 0;
										cord_in_ref[2] = 0;
									}
									if(ino == 2)
									{
										y -= 6;
										z += 6;
										cord_in_ref[0] = 1;
										cord_in_ref[1] = 0;
										cord_in_ref[2] = 0;
									}
									if(ino == 3)
									{
										y -= 6;
										cord_in_ref[0] = 1;
										cord_in_ref[1] = 0;
										cord_in_ref[2] = 1;
									}
								}
								if(aux[isurf] == 4)
								{//5 4 7 6
									if(ino == 0)
									{
										x += 6;
										z += 6;
										cord_in_ref[0] = 0;
										cord_in_ref[1] = -1;
										cord_in_ref[2] = 0;
									}
									if(ino == 1)
									{
										z += 6;
										cord_in_ref[0] = 1;
										cord_in_ref[1] = -1;
										cord_in_ref[2] = 0;
									}
									if(ino == 2)
									{
										y -= 6;
										z += 6;
										cord_in_ref[0] = 1;
										cord_in_ref[1] = 0;
										cord_in_ref[2] = 0;
									}
									if(ino == 3)
									{
										x += 6;
										y -= 6;
										z += 6;
										cord_in_ref[0] = 0;
										cord_in_ref[1] = 0;
										cord_in_ref[2] = 0;
									}
								}
							}
							//corner 6
							if(aux[0] == 0 && aux[1] == 2 && aux[2] == 4 && false)
							{
								if(aux[isurf] == 0)
								{//1 5 6 2
									if(ino == 0)
									{
										x += 6;
										cord_in_ref[0] = 0;
										cord_in_ref[1] = 1;
										cord_in_ref[2] = 1;
									}
									if(ino == 1)
									{
										x += 6;
										z += 6;
										cord_in_ref[0] = 0;
										cord_in_ref[1] = 1;
										cord_in_ref[2] = 0;
									}
									if(ino == 2)
									{
										x += 6;
										y += 6;
										z += 6;
										cord_in_ref[0] = 0;
										cord_in_ref[1] = 0;
										cord_in_ref[2] = 0;
									}
									if(ino == 3)
									{
										x += 6;
										y += 6;
										cord_in_ref[0] = 0;
										cord_in_ref[1] = 0;
										cord_in_ref[2] = 1;
									}
								}
								if(aux[isurf] == 2)
								{//2 6 7 3
									if(ino == 0)
									{
										x += 6;
										y += 6;
										cord_in_ref[0] = 0;
										cord_in_ref[1] = 0;
										cord_in_ref[2] = 1;
									}
									if(ino == 1)
									{
										x += 6;
										y += 6;
										z += 6;
										cord_in_ref[0] = 0;
										cord_in_ref[1] = 0;
										cord_in_ref[2] = 0;
									}
									if(ino == 2)
									{
										y += 6;
										z += 6;
										cord_in_ref[0] = 1;
										cord_in_ref[1] = 0;
										cord_in_ref[2] = 0;
									}
									if(ino == 3)
									{
										y += 6;
										cord_in_ref[0] = 1;
										cord_in_ref[1] = 0;
										cord_in_ref[2] = 1;
									}
								}
								if(aux[isurf] == 4)
								{//5 4 7 6
									if(ino == 0)
									{
										x += 6;
										z += 6;
										cord_in_ref[0] = 0;
										cord_in_ref[1] = 1;
										cord_in_ref[2] = 0;
									}
									if(ino == 1)
									{
										z += 6;
										cord_in_ref[0] = 1;
										cord_in_ref[1] = 1;
										cord_in_ref[2] = 0;
									}
									if(ino == 2)
									{
										y += 6;
										z += 6;
										cord_in_ref[0] = 1;
										cord_in_ref[1] = 0;
										cord_in_ref[2] = 0;
									}
									if(ino == 3)
									{
										x += 6;
										y += 6;
										z += 6;
										cord_in_ref[0] = 0;
										cord_in_ref[1] = 0;
										cord_in_ref[2] = 0;
									}
								}
							}
							//corner 7
							if(aux[0] == 1 && aux[1] == 2 && aux[2] == 4 && false)
							{
								if(aux[isurf] == 1)
								{//1 5 6 2
									if(ino == 0)
									{
										x -= 6;
										cord_in_ref[0] = 0;
										cord_in_ref[1] = 1;
										cord_in_ref[2] = 1;
									}
									if(ino == 1)
									{
										x -= 6;
										z += 6;
										cord_in_ref[0] = 0;
										cord_in_ref[1] = 1;
										cord_in_ref[2] = 0;
									}
									if(ino == 2)
									{
										x -= 6;
										y += 6;
										z += 6;
										cord_in_ref[0] = 0;
										cord_in_ref[1] = 0;
										cord_in_ref[2] = 0;
									}
									if(ino == 3)
									{
										x -= 6;
										y += 6;
										cord_in_ref[0] = 0;
										cord_in_ref[1] = 0;
										cord_in_ref[2] = 1;
									}
								}
								if(aux[isurf] == 2)
								{//2 6 7 3
									if(ino == 0)
									{
										x -= 6;
										y += 6;
										cord_in_ref[0] = 0;
										cord_in_ref[1] = 0;
										cord_in_ref[2] = 1;
									}
									if(ino == 1)
									{
										x -= 6;
										y += 6;
										z += 6;
										cord_in_ref[0] = 0;
										cord_in_ref[1] = 0;
										cord_in_ref[2] = 0;
									}
									if(ino == 2)
									{
										y += 6;
										z += 6;
										cord_in_ref[0] = -1;
										cord_in_ref[1] = 0;
										cord_in_ref[2] = 0;
									}
									if(ino == 3)
									{
										y += 6;
										cord_in_ref[0] = -1;
										cord_in_ref[1] = 0;
										cord_in_ref[2] = 1;
									}
								}
								if(aux[isurf] == 4)
								{//5 4 7 6
									if(ino == 0)
									{
										x -= 6;
										z += 6;
										cord_in_ref[0] = 0;
										cord_in_ref[1] = 1;
										cord_in_ref[2] = 0;
									}
									if(ino == 1)
									{
										z += 6;
										cord_in_ref[0] = -1;
										cord_in_ref[1] = 1;
										cord_in_ref[2] = 0;
									}
									if(ino == 2)
									{
										y += 6;
										z += 6;
										cord_in_ref[0] = -1;
										cord_in_ref[1] = 0;
										cord_in_ref[2] = 0;
									}
									if(ino == 3)
									{
										x -= 6;
										y += 6;
										z += 6;
										cord_in_ref[0] = 0;
										cord_in_ref[1] = 0;
										cord_in_ref[2] = 0;
									}
								}
							}
						}

						//changer le noeud si on besion
						//nao preciso me preocupar com cord_in_ref
						//teoricamente o no ja existe, entao o ref_point
						//nao importa
						if(!exel && true)
						{
							aux_teste_t key;
							size_t position;
							key.orig.x = elem->nodes[FaceNodesMap[aux[isurf]][ino]].x;
							key.orig.y = elem->nodes[FaceNodesMap[aux[isurf]][ino]].y;
							key.orig.z = elem->nodes[FaceNodesMap[aux[isurf]][ino]].z;

							bool lnode = sc_hash_array_lookup(pauxhash, &key, &position);

							if(lnode)
							{
								aux_teste_t* p = (aux_teste_t*) sc_array_index(&pauxhash->a,position);
								x = p->pill.x;
								y = p->pill.y;
								z = p->pill.z;
								new_id[ino] = p->pill.id;

								cord_in_ref[0] = 0;
								cord_in_ref[1] = 0;
								cord_in_ref[2] = 0;
							}
						}

						//ajouter sur la hash de noeuds
						size_t position;
						octant_node_t key;
						key.x = x;
						key.y = y;
						key.z = z;

						GtsPoint* ref_point = LinearMapHex(cord_in_ref, ref_in_x, ref_in_y, ref_in_z);
						int nodeID;
						octant_node_t* r = (octant_node_t*) sc_hash_array_insert_unique(hash_nodes, &key, &position);
						if(r!=NULL){
							r->x = x;
							r->y = y;
							r->z = z;
							r->id = hash_nodes->a.elem_count-1;
							nodeID = r->id;
							coords.push_back(ref_point->x);
							coords.push_back(ref_point->y);
							coords.push_back(ref_point->z);
						}else{
							octant_node_t* r = (octant_node_t*) sc_array_index(&hash_nodes->a,position);
							nodeID = r->id;
						}
						new_id[ino] = nodeID;
					}

					//fazendo a criacao do novo elemento
					//alterando a conectividade do elemento base

					//créer l'élément
					octant_t* pelem;
					pelem = (octant_t*) sc_array_push(&mesh->elements);
					pelem->id = mesh->elements.elem_count-1;
					pelem->n_mat = elem->n_mat;
					pelem->pad = elem->pad;
					pelem->x = elem->x;
					pelem->y = elem->y;
					pelem->z = elem->z;

					//printf("  Face:%d\n",aux[isurf]);
					for(int ino = 0; ino < 4; ino++)
					{
						pelem->nodes[FaceNodesMap[aux[isurf]][ino]].id = connec[FaceNodesMap[aux[isurf]][ino]];
						pelem->nodes[FaceNodesMap[aux[isurf]][ino]].fixed = elemcp->nodes[FaceNodesMap[aux[isurf]][ino]].fixed;
						pelem->nodes[FaceNodesMap_inv[aux[isurf]][ino]].id= new_id[ino];

						elem->nodes[FaceNodesMap[aux[isurf]][ino]].id = new_id[ino];

						pelem->nodes[FaceNodesMap[aux[isurf]][ino]].fixed= 0;
						octant_node_t* node = (octant_node_t*) sc_array_index(&hash_nodes->a,connec[FaceNodesMap[aux[isurf]][ino]]);
						pelem->nodes[FaceNodesMap[aux[isurf]][ino]].x = node->x;
						pelem->nodes[FaceNodesMap[aux[isurf]][ino]].y = node->y;
						pelem->nodes[FaceNodesMap[aux[isurf]][ino]].z = node->z;
						node = (octant_node_t*) sc_array_index(&hash_nodes->a,new_id[ino]);
						pelem->nodes[FaceNodesMap_inv[aux[isurf]][ino]].x = node->x;
						pelem->nodes[FaceNodesMap_inv[aux[isurf]][ino]].y = node->y;
						pelem->nodes[FaceNodesMap_inv[aux[isurf]][ino]].z = node->z;
					}

				}
			}

			sc_array_reset(&toto);
		}
		sc_hash_array_destroy(pillow);
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
	//Pillowing(mesh, coords,nodes_b_mat);
	Pillowing_oct(mesh, coords,nodes_b_mat);
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
