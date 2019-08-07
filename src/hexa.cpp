
#include <math.h>
#include "hexa.h"
#include "pml.h"

#define MAX(a,b) (a > b?a:b)


// Compute node id based on cartesian coordinates.
inline int get_node_id(int nx, int ny, int i, int j, int k)
{
	return (k*((nx+1)*(ny+1))+j*(nx+1)+i+1);
}


inline int get_hexa_id(int nx, int ny, int i, int j, int k)
{
	return (k*((nx)*(ny))+j*(nx)+i+1);
}


void copy_octant(octant_t *orig, octant_t* dest)
{
	memcpy(dest, orig, sizeof(octant_t));
}


void hexa_tree_init(hexa_tree_t* mesh, int max_levels)
{
	sc_array_init(&mesh->elements, sizeof(octant_t));
	mesh->ncellx = 2*(int32_t) pow(3, max_levels);
	mesh->ncelly = mesh->ncellx;
	mesh->ncellz = mesh->ncellx;
	mesh->max_levels = max_levels;
}

void hexa_element_init(octant_t *elem)
{
	elem->level = 0;
	elem->x     = elem->y = elem->z;
	elem->pad   = 0;
	elem->n_mat = 0;
	elem->father = -1;
	elem->tem = -1;
	elem->boundary = false;
	elem->pml_id= PML_NULL;
}

//octant_t * hexa_element_copy(hexa_tree_t* mesh,int iel)
void hexa_element_copy(octant_t* elem, octant_t* elemCopy)
{
	//octant_t * elemCopy;
	//octant_t* elem = (octant_t*) sc_array_index(&mesh->elements, iel);

	elemCopy->level = elem->level;
	elemCopy->x = elem->x;
	elemCopy->y = elem->y;
	elemCopy->z = elem->z;
	elemCopy->pad = elem->pad;
	elemCopy->pml_id = elem->pml_id;
	elemCopy->n_mat = elem->n_mat;
	elemCopy->id = elem->id;
	elemCopy->boundary = elem->boundary;
	elemCopy->father = elem->id;

	for(int ino = 0; ino <8; ino++){
		elemCopy->nodes[ino].x = elem->nodes[ino].x;
		elemCopy->nodes[ino].y = elem->nodes[ino].y;
		elemCopy->nodes[ino].z = elem->nodes[ino].z;
		elemCopy->nodes[ino].fixed = elem->nodes[ino].fixed;
		elemCopy->nodes[ino].id = elem->nodes[ino].id;
		elemCopy->nodes[ino].color = elem->nodes[ino].color;
	}

	for(int iedge = 0; iedge<12; iedge++){
		elemCopy->edge[iedge].ref = 	elem->edge[iedge].ref;
		elemCopy->edge[iedge].id = elem->edge[iedge].id;
		elemCopy->edge[iedge].coord[0] = elem->edge[iedge].coord[0];
		elemCopy->edge[iedge].coord[1] = elem->edge[iedge].coord[1];
	}

	for(int isurf = 0; isurf < 6; isurf++){
		elemCopy->surf[isurf].ext = elem->surf[isurf].ext;
	}
	//return elemCopy;
}

void hexa_element_conn(octant_t* h, int i, int j, int k, int step,  int level)
{
	octant_node_t *node;
	h->x = i;
	h->y = j;
	h->z = k;
	h->level = level;

	//node 1
	node = &h->nodes[0];
	node->x = i;
	node->y = j;
	node->z = k;

	//node 2
	node = &h->nodes[1];
	node->x = i+step;
	node->y = j;
	node->z = k;

	//node 3
	node = &h->nodes[2];
	node->x = i+step;
	node->y = j+step;
	node->z = k;

	//node 4
	node = &h->nodes[3];
	node->x = i      ;
	node->y = j+step;
	node->z = k;

	//node 5
	node = &h->nodes[4];
	node->x = i;
	node->y = j;
	node->z = k+step;

	//node 6
	node = &h->nodes[5];
	node->x = i+step;
	node->y = j;
	node->z = k+step;

	//node 7
	node = &h->nodes[6];
	node->x = i+step;
	node->y = j+step;
	node->z = k+step;

	//node 8
	node = &h->nodes[7];
	node->x = i     ;
	node->y = j+step;
	node->z = k+step;

}

void hexa_get_27tree(hexa_tree_t* mesh, int nx, int ny, int nz, int step, int level)
{
	//int x,y,z;
	for(int k=0, z=nz; k < 3; k++, z+=level)
		for(int j=0, y=ny; j < 3;j++, y+=level)
			for(int i=0,x=nx ; i < 3; i++, x+=level)
			{
				octant_t * elem = (octant_t*) sc_array_push(&mesh->elements);
				elem->id=mesh->elements.elem_count-1;
				hexa_element_init(elem);
				hexa_element_conn(elem,x,y,z, step, level);
			}
}


void hexa_refinement_layer(hexa_tree_t* mesh, int nz, int coarse_step, int internal_step, int level)
{
	for(int ny=0; ny < mesh->ncelly; ny+=coarse_step) {
		for(int nx=0; nx < mesh->ncellx; nx+=coarse_step)
		{
			hexa_get_27tree(mesh,nx,ny,nz,internal_step,level);
			mesh->max_step = MAX(mesh->max_step,internal_step);
		}
	}

}


void hexa_uniform_layer(hexa_tree_t* mesh, int nz, int coarse_step, int internal_step, int level, bool nz_ext)
{
	for(int ny=mesh->y_start; ny < mesh->y_end; ny +=coarse_step) {
		for(int nx=mesh->x_start; nx < mesh->x_end; nx+=coarse_step)
		{
			octant_t * elem = (octant_t*) sc_array_push(&mesh->elements);
			elem->id = mesh->elements.elem_count;
			hexa_element_init(elem);
			hexa_element_conn(elem,nx,ny,nz, coarse_step, level);
			if(ny == mesh->y_start || !((ny + coarse_step) < mesh->y_end) ||
					nx == mesh->x_start || !((nx + coarse_step) < mesh->x_end) || nz_ext) elem->boundary = true;
			mesh->max_step = MAX(mesh->max_step,coarse_step);
			mesh->max_z = nz+mesh->max_step;
		}
	}
}


void hexa_transient_layer(hexa_tree_t* mesh, int nz, int coarse_step, int internal_step, int level)
{
	for(int ny=mesh->y_start; ny < mesh->y_end; ny +=coarse_step) {
		for(int nx=mesh->x_start; nx < mesh->x_end; nx+=coarse_step)
		{
			int ext = 0;
			if(ny == mesh->y_start) ext = 1;
			if(!((ny + coarse_step) < mesh->y_end)) ext = 2;
			if(nx == mesh->x_start) ext = 3;
			if(!((nx + coarse_step) < mesh->x_end)) ext = 4;
			if(ny == mesh->y_start && nx == mesh->x_start) ext = 5;
			if(ny == mesh->y_start && !((nx + coarse_step) < mesh->x_end)) ext = 6;
			if(!((nx + coarse_step) < mesh->x_end) && !((ny + coarse_step) < mesh->y_end)) ext = 7;
			if(!((ny + coarse_step) < mesh->y_end) && nx == mesh->x_start) ext = 8;
			hexa_transition_element(mesh,nx,ny,nz,internal_step,level,ext);
			mesh->max_step = MAX(mesh->max_step,internal_step);
		}
	}
}


void hexa_tree_cube(hexa_tree_t* mesh)
{

	int32_t nx, ny, nz, coarse_step, internal_step;



	if(mesh->mpi_rank == 0)
	{
		printf(" Max cells \n");
		printf("   x-direction: %d\n", mesh->ncellx);
		printf("   y-direction: %d\n", mesh->ncelly);
		printf("   z-direction: %d\n", mesh->ncellz);
	}

	internal_step = 1;
	coarse_step   = 1;
	int level     = 1;
	int layer     = 0;
	int nlayer    = 0;
	mesh->max_step = 0;

	hexa_processors_interval(mesh);
	//TODO
	//preciso achar aqui o numero para dividir esse negocio... assim eu consigo ajustar o numero de camadas e tal...
	nz = 0;
	//int nz_test = nz+internal_step;
	while( (nz+internal_step) <= mesh->ncellz)
	{
		if((nlayer+1)%60 == 0) {
			coarse_step*=3;
			hexa_transient_layer(mesh,nz,coarse_step, internal_step, level);
			internal_step*=3;
		} else {
			bool nz_ext= false;
			if(nz == 0 || !((nz+2*internal_step) <= mesh->ncellz)) nz_ext = true;
			hexa_uniform_layer(mesh,nz,coarse_step, internal_step,level,nz_ext);
		}

		nz+=internal_step;
		nlayer++;
	}
}

void hexa_tree_destroy(hexa_tree_t* mesh)
{

	if(mesh->global_id!=NULL) free(mesh->global_id);
	if(mesh->part_nodes!=NULL) free(mesh->part_nodes);
	sc_array_reset(&mesh->elements);
	sc_array_reset(&mesh->comm_map.RecvFrom);
	sc_array_reset(&mesh->comm_map.SendTo);
	sc_array_reset(&mesh->nodes);
	//sc_array_reset(&mesh->vertex);
	sc_array_reset(&mesh->oct);
	sc_array_reset(&mesh->shared_nodes);

	if(mesh->gdata.bbt!=NULL) gts_bb_tree_destroy(mesh->gdata.bbt, TRUE);

}
