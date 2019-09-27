#include <stdlib.h>
#include <assert.h>
#include <sc.h>
#include <vector>
#include <algorithm>

#include "hexa.h"

void communicate_global_ids(hexa_tree_t* mesh);
void communicate_global_edge_ids(hexa_tree_t* mesh);

unsigned processors_hash_fn (const void *v, const void *u)
{
	const uint32_t *a_ptr = (uint32_t*) v;
	uint32_t        a = *a_ptr;
	uint32_t        b=1, c=1;

	sc_hash_mix (a, b, c);
	sc_hash_final (a, b, c);
	return (unsigned) c;
}

int processors_equal_fn (const void *v1, const void *v2, const void *u)
{
	const uint32_t *q1 = (const uint32_t*) v1;
	const uint32_t *q2 = (const uint32_t*) v2;
	return (*q1 == *q2);
}

int node_comp (const void *v, const void *u)
{

	const octant_node_t *q = (const octant_node_t*) v;
	const octant_node_t *p = (const octant_node_t*) u;

	if(q->x > p->x){
		return 1;
	}else if(q->x == p->x){
		if(q->y > p->y){
			return 1;
		}else if(q->y == p->y){
			if(q->z > p->z){
				return 1;
			}else if(q->y == p->y){
				return 1;
			}else{
				return -1;
			}
		}else{
			return -1;
		}

	}else{
		return -1;
	}
}

uint64_t cantor_function_id(uint32_t c, uint32_t d)
{
	uint64_t e;
	uint32_t a;
	uint32_t b;

	if(c>d){
		a = c;
		b = d;
	}else{
		a=d;
		b=c;
	}

	if(a>=b){
		return e = a*a+a+b;
	}else{
		return e = a+b*b;
	}
}

unsigned edge_hash_id(const void *v, const void *u)
{
	const octant_edge_t *q = (const octant_edge_t*) v;
	uint64_t a, b, c;

	a = (uint32_t) 1;
	b = (uint32_t) q->coord[0];
	c = (uint32_t) q->coord[1];
	sc_hash_mix(a, b, c);
	sc_hash_final(a, b, c);
	return (unsigned) c;
}

int edge_equal_id(const void *v, const void *u, const void *w)
{
	const octant_edge_t *e1 = (const octant_edge_t*) v;
	const octant_edge_t *e2 = (const octant_edge_t*) u;

	return (unsigned) (e1->coord[0] == e2->coord[0] &&
			e1->coord[1] == e2->coord[1]);

}

void hexa_insert_shared_edge(sc_hash_array_t *shared_edges, shared_edge_t* edge, int processor)
{

	size_t position;
	shared_edge_t *sn;
	int i;

	if( processor < 0) return;

	sn = (shared_edge_t*) sc_hash_array_insert_unique (shared_edges, edge, &position);
	if(sn != NULL){
		sn->id = edge->id;
		sn->listSz = 1;
		sn->rankList[0] = processor;
	} else{
		sn = (shared_edge_t*) sc_array_index(&shared_edges->a, position);
		for(i=0; i < sn->listSz; ++i)
			if(sn->rankList[i] == processor) break;
		if(i == sn->listSz){
			sn->rankList[sn->listSz] = processor;
			sn->listSz++;
		}
	}
}

unsigned node_hash_fn (const void *v, const void *u)
{
	const octant_node_t *q = (const octant_node_t*) v;
	uint32_t            a, b, c;

	a = (uint32_t) q->x;
	b = (uint32_t) q->y;
	c = (uint32_t) q->z;
	sc_hash_mix (a, b, c);
	sc_hash_final (a, b, c);
	return (unsigned) c;
}

int node_equal_fn (const void *v1, const void *v2, const void *u)
{
	const octant_node_t *q1 = (const octant_node_t*) v1;
	const octant_node_t *q2 = (const octant_node_t*) v2;
	return (q1->x == q2->x && q1->y == q2->y && q1->z == q2->z);
}

void hexa_insert_shared_node(sc_hash_array_t    *shared_nodes, octant_node_t* node, int processor)
{
	size_t position;
	shared_node_t *sn;
	shared_node_t key;
	int i;

	key.id = node->id;
	key.x = node->x;
	key.y = node->y;
	key.z = node->z;


	if( processor < 0) return;

	sn = (shared_node_t*) sc_hash_array_insert_unique (shared_nodes, &key, &position);
	if(sn != NULL)
	{
		sn->x  = node->x;
		sn->y  = node->y;
		sn->z  = node->z;
		sn->id = node->id;
		sn->listSz = 1;
		sn->rankList[0] = processor;
	} else
	{
		sn = (shared_node_t*) sc_array_index(&shared_nodes->a, position);
		for(i=0; i < sn->listSz; ++i)
			if(sn->rankList[i] == processor) break;
		if(i == sn->listSz){
			sn->rankList[sn->listSz] = processor;
			sn->listSz++;
		}
	}
}

void hexaNodes(hexa_tree_t* mesh,sc_hash_array_t* indep_nodes,sc_hash_array_t* shared_nodes,int64_t *local,int64_t *global)
{

	bool clamped = true;
	int num_indep_nodes = 0;
	//insert internal nodes in the hash_array
	for(int iel = 0; iel < mesh->elements.elem_count; iel++)
	{
		octant_t *h  = (octant_t*) sc_array_index(&mesh->elements, iel);
		h->id = iel;
		for(int ino = 0; ino < 8; ino++)
		{
			size_t position;
			octant_node_t* node = &h->nodes[ino];
			octant_node_t *r = (octant_node_t*) sc_hash_array_insert_unique (indep_nodes, node, &position);
			if(r != NULL)
			{
				r->x     = node->x;
				//if(node->x == (mesh->ncellx+1) ) r->x = -mesh->max_step;
				if(node->x == (mesh->ncellx+2) ) r->x = mesh->ncellx+mesh->max_step;

				r->y     = node->y;
				//if(node->y == (mesh->ncelly+1) ) r->y = -mesh->max_step;
				if(node->y == (mesh->ncelly+2) ) r->y = mesh->ncelly+mesh->max_step;

				r->z     = node->z;
				//if(node->z == (mesh->ncellz+1) ) r->z = -mesh->max_step;
				if(node->z == (mesh->ncellz+2) ) r->z = mesh->ncellz+mesh->max_step;

				r->id    = num_indep_nodes;
				node->id = num_indep_nodes;
				assert(position == num_indep_nodes);
				++num_indep_nodes;
			} else
			{
				node->id = position;
			}
		}
	}
	//extract the nodes from indep_nodes
	sc_hash_array_rip (indep_nodes,  &mesh->nodes);

	//insert the shared nodes in the hash_array
	for(int ino = 0; ino < mesh->nodes.elem_count; ino++)
	{
		octant_node_t* node = (octant_node_t*) sc_array_index (&mesh->nodes, ino);
		if(node->y == mesh->y_start)
		{
			if(node->x == mesh->x_start)
				hexa_insert_shared_node(shared_nodes,node,mesh->neighbors[0]);
			if(node->x == mesh->x_end)
				hexa_insert_shared_node(shared_nodes,node,mesh->neighbors[2]);

			hexa_insert_shared_node(shared_nodes,node,mesh->neighbors[1]);
			//continue;
		}

		if(node->y == mesh->y_end) {
			if(node->x == mesh->x_start)
				hexa_insert_shared_node(shared_nodes,node,mesh->neighbors[6]);
			if(node->x == mesh->x_end)
				hexa_insert_shared_node(shared_nodes,node,mesh->neighbors[8]);

			hexa_insert_shared_node(shared_nodes,node,mesh->neighbors[7]);
			//continue;
		}

		if(node->x == mesh->x_start)
			hexa_insert_shared_node(shared_nodes,node,mesh->neighbors[3]);
		if(node->x == mesh->x_end)
			hexa_insert_shared_node(shared_nodes,node,mesh->neighbors[5]);
	}

	local[0] = mesh->local_n_nodes    = mesh->nodes.elem_count;
	local[1] = mesh->local_n_elements = mesh->elements.elem_count;

#ifdef HEXA_DEBUG_
	fprintf(mesh->fdbg, "Nodes: \n");
	fprintf(mesh->fdbg, "Total: %d\n",mesh->nodes.elem_count);
	for(int i = 0; i < mesh->nodes.elem_count; ++i)
	{
		octant_node_t* n = (octant_node_t*) sc_array_index(&mesh->nodes,i);
		fprintf(mesh->fdbg, "(%d): %d %d %d\n", n->id, n->x, n->y, n->z);
	}

#endif

	// node map
	int not_my_nodes    = 0;
	int my_own_nodes    = 0;
	mesh->global_id     = (int64_t*)malloc(sizeof(int64_t)*mesh->local_n_nodes);
	memset(mesh->global_id,-2,mesh->local_n_nodes*sizeof(int64_t));
	sc_hash_array_t* SendTo   = (sc_hash_array_t *) sc_hash_array_new(sizeof(message_t), processors_hash_fn, processors_equal_fn, &clamped);
	sc_hash_array_t* RecvFrom = (sc_hash_array_t *) sc_hash_array_new(sizeof(message_t), processors_hash_fn, processors_equal_fn, &clamped);

	for(int ino = 0; ino < shared_nodes->a.elem_count; ++ino)
	{
		size_t position;
		shared_node_t* sn = (shared_node_t*) sc_array_index(&shared_nodes->a,ino);
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

	communicate_global_ids(mesh);

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


#ifdef HEXA_DEBUG_
	fprintf(mesh->fdbg, "\nOffset: %d\n", offset);

	fprintf(mesh->fdbg, "Nodes (local) (proc) (global): \n");
	for(int i = 0; i < mesh->nodes.elem_count; ++i)
	{
		octant_node_t* n = (octant_node_t*) sc_array_index(&mesh->nodes,i);
		fprintf(mesh->fdbg, "(%d) (%d) (%ld): %d %d %d\n", n->id, mesh->part_nodes[i], mesh->global_id[i], n->x, n->y, n->z);
	}
#endif   

	/////////////////
	local[0] = mesh->local_n_nodes    = mesh->nodes.elem_count;
	local[1] = mesh->local_n_elements = mesh->elements.elem_count;

	MPI_Allreduce(local, global,2,MPI_LONG_LONG_INT,MPI_SUM,MPI_COMM_WORLD);

	mesh->total_n_nodes    = global[0];
	mesh->total_n_elements = global[1];

}

void hexaEdges(hexa_tree_t* mesh,sc_hash_array_t* indep_nodes,sc_hash_array_t* shared_nodes,int64_t *local,int64_t *global,sc_hash_array_t* indep_edges,sc_hash_array_t* shared_edges)
{
	bool clamped = true;
	int num_indep_edges = 0;

	for (int iel = 0; iel < mesh->elements.elem_count; ++iel)
	{

		octant_t *elem = (octant_t*) sc_array_index(&mesh->elements, iel);

		for (int edge = 0; edge < 12; ++edge)
		{
			int Edge2GNode[2];

			//look for the global_node_id to create the edge_id
			int node1 = mesh->global_id[elem->nodes[EdgeVerticesMap[edge][0]].id];
			int node2 = mesh->global_id[elem->nodes[EdgeVerticesMap[edge][1]].id];
			assert(node1 >= 0);
			assert(node2 >= 0);

			Edge2GNode[0] = node1 <= node2 ? node1 : node2;
			Edge2GNode[1] = node1 >= node2 ? node1 : node2;

			//Cantor function to avoid MPI comm...
			elem->edge[edge].id = cantor_function_id(node1,node2);
			elem->edge[edge].ref = false;
			elem->edge[edge].coord[0] = Edge2GNode[0];
			elem->edge[edge].coord[1] = Edge2GNode[1];

			size_t position;
			octant_edge_t *r;
			octant_edge_t key;
			key.id = elem->edge[edge].id;
			key.coord[0] = node1;
			key.coord[1] = node2;

			r = (octant_edge_t*) sc_hash_array_insert_unique(indep_edges, &key, &position);
			if (r != NULL) {
				r->coord[0] = key.coord[0];
				r->coord[1] = key.coord[1];
				r->id = elem->edge[edge].id;
				r->ref = false;
			}
		}
	}
#ifdef HEXA_DEBUG_
	if(0){
		for (int iel = 0; iel < mesh->elements.elem_count; ++iel)
		{
			octant_t *elem = (octant_t*) sc_array_index(&mesh->elements, iel);
			for (int edge = 0; edge < 12; ++edge)
			{
				fprintf(mesh->fdbg,"El:%lld, edge:%lld, coord[0]:%lld, coord[1]:%lld\n",iel,elem->edge[edge].id,elem->edge[edge].coord[0],elem->edge[edge].coord[1]);
			}
		}
	}
#endif

	// create the shared edges hash
	for (int iel = 0; iel < mesh->elements.elem_count; ++iel)
	{

		octant_t *elem = (octant_t*) sc_array_index(&mesh->elements, iel);
		size_t position0;
		size_t position1;

		for (int edge = 0; edge < 12; ++edge)
		{

			bool out0 = false;
			bool out1 = false;
			octant_node_t* node0 = (octant_node_t*) sc_array_index (&mesh->nodes, elem->nodes[EdgeVerticesMap[edge][0]].id);
			octant_node_t* node1 = (octant_node_t*) sc_array_index (&mesh->nodes, elem->nodes[EdgeVerticesMap[edge][1]].id);

			out0 =  sc_hash_array_lookup(shared_nodes, node0, &position0);
			out1 =  sc_hash_array_lookup(shared_nodes, node1, &position1);

			if(out0 && out1){
				shared_node_t* sn0 = (shared_node_t*) sc_array_index(&shared_nodes->a, position0);
				shared_node_t* sn1 = (shared_node_t*) sc_array_index(&shared_nodes->a, position1);

				for(int i = 0; i < sn0->listSz; i++)
				{
					for(int j = 0; j < sn1->listSz; j++)
					{
						if(sn0->rankList[i] == sn1->rankList[j])
						{
							shared_edge_t      se;
							se.id = elem->edge[edge].id;
							hexa_insert_shared_edge(shared_edges, &se , sn0->rankList[i]);
						}
					}
				}
			}
		}
	}

#ifdef HEXA_DEBUG_
	if(0){
		fprintf(mesh->fdbg,"Shared nodes in global ids:\n");
		fprintf(mesh->fdbg,"Total:%d\n",mesh->shared_nodes.elem_count);
		for(int i = 0; i < mesh->shared_nodes.elem_count; ++i)
		{
			shared_node_t* sn = (shared_node_t*) sc_array_index(&mesh->shared_nodes,i);
			fprintf(mesh->fdbg, "(%ld): %d %d %d\n", sn->id, sn->x, sn->y, sn->z);
			fprintf(mesh->fdbg, "     shared with processors: ");
			for(int j = 0; j < sn->listSz; j++)
			{
				fprintf(mesh->fdbg, "%d ", sn->rankList[j]);
			}
			fprintf(mesh->fdbg, "\n");
		}
	}
#endif
	////////////////////////////////////
	//extract the edges from indep_edges
	sc_hash_array_rip (indep_edges,  &mesh->edges);
	//extract the shared edges from shared_edges
	sc_hash_array_rip (shared_edges, &mesh->shared_edges);
	//sc_array_sort(&mesh->shared_edges,edge_comp);

#ifdef HEXA_DEBUG_
	if (0){
		fprintf(mesh->fdbg,"Edge ids:\n");
		for(int i = 0; i < mesh->edges.elem_count; ++i)
		{
			octant_edge_t* sn = (octant_edge_t*) sc_array_index(&mesh->edges,i);
			fprintf(mesh->fdbg, "Edge id:%lld status:%d\n", sn->id, sn->ref);
		}

		fprintf(mesh->fdbg,"Shared Edge ids:\n");
		for(int i = 0; i < mesh->shared_edges.elem_count; ++i)
		{
			shared_edge_t* sn = (shared_edge_t*) sc_array_index(&mesh->shared_edges,i);
			fprintf(mesh->fdbg, "Edge id:%lld\n", sn->id);
		}
	}
#endif

	//edge map
	int not_my_edges    = 0;
	int my_own_edges    = 0;
	mesh->global_edge_id     = (int64_t*)malloc(sizeof(int64_t)*mesh->edges.elem_count);
	memset(mesh->global_edge_id,-2,mesh->edges.elem_count*sizeof(int64_t));
	sc_hash_array_t* SendTo   = (sc_hash_array_t *) sc_hash_array_new(sizeof(message_t), processors_hash_fn, processors_equal_fn, &clamped);
	sc_hash_array_t* RecvFrom = (sc_hash_array_t *) sc_hash_array_new(sizeof(message_t), processors_hash_fn, processors_equal_fn, &clamped);

	for(int i = 0; i < mesh->shared_edges.elem_count; ++i)
	{
		shared_edge_t* sn = (shared_edge_t*) sc_array_index(&mesh->shared_edges,i);
		size_t position;
		for(int j = 0; j < sn->listSz; j++)
		{
			if(sn->rankList[j] < mesh->mpi_rank)
			{
				message_t* m = (message_t*)sc_hash_array_insert_unique(SendTo,&sn->rankList[j],&position);
				if(m!=NULL)
				{
					m->rank  = sn->rankList[j];
					sc_array_init(&m->idxs, sizeof(uint64_t));
					uint64_t* p = (uint64_t*) sc_array_push(&m->idxs);
					*p = sn->id;
				}
				else
				{
					message_t* m = (message_t*)sc_array_index(&SendTo->a, position);
					uint64_t* p = (uint64_t*) sc_array_push(&m->idxs);
					*p = sn->id;
				}
				my_own_edges++;
				//       mesh->global_edge_id[sn->id] = -3;
			}
			else if (sn->rankList[j] > mesh->mpi_rank)
			{
				message_t *m = (message_t*)sc_hash_array_insert_unique(RecvFrom,&sn->rankList[j],&position);
				if(m!=NULL)
				{
					m->rank  = sn->rankList[j];
					sc_array_init(&m->idxs, sizeof(uint64_t));
					uint64_t* p = (uint64_t*) sc_array_push(&m->idxs);
					*p = sn->id;
				}
				else
				{
					message_t* m = (message_t*)sc_array_index(&RecvFrom->a, position);
					uint64_t* p = (uint64_t*) sc_array_push(&m->idxs);
					*p = sn->id;
				}
				not_my_edges++;
				//mesh->global_edge_id[sn->id] = -1;
			}
		}
		// if(mesh->global_edge_id[sn->id] == -1) not_my_edges++;
		// if(mesh->global_edge_id[sn->id] == -3) my_own_edges++;
	}

	local[2] = mesh->local_n_edges = mesh->edges.elem_count;
	local[2] -= not_my_edges;

	sc_hash_array_rip(RecvFrom, &mesh->comm_map_edge.RecvFrom );
	sc_hash_array_rip(SendTo  , &mesh->comm_map_edge.SendTo   );

#ifdef HEXA_DEBUG_
	if(0){
		fprintf(mesh->fdbg, "Shared Edges: \n");
		fprintf(mesh->fdbg, "Total: %d\n",mesh->shared_edges.elem_count);
		for(int i = 0; i < mesh->shared_edges.elem_count; ++i)
		{
			shared_edge_t* sn = (shared_edge_t*) sc_array_index(&mesh->shared_edges,i);
			fprintf(mesh->fdbg, "%lld\n", sn->id);
			fprintf(mesh->fdbg, "     shared with processors: ");
			for(int j = 0; j < sn->listSz; j++)
			{
				fprintf(mesh->fdbg, "%d ", sn->rankList[j]);
			}
			fprintf(mesh->fdbg, "\n");
		}
		fprintf(mesh->fdbg,"Edges:\n");
		fprintf(mesh->fdbg, "Recv from %ld processors\n", mesh->comm_map_edge.RecvFrom.elem_count);
		for(int i=0; i < mesh->comm_map_edge.RecvFrom.elem_count; i++)
		{
			message_t* m = (message_t*) sc_array_index(&mesh->comm_map_edge.RecvFrom, i);
			fprintf(mesh->fdbg, "  \n Recv %ld edges from %d\n", m->idxs.elem_count, m->rank);
			for(int j =0; j < m->idxs.elem_count; j++)
			{
				int64_t *id = (int64_t*) sc_array_index(&m->idxs,j);
				fprintf(mesh->fdbg, "%lld ", *id);
				if( (j+1) % 5 == 0 )fprintf(mesh->fdbg, "\n");
			}
		}
		fprintf(mesh->fdbg,"\n");
		fprintf(mesh->fdbg, "Send to %ld processors\n", mesh->comm_map_edge.SendTo.elem_count);
		for(int i=0; i < mesh->comm_map_edge.SendTo.elem_count; i++)
		{
			message_t* m = (message_t*) sc_array_index(&mesh->comm_map_edge.SendTo, i);
			fprintf(mesh->fdbg, "\n Sending %ld edges from %d\n", m->idxs.elem_count, m->rank);
			for(int j =0; j < m->idxs.elem_count; j++)
			{
				int64_t *id = (int64_t*) sc_array_index(&m->idxs,j);
				fprintf(mesh->fdbg, "%lld ", *id);
				if( (j+1) % 5 == 0 )fprintf(mesh->fdbg, "\n");
			}
		}
	}
#endif

	//size for the edges message
	mesh->comm_map_edge.max_recvbuf_size = 0;
	for(int i = 0; i < mesh->comm_map_edge.RecvFrom.elem_count; i++)
	{
		message_t* m = (message_t*) sc_array_index(&mesh->comm_map_edge.RecvFrom, i);
		mesh->comm_map_edge.max_recvbuf_size +=  m->idxs.elem_count;
	}

	mesh->comm_map_edge.max_sendbuf_size = 0;
	for(int i = 0; i < mesh->comm_map_edge.SendTo.elem_count; i++)
	{
		message_t* m = (message_t*) sc_array_index(&mesh->comm_map_edge.SendTo, i);
		mesh->comm_map_edge.max_sendbuf_size +=  m->idxs.elem_count;
	}

	mesh->comm_map_edge.nrequests = mesh->comm_map_edge.RecvFrom.elem_count +
			mesh->comm_map_edge.SendTo.elem_count;

	for(int i=0; i < mesh->edges.elem_count; ++i){
		octant_edge_t * edge = (octant_edge_t*) sc_array_index(&mesh->edges,i);
		mesh->global_edge_id[i] = edge->id;
	}

	/*
	    //find the global edge_id
	    offset = 0;
	    MPI_Scan(&my_own_edges, &offset, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
	    count = offset;
	    for(int i=0; i < mesh->edges.elem_count; ++i)
	        if(mesh->global_edge_id[i] != -1) mesh->global_edge_id[i] = count++;
	    communicate_global_edge_ids(mesh);
	    //redo the edge_comm_map
	    for(int ied = 0; ied < mesh->shared_edges.elem_count; ++ied){
	        shared_edge_t *ed = (shared_edge_t*) sc_array_index(&mesh->shared_edges, ied);
	        int id = mesh->global_edge_id[ed->id];
	        ed->id= id;
	    }
	    SendTo   = (sc_hash_array_t *) sc_hash_array_new(sizeof(message_edge_t), processors_hash_fn, processors_equal_fn, &clamped);
	    RecvFrom = (sc_hash_array_t *) sc_hash_array_new(sizeof(message_edge_t), processors_hash_fn, processors_equal_fn, &clamped);
	    for(int i = 0; i < mesh->shared_edges.elem_count; ++i){
	        shared_edge_t* sn = (shared_edge_t*) sc_array_index(&mesh->shared_edges,i);
	        for(int j = 0; j < sn->listSz; j++){
	            if(sn->rankList[j] < mesh->mpi_rank) {
	                message_edge_t* m = (message_edge_t*)sc_hash_array_insert_unique(SendTo,&sn->rankList[j],&position);
	                if(m!=NULL){
	                    m->rank  = sn->rankList[j];
	                    sc_array_init(&m->idxs, sizeof(uint32_t));
	                    sc_array_init(&m->ref, sizeof(uint8_t));
	                    uint32_t* p = (uint32_t*) sc_array_push(&m->idxs);
	 *p = sn->id;
	                }else{
	                    message_edge_t* m = (message_edge_t*)sc_array_index(&SendTo->a, position);
	                    uint32_t* p = (uint32_t*) sc_array_push(&m->idxs);
	 *p = sn->id;
	                }
	            }
	            else if (sn->rankList[j] > mesh->mpi_rank){
	                message_edge_t *m = (message_edge_t*)sc_hash_array_insert_unique(RecvFrom,&sn->rankList[j],&position);
	                if(m!=NULL){
	                    m->rank  = sn->rankList[j];
	                    sc_array_init(&m->idxs, sizeof(uint32_t));
	                    uint32_t* p = (uint32_t*) sc_array_push(&m->idxs);
	 *p = sn->id;
	                }else{
	                    message_edge_t* m = (message_edge_t*)sc_array_index(&RecvFrom->a, position);
	                    uint32_t* p = (uint32_t*) sc_array_push(&m->idxs);
	 *p = sn->id;
	                }
	            }
	        }
	    }
	    sc_hash_array_rip(RecvFrom, &mesh->comm_map_edge.RecvFrom );
	    sc_hash_array_rip(SendTo  , &mesh->comm_map_edge.SendTo   );
	#ifdef HEXA_DEBUG_
	    fprintf(mesh->fdbg, "Offset: %d\n", offset);
	    fprintf(mesh->fdbg, "Edges: \n");
	    fprintf(mesh->fdbg, "Edges number:%d \n",mesh->edges.elem_count);
	    for(int i = 0; i < mesh->edges.elem_count; ++i){
	        octant_edge_t* n = (octant_edge_t*) sc_array_index(&mesh->edges,i);
	        fprintf(mesh->fdbg, "id:%d  global:%ld\n", n->id, mesh->global_edge_id[i]);
	    }
	    fprintf(mesh->fdbg, "Shared Edges in global: \n");
	    fprintf(mesh->fdbg, "Total: %d\n",mesh->shared_edges.elem_count);
	    for(int i = 0; i < mesh->shared_edges.elem_count; ++i){
	        shared_edge_t* sn = (shared_edge_t*) sc_array_index(&mesh->shared_edges,i);
	        fprintf(mesh->fdbg, "%d\n", sn->id);
	        fprintf(mesh->fdbg, "     shared with processors: ");
	            for(int j = 0; j < sn->listSz; j++){
	                fprintf(mesh->fdbg, "%d ", sn->rankList[j]);
	            }
	            fprintf(mesh->fdbg, "\n");
	    }
	    fprintf(mesh->fdbg,"Edges:\n");
	    fprintf(mesh->fdbg, "Recv from %ld processors\n", mesh->comm_map_edge.RecvFrom.elem_count);
	    for(int i=0; i < mesh->comm_map_edge.RecvFrom.elem_count; i++){
	        message_t* m = (message_t*) sc_array_index(&mesh->comm_map_edge.RecvFrom, i);
	        fprintf(mesh->fdbg, "  \n Recv %ld edges from %d\n", m->idxs.elem_count, m->rank);
	        for(int j =0; j < m->idxs.elem_count; j++){
	            int32_t *id = (int32_t*) sc_array_index(&m->idxs,j);
	            fprintf(mesh->fdbg, "%d ", *id);
	            if( (j+1) % 5 == 0 )fprintf(mesh->fdbg, "\n");
	        }
	    }
	    fprintf(mesh->fdbg,"\n");
	    fprintf(mesh->fdbg, "Send to %ld processors\n", mesh->comm_map_edge.SendTo.elem_count);
	    for(int i=0; i < mesh->comm_map_edge.SendTo.elem_count; i++){
	        message_t* m = (message_t*) sc_array_index(&mesh->comm_map_edge.SendTo, i);
	        fprintf(mesh->fdbg, "\n Sending %ld edges from %d\n", m->idxs.elem_count, m->rank);
	        for(int j =0; j < m->idxs.elem_count; j++){
	            int32_t *id = (int32_t*) sc_array_index(&m->idxs,j);
	            fprintf(mesh->fdbg, "%d ", *id);
	            if( (j+1) % 5 == 0 )fprintf(mesh->fdbg, "\n");
	        }
	    }
	#endif
	    //update edge_id
	    for(int i = 0; i < mesh->edges.elem_count; ++i){
	        octant_edge_t* n = (octant_edge_t*) sc_array_index(&mesh->edges,i);
	        n->id=mesh->global_edge_id[i];
	    }
	    //fprintf(mesh->fdbg,"local to global in the edges\n");
	    for(int iel = 0; iel < mesh->elements.elem_count;iel++){
	        octant_t *elem = (octant_t*) sc_array_index(&mesh->elements, iel);
	        for (int edge = 0; edge < 12; ++edge) {
	            int Edge2GNode[2];
	            int id = mesh->global_edge_id[elem->edge[edge].id];
	            //fprintf(mesh->fdbg, "id:%d  global:%ld\n", elem->edge[edge].id, id);
	            //elem->edge[edge].id = id;
	        }
	    }
	 */
}

void communicate_global_ids(hexa_tree_t* mesh){
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

	int offset = 0;

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

}

void communicate_global_edge_ids(hexa_tree_t* mesh){
	int          n_requests;

	MPI_Request *requests;
	MPI_Status  *statuses;
	long long    *recvbuf;
	long long   *sendbuf;

	n_requests = mesh->comm_map_edge.nrequests;
	recvbuf    = (long long*)malloc(mesh->comm_map_edge.max_recvbuf_size*sizeof(long long));
	sendbuf    = (long long*)malloc(mesh->comm_map_edge.max_sendbuf_size*sizeof(long long));

	requests = (MPI_Request*) malloc (n_requests*sizeof(MPI_Request));
	statuses = (MPI_Status*)  malloc (n_requests*sizeof(MPI_Status));
	int c = 0;

	int offset = 0;

	// post all non-blocking receives
	for(int i = 0; i < mesh->comm_map_edge.RecvFrom.elem_count; ++i)
	{
		message_t *m = (message_t*) sc_array_index(&mesh->comm_map_edge.RecvFrom, i);
		MPI_Irecv(&recvbuf[offset], m->idxs.elem_count, MPI_LONG_LONG, m->rank,0,MPI_COMM_WORLD, &requests[c]);
		offset += m->idxs.elem_count;
		c++;
	}

	assert(offset == mesh->comm_map_edge.max_recvbuf_size);

	offset = 0;
	for(int i = 0; i < mesh->comm_map_edge.SendTo.elem_count; ++i) {
		message_t *m = (message_t*) sc_array_index(&mesh->comm_map_edge.SendTo, i);
		for(int j = 0; j < m->idxs.elem_count; ++j)
		{
			int32_t *id = (int32_t*) sc_array_index(&m->idxs,j);
			sendbuf[offset+j] = (long long) mesh->global_edge_id[*id];
		}
		MPI_Isend(&sendbuf[offset], m->idxs.elem_count, MPI_LONG_LONG, m->rank,0,MPI_COMM_WORLD, &requests[c]);
		offset += m->idxs.elem_count;
		c++;
	}
	assert(offset == mesh->comm_map_edge.max_sendbuf_size);

	assert(c == n_requests);

	MPI_Waitall(n_requests,requests,statuses);

	offset = 0;
	for(int i = 0; i < mesh->comm_map_edge.RecvFrom.elem_count; ++i) {
		message_t *m = (message_t*) sc_array_index(&mesh->comm_map_edge.RecvFrom, i);
		for(int j = 0; j < m->idxs.elem_count; ++j)
		{
			int32_t *id = (int32_t*) sc_array_index(&m->idxs,j);
			mesh->global_edge_id[*id] = recvbuf[offset+j];
		}
		offset += m->idxs.elem_count;
	}

	free(&recvbuf[0]);
	free(&sendbuf[0]);
	free(requests);
	free(statuses);
}

void hexa_mesh(hexa_tree_t* mesh)
{

	bool                clamped = true;
	int64_t    local[3], global[3];

	sc_hash_array_t    *indep_nodes     = (sc_hash_array_t *)sc_hash_array_new(sizeof (octant_node_t), node_hash_fn, node_equal_fn, &clamped);
	sc_hash_array_t    *shared_nodes    = (sc_hash_array_t *)sc_hash_array_new(sizeof (shared_node_t), node_hash_fn, node_equal_fn, &clamped);

	hexaNodes(mesh,indep_nodes,shared_nodes,local,global);

	/*
	sc_hash_array_t    *indep_edges     = (sc_hash_array_t *)sc_hash_array_new(sizeof (octant_edge_t), edge_hash_id, edge_equal_id, &clamped);
	sc_hash_array_t    *shared_edges    = (sc_hash_array_t *)sc_hash_array_new(sizeof (shared_edge_t), edge_hash_id, edge_equal_id, &clamped);

	hexaEdges(mesh,indep_nodes,shared_nodes,local,global,indep_edges,shared_edges);
*/
	////////////////////////////////////
	//extract the share nodes from shared_nodes
	sc_hash_array_rip (shared_nodes, &mesh->shared_nodes);
	sc_array_sort(&mesh->shared_nodes,node_comp);

	//hexaElements(mesh,indep_nodes,shared_nodes,local,global);

	if(mesh->mpi_rank == 0)
	{
		printf("Total number of elements: %lld\n", mesh->total_n_elements);
		printf("Total number of nodes: %lld\n", mesh->total_n_nodes);
	}


}
