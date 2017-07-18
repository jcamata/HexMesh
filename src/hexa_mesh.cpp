#include <stdlib.h>
#include <assert.h>
#include <sc.h>

#include "hexa.h"
#include "refinement.h"
void communicate_global_ids(hexa_tree_t* mesh);
void communicate_global_edge_ids(hexa_tree_t* mesh);

int node_comp (const void *v, const void *u){
  
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

int edge_sort (const void *v, const void *u){
  
    const shared_edge_t *q = (const shared_edge_t*) v;
    const shared_edge_t *p = (const shared_edge_t*) u;
    
      if(q->id < p->id){
          return 1;
      }
}

int edge_comp (const void *v, const void *u){
  
    const shared_edge_t *q = (const shared_edge_t*) v;
    const shared_edge_t *p = (const shared_edge_t*) u;
    
      if(q->id == p->id){
          return 1;
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

void hexa_insert_shared_node(sc_hash_array_t    *shared_nodes, octant_node_t* node, int processor)
{
    size_t position;
    shared_node_t *sn;
    int i;
    
    if( processor < 0) return;
    
    sn = (shared_node_t*) sc_hash_array_insert_unique (shared_nodes, node, &position);
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

unsigned edge_hash_f(const void *v, const void *u) {
	const octant_edge_t *q = (const octant_edge_t*) v;
	uint32_t a, b, c;

	a = (uint32_t) q->coord[0];
	b = (uint32_t) q->coord[1];
	c = (uint32_t) 1;
	sc_hash_mix(a, b, c);
	sc_hash_final(a, b, c);
	return (unsigned) c;
}

int edge_equal_f(const void *v, const void *u, const void *w) {
	const octant_edge_t *e1 = (const octant_edge_t*) v;
	const octant_edge_t *e2 = (const octant_edge_t*) u;

	return (unsigned) ((e1->coord[0] == e2->coord[0]) && (e1->coord[1] == e2->coord[1]));

}

unsigned edge_hash_id_f(const void *v, const void *u) {
	const octant_edge_t *q = (const octant_edge_t*) v;
	uint32_t a, b, c;

	a = (uint32_t) q->id;
	b = (uint32_t) 1;
	c = (uint32_t) 1;
	sc_hash_mix(a, b, c);
	sc_hash_final(a, b, c);
	return (unsigned) c;
}

int edge_equal_id_f(const void *v, const void *u, const void *w) {
	const octant_edge_t *e1 = (const octant_edge_t*) v;
	const octant_edge_t *e2 = (const octant_edge_t*) u;

	return (unsigned) (e1->id == e2->id);

}

unsigned int edge_id(int* nodes, sc_hash_array_t* hash, int &npoints) {
	size_t position;
	octant_edge_t *r;
	octant_edge_t key;
	key.coord[0] = nodes[0];
	key.coord[1] = nodes[1];

	r = (octant_edge_t*) sc_hash_array_insert_unique(hash, &key, &position);
	if (r != NULL) {
		r->coord[0] = key.coord[0];
		r->coord[1] = key.coord[1];
		r->id = npoints;
		npoints++;
		return r->id;
	} else {
		r = (octant_edge_t*) sc_array_index(&hash->a, position);
		return r->id;
	}
}

void hexa_insert_shared_edge(sc_hash_array_t *shared_edges, shared_edge_t* edge, int processor){
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

void hexa_mesh(hexa_tree_t* mesh){
    
    bool                clamped = true;
    sc_hash_array_t    *indep_nodes;
    sc_hash_array_t    *indep_edges;
    sc_hash_array_t    *shared_nodes;
    sc_hash_array_t    *shared_edges;
    sc_hash_array_t    *SendTo;
    sc_hash_array_t    *RecvFrom;
    size_t              position;
    octant_node_t      *r;
    int num_indep_nodes = 0;
    int64_t    local[3], global[3];

    indep_nodes     = (sc_hash_array_t *)sc_hash_array_new(sizeof (octant_node_t), node_hash_fn, node_equal_fn, &clamped);
    indep_edges     = (sc_hash_array_t *)sc_hash_array_new(sizeof (octant_edge_t), edge_hash_f, edge_equal_f, &clamped);
    shared_nodes    = (sc_hash_array_t *)sc_hash_array_new(sizeof (shared_node_t), node_hash_fn, node_equal_fn, &clamped);
    shared_edges    = (sc_hash_array_t *)sc_hash_array_new(sizeof (shared_edge_t), edge_hash_id_f, edge_equal_id_f, &clamped);

    //insert internal nodes in the hash_array
    for(int i = 0; i < mesh->elements.elem_count; i++)
    {
        octant_t *h  = (octant_t*) sc_array_index(&mesh->elements, i);
        for(int j = 0; j < 8; j++) {
            octant_node_t* node = &h->nodes[j];
            //initialization for the node color
            node->color = -1;
            r = (octant_node_t*) sc_hash_array_insert_unique (indep_nodes, node, &position);
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

    
    //extract the share nodes from shared_nodes
    sc_hash_array_rip (shared_nodes, &mesh->shared_nodes);
    sc_array_sort(&mesh->shared_nodes,node_comp);
    
    local[0] = mesh->local_n_nodes    = mesh->nodes.elem_count;
    local[1] = mesh->local_n_elements = mesh->elements.elem_count;
    
    
#ifdef HEXA_DEBUG_
    fprintf(mesh->fdbg, "Nodes: \n");
    for(int i = 0; i < mesh->nodes.elem_count; ++i)
    {
        octant_node_t* n = (octant_node_t*) sc_array_index(&mesh->nodes,i);
        fprintf(mesh->fdbg, "(%d): %d %d %d\n", n->id, n->x, n->y, n->z);
    }

    fprintf(mesh->fdbg, "Shared Nodes: \n");
    fprintf(mesh->fdbg, "Total: %d\n",mesh->shared_nodes.elem_count);
    for(int i = 0; i < mesh->shared_nodes.elem_count; ++i)
    {
        shared_node_t* sn = (shared_node_t*) sc_array_index(&mesh->shared_nodes,i);
        fprintf(mesh->fdbg, "(%d): %d %d %d\n", sn->id, sn->x, sn->y, sn->z);
        fprintf(mesh->fdbg, "     shared with processors: ");
        for(int j = 0; j < sn->listSz; j++){
            fprintf(mesh->fdbg, "%d ", sn->rankList[j]);
        }
        fprintf(mesh->fdbg, "\n");
    }
#endif
    
    // node map
    int not_my_nodes    = 0;
    int my_own_nodes    = 0;
    mesh->global_id     = (int64_t*)malloc(sizeof(int64_t)*mesh->local_n_nodes);
    memset(mesh->global_id,-2,mesh->local_n_nodes*sizeof(int64_t));
    SendTo   = (sc_hash_array_t *) sc_hash_array_new(sizeof(message_t), processors_hash_fn, processors_equal_fn, &clamped);
    RecvFrom = (sc_hash_array_t *) sc_hash_array_new(sizeof(message_t), processors_hash_fn, processors_equal_fn, &clamped);

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
    MPI_Scan(&my_own_nodes, &offset, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

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
    fprintf(mesh->fdbg, "Offset: %d\n", offset);

    fprintf(mesh->fdbg, "Nodes: \n");
    for(int i = 0; i < mesh->nodes.elem_count; ++i)
    {
        octant_node_t* n = (octant_node_t*) sc_array_index(&mesh->nodes,i);
        fprintf(mesh->fdbg, "(%d) (%d) (%ld): %d %d %d\n", n->id, mesh->part_nodes[i], mesh->global_id[i], n->x, n->y, n->z);
    }

#endif
     
    ///////////////////////////////
    ///////////////////////////////
    ///////////////////////////////
    //update the local node_id to global node_id
    for(int i = 0; i < mesh->nodes.elem_count; ++i){
        octant_node_t* n = (octant_node_t*) sc_array_index(&mesh->nodes,i);
        n->id=mesh->global_id[i];
    }
    //insert the shared nodes in the hash_array now with global_id
    shared_nodes    = (sc_hash_array_t *)sc_hash_array_new(sizeof (shared_node_t), node_hash_fn, node_equal_fn, &clamped);
    for(int i = 0; i < mesh->nodes.elem_count; i++){
        octant_node_t* node = (octant_node_t*) sc_array_index (&mesh->nodes, i);
         if(node->y == mesh->y_start){
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

    ///////////////////////////////
    ///////////////////////////////
    ///////////////////////////////
    
    //edges
    //add edge_identification    
    int nedges = 0;    
    //initialization for the edge id
    for (int iel = 0; iel < mesh->elements.elem_count; ++iel) {

        octant_t *elem = (octant_t*) sc_array_index(&mesh->elements, iel);

        for (int edge = 0; edge < 12; ++edge) {
            int Edge2GNode[2];

            int node1 = mesh->global_id[elem->nodes[EdgeVerticesMap[edge][0]].id];
            int node2 = mesh->global_id[elem->nodes[EdgeVerticesMap[edge][1]].id];
            assert(node1 >= 0);
            assert(node2 >= 0);
            
            Edge2GNode[0] = node1 <= node2 ? node1 : node2;
            Edge2GNode[1] = node1 >= node2 ? node1 : node2;
             
            elem->edge[edge].id = edge_id(Edge2GNode, indep_edges, nedges);
            elem->edge[edge].ref = false;  
            elem->edge[edge].coord[0] = Edge2GNode[0];
            elem->edge[edge].coord[1] = Edge2GNode[1];
        }
    }
    // create the shared edges hash
    for (int iel = 0; iel < mesh->elements.elem_count; ++iel) {

        octant_t *elem = (octant_t*) sc_array_index(&mesh->elements, iel);
        size_t position0;
        size_t position1;

        for (int edge = 0; edge < 12; ++edge) {

            bool out0 = false;
            bool out1 = false;
            octant_node_t* node0 = (octant_node_t*) sc_array_index (&mesh->nodes, elem->nodes[EdgeVerticesMap[edge][0]].id);
            octant_node_t* node1 = (octant_node_t*) sc_array_index (&mesh->nodes, elem->nodes[EdgeVerticesMap[edge][1]].id);

            out0 =  sc_hash_array_lookup(shared_nodes, node0, &position0);
            out1 =  sc_hash_array_lookup(shared_nodes, node1, &position1);	

            if(out0 && out1){              
                shared_node_t* sn0 = (shared_node_t*) sc_array_index(&shared_nodes->a, position0);
                shared_node_t* sn1 = (shared_node_t*) sc_array_index(&shared_nodes->a, position1);

                for(int i = 0; i < sn0->listSz; i++){
                    for(int j = 0; j < sn1->listSz; j++){
                        if(sn0->rankList[i] == sn1->rankList[j]){
                            shared_edge_t      e;
                            e.id = elem->edge[edge].id;
                            hexa_insert_shared_edge(shared_edges, &e , sn0->rankList[i]);
                        }
                    }
                }
            }     
        }
    }
    ////////////////////////////////////
    //extract the share nodes from shared_nodes
    sc_hash_array_rip (shared_nodes, &mesh->shared_nodes);

#ifdef HEXA_DEBUG_    
    fprintf(mesh->fdbg,"Shared nodes in global ids:\n");
    for(int i = 0; i < mesh->shared_nodes.elem_count; ++i){
        shared_node_t* sn = (shared_node_t*) sc_array_index(&mesh->shared_nodes,i);
        fprintf(mesh->fdbg, "(%d): %d %d %d\n", sn->id, sn->x, sn->y, sn->z);
        fprintf(mesh->fdbg, "     shared with processors: ");
        for(int j = 0; j < sn->listSz; j++){
            fprintf(mesh->fdbg, "%d ", sn->rankList[j]);
        }
        fprintf(mesh->fdbg, "\n");
    }
#endif
    ////////////////////////////////////
    //extract the edges from indep_edges
    sc_hash_array_rip (indep_edges,  &mesh->edges);
    //extract the shared edges from shared_edges
    sc_hash_array_rip (shared_edges, &mesh->shared_edges);
  
#ifdef HEXA_DEBUG_ 
    if (0){
    fprintf(mesh->fdbg,"Edge ids:\n");
    for(int i = 0; i < mesh->edges.elem_count; ++i){
        octant_edge_t* sn = (octant_edge_t*) sc_array_index(&mesh->edges,i);
        fprintf(mesh->fdbg, "Edge id:%d status:%d\n", sn->id, sn->ref);
    }
    }
#endif
    
    //edge map
    int not_my_edges    = 0;
    int my_own_edges    = 0;
    mesh->global_edge_id     = (int64_t*)malloc(sizeof(int64_t)*mesh->edges.elem_count);
    memset(mesh->global_edge_id,-2,mesh->edges.elem_count*sizeof(int64_t));
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
                mesh->global_edge_id[sn->id] = -3;
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
                mesh->global_edge_id[sn->id] = -1;
            }
        }   
        if(mesh->global_edge_id[sn->id] == -1) not_my_edges++;
        if(mesh->global_edge_id[sn->id] == -3) my_own_edges++;
    }
    
    local[2] = mesh->local_n_edges = mesh->edges.elem_count;
    local[2] -= not_my_edges;
    
    sc_hash_array_rip(RecvFrom, &mesh->comm_map_edge.RecvFrom );
    sc_hash_array_rip(SendTo  , &mesh->comm_map_edge.SendTo   );

#ifdef HEXA_DEBUG_
    fprintf(mesh->fdbg, "Shared Edges: \n");
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
    
    //size for the edges message
    mesh->comm_map_edge.max_recvbuf_size = 0;
    for(int i = 0; i < mesh->comm_map_edge.RecvFrom.elem_count; i++){
        message_t* m = (message_t*) sc_array_index(&mesh->comm_map_edge.RecvFrom, i);
        mesh->comm_map_edge.max_recvbuf_size +=  m->idxs.elem_count;
    }
    
    mesh->comm_map_edge.max_sendbuf_size = 0;
    for(int i = 0; i < mesh->comm_map_edge.SendTo.elem_count; i++){
        message_t* m = (message_t*) sc_array_index(&mesh->comm_map_edge.SendTo, i);
        mesh->comm_map_edge.max_sendbuf_size +=  m->idxs.elem_count;
    }
 
    mesh->comm_map_edge.nrequests = mesh->comm_map_edge.RecvFrom.elem_count + 
                                     mesh->comm_map_edge.SendTo.elem_count;
    
    //find the global edge_id
    offset = 0;
    MPI_Scan(&my_own_edges, &offset, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

    count = offset;
    for(int i=0; i < mesh->edges.elem_count; ++i)
        if(mesh->global_edge_id[i] != -1) mesh->global_edge_id[i] = count++;

    communicate_global_edge_ids(mesh);
    
    //redo the edge_comm_map
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
            elem->edge[edge].id = id;

        }
    }
    
    
    local[0] = mesh->local_n_nodes    = mesh->nodes.elem_count;
    local[1] = mesh->local_n_elements = mesh->elements.elem_count;
    
    
    MPI_Allreduce(local, global,3,MPI_LONG_LONG_INT,MPI_SUM,MPI_COMM_WORLD);
    
    mesh->total_n_nodes    = global[0];
    mesh->total_n_elements = global[1];
    mesh->total_n_edges = global[2];
   
    if(mesh->mpi_rank == 0)
    {
        printf("Total number of elements: %lld\n", mesh->total_n_elements);
        printf("Total number of edges: %lld\n", mesh->total_n_edges);
        printf("Total number of nodes: %lld\n", mesh->total_n_nodes);
    }
    
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
        
    free(&recvbuf[0]);
    free(&sendbuf[0]);
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
    for(int i = 0; i < mesh->comm_map_edge.RecvFrom.elem_count; ++i) {
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