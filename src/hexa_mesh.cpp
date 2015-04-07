#include <stdlib.h>
#include <assert.h>

#include "hexa.h"

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
  const int32_t *a_ptr = (int32_t*) v;
  uint32_t        a = *a_ptr;
  uint32_t        b=1, c=1;

  sc_hash_mix (a, b, c);
  sc_hash_final (a, b, c);
  return (unsigned) c;
}

int processors_equal_fn (const void *v1, const void *v2, const void *u)
{
  const int32_t *q1 = (const int32_t*) v1;
  const int32_t *q2 = (const int32_t*) v2;
  return (*q1 == *q2);
}


void hexa_insert_shared_node(sc_hash_array_t    *shared_nodes, octant_node_t* node, int processor)
{
    size_t position;
    shared_node_t *sn;
    int i;
    
    if( processor == -1) return;
    
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



void hexa_mesh(hexa_tree_t* tree)
{
    
    bool                clamped = true;
    sc_hash_array_t    *indep_nodes;
    sc_hash_array_t    *shared_nodes;
    sc_hash_array_t    *SendTo;
    sc_hash_array_t    *RecvFrom;
    size_t              position;
    octant_node_t      *r;
    int num_indep_nodes = 0;
    int64_t    local[2], global[2];
    
    indep_nodes     = sc_hash_array_new(sizeof (octant_node_t), node_hash_fn, node_equal_fn, &clamped);
    shared_nodes    = sc_hash_array_new(sizeof (shared_node_t), node_hash_fn, node_equal_fn, &clamped);
    
    for(int i = 0; i < tree->elements.elem_count; i++)
    {
        octant_t *h  = (octant_t*) sc_array_index(&tree->elements, i);
        for(int j = 0; j < 8; j++) {
            octant_node_t* node = (octant_node_t*) sc_array_index(&h->nodes, j);
            r = (octant_node_t*) sc_hash_array_insert_unique (indep_nodes, node, &position);
            if(r != NULL)
            {
               r->x     = node->x;
               r->y     = node->y;
               r->z     = node->z;
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
    sc_hash_array_rip (indep_nodes,  &tree->nodes);
    
    for(int i = 0; i < tree->nodes.elem_count; i++)
    {
        octant_node_t* node = (octant_node_t*) sc_array_index (&tree->nodes, i);
         if(node->y == tree->y_start) 
         {
                if(node->x == tree->x_start)
                    hexa_insert_shared_node(shared_nodes,node,tree->neighbors[0]);
                if(node->x == tree->x_end)
                    hexa_insert_shared_node(shared_nodes,node,tree->neighbors[2]);
               
                hexa_insert_shared_node(shared_nodes,node,tree->neighbors[1]);
                continue;
         }
            
         if(node->y == tree->y_end) {
                if(node->x == tree->x_start)
                    hexa_insert_shared_node(shared_nodes,node,tree->neighbors[6]);
                if(node->x == tree->x_end)
                    hexa_insert_shared_node(shared_nodes,node,tree->neighbors[8]);
             
                hexa_insert_shared_node(shared_nodes,node,tree->neighbors[7]);
                continue;
          }
            
            if(node->x == tree->x_start)
                 hexa_insert_shared_node(shared_nodes,node,tree->neighbors[3]);
            if(node->x == tree->x_end)
                 hexa_insert_shared_node(shared_nodes,node,tree->neighbors[5]);
    }
    
    
    //
    // getting global index
    // building communication map
    //
    
    sc_hash_array_rip (shared_nodes, &tree->shared_nodes);
    
    local[0] = tree->local_n_nodes    = tree->nodes.elem_count;
    local[1] = tree->local_n_elements = tree->elements.elem_count;
    
    MPI_Allreduce(local, global,2,MPI_LONG_LONG_INT,MPI_SUM,MPI_COMM_WORLD);
    
    tree->total_n_nodes    = global[0];
    tree->total_n_elements = global[1];
   
    
#ifdef HEXA_DEBUG_
#if  0
    fprintf(tree->fdbg, "Nodes: \n");
    for(int i = 0; i < tree->nodes.elem_count; ++i)
    {
        octant_node_t* n = (octant_node_t*) sc_array_index(&tree->nodes,i);
        fprintf(tree->fdbg, "(%d): %d %d %d\n", n->id, n->x, n->y, n->z);
    }
    
    fprintf(tree->fdbg, "Shared Nodes: \n");
    for(int i = 0; i < tree->shared_nodes.elem_count; ++i)
    {
        shared_node_t* sn = (shared_node_t*) sc_array_index(&tree->shared_nodes,i);
        fprintf(tree->fdbg, "(%d): %d %d %d\n", sn->id, sn->x, sn->y, sn->z);
        fprintf(tree->fdbg, "     shared with processors: ");
        for(int j = 0; j < sn->listSz; j++)
            fprintf(tree->fdbg, "%d ", sn->rankList[j]);
        fprintf(tree->fdbg, "\n");
    }
#endif
#endif
    

    int not_my_nodes    = 0;
    tree->global_id = (int64_t*)malloc(sizeof(int64_t)*tree->local_n_nodes);
    memset(tree->global_id,-2,tree->local_n_nodes*sizeof(int64_t));
    SendTo   = sc_hash_array_new(sizeof(message_t), processors_hash_fn, processors_equal_fn, &clamped);
    RecvFrom = sc_hash_array_new(sizeof(message_t), processors_hash_fn, processors_equal_fn, &clamped);
    
    for(int i = 0; i < tree->shared_nodes.elem_count; ++i)
    { 
        shared_node_t* sn = (shared_node_t*) sc_array_index(&tree->shared_nodes,i);
        for(int j = 0; j < sn->listSz; j++)
            
            if(sn->rankList[j] < tree->mpi_rank) {
                //my_nodes++;
               
                message_t* m = (message_t*)sc_hash_array_insert_unique(SendTo,&sn->rankList[j],&position);
                if(m!=NULL)
                {
                    m->rank  = sn->rankList[j];
                    sc_array_init(&m->idxs, sizeof(int32_t));
                    int32_t* p = (int32_t*) sc_array_push(&m->idxs);
                    *p = sn->id;
                } else
                {
                    message_t* m = (message_t*)sc_array_index(&SendTo->a, position);
                    int32_t* p = (int32_t*) sc_array_push(&m->idxs);
                    *p = sn->id;   
                }
            } 
            else if (sn->rankList[j] > tree->mpi_rank)
            {
                message_t *m = (message_t*)sc_hash_array_insert_unique(RecvFrom,&sn->rankList[j],&position);
                if(m!=NULL)
                {
                    m->rank  = sn->rankList[j];
                    sc_array_init(&m->idxs, sizeof(int32_t));
                    int32_t* p = (int32_t*) sc_array_push(&m->idxs);
                    *p = sn->id;
                } else
                {
                    message_t* m = (message_t*)sc_array_index(&RecvFrom->a, position);
                    int32_t* p = (int32_t*) sc_array_push(&m->idxs);
                    *p = sn->id;   
                }
                tree->global_id[sn->id] = -1;
            }
          
        if(tree->global_id[sn->id] == -1) not_my_nodes++;
         
    }
    
    sc_hash_array_rip(RecvFrom, &tree->comm_map.RecvFrom );
    sc_hash_array_rip(SendTo  , &tree->comm_map.SendTo   );
    
#ifdef HEXA_DEBUG_
#if 0
    fprintf(tree->fdbg, "Recv from %ld processors\n", tree->comm_map.RecvFrom.elem_count);
    for(int i=0; i < tree->comm_map.RecvFrom.elem_count; i++)
    {
        message_t* m = (message_t*) sc_array_index(&tree->comm_map.RecvFrom, i);
        fprintf(tree->fdbg, "  Recv %ld nodes from %d\n", m->idxs.elem_count, m->rank);
    }
    fprintf(tree->fdbg, "Send to %ld processors\n", tree->comm_map.SendTo.elem_count);
    for(int i=0; i < tree->comm_map.SendTo.elem_count; i++)
    {
        message_t* m = (message_t*) sc_array_index(&tree->comm_map.SendTo, i);
        fprintf(tree->fdbg, "  Sending %ld nodes from %d\n", m->idxs.elem_count, m->rank);
    }
#endif
#endif
    
#if 1
    
    tree->comm_map.max_recvbuf_size = 0;
    for(int i = 0; i < tree->comm_map.RecvFrom.elem_count; i++)
    {
        message_t* m = (message_t*) sc_array_index(&tree->comm_map.RecvFrom, i);
        tree->comm_map.max_recvbuf_size +=  m->idxs.elem_count;
    }
    
    tree->comm_map.max_sendbuf_size = 0;
    for(int i = 0; i < tree->comm_map.SendTo.elem_count; i++)
    {
        message_t* m = (message_t*) sc_array_index(&tree->comm_map.SendTo, i);
        tree->comm_map.max_sendbuf_size +=  m->idxs.elem_count;
    }
    

    
    tree->comm_map.nrequests = tree->comm_map.RecvFrom.elem_count + 
                                     tree->comm_map.SendTo.elem_count;
    
    int offset = 0;
    int my_nodes = tree->local_n_nodes - not_my_nodes;
    MPI_Scan(&my_nodes, &offset, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    
    int count = offset;
    for(int i=0; i < tree->local_n_nodes; ++i)
        if(tree->global_id[i] != -1) tree->global_id[i] = count++;

    void communicate_global_ids(hexa_tree_t* mesh);
#endif   
    
#ifdef HEXA_DEBUG_
#if  1
    fprintf(tree->fdbg, "Offset: %d\n", offset);
    
    fprintf(tree->fdbg, "Nodes: \n");
    for(int i = 0; i < tree->nodes.elem_count; ++i)
    {
        octant_node_t* n = (octant_node_t*) sc_array_index(&tree->nodes,i);
        fprintf(tree->fdbg, "(%d) (%ld): %d %d %d\n", n->id, tree->global_id[i], n->x, n->y, n->z);
    }
#endif
#endif
    
}


void communicate_global_ids(hexa_tree_t* mesh)
{
    int          n_requests;
    
    MPI_Request *requests;
    MPI_Status  *statuses;
    uint32_t    *recvbuf;
    uint32_t    *sendbuf;
    
    n_requests = mesh->comm_map.RecvFrom.elem_count + mesh->comm_map.SendTo.elem_count;
    recvbuf    = (uint32_t*)malloc(mesh->comm_map.max_recvbuf_size*sizeof(uint32_t));
    sendbuf    = (uint32_t*)malloc(mesh->comm_map.max_sendbuf_size*sizeof(uint32_t));
    
    requests = (MPI_Request*) malloc (n_requests*sizeof(MPI_Request));
    statuses = (MPI_Status*)  malloc (n_requests*sizeof(MPI_Status));
    int c = 0;
    
    int offset = 0;
    // post all non-blocking receives
    for(int i = 0; i < mesh->comm_map.RecvFrom.elem_count; ++i) {
        message_t *m = (message_t*) sc_array_index(&mesh->comm_map.RecvFrom, i);
        MPI_Irecv(&recvbuf[offset], m->idxs.elem_count, MPI_LONG, m->rank,0,MPI_COMM_WORLD, &requests[c]);
        offset+=m->idxs.elem_count;
        c++;
    }
    
    offset = 0;
    for(int i = 0; i < mesh->comm_map.SendTo.elem_count; ++i) {
         message_t *m = (message_t*) sc_array_index(&mesh->comm_map.SendTo, i);
         for(int j = 0; j < m->idxs.elem_count; ++j)
         {
             int32_t *id = (int32_t*) sc_array_index(&m->idxs,j);
             sendbuf[offset+j] = mesh->global_id[*id];
         }
         MPI_Isend(&sendbuf[offset], m->idxs.elem_count, MPI_LONG, m->rank,0,MPI_COMM_WORLD, &requests[c]);
         offset += m->idxs.elem_count;
         c++;
    }
    
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
    
}

