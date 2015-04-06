
#include <mpi.h>

#include "hexa.h"

void hexa_init(int argc, char* argv[], hexa_tree_t* mesh)
{
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &mesh->mpi_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &mesh->mpi_rank);
#ifdef HEXA_DEBUG_
    char fdname[80];
    sprintf(fdname,"debug_%04d_%04d.dbg", mesh->mpi_size, mesh->mpi_rank);
    mesh->fdbg = fopen(fdname,"w");
    
    fprintf(mesh->fdbg,"MESH INIT \n");
    fprintf(mesh->fdbg," Number of processors: %d\n Processor rank: %d\n", mesh->mpi_size, mesh->mpi_rank);
#endif
}


void hexa_processors_interval(hexa_tree_t* mesh)
{
    int x,y,i,j,m, n;
    int M = mesh->ncellx;
    int N = mesh->ncelly;
    int *lx;
    int *ly;
    
    /* try for squarish distribution */
    m = (int)(0.5 + sqrt(((double)M)*((double)mesh->mpi_size)/((double)N)));
    if (!m) m = 1;
    while (m > 0) {
        n = mesh->mpi_size/m;
        if (m*n == mesh->mpi_size) break;
        m--;
      }
      if (M > N && m < n) {int _m = m; m = n; n = _m;}
    
     if (m*n != mesh->mpi_size) printf("Given Bad partition\n");

     if (M < m) printf("Partition in x direction is too fine! %d %d\n",M,m);
     if (N < n) printf("Partition in y direction is too fine! %d %d\n",N,n);
    
    
    
    lx = (int*) malloc(m*sizeof(int));
    for (i=0; i<m; i++) {
      lx[i] = M/m + ((M % m) > i);
    }
    x  = lx[mesh->mpi_rank % m];
    int xs = 0;
    for (i=0; i<(mesh->mpi_rank % m); i++) {
        xs += lx[i];
    }
    
    
    if(mesh->mpi_rank == 0)
    {
        printf(" Processors grid: \n");
        printf("   x-direction: %d processors\n",  m);
        printf("   y-direction: %d processors\n", n );
    }
    
    /*
     Determine locally owned region
     ys is the first local node number, y is the number of local nodes
    */
    ly = (int*) malloc(n*sizeof(int));
    for (i=0; i<n; i++) {
      ly[i] = N/n + ((N % n) > i);
    }
  
    y  = ly[mesh->mpi_rank/m];
    int  ys = 0;
    for (i=0; i<(mesh->mpi_rank/m); i++) {
            ys += ly[i];
    }
    
    mesh->x_start = xs;
    mesh->x_end   = xs + x;
        
    mesh->y_start = ys;
    mesh->y_end   = ys + y;
    
        /* Assume the Non-Periodic Case */
    int n0, n1, n2, n3, n4, n5, n6, n7, n8;
    
    n1 = mesh->mpi_rank - m;
    if (mesh->mpi_rank % m) {
        n0 = n1 - 1;
    } else {
        n0 = -1;
    }
    if ((mesh->mpi_rank+1) % m) 
    {
        n2 = n1 + 1;
        n5 = mesh->mpi_rank + 1;
        n8 = mesh->mpi_rank + m + 1; if (n8 >= m*n) n8 = -1;
    } else {
        n2 = -1; n5 = -1; n8 = -1;
    }
    if (mesh->mpi_rank % m) {
        n3 = mesh->mpi_rank - 1;
        n6 = n3 + m; if (n6 >= m*n) n6 = -1;
    } else 
    {
        n3 = -1; n6 = -1;
    }
    n7 = mesh->mpi_rank + m; 
    if (n7 >= m*n) n7 = -1;
    
    
/* determine who lies on each side of us stored in      n6 n7 n8
                                                        n3    n5
                                                        n0 n1 n2
*/
    
    mesh->neighbors[0] = n0; 
    mesh->neighbors[1] = n1;
    mesh->neighbors[2] = n2;
    mesh->neighbors[3] = n3;
    mesh->neighbors[4] = mesh->mpi_rank;
    mesh->neighbors[5] = n5;
    mesh->neighbors[6] = n6;
    mesh->neighbors[7] = n7;
    mesh->neighbors[8] = n8;
  
#ifdef HEXA_DEBUG_
    fprintf(mesh->fdbg," Locally owned region\n");
    fprintf(mesh->fdbg,"  x ranges from %d to %d (len = %d) \n", mesh->x_start, mesh->x_end, x);
    fprintf(mesh->fdbg,"  y ranges from %d to %d (len = %d) \n", mesh->y_start, mesh->y_end, y);
    fprintf(mesh->fdbg,"Neighbors processors:   \n");
    fprintf(mesh->fdbg,"  Corner left-bottom  : %d\n", mesh->neighbors[0]);
    fprintf(mesh->fdbg,"  Bottom              : %d\n", mesh->neighbors[1]);
    fprintf(mesh->fdbg,"  Corner right-bottom : %d\n", mesh->neighbors[2]);
    fprintf(mesh->fdbg,"  Left                : %d\n", mesh->neighbors[3]);
    fprintf(mesh->fdbg,"  Right               : %d\n", mesh->neighbors[5]);
    fprintf(mesh->fdbg,"  Corner left-top     : %d\n", mesh->neighbors[6]);
    fprintf(mesh->fdbg,"  top                 : %d\n", mesh->neighbors[7]);
    fprintf(mesh->fdbg,"  Corner right-top    : %d\n", mesh->neighbors[8]);
#endif
    

    
    free(lx);
    free(ly);
}


void hexa_finalize(hexa_tree_t* mesh)
{
#ifdef HEXA_DEBUG_
    fclose(mesh->fdbg);
#endif
    MPI_Finalize();
}