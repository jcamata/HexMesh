
#include <stdlib.h>
#include <stdio.h>
#include <sc_io.h>
#include <sc_containers.h>


#include "hexa.h"


#define VTK_FLOAT_NAME "Float32"
#define VTK_FLOAT_TYPE float
#define VTK_LOCIDX     "Int32"

#ifndef VTK_BINARY
#define VTK_ASCII 1
#define VTK_FORMAT_STRING "ascii"
#else
#define VTK_FORMAT_STRING "binary"

static int vtk_write_binary (FILE * vtkfile, char *numeric_data, size_t byte_length)
{
#ifndef VTK_COMPRESSION
  return sc_vtk_write_binary (vtkfile, numeric_data, byte_length);
#else
  return sc_vtk_write_compressed (vtkfile, numeric_data, byte_length);
#endif 
}

#endif


int hexa_tree_write_vtk(hexa_tree_t* mesh,  const char *filename)
{
 
  const int32_t Ncells = mesh->elements.elem_count;
  const int32_t Ntotal = 8 * Ncells;        /* type ok */
  

#ifdef VTK_ASCII
  double              wx, wy, wz;
  int32_t      sk;
#else
  int                 retval;
  uint8_t            *uint8_data;
  int32_t             *locidx_data;
#endif
  int                 xi, yi, zi, j, k;
  double              h2, eta_x, eta_y, eta_z;
  size_t              count, zz;
  int32_t             jt;
  int32_t             vt[8];
  int32_t             quad_count = 0;
  int32_t         il;
  VTK_FLOAT_TYPE *float_data;

  octant_t   *h;
  char                vtufilename[BUFSIZ];
  FILE               *vtufile;
 

 
  /* Have each proc write to its own file */
  snprintf (vtufilename, BUFSIZ, "%s_%d_%d.vtu", filename, mesh->mpi_size, mesh->mpi_rank);
  vtufile = fopen (vtufilename, "w");
  if (vtufile == NULL) {
    printf("Could not open %s for output!\n", vtufilename);
    return -1;
  }

  fprintf (vtufile, "<?xml version=\"1.0\"?>\n");
  fprintf (vtufile, "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\"  byte_order=\"LittleEndian\">\n");
  fprintf (vtufile, "  <UnstructuredGrid>\n");
  fprintf (vtufile,"    <Piece NumberOfPoints=\"%lld\" NumberOfCells=\"%lld\">\n", (long long) Ntotal, (long long) Ncells);
  fprintf (vtufile, "      <Points>\n");
  float_data = (VTK_FLOAT_TYPE*)malloc(sizeof(VTK_FLOAT_TYPE)*3* Ntotal);

  /* write point position data */
  fprintf (vtufile, "        <DataArray type=\"%s\" Name=\"Position\"" " NumberOfComponents=\"3\" format=\"%s\">\n", VTK_FLOAT_NAME, VTK_FORMAT_STRING);


    /* loop over the elements in the tree and calculated vertex coordinates */
    count = 0;
    for (zz = 0; zz < Ncells; ++zz) {
      octant_t* h = (octant_t*) sc_array_index (&mesh->elements, zz);
      for(jt=0; jt < 8; ++jt){
          octant_node_t* node = &h->nodes[jt];
          float_data[count*3  ] = (VTK_FLOAT_TYPE) node->x;
          float_data[count*3+1] = (VTK_FLOAT_TYPE) node->y;
          float_data[count*3+2] = (VTK_FLOAT_TYPE) node->z;
          count++;
       }
    }
   
  for (il = 0; il < Ntotal; ++il) {
    wx = float_data[3 * il + 0];
    wy = float_data[3 * il + 1];
    wz = float_data[3 * il + 2];
#ifdef VTK_DOUBLES
    fprintf (vtufile, "     %24.16e %24.16e %24.16e\n", wx, wy, wz);
#else
    fprintf (vtufile, "          %16.8e %16.8e %16.8e\n", wx, wy, wz);
#endif
  }
  free(float_data);
  fprintf (vtufile, "        </DataArray>\n");
  fprintf (vtufile, "      </Points>\n");
  fprintf (vtufile, "      <Cells>\n");

  /* write connectivity data */
  fprintf (vtufile, "        <DataArray type=\"%s\" Name=\"connectivity\" format=\"%s\">\n", VTK_LOCIDX, VTK_FORMAT_STRING);
  for (il = 0; il < Ncells; ++il) {
    sk = 8 * il;   /* type ok */
    fprintf (vtufile, "          %lld %lld %lld %lld %lld %lld %lld %lld\n",
             (long long) sk + 0, (long long) sk + 1,
             (long long) sk + 2, (long long) sk + 3,
             (long long) sk + 4, (long long) sk + 5,
             (long long) sk + 6, (long long) sk + 7);
  }
  fprintf (vtufile, "        </DataArray>\n");

  /* write offset data */
  fprintf (vtufile, "        <DataArray type=\"%s\" Name=\"offsets\" format=\"%s\">\n", VTK_LOCIDX, VTK_FORMAT_STRING);

  fprintf (vtufile, "         ");
  for (il = 1, sk = 1; il <= Ncells; ++il, ++sk) {
    fprintf (vtufile, " %lld", (long long) (8 * il));
    if (!(sk % 8) && il != Ncells)
      fprintf (vtufile, "\n         ");
  }
  fprintf (vtufile, "\n");
  fprintf (vtufile, "        </DataArray>\n");

  /* write type data */
  fprintf (vtufile, "        <DataArray type=\"UInt8\" Name=\"types\""
           " format=\"%s\">\n", VTK_FORMAT_STRING);

  fprintf (vtufile, "         ");
  for (il = 0, sk = 1; il < Ncells; ++il, ++sk) {
    fprintf (vtufile, " 12");
    if (!(sk % 20) && il != (Ncells - 1))
      fprintf (vtufile, "\n         ");
  }
  fprintf (vtufile, "\n");
  fprintf (vtufile, "        </DataArray>\n");
  fprintf (vtufile, "      </Cells>\n");

  fprintf (vtufile, "    </Piece>\n");
  fprintf (vtufile, "  </UnstructuredGrid>\n");
  fprintf (vtufile, "</VTKFile>\n");

  if (ferror (vtufile)) {
    printf ("VTKIO: Error writing footer\n");
    fclose (vtufile);
    return -1;
  }
  if (fclose (vtufile)) {
    printf ("VTKIO: Error closing footer\n");
    return -1;
  }
  
  if(mesh->mpi_rank == 0)
  {
        char                pvtufilename[BUFSIZ];
        FILE               *pvtufile;
        snprintf (pvtufilename, BUFSIZ, "%s.pvtu", filename);

        pvtufile = fopen (pvtufilename, "w");
        if (!pvtufile) {
            printf ("Could not open %s for output!\n", vtufilename);
            return -1;
        }

        fprintf (pvtufile, "<?xml version=\"1.0\"?>\n");
        fprintf (pvtufile, "<VTKFile type=\"PUnstructuredGrid\" version=\"0.1\"");
#if defined(VTK_BINARY) && defined(VTK_COMPRESSION)
        fprintf (pvtufile, " compressor=\"vtkZLibDataCompressor\"");
#endif
#ifdef SC_WORDS_BIGENDIAN
        fprintf (pvtufile, " byte_order=\"BigEndian\">\n");
#else
        fprintf (pvtufile, " byte_order=\"LittleEndian\">\n");
#endif

        fprintf (pvtufile, "  <PUnstructuredGrid GhostLevel=\"0\">\n");
        fprintf (pvtufile, "    <PPoints>\n");
        fprintf (pvtufile, "      <PDataArray type=\"%s\" Name=\"Position\""
                 " NumberOfComponents=\"3\" format=\"%s\"/>\n",
                    VTK_FLOAT_NAME, VTK_FORMAT_STRING);
        fprintf (pvtufile, "    </PPoints>\n");
        for (int p = 0; p < mesh->mpi_size; ++p) {
            fprintf (pvtufile, "    <Piece Source=\"%s_%d_%d.vtu\"/>\n", filename, mesh->mpi_size, p);
        }
        fprintf (pvtufile, "  </PUnstructuredGrid>\n");
        fprintf (pvtufile, "</VTKFile>\n");

        if (ferror (pvtufile)) {
            printf ("vtk: Error writing parallel footer\n");
            fclose (pvtufile);
            return -1;
        }
        if (fclose (pvtufile)) {
            printf ("vtk: Error closing parallel footer\n");
            return -1;
        }
        
    }
 
  return 0;
    
}

int hexa_mesh_write_vtk(hexa_tree_t* mesh,  const char *filename, std::vector<double> *coords = NULL)
{
 
  const int32_t Ncells = mesh->elements.elem_count;
  const int32_t Ntotal = mesh->nodes.elem_count; 
  

#ifdef VTK_ASCII
  double              wx, wy, wz;
  int32_t      sk;
#else
  int                 retval;
  uint8_t            *uint8_data;
  int32_t             *locidx_data;
#endif
  int                 xi, yi, zi, j, k;
  double              h2, eta_x, eta_y, eta_z;
  size_t              count, zz;
  int32_t             jt;
  int32_t             vt[8];
  int32_t             quad_count = 0;
  int32_t         il;

  octant_t   *h;
  char                vtufilename[BUFSIZ];
  FILE               *vtufile;
 

 
  /* Have each proc write to its own file */
  snprintf (vtufilename, BUFSIZ, "%s_%d_%d.vtu", filename, mesh->mpi_size, mesh->mpi_rank);
  vtufile = fopen (vtufilename, "w");
  if (vtufile == NULL) {
    printf("Could not open %s for output!\n", vtufilename);
    return -1;
  }

  fprintf (vtufile, "<?xml version=\"1.0\"?>\n");
  fprintf (vtufile, "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\"  byte_order=\"LittleEndian\">\n");
  fprintf (vtufile, "  <UnstructuredGrid>\n");
  fprintf (vtufile,"    <Piece NumberOfPoints=\"%lld\" NumberOfCells=\"%lld\">\n", (long long) Ntotal, (long long) Ncells);
  fprintf (vtufile, "      <Points>\n");
  

  /* write point position data */
  fprintf (vtufile, "        <DataArray type=\"%s\" Name=\"Position\"" " NumberOfComponents=\"3\" format=\"%s\">\n", VTK_FLOAT_NAME, VTK_FORMAT_STRING);


    /* loop over the elements in the tree and calculated vertex coordinates */
    count = 0;
    if(coords == NULL) {
        for(zz=0; zz < Ntotal; ++zz)
        {
            octant_node_t* node   = (octant_node_t*) sc_array_index (&mesh->nodes, zz);
            VTK_FLOAT_TYPE wx  = (VTK_FLOAT_TYPE) node->x;
            VTK_FLOAT_TYPE wy  = (VTK_FLOAT_TYPE) node->y;
            VTK_FLOAT_TYPE wz  = (VTK_FLOAT_TYPE) node->z;
         
#ifdef VTK_DOUBLES
            fprintf (vtufile, "     %24.16e %24.16e %24.16e\n", wx, wy, wz);
#else
            fprintf (vtufile, "     %16.8e %16.8e %16.8e\n", wx, wy, wz);
#endif
        }
    } else
    {
        for(zz=0; zz < Ntotal; ++zz)
        {
            
            VTK_FLOAT_TYPE wx  = (VTK_FLOAT_TYPE) (*coords)[zz*3  ];
            VTK_FLOAT_TYPE wy  = (VTK_FLOAT_TYPE) (*coords)[zz*3+1];
            VTK_FLOAT_TYPE wz  = (VTK_FLOAT_TYPE) (*coords)[zz*3+2];
         
#ifdef VTK_DOUBLES
            fprintf (vtufile, "     %24.16e %24.16e %24.16e\n", wx, wy, wz);
#else
            fprintf (vtufile, "     %16.8e %16.8e %16.8e\n", wx, wy, wz);
#endif 
        }
    }
        fprintf (vtufile, "        </DataArray>\n");
        fprintf (vtufile, "      </Points>\n");
    
    fprintf (vtufile, "      <Cells>\n");

    /* write connectivity data */
    fprintf (vtufile, "        <DataArray type=\"%s\" Name=\"connectivity\" format=\"%s\">\n", VTK_LOCIDX, VTK_FORMAT_STRING);
    for (zz = 0; zz < Ncells; ++zz) {
        octant_t *h     = (octant_t*) sc_array_index(&mesh->elements, zz);
        fprintf (vtufile, "          ");
        for(il=0; il < 8; il++)
        {
            octant_node_t* node = &h->nodes[il];
            fprintf (vtufile, "%ld ", node->id);     
        }
        fprintf (vtufile, "\n");
  }
    
  fprintf (vtufile, "        </DataArray>\n");
  /* write offset data */
  fprintf (vtufile, "        <DataArray type=\"%s\" Name=\"offsets\" format=\"%s\">\n", VTK_LOCIDX, VTK_FORMAT_STRING);
  fprintf (vtufile, "         ");
  for (il = 1, sk = 1; il <= Ncells; ++il, ++sk) {
    fprintf (vtufile, " %lld", (long long) (8 * il));
    if (!(sk % 8) && il != Ncells)
      fprintf (vtufile, "\n         ");
  }
  fprintf (vtufile, "\n");
  fprintf (vtufile, "        </DataArray>\n");

  /* write type data */
  fprintf (vtufile, "        <DataArray type=\"UInt8\" Name=\"types\""
           " format=\"%s\">\n", VTK_FORMAT_STRING);

  fprintf (vtufile, "         ");
  for (il = 0, sk = 1; il < Ncells; ++il, ++sk) {
    fprintf (vtufile, " 12");
    if (!(sk % 20) && il != (Ncells - 1))
      fprintf (vtufile, "\n         ");
  }
  fprintf (vtufile, "\n");
  fprintf (vtufile, "        </DataArray>\n");
  fprintf (vtufile, "      </Cells>\n");
  fprintf (vtufile, "      <CellData Scalars=\"ElemType\" >\n");
  /* write connectivity data */
  fprintf (vtufile, "        <DataArray type=\"%s\" Name=\"ElemType\" format=\"%s\">\n", VTK_LOCIDX, VTK_FORMAT_STRING);
  for (il = 0, sk = 1; il < Ncells; ++il, ++sk) {
        octant_t *h     = (octant_t*) sc_array_index(&mesh->elements, il);
        fprintf (vtufile, " %d", h->pad);
        if (!(sk % 20) && il != (Ncells - 1))
            fprintf (vtufile, "\n         ");
  }
  fprintf (vtufile, "\n");
  fprintf (vtufile, "        </DataArray>\n");
  fprintf (vtufile, "      </CellData>\n"); 
  fprintf (vtufile, "      <PointData Scalars=\"NodePart\" >\n");
  /* write connectivity data */
  fprintf (vtufile, "        <DataArray type=\"%s\" Name=\"NodePart\" format=\"%s\">\n", VTK_LOCIDX, VTK_FORMAT_STRING);
  for (il = 0, sk = 1; il < Ntotal; ++il, ++sk) {
        fprintf (vtufile, " %d", mesh->part_nodes[il]);
        if (!(sk % 20) && il != (Ntotal - 1))
            fprintf (vtufile, "\n         ");
  }
  fprintf (vtufile, "\n");
  fprintf (vtufile, "        </DataArray>\n");
  fprintf (vtufile, "      </PointData>\n"); 
  fprintf (vtufile, "    </Piece>\n");
  fprintf (vtufile, "  </UnstructuredGrid>\n");
  fprintf (vtufile, "</VTKFile>\n");

  if (ferror (vtufile)) {
    printf ("VTKIO: Error writing footer\n");
    fclose (vtufile);
    return -1;
  }
  if (fclose (vtufile)) {
    printf ("VTKIO: Error closing footer\n");
    return -1;
  }
  
  if(mesh->mpi_rank == 0)
  {
        char                pvtufilename[BUFSIZ];
        FILE               *pvtufile;
        snprintf (pvtufilename, BUFSIZ, "%s.pvtu", filename);

        pvtufile = fopen (pvtufilename, "w");
        if (!pvtufile) {
            printf ("Could not open %s for output!\n", vtufilename);
            return -1;
        }

        fprintf (pvtufile, "<?xml version=\"1.0\"?>\n");
        fprintf (pvtufile, "<VTKFile type=\"PUnstructuredGrid\" version=\"0.1\"");
#if defined(VTK_BINARY) && defined(VTK_COMPRESSION)
        fprintf (pvtufile, " compressor=\"vtkZLibDataCompressor\"");
#endif
#ifdef SC_WORDS_BIGENDIAN
        fprintf (pvtufile, " byte_order=\"BigEndian\">\n");
#else
        fprintf (pvtufile, " byte_order=\"LittleEndian\">\n");
#endif

        fprintf (pvtufile, "  <PUnstructuredGrid GhostLevel=\"0\">\n");
        fprintf (pvtufile, "    <PPoints>\n");
        fprintf (pvtufile, "      <PDataArray type=\"%s\" Name=\"Position\" NumberOfComponents=\"3\" format=\"%s\"/>\n",VTK_FLOAT_NAME, VTK_FORMAT_STRING);
        fprintf (pvtufile, "    </PPoints>\n");
        fprintf (vtufile, "     <PCellData Scalars=\"ElemType\" >\n");
        fprintf (vtufile, "        <PDataArray type=\"%s\" Name=\"ElemType\" format=\"%s\" />\n", VTK_LOCIDX, VTK_FORMAT_STRING);
        fprintf (vtufile, "      </PCellData>\n"); 
        fprintf (vtufile, "     <PPointData Scalars=\"NodePart\" >\n");
        fprintf (vtufile, "        <PDataArray type=\"%s\" Name=\"NodePart\" format=\"%s\" />\n", VTK_LOCIDX, VTK_FORMAT_STRING);
        fprintf (vtufile, "      </PPointData>\n");
        for (int p = 0; p < mesh->mpi_size; ++p) {
            fprintf (pvtufile, "    <Piece Source=\"%s_%d_%d.vtu\"/>\n", filename, mesh->mpi_size, p);
        }
        fprintf (pvtufile, "  </PUnstructuredGrid>\n");
        fprintf (pvtufile, "</VTKFile>\n");

        if (ferror (pvtufile)) {
            printf ("vtk: Error writing parallel footer\n");
            fclose (pvtufile);
            return -1;
        }
        if (fclose (pvtufile)) {
            printf ("vtk: Error closing parallel footer\n");
            return -1;
        }
        
    }
 
  return 0;
    
}


