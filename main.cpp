/*
 * File:   main.cpp
 * Author: camata
 *
 * Created on March 19, 2015, 1:42 PM
 */

#include <cstdlib>
#include <mpi.h>
#include <vector>

#include "hexa.h"
#include "hilbert.h"
#include "verify_mesh.h"
#include <ctime>
#include <sc.h>
#include <sc_containers.h>
#include <sc_io.h>

#include <chrono>
#include <iostream>
/*
 *
 */

int main(int argc, char **argv) {

  hexa_tree_t mesh{};

  std::vector<double> coords;
  std::vector<int> element_ids;
  std::vector<int> nodes_b_mat;
  auto start = std::chrono::steady_clock::now();

  // read input file
  inpreader(&mesh);
  int l = mesh.input.ref;
  // mpi init
  hexa_init(l, argv, &mesh);
  // set the initial number of elements in x,y,z
  // hexa_tree_init(&mesh, l);
  hexa_tree_init(&mesh, mesh.input.ref);
  // build the reference mesh
  hexa_tree_cube(&mesh);

  // deal with the mpi com
  hexa_mesh(&mesh);
  auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(
      std::chrono::steady_clock::now() - start);
  fprintf(mesh.profile, "Time in the initialization %lld millisecond(s).\n",
          elapsed.count());
  std::cout << "Time in the initialization " << elapsed.count()
            << " millisecond(s)." << std::endl;

  const char *bathy;
  const char *topo;
  printf("GetMeshFromSurface\n");
  topo = mesh.input.topo.c_str();
  bathy = mesh.input.inter.c_str();
  printf("GetMeshFromSurface\n");
  printf("Loading files:\n \t %s \n", topo);
  if (mesh.input.inter_files.empty()) {
    printf(" \t %s \n", bathy);
  } else {
    for (const auto &f : mesh.input.inter_files) printf(" \t %s \n", f.c_str());
  }
  start = std::chrono::steady_clock::now();
  // Note that here we use a gts file.
  // There is a tool called stl2gts that convert STL files to GTS.
  // It is installed together with the gts library.
  // create the geometrical mesh
  GetMeshFromSurface(&mesh, topo, coords);
  elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(
      std::chrono::steady_clock::now() - start);
  fprintf(mesh.profile, "Time in the GetMeshFromSurface %lld millisecond(s).\n",
          elapsed.count());
  std::cout << "Time in the GetMeshFromSurface " << elapsed.count()
            << " millisecond(s)." << std::endl;

  if (mesh.input.interfaceNumber == 0) {

  } else {
    // find the elements intercepted by the bathy
    start = std::chrono::steady_clock::now();
    GetInterceptedElements(&mesh, coords, element_ids, bathy);
    printf(" Elements intercepted: %zu\n\n", element_ids.size());
    elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(
        std::chrono::steady_clock::now() - start);
    fprintf(mesh.profile,
            "Time in the GetInterceptedElements %lld millisecond(s).\n",
            elapsed.count());
    std::cout << "Time in GetInterceptedElements " << elapsed.count()
              << " millisecond(s)." << std::endl;

    if (mesh.input.movingNodes == 1) {
      // apply a deformation in the mesh to fit the bathy
      start = std::chrono::steady_clock::now();
      printf(" Project nodes to the surface\n\n");
      MovingNodes(&mesh, coords, nodes_b_mat);
      elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(
          std::chrono::steady_clock::now() - start);
      fprintf(mesh.profile, "Time in the MovingNodes %lld millisecond(s).\n",
              elapsed.count());
      std::cout << "Time in MovingNodes " << elapsed.count()
                << " millisecond(s)." << std::endl;
    }
  }
  // apply material
  start = std::chrono::steady_clock::now();
  printf(" Applying material \n\n");
  element_ids.clear();
  Apply_material(&mesh, coords);
  elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(
      std::chrono::steady_clock::now() - start);
  fprintf(mesh.profile, "Time in the Apply_material %lld millisecond(s).\n",
          elapsed.count());
  std::cout << "Time in Apply_material " << elapsed.count()
            << " millisecond(s)." << std::endl;

  if (mesh.input.movingNodes == 0) {
    // do nothing
  } else {
    // do the pillow
    start = std::chrono::steady_clock::now();
    ApplyDoublePillowing(&mesh, coords, nodes_b_mat);
    elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(
        std::chrono::steady_clock::now() - start);
    fprintf(mesh.profile,
            "Time in the PillowingInterface %lld millisecond(s).\n",
            elapsed.count());
    std::cout << "Time in PillowingInterface " << elapsed.count()
              << " millisecond(s)." << std::endl;
  }

  if (mesh.input.meshOpt) {
    // Untangler runs; the size (time-step) optimization inside is switched off
    // for now -- see optimize_size in optimize_mesh.cpp (defined, currently
    // unused)
    start = std::chrono::steady_clock::now();
    MeshOptimization(&mesh, coords, nodes_b_mat);
    elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(
        std::chrono::steady_clock::now() - start);
    fprintf(mesh.profile, "Time in the MeshOptimization %lld millisecond(s).\n",
            elapsed.count());
    std::cout << "Time in MeshOptimization " << elapsed.count()
              << " millisecond(s)." << std::endl;
  }

  // final checks: orientation, then face planarity
  VerifyMeshInversion(&mesh, &coords);
  VerifyFacePlanarity(&mesh, coords);

  if (mesh.input.PML == 0) {
    // do nothing
  } else {
    // add pml
    // SurfaceIdentification populates mesh->outsurf (ExtrudePMLElements' only source of
    // elements to extrude) and tags domain-boundary node colors ExtrudePMLElements reads
    // directly -- required regardless of movingNodes. RedoNodeMapping's integer-lattice
    // rescale is legacy setup for the old PillowingInterface() pillow path (now dead code,
    // ApplyDoublePillowing replaced it for movingNodes==1) and SurfaceIdentification's own
    // domain-boundary check (node.x == 3*mesh->ncellx) is scale-invariant under that rescale,
    // so it is not needed here.
    if (mesh.input.movingNodes == 0) {
      RedoNodeMapping(&mesh);
    }
    SurfaceIdentification(&mesh, coords);
    start = std::chrono::steady_clock::now();
    printf(" Extrude elements\n\n");
    ExtrudePMLElements(&mesh, coords);
    elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(
        std::chrono::steady_clock::now() - start);
    fprintf(mesh.profile,
            "Time in the ExtrudePMLElements %lld millisecond(s).\n",
            elapsed.count());
    std::cout << "Time in ExtrudePMLElements " << elapsed.count()
              << " millisecond(s)." << std::endl;
  }
  // clean vectors
  // std::vector<int>().swap(element_ids);
  // std::vector<int>().swap(nodes_b_mat);

  start = std::chrono::steady_clock::now();
  printf(" Writing output files \n\n");
  hexa_mesh_write_vtk(&mesh, "mesh", &coords);
  // hexa_mesh_write_msh(&mesh, "mesh", &coords);
  // hexa_mesh_write_h5(&mesh, "mesh", coords);
  elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(
      std::chrono::steady_clock::now() - start);
  fprintf(mesh.profile, "Time in Writing output files %lld millisecond(s).\n",
          elapsed.count());
  std::cout << "Time in Writing output files " << elapsed.count()
            << " millisecond(s)." << std::endl;

  hexa_mesh_write_vtk(&mesh, "test", NULL);
  start = std::chrono::steady_clock::now();

  printf(" Cleaning variables \n\n");

  hexa_tree_destroy(&mesh);
  hexa_finalize(&mesh);
  std::vector<double>().swap(coords);
  elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(
      std::chrono::steady_clock::now() - start);
  fprintf(mesh.profile, "Time in the Cleaning variables %lld millisecond(s).\n",
          elapsed.count());
  std::cout << "Time in Cleaning variables " << elapsed.count()
            << " millisecond(s)." << std::endl;

  return 0;
}
