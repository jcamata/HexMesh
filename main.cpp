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
#include "mesh_geom.h"
#include "stability.h"
#include <ctime>
#include <sc.h>
#include <sc_containers.h>
#include <sc_io.h>

#include <chrono>
#include <iostream>
#include <array>
#include <map>
#include <algorithm>
/*
 *
 */

int main(int argc, char **argv) {
  MPI_Init(&argc, &argv);

  hexa_tree_t mesh{};

  std::vector<double> coords;
  std::vector<int> element_ids;
  std::vector<int> nodes_b_mat;
  auto start = std::chrono::steady_clock::now();

  // read input file
  const char *input_path = (argc > 1) ? argv[1] : nullptr;
  inpreader(&mesh, input_path);
  // mpi init
  hexa_init(argc, argv, &mesh);
  // set the initial number of elements in x,y,z
  // hexa_tree_init(&mesh, l);
  hexa_tree_init(&mesh, mesh.input.ref);
  // build the reference mesh
  hexa_tree_cube(&mesh);

  // deal with the mpi com
  hexa_mesh(&mesh);
  auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(
      std::chrono::steady_clock::now() - start);
  // fprintf(mesh.profile, "Time in the initialization %lld millisecond(s).\n",
  //         elapsed.count());
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
  // fprintf(mesh.profile, "Time in the GetMeshFromSurface %lld millisecond(s).\n",
  //         elapsed.count());
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
    // fprintf(mesh.profile,
    //         "Time in the GetInterceptedElements %lld millisecond(s).\n",
    //         elapsed.count());
    std::cout << "Time in GetInterceptedElements " << elapsed.count()
              << " millisecond(s)." << std::endl;

    if (mesh.input.movingNodes == 1) {
      // apply a deformation in the mesh to fit the bathy
      start = std::chrono::steady_clock::now();
      printf(" Project nodes to the surface\n\n");
      MovingNodes(&mesh, coords, nodes_b_mat);
      elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(
          std::chrono::steady_clock::now() - start);
      // fprintf(mesh.profile, "Time in the MovingNodes %lld millisecond(s).\n",
      //         elapsed.count());
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
  // fprintf(mesh.profile, "Time in the Apply_material %lld millisecond(s).\n",
  //         elapsed.count());
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
    // fprintf(mesh.profile,
    //         "Time in the PillowingInterface %lld millisecond(s).\n",
    //         elapsed.count());
    std::cout << "Time in PillowingInterface " << elapsed.count()
              << " millisecond(s)." << std::endl;
    DumpInversionMap(&mesh, coords, std::vector<double>(), "postpillow");
  }

  if (mesh.input.meshOpt) {
    // Untangler runs; the size (time-step) optimization inside is switched off
    // for now -- see optimize_size in optimize_mesh.cpp (defined, currently
    // unused)
    start = std::chrono::steady_clock::now();
    MeshOptimization(&mesh, coords, nodes_b_mat);
    elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(
        std::chrono::steady_clock::now() - start);
    // fprintf(mesh.profile, "Time in the MeshOptimization %lld millisecond(s).\n",
    //         elapsed.count());
    std::cout << "Time in MeshOptimization " << elapsed.count()
              << " millisecond(s)." << std::endl;
  }

  // final checks: orientation, then face planarity
  VerifyMeshInversion(&mesh, &coords);
  VerifyFacePlanarity(&mesh, coords);

  if (!mesh.input.meshOpt) {
    MeshStabilityAnalysis stab = analyze_mesh_stability(&mesh, coords, mesh.input.gll_order);
    print_stability_report("FINAL MESH", stab, mesh.input.gll_order);
  }

  { // ponytail: bug3 evidence -- dump all inverted element ids for A/B comparison, remove after
    const char *dbgpath = getenv("HEXMESH_INVDUMP");
    if (dbgpath) {
      MeshAnalysis a = analyze_mesh(&mesh, coords);
      FILE *f = fopen(dbgpath, "w");
      for (int id : a.inverted_ids) fprintf(f, "%d\n", id);
      fclose(f);
    }
  }

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
    // fprintf(mesh.profile,
    //         "Time in the ExtrudePMLElements %lld millisecond(s).\n",
    //         elapsed.count());
    std::cout << "Time in ExtrudePMLElements " << elapsed.count()
              << " millisecond(s)." << std::endl;
  }
  // clean vectors
  // std::vector<int>().swap(element_ids);
  // std::vector<int>().swap(nodes_b_mat);

  // Per-element diagnostic tag for visual inspection in ParaView: 2 = inverted element,
  // 1 = face-neighbour of an inverted element, 0 = everything else. Computed on the FINAL
  // mesh (after PML, if any) from the real physical coords -- the authoritative geometry --
  // so it applies unchanged to both mesh.pvtu (real coords) and test.pvtu (integer lattice
  // coords): the tag is about which elements are geometrically bad, not about which
  // coordinate system is being rendered.
  std::vector<int> invtag(mesh.elements.elem_count, 0);
  std::vector<int> foldtag(mesh.elements.elem_count, 0);
  {
    MeshAnalysis a2 = analyze_mesh(&mesh, coords);
    for (int id : a2.inverted_ids) invtag[id] = 2;

    // Face-adjacency scan: pairs every element's 6 faces against every other element's, once.
    // Reused for two independent checks:
    //  (a) InvertedTag neighbour marking (2=inverted, 1=neighbour of one, existing behaviour).
    //  (b) NeighborFold: two face-adjacent elements should extend AWAY from their shared face,
    //      one to each side. Project each element's own OTHER 4 corners' centroid onto the
    //      shared face's normal; a normal (non-folded) pair lands on opposite sides. Landing on
    //      the SAME side means one of the pair has folded back over the shared face into the
    //      other's own volume -- real 3D interpenetration, not just a bad-quality-but-present
    //      element (which is what InvertedTag / the corner-Jacobian check alone measures, and
    //      does NOT catch this: two elements can each be locally valid, non-inverted, and still
    //      overlap each other -- confirmed on 2026-08-19 with elems 13011/527802, 527818/527819/
    //      527810, all invtag==0 but visibly folded over one another in ParaView).
    std::map<std::array<int,4>, std::pair<int,int>> face_owner; // key -> (elem id, face index)
    for (size_t iel = 0; iel < mesh.elements.elem_count; iel++) {
      octant_t *e = (octant_t *)sc_array_index(&mesh.elements, iel);
      for (int f = 0; f < 6; f++) {
        std::array<int,4> key;
        for (int k = 0; k < 4; k++) key[k] = e->nodes[FaceNodesMap[f][k]].id;
        std::array<int,4> sorted_key = key;
        std::sort(sorted_key.begin(), sorted_key.end());
        auto it = face_owner.find(sorted_key);
        if (it == face_owner.end()) {
          face_owner[sorted_key] = {(int)iel, f};
        } else {
          int other = it->second.first;

          if (invtag[iel] == 2 && invtag[other] == 0) invtag[other] = 1;
          if (invtag[other] == 2 && invtag[iel] == 0) invtag[iel] = 1;

          // Shared-face corners in THIS element's own traversal order (arbitrary but
          // consistent -- the fold test only needs the plane, and a plane doesn't care
          // which side supplied it).
          double q[4][3];
          for (int k = 0; k < 4; k++)
            for (int d = 0; d < 3; d++) q[k][d] = coords[3*key[k]+d];

          auto complement_centroid = [&](octant_t *elem, int face) -> std::array<double,3> {
            bool on_face[8] = {false};
            for (int k = 0; k < 4; k++) on_face[FaceNodesMap[face][k]] = true;
            std::array<double,3> c = {0,0,0};
            int cnt = 0;
            for (int k = 0; k < 8; k++) {
              if (on_face[k]) continue;
              int nid = elem->nodes[k].id;
              for (int d = 0; d < 3; d++) c[d] += coords[3*nid+d];
              cnt++;
            }
            for (int d = 0; d < 3; d++) c[d] /= cnt; // cnt is always 4
            return c;
          };
          octant_t *eo = (octant_t *)sc_array_index(&mesh.elements, other);
          auto cA = complement_centroid(e, f);
          auto cB = complement_centroid(eo, it->second.second);
          if (mgeom::faces_folded(q, cA.data(), cB.data())) {
            foldtag[iel] = 1;
            foldtag[other] = 1;
          }

          face_owner.erase(it);
        }
      }
    }
    int n_folded = 0, n_both = 0;
    for (size_t i = 0; i < foldtag.size(); i++) {
      if (foldtag[i]) n_folded++;
      if (foldtag[i] && invtag[i] == 2) n_both++;
    }
    printf("    Neighbour-fold elements (real 3D overlap with a face-adjacent neighbour, "
           "not caught by the corner-Jacobian inversion check): %d / %zu (%d also inverted)\n",
           n_folded, mesh.elements.elem_count, n_both);

  }

  MeshStabilityAnalysis stab_final = analyze_mesh_stability(&mesh, coords, mesh.input.gll_order);
  if (!mesh.input.meshOpt) {
    print_stability_report("FINAL MESH", stab_final, mesh.input.gll_order);
  }
  // Duas quantidades, exportadas lado a lado: dt_crit e' o limite de Irons
  // (2/sqrt(lambda_max), iteracao de potencia sobre M_e^-1 K_e) e dt_cfl e' a CFL
  // classica (C*h_min/vp). Em malha regular diferem ~1.4x; em hexaedros
  // deformados a CFL classica sobrestima o passo por ordens de grandeza, porque
  // mede distancias e nao o operador. Exportar as duas torna a diferenca
  // mensuravel na malha em vez de suposta.
  std::vector<double> dtcrit_elem(mesh.elements.elem_count, 0.0);
  std::vector<double> dtcfl_elem(mesh.elements.elem_count, 0.0);
  std::vector<double> jacratio_elem(mesh.elements.elem_count, 0.0);
  for (size_t i = 0; i < stab_final.elem_stability.size(); i++) {
    dtcrit_elem[i]   = stab_final.elem_stability[i].dt_crit;
    dtcfl_elem[i]    = stab_final.elem_stability[i].dt_cfl;
    jacratio_elem[i] = stab_final.elem_stability[i].jac_ratio;
  }

  start = std::chrono::steady_clock::now();
  std::string out_prefix = mesh.input.output_prefix.empty() ? "mesh" : mesh.input.output_prefix;
  printf(" Writing output files (%s.pvtu)\n\n", out_prefix.c_str());
  hexa_mesh_write_vtk(&mesh, out_prefix.c_str(), &coords, &invtag, &foldtag, &dtcrit_elem, &dtcfl_elem, &jacratio_elem);
  // hexa_mesh_write_msh(&mesh, out_prefix.c_str(), &coords);
  // Saida HDF5: desligada por omissao (um sweep completo do run_cases.sh enche o
  // disco), ligada por caso com "writeH5 = 1" no .input. E' o formato que o
  // MeshClass do FEM le' diretamente, ao contrario do VTU.
  if (mesh.input.writeH5) {
    printf(" Writing HDF5 output (%s_%04d_%04d.h5)\n", out_prefix.c_str(),
           mesh.mpi_size, mesh.mpi_rank);
    hexa_mesh_write_h5(&mesh, out_prefix.c_str(), coords, &dtcrit_elem, &invtag, &foldtag, &dtcfl_elem, &jacratio_elem);
  }
  elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(
      std::chrono::steady_clock::now() - start);
  // fprintf(mesh.profile, "Time in Writing output files %lld millisecond(s).\n",
  //         elapsed.count());
  std::cout << "Time in Writing output files " << elapsed.count()
            << " millisecond(s)." << std::endl;

  std::string test_prefix = out_prefix + "_lattice";
  hexa_mesh_write_vtk(&mesh, test_prefix.c_str(), NULL, &invtag, &foldtag);
  start = std::chrono::steady_clock::now();

  printf(" Cleaning variables \n\n");

  hexa_tree_destroy(&mesh);
  hexa_finalize(&mesh);
  std::vector<double>().swap(coords);
  elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(
      std::chrono::steady_clock::now() - start);
  // fprintf(mesh.profile, "Time in the Cleaning variables %lld millisecond(s).\n",
  //         elapsed.count());
  std::cout << "Time in Cleaning variables " << elapsed.count()
            << " millisecond(s)." << std::endl;

  return 0;
}
