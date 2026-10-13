#ifndef CONFORMING_REFINEMENT_H
#define CONFORMING_REFINEMENT_H

#include "refinement_templates.h"
#include <cstddef>
#include <utility>

namespace refinement {
// One serial, transactional pass. Points retain VTK/H5 ordering; identities
// are topological interpolation weights, never rounded physical coordinates.
struct MeshData {
    std::vector<Point3> points;
    std::vector<HexConnectivity> cells;
};
struct RefinedMesh {
    MeshData mesh;
    std::vector<int> parent; // original parent for every output cell
    std::vector<int> pads;   // chosen template per original parent; 0=unchanged
    // Weights of ORIGINAL nodes, integer numerator / 24^3.
    std::vector<std::vector<std::pair<int,int>>> node_weights;
    std::size_t requested=0, replaced=0;
};
RefinedMesh refine_conforming(const MeshData &mesh,
                             const std::vector<int> &requests,
                             std::size_t max_elements,
                             std::size_t max_search_steps=20000);
} // namespace refinement
#endif
