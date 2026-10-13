#ifndef ADAPTIVE_REFINEMENT_H
#define ADAPTIVE_REFINEMENT_H
#include "hexa.h"
void ApplyAdaptiveRefinement(hexa_tree_t *mesh, std::vector<double> &coords,
                             std::vector<int> &nodes_b_mat);
#endif
