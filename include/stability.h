#ifndef STABILITY_H
#define STABILITY_H

#include <vector>
#include <string>
#include "hexa.h"

struct ElementStability {
	int elem_id;
	int mat_id;
	double dt_crit;      // Critical time step: 2 / sqrt(lambda_max)
	double lambda_max;   // Highest generalized eigenvalue
	double omega_max;    // Maximum natural frequency sqrt(lambda_max)
	bool valid;          // True if element is geometrically valid and positive-definite
};

struct MeshStabilityAnalysis {
	double dt_crit_min;           // Minimum dt across all valid elements (Irons bound)
	double dt_crit_mean;          // Mean dt
	double dt_crit_p01;           // 1st percentile dt
	double dt_crit_p05;           // 5th percentile dt
	double dt_crit_p50;           // Median dt
	double dt_crit_max;           // Maximum dt
	double max_freq;              // Maximum natural frequency (rad/s)
	int worst_elem_id;            // Element ID that restricts dt_crit_min
	int worst_mat_id;             // Material of the worst element
	int n_unstable_elems;         // Number of invalid/inverted elements with dt=0
	std::vector<ElementStability> elem_stability;
};

// Generate Gauss-Lobatto-Legendre nodes and weights of arbitrary order N in [1, 20]
void compute_gll_nodes_and_weights(int order, std::vector<double> &nodes, std::vector<double> &weights);

// Compute element stability for a single 8-node hex using GLL order N
ElementStability compute_hex_gll_stability(int elem_id, int mat_id, const Material &mat,
                                           const double X[8], const double Y[8], const double Z[8],
                                           int gll_order,
                                           const std::vector<double> &gll_nodes,
                                           const std::vector<double> &gll_weights);

// Compute critical time step for every element in the mesh using the element-wise
// Spectral Element GLL quadrature of polynomial order N in [1, 20] and generalized eigenproblem.
MeshStabilityAnalysis analyze_mesh_stability(hexa_tree_t *mesh, const std::vector<double> &coords, int gll_order = 4);

// Print a formatted report of the stability analysis.
void print_stability_report(const char *label, const MeshStabilityAnalysis &analysis, int gll_order = 4);

#endif // STABILITY_H
