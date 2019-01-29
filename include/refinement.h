#include <stdio.h>
#include <stdlib.h>
#include <sc.h>
#include <sc_containers.h>
#include <gts.h>
#include <vector>

void IdentifyTemplate(hexa_tree_t* mesh, const std::vector<double>& coords, std::vector<int>& elements_ids);

void ApplyOctreeTemplate(hexa_tree_t* mesh, std::vector<double>& coords, std::vector<int>& elements_ids);

int AddPoint(hexa_tree_t* mesh, sc_hash_array_t* hash, GtsPoint *p, std::vector<double> &coords);

void IdentifyTemplate(hexa_tree_t* mesh, const std::vector<double>& coords, std::vector<int>& elements_ids);






