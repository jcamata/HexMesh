// Implementation of reference templates and orientation data.
// Recovered from origin/template_ref pinned at commit 1f146f1e7ddd85ff5f63700dc10f197387cba6d1.

#include "refinement_templates.h"
#include <cmath>
#include <algorithm>
#include <stdexcept>
#include <map>

namespace refinement {

// Standard Hex Face Definitions: outward-oriented quad loops in VTK/H5 order
static constexpr int HEX_FACES[6][4] = {
    {0, 3, 2, 1}, // Face 0: bottom (-z)
    {4, 5, 6, 7}, // Face 1: top (+z)
    {0, 1, 5, 4}, // Face 2: front (-y)
    {1, 2, 6, 5}, // Face 3: right (+x)
    {2, 3, 7, 6}, // Face 4: back (+y)
    {3, 0, 4, 7}  // Face 5: left (-x)
};

// Child corner reversal permutation when reflection parity is present
// Swaps (0,3)(1,2)(4,7)(5,6) to invert orientation while preserving hex topology
static constexpr int FLIP_CORNER_PERM[8] = {3, 2, 1, 0, 7, 6, 5, 4};

// Template 1: 5 children, 16 unique vertices
static const RationalCoord3 T1_VERTICES[16] = {
    {0, 0, 0},
    {24, 0, 0},
    {24, 24, 0},
    {0, 24, 0},
    {0, 0, 24},
    {24, 0, 24},
    {24, 24, 24},
    {0, 24, 24},
    {8, 0, 0},
    {8, 16, 0},
    {8, 0, 16},
    {8, 16, 16},
    {16, 0, 0},
    {16, 16, 0},
    {16, 0, 16},
    {16, 16, 16},
};

static const HexConnectivity T1_CHILDREN[5] = {
    {0, 8, 9, 3, 4, 10, 11, 7},
    {8, 12, 13, 9, 10, 14, 15, 11},
    {12, 1, 2, 13, 14, 5, 6, 15},
    {9, 13, 2, 3, 11, 15, 6, 7},
    {10, 14, 15, 11, 4, 5, 6, 7},
};

// Template 2: 3 children, 16 unique vertices
static const RationalCoord3 T2_VERTICES[16] = {
    {0, 0, 0},
    {24, 0, 0},
    {24, 24, 0},
    {0, 24, 0},
    {0, 0, 24},
    {24, 0, 24},
    {24, 24, 24},
    {0, 24, 24},
    {0, 0, 8},
    {24, 0, 8},
    {24, 24, 8},
    {0, 24, 8},
    {0, 0, 16},
    {24, 0, 16},
    {24, 24, 16},
    {0, 24, 16},
};

static const HexConnectivity T2_CHILDREN[3] = {
    {0, 1, 2, 3, 8, 9, 10, 11},
    {8, 9, 10, 11, 12, 13, 14, 15},
    {12, 13, 14, 15, 4, 5, 6, 7},
};

// Template 3: 4 children, 16 unique vertices
static const RationalCoord3 T3_VERTICES[16] = {
    {0, 0, 0},
    {24, 0, 0},
    {24, 24, 0},
    {0, 24, 0},
    {0, 0, 24},
    {24, 0, 24},
    {24, 24, 24},
    {0, 24, 24},
    {24, 8, 0},
    {0, 8, 0},
    {24, 8, 16},
    {0, 8, 16},
    {24, 16, 0},
    {0, 16, 0},
    {24, 16, 16},
    {0, 16, 16},
};

static const HexConnectivity T3_CHILDREN[4] = {
    {0, 1, 8, 9, 4, 5, 10, 11},
    {9, 8, 12, 13, 11, 10, 14, 15},
    {13, 12, 2, 3, 15, 14, 6, 7},
    {11, 10, 14, 15, 4, 5, 6, 7},
};

// Template 4: 5 children, 20 unique vertices
static const RationalCoord3 T4_VERTICES[20] = {
    {0, 0, 0},
    {24, 0, 0},
    {24, 24, 0},
    {0, 24, 0},
    {0, 0, 24},
    {24, 0, 24},
    {24, 24, 24},
    {0, 24, 24},
    {8, 0, 0},
    {8, 8, 0},
    {0, 8, 0},
    {8, 0, 24},
    {8, 8, 24},
    {0, 8, 24},
    {16, 0, 0},
    {16, 16, 0},
    {16, 0, 24},
    {16, 16, 24},
    {0, 16, 0},
    {0, 16, 24},
};

static const HexConnectivity T4_CHILDREN[5] = {
    {0, 8, 9, 10, 4, 11, 12, 13},
    {8, 14, 15, 9, 11, 16, 17, 12},
    {10, 9, 15, 18, 13, 12, 17, 19},
    {14, 1, 2, 15, 16, 5, 6, 17},
    {18, 15, 2, 3, 19, 17, 6, 7},
};

// Template 5: 13 children, 32 unique vertices
static const RationalCoord3 T5_VERTICES[32] = {
    {0, 0, 0},
    {24, 0, 0},
    {24, 24, 0},
    {0, 24, 0},
    {0, 0, 24},
    {24, 0, 24},
    {24, 24, 24},
    {0, 24, 24},
    {0, 8, 8},
    {24, 8, 8},
    {24, 16, 8},
    {0, 16, 8},
    {8, 8, 16},
    {16, 8, 16},
    {16, 16, 16},
    {8, 16, 16},
    {8, 8, 24},
    {16, 8, 24},
    {16, 16, 24},
    {8, 16, 24},
    {8, 0, 8},
    {16, 0, 8},
    {16, 24, 8},
    {8, 24, 8},
    {8, 0, 24},
    {0, 8, 24},
    {16, 0, 24},
    {24, 8, 24},
    {24, 16, 24},
    {16, 24, 24},
    {8, 24, 24},
    {0, 16, 24},
};

static const HexConnectivity T5_CHILDREN[13] = {
    {0, 1, 2, 3, 8, 9, 10, 11},
    {8, 9, 10, 11, 12, 13, 14, 15},
    {12, 13, 14, 15, 16, 17, 18, 19},
    {0, 1, 9, 8, 20, 21, 13, 12},
    {11, 10, 2, 3, 15, 14, 22, 23},
    {0, 20, 12, 8, 4, 24, 16, 25},
    {20, 21, 13, 12, 24, 26, 17, 16},
    {21, 1, 9, 13, 26, 5, 27, 17},
    {13, 9, 10, 14, 17, 27, 28, 18},
    {14, 10, 2, 22, 18, 28, 6, 29},
    {15, 14, 22, 23, 19, 18, 29, 30},
    {11, 15, 23, 3, 31, 19, 30, 7},
    {8, 12, 15, 11, 25, 16, 19, 31},
};

// Template 6: 10 children, 32 unique vertices
static const RationalCoord3 T6_VERTICES[32] = {
    {0, 0, 0},
    {24, 0, 0},
    {24, 24, 0},
    {0, 24, 0},
    {0, 0, 24},
    {24, 0, 24},
    {24, 24, 24},
    {0, 24, 24},
    {8, 0, 0},
    {12, 8, 0},
    {6, 8, 0},
    {8, 0, 24},
    {12, 8, 24},
    {6, 8, 24},
    {16, 0, 0},
    {18, 8, 0},
    {16, 0, 24},
    {18, 8, 24},
    {24, 8, 0},
    {24, 8, 24},
    {12, 16, 0},
    {6, 16, 0},
    {12, 16, 24},
    {6, 16, 24},
    {18, 16, 0},
    {18, 16, 24},
    {24, 16, 0},
    {24, 16, 24},
    {8, 24, 0},
    {8, 24, 24},
    {16, 24, 0},
    {16, 24, 24},
};

static const HexConnectivity T6_CHILDREN[10] = {
    {0, 8, 9, 10, 4, 11, 12, 13},
    {8, 14, 15, 9, 11, 16, 17, 12},
    {14, 1, 18, 15, 16, 5, 19, 17},
    {10, 9, 20, 21, 13, 12, 22, 23},
    {9, 15, 24, 20, 12, 17, 25, 22},
    {15, 18, 26, 24, 17, 19, 27, 25},
    {21, 20, 28, 3, 23, 22, 29, 7},
    {20, 24, 30, 28, 22, 25, 31, 29},
    {24, 26, 2, 30, 25, 27, 6, 31},
    {0, 10, 21, 3, 4, 13, 23, 7},
};

// Template 7: 7 children, 22 unique vertices
static const RationalCoord3 T7_VERTICES[22] = {
    {0, 0, 0},
    {24, 0, 0},
    {24, 24, 0},
    {0, 24, 0},
    {0, 0, 24},
    {24, 0, 24},
    {24, 24, 24},
    {0, 24, 24},
    {8, 0, 0},
    {8, 8, 0},
    {0, 8, 0},
    {0, 0, 8},
    {8, 0, 8},
    {8, 8, 8},
    {0, 8, 8},
    {16, 0, 0},
    {16, 16, 0},
    {16, 0, 16},
    {16, 16, 16},
    {0, 16, 0},
    {0, 16, 16},
    {0, 0, 16},
};

static const HexConnectivity T7_CHILDREN[7] = {
    {0, 8, 9, 10, 11, 12, 13, 14},
    {8, 15, 16, 9, 12, 17, 18, 13},
    {15, 1, 2, 16, 17, 5, 6, 18},
    {10, 9, 16, 19, 14, 13, 18, 20},
    {11, 12, 13, 14, 21, 17, 18, 20},
    {21, 17, 18, 20, 4, 5, 6, 7},
    {19, 16, 2, 3, 20, 18, 6, 7},
};

// Template 8: 9 children, 32 unique vertices
static const RationalCoord3 T8_VERTICES[32] = {
    {0, 0, 0},
    {24, 0, 0},
    {24, 24, 0},
    {0, 24, 0},
    {0, 0, 24},
    {24, 0, 24},
    {24, 24, 24},
    {0, 24, 24},
    {8, 0, 0},
    {8, 8, 0},
    {0, 8, 0},
    {8, 0, 24},
    {8, 8, 24},
    {0, 8, 24},
    {16, 0, 0},
    {16, 8, 0},
    {16, 0, 24},
    {16, 8, 24},
    {24, 8, 0},
    {24, 8, 24},
    {8, 16, 0},
    {0, 16, 0},
    {8, 16, 24},
    {0, 16, 24},
    {16, 16, 0},
    {16, 16, 24},
    {24, 16, 0},
    {24, 16, 24},
    {8, 24, 0},
    {8, 24, 24},
    {16, 24, 0},
    {16, 24, 24},
};

static const HexConnectivity T8_CHILDREN[9] = {
    {0, 8, 9, 10, 4, 11, 12, 13},
    {8, 14, 15, 9, 11, 16, 17, 12},
    {14, 1, 18, 15, 16, 5, 19, 17},
    {10, 9, 20, 21, 13, 12, 22, 23},
    {9, 15, 24, 20, 12, 17, 25, 22},
    {15, 18, 26, 24, 17, 19, 27, 25},
    {21, 20, 28, 3, 23, 22, 29, 7},
    {20, 24, 30, 28, 22, 25, 31, 29},
    {24, 26, 2, 30, 25, 27, 6, 31},
};

// Template 9: 31 children, 64 unique vertices
static const RationalCoord3 T9_VERTICES[64] = {
    {0, 0, 0},
    {24, 0, 0},
    {24, 24, 0},
    {0, 24, 0},
    {0, 0, 24},
    {24, 0, 24},
    {24, 24, 24},
    {0, 24, 24},
    {0, 8, 6},
    {24, 8, 6},
    {24, 16, 6},
    {0, 16, 6},
    {8, 8, 9},
    {16, 8, 9},
    {16, 16, 9},
    {8, 16, 9},
    {8, 0, 6},
    {16, 0, 6},
    {16, 24, 6},
    {8, 24, 6},
    {8, 8, 12},
    {16, 8, 12},
    {16, 16, 12},
    {8, 16, 12},
    {8, 8, 18},
    {16, 8, 18},
    {16, 16, 18},
    {8, 16, 18},
    {8, 8, 24},
    {16, 8, 24},
    {16, 16, 24},
    {8, 16, 24},
    {0, 0, 8},
    {8, 0, 12},
    {0, 8, 12},
    {16, 0, 12},
    {24, 0, 8},
    {24, 8, 12},
    {24, 16, 12},
    {24, 24, 8},
    {16, 24, 12},
    {8, 24, 12},
    {0, 16, 12},
    {0, 24, 8},
    {0, 0, 16},
    {8, 0, 18},
    {0, 8, 18},
    {16, 0, 18},
    {24, 0, 16},
    {24, 8, 18},
    {24, 16, 18},
    {24, 24, 16},
    {16, 24, 18},
    {8, 24, 18},
    {0, 16, 18},
    {0, 24, 16},
    {8, 0, 24},
    {0, 8, 24},
    {16, 0, 24},
    {24, 8, 24},
    {24, 16, 24},
    {16, 24, 24},
    {8, 24, 24},
    {0, 16, 24},
};

static const HexConnectivity T9_CHILDREN[31] = {
    {0, 1, 2, 3, 8, 9, 10, 11},
    {8, 9, 10, 11, 12, 13, 14, 15},
    {0, 1, 9, 8, 16, 17, 13, 12},
    {11, 10, 2, 3, 15, 14, 18, 19},
    {12, 13, 14, 15, 20, 21, 22, 23},
    {20, 21, 22, 23, 24, 25, 26, 27},
    {24, 25, 26, 27, 28, 29, 30, 31},
    {0, 16, 12, 8, 32, 33, 20, 34},
    {16, 17, 13, 12, 33, 35, 21, 20},
    {17, 1, 9, 13, 35, 36, 37, 21},
    {13, 9, 10, 14, 21, 37, 38, 22},
    {14, 10, 2, 18, 22, 38, 39, 40},
    {15, 14, 18, 19, 23, 22, 40, 41},
    {11, 15, 19, 3, 42, 23, 41, 43},
    {8, 12, 15, 11, 34, 20, 23, 42},
    {32, 33, 20, 34, 44, 45, 24, 46},
    {33, 35, 21, 20, 45, 47, 25, 24},
    {35, 36, 37, 21, 47, 48, 49, 25},
    {21, 37, 38, 22, 25, 49, 50, 26},
    {22, 38, 39, 40, 26, 50, 51, 52},
    {23, 22, 40, 41, 27, 26, 52, 53},
    {42, 23, 41, 43, 54, 27, 53, 55},
    {34, 20, 23, 42, 46, 24, 27, 54},
    {44, 45, 24, 46, 4, 56, 28, 57},
    {45, 47, 25, 24, 56, 58, 29, 28},
    {47, 48, 49, 25, 58, 5, 59, 29},
    {25, 49, 50, 26, 29, 59, 60, 30},
    {26, 50, 51, 52, 30, 60, 6, 61},
    {27, 26, 52, 53, 31, 30, 61, 62},
    {54, 27, 53, 55, 63, 31, 62, 7},
    {46, 24, 27, 54, 57, 28, 31, 63},
};

// Template 10: 17 children, 40 unique vertices
static const RationalCoord3 T10_VERTICES[40] = {
    {0, 0, 0},
    {24, 0, 0},
    {24, 24, 0},
    {0, 24, 0},
    {0, 0, 24},
    {24, 0, 24},
    {24, 24, 24},
    {0, 24, 24},
    {8, 0, 0},
    {8, 6, 0},
    {0, 8, 0},
    {0, 0, 8},
    {8, 0, 6},
    {8, 6, 6},
    {0, 8, 8},
    {16, 0, 0},
    {16, 6, 0},
    {16, 0, 6},
    {16, 6, 6},
    {24, 8, 0},
    {24, 0, 8},
    {24, 8, 8},
    {0, 0, 16},
    {8, 0, 12},
    {8, 12, 12},
    {0, 16, 16},
    {16, 0, 12},
    {16, 12, 12},
    {24, 0, 16},
    {24, 16, 16},
    {8, 0, 18},
    {8, 18, 18},
    {16, 0, 18},
    {16, 18, 18},
    {8, 12, 0},
    {0, 16, 0},
    {16, 12, 0},
    {24, 16, 0},
    {8, 18, 0},
    {16, 18, 0},
};

static const HexConnectivity T10_CHILDREN[17] = {
    {0, 8, 9, 10, 11, 12, 13, 14},
    {8, 15, 16, 9, 12, 17, 18, 13},
    {15, 1, 19, 16, 17, 20, 21, 18},
    {11, 12, 13, 14, 22, 23, 24, 25},
    {12, 17, 18, 13, 23, 26, 27, 24},
    {17, 20, 21, 18, 26, 28, 29, 27},
    {22, 23, 24, 25, 4, 30, 31, 7},
    {23, 26, 27, 24, 30, 32, 33, 31},
    {26, 28, 29, 27, 32, 5, 6, 33},
    {10, 9, 34, 35, 14, 13, 24, 25},
    {9, 16, 36, 34, 13, 18, 27, 24},
    {16, 19, 37, 36, 18, 21, 29, 27},
    {35, 34, 38, 3, 25, 24, 31, 7},
    {34, 36, 39, 38, 24, 27, 33, 31},
    {36, 37, 2, 39, 27, 29, 6, 33},
    {30, 32, 33, 31, 4, 5, 6, 7},
    {38, 39, 2, 3, 31, 33, 6, 7},
};

// Template 11: 27 children, 64 unique vertices
static const RationalCoord3 T11_VERTICES[64] = {
    {0, 0, 0},
    {24, 0, 0},
    {24, 24, 0},
    {0, 24, 0},
    {0, 0, 24},
    {24, 0, 24},
    {24, 24, 24},
    {0, 24, 24},
    {8, 0, 0},
    {8, 8, 0},
    {0, 8, 0},
    {0, 0, 8},
    {8, 0, 8},
    {8, 8, 8},
    {0, 8, 8},
    {0, 0, 16},
    {8, 0, 16},
    {8, 8, 16},
    {0, 8, 16},
    {8, 0, 24},
    {8, 8, 24},
    {0, 8, 24},
    {8, 16, 0},
    {0, 16, 0},
    {8, 16, 8},
    {0, 16, 8},
    {8, 16, 16},
    {0, 16, 16},
    {8, 16, 24},
    {0, 16, 24},
    {8, 24, 0},
    {8, 24, 8},
    {0, 24, 8},
    {8, 24, 16},
    {0, 24, 16},
    {8, 24, 24},
    {16, 0, 0},
    {16, 8, 0},
    {16, 0, 8},
    {16, 8, 8},
    {16, 0, 16},
    {16, 8, 16},
    {16, 0, 24},
    {16, 8, 24},
    {16, 16, 0},
    {16, 16, 8},
    {16, 16, 16},
    {16, 16, 24},
    {16, 24, 0},
    {16, 24, 8},
    {16, 24, 16},
    {16, 24, 24},
    {24, 8, 0},
    {24, 0, 8},
    {24, 8, 8},
    {24, 0, 16},
    {24, 8, 16},
    {24, 8, 24},
    {24, 16, 0},
    {24, 16, 8},
    {24, 16, 16},
    {24, 16, 24},
    {24, 24, 8},
    {24, 24, 16},
};

static const HexConnectivity T11_CHILDREN[27] = {
    {0, 8, 9, 10, 11, 12, 13, 14},
    {11, 12, 13, 14, 15, 16, 17, 18},
    {15, 16, 17, 18, 4, 19, 20, 21},
    {10, 9, 22, 23, 14, 13, 24, 25},
    {14, 13, 24, 25, 18, 17, 26, 27},
    {18, 17, 26, 27, 21, 20, 28, 29},
    {23, 22, 30, 3, 25, 24, 31, 32},
    {25, 24, 31, 32, 27, 26, 33, 34},
    {27, 26, 33, 34, 29, 28, 35, 7},
    {8, 36, 37, 9, 12, 38, 39, 13},
    {12, 38, 39, 13, 16, 40, 41, 17},
    {16, 40, 41, 17, 19, 42, 43, 20},
    {9, 37, 44, 22, 13, 39, 45, 24},
    {13, 39, 45, 24, 17, 41, 46, 26},
    {17, 41, 46, 26, 20, 43, 47, 28},
    {22, 44, 48, 30, 24, 45, 49, 31},
    {24, 45, 49, 31, 26, 46, 50, 33},
    {26, 46, 50, 33, 28, 47, 51, 35},
    {36, 1, 52, 37, 38, 53, 54, 39},
    {38, 53, 54, 39, 40, 55, 56, 41},
    {40, 55, 56, 41, 42, 5, 57, 43},
    {37, 52, 58, 44, 39, 54, 59, 45},
    {39, 54, 59, 45, 41, 56, 60, 46},
    {41, 56, 60, 46, 43, 57, 61, 47},
    {44, 58, 2, 48, 45, 59, 62, 49},
    {45, 59, 62, 49, 46, 60, 63, 50},
    {46, 60, 63, 50, 47, 61, 6, 51},
};

static const CanonicalTemplate CANONICAL_TEMPLATES[11] = {
    {1, 5, 16, T1_VERTICES, T1_CHILDREN},
    {2, 3, 16, T2_VERTICES, T2_CHILDREN},
    {3, 4, 16, T3_VERTICES, T3_CHILDREN},
    {4, 5, 20, T4_VERTICES, T4_CHILDREN},
    {5, 13, 32, T5_VERTICES, T5_CHILDREN},
    {6, 10, 32, T6_VERTICES, T6_CHILDREN},
    {7, 7, 22, T7_VERTICES, T7_CHILDREN},
    {8, 9, 32, T8_VERTICES, T8_CHILDREN},
    {9, 31, 64, T9_VERTICES, T9_CHILDREN},
    {10, 17, 40, T10_VERTICES, T10_CHILDREN},
    {11, 27, 64, T11_VERTICES, T11_CHILDREN},
};

static const PadOrientation PAD_TABLE[135] = {
    {10, 1, {0, 0, 0}, {0, 0, 0}, {0, 1, 2, 3, 4, 5, 6, 7}, false, 0},
    {11, 1, {0, 0, -1}, {0, 0, 0}, {1, 2, 3, 0, 5, 6, 7, 4}, false, 2},
    {12, 1, {0, 0, 0}, {0, 1, 0}, {3, 2, 1, 0, 7, 6, 5, 4}, true, 4},
    {13, 1, {0, 0, 1}, {0, 0, 0}, {3, 0, 1, 2, 7, 4, 5, 6}, false, 8},
    {14, 1, {0, -1, 0}, {0, 0, 0}, {4, 0, 3, 7, 5, 1, 2, 6}, false, 16},
    {15, 1, {0, 1, 0}, {0, 0, 0}, {1, 5, 6, 2, 0, 4, 7, 3}, false, 32},
    {16, 1, {0, 1, 0}, {0, 1, 0}, {2, 6, 5, 1, 3, 7, 4, 0}, true, 64},
    {17, 1, {0, -1, 0}, {0, 1, 0}, {7, 3, 0, 4, 6, 2, 1, 5}, true, 128},
    {18, 1, {0, 0, 0}, {0, 0, 1}, {4, 5, 6, 7, 0, 1, 2, 3}, true, 256},
    {19, 1, {0, 0, -1}, {0, 0, 1}, {5, 6, 7, 4, 1, 2, 3, 0}, true, 512},
    {20, 1, {-1, 0, 0}, {0, 1, 0}, {7, 6, 2, 3, 4, 5, 1, 0}, true, 1024},
    {21, 1, {0, 0, 1}, {0, 0, 1}, {7, 4, 5, 6, 3, 0, 1, 2}, true, 2048},
    {22, 2, {0, 1, 0}, {0, 0, 0}, {1, 5, 6, 2, 0, 4, 7, 3}, false, 0},
    {23, 2, {1, 0, 0}, {0, 0, 0}, {4, 5, 1, 0, 7, 6, 2, 3}, false, 2570},
    {24, 2, {0, 0, 0}, {0, 0, 0}, {0, 1, 2, 3, 4, 5, 6, 7}, false, 240},
    {25, 3, {0, 0, 1}, {0, 0, 0}, {3, 0, 1, 2, 7, 4, 5, 6}, false, 0},
    {26, 3, {0, 0, 0}, {0, 0, 0}, {0, 1, 2, 3, 4, 5, 6, 7}, false, 10},
    {27, 3, {0, 0, 0}, {0, 0, 1}, {4, 5, 6, 7, 0, 1, 2, 3}, true, 1280},
    {28, 3, {0, 0, 1}, {0, 0, 1}, {7, 4, 5, 6, 3, 0, 1, 2}, true, 2560},
    {29, 3, {1, 0, 0}, {0, 0, 0}, {4, 5, 1, 0, 7, 6, 2, 3}, false, 48},
    {30, 3, {-1, 1, 0}, {0, 0, 0}, {2, 1, 5, 6, 3, 0, 4, 7}, false, 96},
    {31, 3, {-1, 0, 0}, {0, 0, 0}, {3, 2, 6, 7, 0, 1, 5, 4}, false, 192},
    {32, 3, {-1, -1, 0}, {0, 0, 0}, {0, 3, 7, 4, 1, 2, 6, 5}, false, 144},
    {33, 3, {1, 0, 1}, {0, 0, 0}, {0, 4, 5, 1, 3, 7, 6, 2}, false, 257},
    {34, 3, {0, 1, 0}, {0, 0, 0}, {1, 5, 6, 2, 0, 4, 7, 3}, false, 514},
    {35, 3, {-1, 0, 1}, {0, 0, 0}, {7, 3, 2, 6, 4, 0, 1, 5}, false, 1028},
    {36, 3, {0, -1, 0}, {0, 0, 0}, {4, 0, 3, 7, 5, 1, 2, 6}, false, 2056},
    {37, 4, {0, 0, 0}, {1, 0, 0}, {1, 0, 3, 2, 5, 4, 7, 6}, true, 0},
    {38, 4, {0, 0, -1}, {1, 0, 0}, {2, 1, 0, 3, 6, 5, 4, 7}, true, 1542},
    {39, 4, {0, 0, 1}, {0, 0, 0}, {3, 0, 1, 2, 7, 4, 5, 6}, false, 3084},
    {40, 4, {0, 0, 0}, {0, 0, 0}, {0, 1, 2, 3, 4, 5, 6, 7}, false, 2313},
    {41, 4, {1, 0, 0}, {0, 1, 0}, {0, 1, 5, 4, 3, 2, 6, 7}, true, 149},
    {42, 4, {1, 0, 0}, {1, 1, 0}, {1, 0, 4, 5, 2, 3, 7, 6}, false, 101},
    {43, 4, {1, 0, 0}, {1, 0, 0}, {5, 4, 0, 1, 6, 7, 3, 2}, true, 1376},
    {44, 4, {1, 0, 0}, {0, 0, 0}, {4, 5, 1, 0, 7, 6, 2, 3}, false, 1424},
    {45, 4, {0, 1, 0}, {0, 0, 0}, {1, 5, 6, 2, 0, 4, 7, 3}, false, 58},
    {46, 4, {0, 1, 0}, {0, 1, 0}, {2, 6, 5, 1, 3, 7, 4, 0}, true, 202},
    {47, 4, {0, 1, 0}, {1, 1, 0}, {6, 2, 1, 5, 7, 3, 0, 4}, false, 2752},
    {48, 4, {0, 1, 0}, {1, 0, 0}, {5, 1, 2, 6, 4, 0, 3, 7}, true, 2608},
    {49, 5, {0, 0, 0}, {0, 0, 1}, {4, 5, 6, 7, 0, 1, 2, 3}, true, 0},
    {50, 5, {0, 0, 0}, {0, 0, 0}, {0, 1, 2, 3, 4, 5, 6, 7}, false, 3840},
    {51, 5, {-1, 0, 0}, {0, 0, 0}, {3, 2, 6, 7, 0, 1, 5, 4}, false, 305},
    {52, 5, {0, -1, 0}, {0, 0, 0}, {4, 0, 3, 7, 5, 1, 2, 6}, false, 610},
    {53, 5, {1, 0, 0}, {0, 0, 0}, {4, 5, 1, 0, 7, 6, 2, 3}, false, 1220},
    {54, 5, {0, 1, 0}, {0, 0, 0}, {1, 5, 6, 2, 0, 4, 7, 3}, false, 2200},
    {55, 5, {0, 0, 0}, {0, 0, 1}, {4, 5, 6, 7, 0, 1, 2, 3}, true, 11},
    {56, 5, {0, 0, 0}, {0, 0, 1}, {4, 5, 6, 7, 0, 1, 2, 3}, true, 7},
    {57, 5, {0, 0, 0}, {0, 0, 1}, {4, 5, 6, 7, 0, 1, 2, 3}, true, 14},
    {58, 5, {0, 0, 0}, {0, 0, 1}, {4, 5, 6, 7, 0, 1, 2, 3}, true, 13},
    {59, 5, {0, 0, 0}, {0, 0, 0}, {0, 1, 2, 3, 4, 5, 6, 7}, false, 1792},
    {60, 5, {0, 0, 0}, {0, 0, 0}, {0, 1, 2, 3, 4, 5, 6, 7}, false, 2816},
    {61, 5, {0, 0, 0}, {0, 0, 0}, {0, 1, 2, 3, 4, 5, 6, 7}, false, 3584},
    {62, 5, {0, 0, 0}, {0, 0, 0}, {0, 1, 2, 3, 4, 5, 6, 7}, false, 3328},
    {63, 5, {-1, 0, 0}, {0, 0, 0}, {3, 2, 6, 7, 0, 1, 5, 4}, false, 49},
    {64, 5, {-1, 0, 0}, {0, 0, 0}, {3, 2, 6, 7, 0, 1, 5, 4}, false, 289},
    {65, 5, {-1, 0, 0}, {0, 0, 0}, {3, 2, 6, 7, 0, 1, 5, 4}, false, 304},
    {66, 5, {-1, 0, 0}, {0, 0, 0}, {3, 2, 6, 7, 0, 1, 5, 4}, false, 273},
    {67, 5, {0, -1, 0}, {0, 0, 0}, {4, 0, 3, 7, 5, 1, 2, 6}, false, 98},
    {68, 5, {0, -1, 0}, {0, 0, 0}, {4, 0, 3, 7, 5, 1, 2, 6}, false, 546},
    {69, 5, {0, -1, 0}, {0, 0, 0}, {4, 0, 3, 7, 5, 1, 2, 6}, false, 608},
    {70, 5, {0, -1, 0}, {0, 0, 0}, {4, 0, 3, 7, 5, 1, 2, 6}, false, 578},
    {71, 5, {1, 0, 0}, {0, 0, 0}, {4, 5, 1, 0, 7, 6, 2, 3}, false, 196},
    {72, 5, {1, 0, 0}, {0, 0, 0}, {4, 5, 1, 0, 7, 6, 2, 3}, false, 1092},
    {73, 5, {1, 0, 0}, {0, 0, 0}, {4, 5, 1, 0, 7, 6, 2, 3}, false, 1216},
    {74, 5, {1, 0, 0}, {0, 0, 0}, {4, 5, 1, 0, 7, 6, 2, 3}, false, 1156},
    {75, 5, {0, 1, 0}, {0, 0, 0}, {1, 5, 6, 2, 0, 4, 7, 3}, false, 152},
    {76, 5, {0, 1, 0}, {0, 0, 0}, {1, 5, 6, 2, 0, 4, 7, 3}, false, 2072},
    {77, 5, {0, 1, 0}, {0, 0, 0}, {1, 5, 6, 2, 0, 4, 7, 3}, false, 2192},
    {78, 5, {0, 1, 0}, {0, 0, 0}, {1, 5, 6, 2, 0, 4, 7, 3}, false, 2184},
    {79, 5, {0, 0, 0}, {0, 0, 1}, {4, 5, 6, 7, 0, 1, 2, 3}, true, 3},
    {80, 5, {0, 0, 0}, {0, 0, 1}, {4, 5, 6, 7, 0, 1, 2, 3}, true, 6},
    {81, 5, {0, 0, 0}, {0, 0, 1}, {4, 5, 6, 7, 0, 1, 2, 3}, true, 12},
    {82, 5, {0, 0, 0}, {0, 0, 1}, {4, 5, 6, 7, 0, 1, 2, 3}, true, 9},
    {83, 5, {0, 0, 0}, {0, 0, 0}, {0, 1, 2, 3, 4, 5, 6, 7}, false, 768},
    {84, 5, {0, 0, 0}, {0, 0, 0}, {0, 1, 2, 3, 4, 5, 6, 7}, false, 1536},
    {85, 5, {0, 0, 0}, {0, 0, 0}, {0, 1, 2, 3, 4, 5, 6, 7}, false, 3072},
    {86, 5, {0, 0, 0}, {0, 0, 0}, {0, 1, 2, 3, 4, 5, 6, 7}, false, 2304},
    {87, 5, {-1, 0, 0}, {0, 0, 0}, {3, 2, 6, 7, 0, 1, 5, 4}, false, 17},
    {88, 5, {-1, 0, 0}, {0, 0, 0}, {3, 2, 6, 7, 0, 1, 5, 4}, false, 33},
    {89, 5, {-1, 0, 0}, {0, 0, 0}, {3, 2, 6, 7, 0, 1, 5, 4}, false, 288},
    {90, 5, {-1, 0, 0}, {0, 0, 0}, {3, 2, 6, 7, 0, 1, 5, 4}, false, 272},
    {91, 5, {0, -1, 0}, {0, 0, 0}, {4, 0, 3, 7, 5, 1, 2, 6}, false, 544},
    {92, 5, {0, -1, 0}, {0, 0, 0}, {4, 0, 3, 7, 5, 1, 2, 6}, false, 576},
    {93, 5, {0, -1, 0}, {0, 0, 0}, {4, 0, 3, 7, 5, 1, 2, 6}, false, 66},
    {94, 5, {0, -1, 0}, {0, 0, 0}, {4, 0, 3, 7, 5, 1, 2, 6}, false, 34},
    {95, 5, {1, 0, 0}, {0, 0, 0}, {4, 5, 1, 0, 7, 6, 2, 3}, false, 68},
    {96, 5, {1, 0, 0}, {0, 0, 0}, {4, 5, 1, 0, 7, 6, 2, 3}, false, 1088},
    {97, 5, {1, 0, 0}, {0, 0, 0}, {4, 5, 1, 0, 7, 6, 2, 3}, false, 1152},
    {98, 5, {1, 0, 0}, {0, 0, 0}, {4, 5, 1, 0, 7, 6, 2, 3}, false, 132},
    {99, 5, {0, 1, 0}, {0, 0, 0}, {1, 5, 6, 2, 0, 4, 7, 3}, false, 24},
    {100, 5, {0, 1, 0}, {0, 0, 0}, {1, 5, 6, 2, 0, 4, 7, 3}, false, 80},
    {101, 5, {0, 1, 0}, {0, 0, 0}, {1, 5, 6, 2, 0, 4, 7, 3}, false, 192},
    {102, 5, {0, 1, 0}, {0, 0, 0}, {1, 5, 6, 2, 0, 4, 7, 3}, false, 136},
    {103, 6, {0, 0, 1}, {0, 0, 0}, {3, 0, 1, 2, 7, 4, 5, 6}, false, 0},
    {104, 6, {0, 0, 0}, {0, 0, 0}, {0, 1, 2, 3, 4, 5, 6, 7}, false, 1799},
    {105, 6, {0, 0, -1}, {0, 0, 0}, {1, 2, 3, 0, 5, 6, 7, 4}, false, 3598},
    {106, 6, {0, 0, 0}, {1, 0, 0}, {1, 0, 3, 2, 5, 4, 7, 6}, true, 3341},
    {107, 6, {1, 0, 0}, {0, 0, 0}, {4, 5, 1, 0, 7, 6, 2, 3}, false, 1381},
    {108, 6, {1, 0, 1}, {0, 0, 0}, {0, 4, 5, 1, 3, 7, 6, 2}, false, 1520},
    {109, 6, {1, 0, 0}, {1, 0, 0}, {5, 4, 0, 1, 6, 7, 3, 2}, true, 1429},
    {110, 6, {1, 0, -1}, {0, 0, 0}, {5, 1, 0, 4, 6, 2, 3, 7}, false, 245},
    {111, 6, {0, 1, 0}, {1, 0, 0}, {5, 1, 2, 6, 4, 0, 3, 7}, true, 250},
    {112, 6, {1, 1, 0}, {1, 0, 0}, {6, 5, 1, 2, 7, 4, 0, 3}, true, 2618},
    {113, 6, {0, 1, 0}, {0, 0, 0}, {1, 5, 6, 2, 0, 4, 7, 3}, false, 2800},
    {114, 6, {1, 1, 0}, {0, 0, 0}, {5, 6, 2, 1, 4, 7, 3, 0}, false, 2762},
    {115, 7, {0, 0, 0}, {0, 0, 0}, {0, 1, 2, 3, 4, 5, 6, 7}, false, 0},
    {116, 7, {0, 0, -1}, {0, 0, 0}, {1, 2, 3, 0, 5, 6, 7, 4}, false, 35},
    {117, 7, {0, 0, -1}, {1, 0, 0}, {2, 1, 0, 3, 6, 5, 4, 7}, true, 70},
    {118, 7, {0, 0, 1}, {0, 0, 0}, {3, 0, 1, 2, 7, 4, 5, 6}, false, 140},
    {119, 7, {0, 0, 0}, {0, 0, 1}, {4, 5, 6, 7, 0, 1, 2, 3}, true, 2320},
    {120, 7, {0, 0, -1}, {0, 0, 1}, {5, 6, 7, 4, 1, 2, 3, 0}, true, 800},
    {121, 7, {0, 0, -1}, {1, 0, 1}, {6, 5, 4, 7, 2, 1, 0, 3}, false, 1600},
    {122, 7, {0, 0, 1}, {0, 0, 1}, {7, 4, 5, 6, 3, 0, 1, 2}, true, 3200},
    {123, 8, {0, 0, 0}, {0, 0, 0}, {0, 1, 2, 3, 4, 5, 6, 7}, false, 0},
    {124, 8, {1, 0, 0}, {0, 0, 0}, {4, 5, 1, 0, 7, 6, 2, 3}, false, 1525},
    {125, 8, {0, 1, 0}, {0, 0, 0}, {1, 5, 6, 2, 0, 4, 7, 3}, false, 2810},
    {126, 9, {0, 0, 0}, {0, 0, 1}, {4, 5, 6, 7, 0, 1, 2, 3}, true, 0},
    {127, 9, {0, 0, 0}, {0, 0, 0}, {0, 1, 2, 3, 4, 5, 6, 7}, false, 4080},
    {128, 9, {-1, 0, 0}, {0, 0, 0}, {3, 2, 6, 7, 0, 1, 5, 4}, false, 2875},
    {129, 9, {1, 0, 0}, {0, 0, 0}, {4, 5, 1, 0, 7, 6, 2, 3}, false, 3790},
    {130, 9, {0, -1, 0}, {0, 0, 0}, {4, 0, 3, 7, 5, 1, 2, 6}, false, 1895},
    {131, 9, {0, 1, 0}, {0, 0, 0}, {1, 5, 6, 2, 0, 4, 7, 3}, false, 3485},
    {132, 10, {0, 0, 0}, {0, 0, 0}, {0, 1, 2, 3, 4, 5, 6, 7}, false, 0},
    {133, 10, {0, 0, -1}, {0, 0, 0}, {1, 2, 3, 0, 5, 6, 7, 4}, false, 103},
    {134, 10, {0, 0, 0}, {0, 1, 0}, {3, 2, 1, 0, 7, 6, 5, 4}, true, 206},
    {135, 10, {0, 0, 1}, {0, 0, 0}, {3, 0, 1, 2, 7, 4, 5, 6}, false, 157},
    {136, 10, {0, -1, 0}, {0, 0, 0}, {4, 0, 3, 7, 5, 1, 2, 6}, false, 2329},
    {137, 10, {0, 1, 0}, {0, 0, 0}, {1, 5, 6, 2, 0, 4, 7, 3}, false, 803},
    {138, 10, {0, 1, 0}, {0, 1, 0}, {2, 6, 5, 1, 3, 7, 4, 0}, true, 1606},
    {139, 10, {0, -1, 0}, {0, 1, 0}, {7, 3, 0, 4, 6, 2, 1, 5}, true, 3212},
    {140, 10, {0, 0, 0}, {0, 0, 1}, {4, 5, 6, 7, 0, 1, 2, 3}, true, 2864},
    {141, 10, {0, 0, -1}, {0, 0, 1}, {5, 6, 7, 4, 1, 2, 3, 0}, true, 1888},
    {142, 10, {0, 0, 0}, {0, 1, 1}, {7, 6, 5, 4, 3, 2, 1, 0}, false, 3776},
    {143, 10, {0, 0, 1}, {0, 0, 1}, {7, 4, 5, 6, 3, 0, 1, 2}, true, 3472},
    {144, 11, {0, 0, 0}, {0, 0, 0}, {0, 1, 2, 3, 4, 5, 6, 7}, false, 4095},
};


const CanonicalTemplate& get_canonical_template(int template_id) {
    if (template_id < 1 || template_id > 11) {
        throw std::out_of_range("Invalid template_id: must be 1..11");
    }
    return CANONICAL_TEMPLATES[template_id - 1];
}

const PadOrientation* get_pad_orientation(int pad) {
    // Binary search in sorted PAD_TABLE
    int low = 0, high = static_cast<int>(sizeof(PAD_TABLE) / sizeof(PAD_TABLE[0])) - 1;
    while (low <= high) {
        int mid = low + (high - low) / 2;
        if (PAD_TABLE[mid].pad == pad) return &PAD_TABLE[mid];
        if (PAD_TABLE[mid].pad < pad) low = mid + 1;
        else high = mid - 1;
    }
    return nullptr;
}

const std::vector<int>& get_all_pad_codes() {
    static std::vector<int> codes;
    if (codes.empty()) {
        codes.reserve(sizeof(PAD_TABLE) / sizeof(PAD_TABLE[0]));
        for (const auto &p : PAD_TABLE) codes.push_back(p.pad);
    }
    return codes;
}

std::array<int, 8> rotate_hex(const int rot[3], const int sym[3]) {
    std::array<int, 8> order = {0, 1, 2, 3, 4, 5, 6, 7};
    std::array<int, 8> aux;

    if (rot[0] == -1) {
        aux = order;
        order = {aux[3], aux[2], aux[6], aux[7], aux[0], aux[1], aux[5], aux[4]};
    }
    if (rot[0] == 1) {
        aux = order;
        order = {aux[4], aux[5], aux[1], aux[0], aux[7], aux[6], aux[2], aux[3]};
    }
    if (rot[1] == 1) {
        aux = order;
        order = {aux[1], aux[5], aux[6], aux[2], aux[0], aux[4], aux[7], aux[3]};
    }
    if (rot[1] == -1) {
        aux = order;
        order = {aux[4], aux[0], aux[3], aux[7], aux[5], aux[1], aux[2], aux[6]};
    }
    if (rot[2] == 1) {
        aux = order;
        order = {aux[3], aux[0], aux[1], aux[2], aux[7], aux[4], aux[5], aux[6]};
    }
    if (rot[2] == -1) {
        aux = order;
        order = {aux[1], aux[2], aux[3], aux[0], aux[5], aux[6], aux[7], aux[4]};
    }
    if (sym[0] == 1) {
        aux = order;
        order = {aux[1], aux[0], aux[3], aux[2], aux[5], aux[4], aux[7], aux[6]};
    }
    if (sym[1] == 1) {
        aux = order;
        order = {aux[3], aux[2], aux[1], aux[0], aux[7], aux[6], aux[5], aux[4]};
    }
    if (sym[2] == 1) {
        aux = order;
        order = {aux[4], aux[5], aux[6], aux[7], aux[0], aux[1], aux[2], aux[3]};
    }
    return order;
}

Point3 trilinear_map(const Point3 &ref_pt, const HexCorners &parent) {
    const double u = ref_pt.x, v = ref_pt.y, w = ref_pt.z;
    const double N[8] = {
        (1.0 - u) * (1.0 - v) * (1.0 - w),
        u         * (1.0 - v) * (1.0 - w),
        u         * v         * (1.0 - w),
        (1.0 - u) * v         * (1.0 - w),
        (1.0 - u) * (1.0 - v) * w,
        u         * (1.0 - v) * w,
        u         * v         * w,
        (1.0 - u) * v         * w
    };
    Point3 pt = {0.0, 0.0, 0.0};
    for (int i = 0; i < 8; ++i) {
        pt.x += N[i] * parent[i].x;
        pt.y += N[i] * parent[i].y;
        pt.z += N[i] * parent[i].z;
    }
    return pt;
}

Point3 trilinear_map(const RationalCoord3 &ref_pt, const HexCorners &parent) {
    Point3 p_double = {ref_pt.x_double(), ref_pt.y_double(), ref_pt.z_double()};
    return trilinear_map(p_double, parent);
}

double hex_jacobian_at(const HexCorners &hex, double u, double v, double w) {
    const double dNu[8] = {
        -(1.0 - v) * (1.0 - w),  (1.0 - v) * (1.0 - w),  v * (1.0 - w), -v * (1.0 - w),
        -(1.0 - v) * w,          (1.0 - v) * w,          v * w,         -v * w
    };
    const double dNv[8] = {
        -(1.0 - u) * (1.0 - w), -u * (1.0 - w),          u * (1.0 - w),  (1.0 - u) * (1.0 - w),
        -(1.0 - u) * w,         -u * w,                  u * w,          (1.0 - u) * w
    };
    const double dNw[8] = {
        -(1.0 - u) * (1.0 - v), -u * (1.0 - v),         -u * v,         -(1.0 - u) * v,
         (1.0 - u) * (1.0 - v),  u * (1.0 - v),          u * v,          (1.0 - u) * v
    };

    double J[3][3] = {{0}};
    for (int i = 0; i < 8; ++i) {
        J[0][0] += dNu[i] * hex[i].x; J[0][1] += dNu[i] * hex[i].y; J[0][2] += dNu[i] * hex[i].z;
        J[1][0] += dNv[i] * hex[i].x; J[1][1] += dNv[i] * hex[i].y; J[1][2] += dNv[i] * hex[i].z;
        J[2][0] += dNw[i] * hex[i].x; J[2][1] += dNw[i] * hex[i].y; J[2][2] += dNw[i] * hex[i].z;
    }

    const double detJ = J[0][0] * (J[1][1] * J[2][2] - J[1][2] * J[2][1])
                      - J[0][1] * (J[1][0] * J[2][2] - J[1][2] * J[2][0])
                      + J[0][2] * (J[1][0] * J[2][1] - J[1][1] * J[2][0]);
    return detJ;
}

void hex_volume_and_jacobian_bounds(const HexCorners &hex, double &volume, double &min_jac) {
    const double a = 0.5 - 0.5 / std::sqrt(3.0), b = 1.0 - a;
    const double gp[2] = {a, b};
    volume = 0.0;
    min_jac = 1e18;
    for (double u : gp) {
        for (double v : gp) {
            for (double w : gp) {
                const double detJ = hex_jacobian_at(hex, u, v, w);
                volume += detJ / 8.0;
                if (detJ < min_jac) min_jac = detJ;
            }
        }
    }
}

std::vector<HexCorners> generate_children_for_pad(int pad, const HexCorners &parent) {
    const PadOrientation *po = get_pad_orientation(pad);
    if (!po) {
        throw std::invalid_argument("Unknown pad code: " + std::to_string(pad));
    }
    const CanonicalTemplate &ct = get_canonical_template(po->template_id);

    // Permute parent corners
    HexCorners permuted_parent;
    for (int i = 0; i < 8; ++i) {
        permuted_parent[i] = parent[po->parent_order[i]];
    }

    // Evaluate mapped vertices
    std::vector<Point3> mapped_verts(ct.num_unique_vertices);
    for (int i = 0; i < ct.num_unique_vertices; ++i) {
        mapped_verts[i] = trilinear_map(ct.vertices[i], permuted_parent);
    }

    // Build children
    std::vector<HexCorners> result(ct.num_children);
    for (int c = 0; c < ct.num_children; ++c) {
        const auto &conn = ct.children[c];
        HexCorners child;
        for (int v = 0; v < 8; ++v) {
            child[v] = mapped_verts[conn[v]];
        }
        if (po->flip_parity) {
            HexCorners flipped;
            for (int v = 0; v < 8; ++v) {
                flipped[v] = child[FLIP_CORNER_PERM[v]];
            }
            child = flipped;
        }
        result[c] = child;
    }
    return result;
}

static std::array<RationalCoord2, 4> extract_face_quad_2d(const std::array<RationalCoord3, 4> &quad3d, int face_idx) {
    std::array<RationalCoord2, 4> q2d;
    for (int i = 0; i < 4; ++i) {
        const auto &p = quad3d[i];
        switch (face_idx) {
            case 0: // bottom: -z (z=0), coordinates (u, v) = (x, y)
                q2d[i] = {p.x, p.y}; break;
            case 1: // top: +z (z=24), coordinates (u, v) = (x, y)
                q2d[i] = {p.x, p.y}; break;
            case 2: // front: -y (y=0), coordinates (u, v) = (x, z)
                q2d[i] = {p.x, p.z}; break;
            case 3: // right: +x (x=24), coordinates (u, v) = (y, z)
                q2d[i] = {p.y, p.z}; break;
            case 4: // back: +y (y=24), coordinates (u, v) = (x, z)
                q2d[i] = {p.x, p.z}; break;
            case 5: // left: -x (x=0), coordinates (u, v) = (y, z)
                q2d[i] = {p.y, p.z}; break;
            default:
                throw std::out_of_range("Invalid face_idx");
        }
    }
    // Canonicalize quad: cyclically shift so smallest vertex is first
    auto min_it = std::min_element(q2d.begin(), q2d.end());
    int shift = std::distance(q2d.begin(), min_it);
    std::array<RationalCoord2, 4> canonical;
    for (int i = 0; i < 4; ++i) canonical[i] = q2d[(i + shift) % 4];
    return canonical;
}

FaceSignature compute_face_signature(int template_id, int face_idx) {
    const CanonicalTemplate &ct = get_canonical_template(template_id);
    FaceSignature sig;
    sig.face_idx = face_idx;

    // Collect all boundary quads lying on face_idx
    // First map all internal and boundary faces of the canonical template
    using QuadKey = std::array<RationalCoord3, 4>;
    auto canonical_key = [](QuadKey q) {
        // cycle to minimum
        auto min_it = std::min_element(q.begin(), q.end());
        int shift = std::distance(q.begin(), min_it);
        QuadKey c1;
        for (int i = 0; i < 4; ++i) c1[i] = q[(i + shift) % 4];
        // reversed cycle
        QuadKey rev = {q[3], q[2], q[1], q[0]};
        auto rev_min = std::min_element(rev.begin(), rev.end());
        int rev_shift = std::distance(rev.begin(), rev_min);
        QuadKey c2;
        for (int i = 0; i < 4; ++i) c2[i] = rev[(i + rev_shift) % 4];
        return std::min(c1, c2);
    };

    std::map<QuadKey, std::vector<QuadKey>> face_counts;
    for (int c = 0; c < ct.num_children; ++c) {
        for (int f = 0; f < 6; ++f) {
            QuadKey q;
            for (int v = 0; v < 4; ++v) {
                q[v] = ct.vertices[ct.children[c][HEX_FACES[f][v]]];
            }
            face_counts[canonical_key(q)].push_back(q);
        }
    }

    for (const auto &pair : face_counts) {
        if (pair.second.size() == 1) { // boundary face
            const auto &q = pair.second[0];
            bool on_face = false;
            switch (face_idx) {
                case 0: on_face = (q[0].z == 0  && q[1].z == 0  && q[2].z == 0  && q[3].z == 0); break;
                case 1: on_face = (q[0].z == 24 && q[1].z == 24 && q[2].z == 24 && q[3].z == 24); break;
                case 2: on_face = (q[0].y == 0  && q[1].y == 0  && q[2].y == 0  && q[3].y == 0); break;
                case 3: on_face = (q[0].x == 24 && q[1].x == 24 && q[2].x == 24 && q[3].x == 24); break;
                case 4: on_face = (q[0].y == 24 && q[1].y == 24 && q[2].y == 24 && q[3].y == 24); break;
                case 5: on_face = (q[0].x == 0  && q[1].x == 0  && q[2].x == 0  && q[3].x == 0); break;
            }
            if (on_face) {
                sig.quads.push_back(extract_face_quad_2d(q, face_idx));
            }
        }
    }
    std::sort(sig.quads.begin(), sig.quads.end());
    return sig;
}

FaceSignature compute_face_signature_for_pad(int pad, int face_idx) {
    const PadOrientation *po = get_pad_orientation(pad);
    if (!po) {
        throw std::invalid_argument("Unknown pad code: " + std::to_string(pad));
    }
    const CanonicalTemplate &ct = get_canonical_template(po->template_id);

    // Build unit cube parent
    HexCorners unit_parent = {
        Point3{0.0, 0.0, 0.0}, Point3{1.0, 0.0, 0.0}, Point3{1.0, 1.0, 0.0}, Point3{0.0, 1.0, 0.0},
        Point3{0.0, 0.0, 1.0}, Point3{1.0, 0.0, 1.0}, Point3{1.0, 1.0, 1.0}, Point3{0.0, 1.0, 1.0}
    };
    auto children = generate_children_for_pad(pad, unit_parent);

    // Convert children corners to RationalCoord3
    std::vector<std::array<RationalCoord3, 8>> rat_children(children.size());
    for (size_t c = 0; c < children.size(); ++c) {
        for (int v = 0; v < 8; ++v) {
            rat_children[c][v] = {
                static_cast<int8_t>(std::round(children[c][v].x * 24.0)),
                static_cast<int8_t>(std::round(children[c][v].y * 24.0)),
                static_cast<int8_t>(std::round(children[c][v].z * 24.0))
            };
        }
    }

    using QuadKey = std::array<RationalCoord3, 4>;
    auto canonical_key = [](QuadKey q) {
        auto min_it = std::min_element(q.begin(), q.end());
        int shift = std::distance(q.begin(), min_it);
        QuadKey c1;
        for (int i = 0; i < 4; ++i) c1[i] = q[(i + shift) % 4];
        QuadKey rev = {q[3], q[2], q[1], q[0]};
        auto rev_min = std::min_element(rev.begin(), rev.end());
        int rev_shift = std::distance(rev.begin(), rev_min);
        QuadKey c2;
        for (int i = 0; i < 4; ++i) c2[i] = rev[(i + rev_shift) % 4];
        return std::min(c1, c2);
    };

    std::map<QuadKey, std::vector<QuadKey>> face_counts;
    for (size_t c = 0; c < rat_children.size(); ++c) {
        for (int f = 0; f < 6; ++f) {
            QuadKey q;
            for (int v = 0; v < 4; ++v) {
                q[v] = rat_children[c][HEX_FACES[f][v]];
            }
            face_counts[canonical_key(q)].push_back(q);
        }
    }

    FaceSignature sig;
    sig.face_idx = face_idx;
    for (const auto &pair : face_counts) {
        if (pair.second.size() == 1) {
            const auto &q = pair.second[0];
            bool on_face = false;
            switch (face_idx) {
                case 0: on_face = (q[0].z == 0  && q[1].z == 0  && q[2].z == 0  && q[3].z == 0); break;
                case 1: on_face = (q[0].z == 24 && q[1].z == 24 && q[2].z == 24 && q[3].z == 24); break;
                case 2: on_face = (q[0].y == 0  && q[1].y == 0  && q[2].y == 0  && q[3].y == 0); break;
                case 3: on_face = (q[0].x == 24 && q[1].x == 24 && q[2].x == 24 && q[3].x == 24); break;
                case 4: on_face = (q[0].y == 24 && q[1].y == 24 && q[2].y == 24 && q[3].y == 24); break;
                case 5: on_face = (q[0].x == 0  && q[1].x == 0  && q[2].x == 0  && q[3].x == 0); break;
            }
            if (on_face) {
                sig.quads.push_back(extract_face_quad_2d(q, face_idx));
            }
        }
    }
    std::sort(sig.quads.begin(), sig.quads.end());
    return sig;
}

TemplateValidationReport validate_template_pad(int pad, const HexCorners &parent) {
    TemplateValidationReport report;
    const PadOrientation *po = get_pad_orientation(pad);
    if (!po) {
        report.error_message = "Unknown pad: " + std::to_string(pad);
        return report;
    }

    std::vector<HexCorners> children;
    try {
        children = generate_children_for_pad(pad, parent);
    } catch (const std::exception &e) {
        report.error_message = e.what();
        return report;
    }

    report.child_count = static_cast<int>(children.size());
    double parent_vol = 0.0, parent_min_j = 0.0;
    hex_volume_and_jacobian_bounds(parent, parent_vol, parent_min_j);

    report.total_volume = 0.0;
    report.min_jacobian = 1e18;

    for (const auto &child : children) {
        double vol = 0.0, min_j = 0.0;
        hex_volume_and_jacobian_bounds(child, vol, min_j);
        report.total_volume += vol;
        if (min_j < report.min_jacobian) report.min_jacobian = min_j;
        if (min_j <= 0.0) {
            report.error_message = "Child has non-positive Jacobian: " + std::to_string(min_j);
            return report;
        }
    }

    if (std::abs(report.total_volume - parent_vol) > 1e-6 * std::abs(parent_vol)) {
        report.error_message = "Volume mismatch: total child vol=" + std::to_string(report.total_volume) +
                               " vs parent vol=" + std::to_string(parent_vol);
        return report;
    }

    report.valid = true;
    return report;
}

} // namespace refinement
