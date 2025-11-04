#include <stdio.h>
#include <stdlib.h>
#include <math.h>

// Function to compute 3x3 matrix determinant
double det3x3(double J[3][3]) {
    return J[0][0] * (J[1][1] * J[2][2] - J[1][2] * J[2][1]) -
           J[0][1] * (J[1][0] * J[2][2] - J[1][2] * J[2][0]) +
           J[0][2] * (J[1][0] * J[2][1] - J[1][1] * J[2][0]);
}

// Function to approximate condition number (max/min singular values) using Frobenius norm
double conditionNumber(double J[3][3]) {
    double frobenius = 0.0;
    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++)
            frobenius += J[i][j] * J[i][j];
    frobenius = sqrt(frobenius);

    double det = det3x3(J);
    if (fabs(det) < 1e-10) {
        printf("Warning: Singular Jacobian matrix.\n");
        return 1e10; // Large value for degenerate case
    }

    double invJ[3][3];
    invJ[0][0] = (J[1][1] * J[2][2] - J[1][2] * J[2][1]) / det;
    invJ[0][1] = -(J[0][1] * J[2][2] - J[0][2] * J[2][1]) / det;
    invJ[0][2] = (J[0][1] * J[1][2] - J[0][2] * J[1][1]) / det;
    invJ[1][0] = -(J[1][0] * J[2][2] - J[1][2] * J[2][0]) / det;
    invJ[1][1] = (J[0][0] * J[2][2] - J[0][2] * J[2][0]) / det;
    invJ[1][2] = -(J[0][0] * J[1][2] - J[0][2] * J[1][0]) / det;
    invJ[2][0] = (J[1][0] * J[2][1] - J[1][1] * J[2][0]) / det;
    invJ[2][1] = -(J[0][0] * J[2][1] - J[0][1] * J[2][0]) / det;
    invJ[2][2] = (J[0][0] * J[1][1] - J[0][1] * J[1][0]) / det;

    double frobeniusInv = 0.0;
    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++)
            frobeniusInv += invJ[i][j] * invJ[i][j];
    frobeniusInv = sqrt(frobeniusInv);

    return frobenius * frobeniusInv / 3.0; // Normalized for 3x3 matrix
}

// Function to compute volume using tetrahedral decomposition
double computeHexVolume(double nodes[8][3]) {
    int tetrahedra[5][4] = {
        {0, 1, 3, 4}, // Tet 1: nodes 0,1,3,4
        {1, 2, 3, 6}, // Tet 2: nodes 1,2,3,6
        {1, 4, 5, 6}, // Tet 3: nodes 1,4,5,6
        {3, 4, 6, 7}, // Tet 4: nodes 3,4,6,7
        {1, 3, 4, 6}  // Tet 5: nodes 1,3,4,6 (central)
    };

    double volume = 0.0;
    for (int i = 0; i < 5; i++) {
        double *n0 = nodes[tetrahedra[i][0]];
        double *n1 = nodes[tetrahedra[i][1]];
        double *n2 = nodes[tetrahedra[i][2]];
        double *n3 = nodes[tetrahedra[i][3]];

        double v1[3] = {n1[0] - n0[0], n1[1] - n0[1], n1[2] - n0[2]};
        double v2[3] = {n2[0] - n0[0], n2[1] - n0[1], n2[2] - n0[2]};
        double v3[3] = {n3[0] - n0[0], n3[1] - n0[1], n3[2] - n0[2]};

        double det = v1[0] * (v2[1] * v3[2] - v2[2] * v3[1]) -
                     v1[1] * (v2[0] * v3[2] - v2[2] * v3[0]) +
                     v1[2] * (v2[0] * v3[1] - v2[1] * v3[0]);

        volume += fabs(det) / 6.0;
    }
    return volume;
}

// Function to compute Jacobian determinant and skewness at centroid
void computeJacobianAndSkew(double nodes[8][3], double *jacobianDet, double *skew) {
    double xi = 0.0, eta = 0.0, zeta = 0.0;

    double dNdxi[8][3] = {
        {-0.125 * (1-eta) * (1-zeta), -0.125 * (1-xi) * (1-zeta), -0.125 * (1-xi) * (1-eta)}, // dN0
        { 0.125 * (1-eta) * (1-zeta), -0.125 * (1+xi) * (1-zeta), -0.125 * (1+xi) * (1-eta)}, // dN1
        { 0.125 * (1+eta) * (1-zeta),  0.125 * (1+xi) * (1-zeta), -0.125 * (1+xi) * (1+eta)}, // dN2
        {-0.125 * (1+eta) * (1-zeta),  0.125 * (1-xi) * (1-zeta), -0.125 * (1-xi) * (1+eta)}, // dN3
        {-0.125 * (1-eta) * (1+zeta), -0.125 * (1-xi) * (1+zeta),  0.125 * (1-xi) * (1-eta)}, // dN4
        { 0.125 * (1-eta) * (1+zeta), -0.125 * (1+xi) * (1+zeta),  0.125 * (1+xi) * (1-eta)}, // dN5
        { 0.125 * (1+eta) * (1+zeta),  0.125 * (1+xi) * (1+zeta),  0.125 * (1+xi) * (1+eta)}, // dN6
        {-0.125 * (1+eta) * (1+zeta),  0.125 * (1-xi) * (1+zeta),  0.125 * (1-xi) * (1+eta)}  // dN7
    };

    double J[3][3] = {{0.0}};
    for (int i = 0; i < 8; i++) {
        J[0][0] += dNdxi[i][0] * nodes[i][0]; // dx/dxi
        J[0][1] += dNdxi[i][1] * nodes[i][0]; // dx/deta
        J[0][2] += dNdxi[i][2] * nodes[i][0]; // dx/dzeta
        J[1][0] += dNdxi[i][0] * nodes[i][1]; // dy/dxi
        J[1][1] += dNdxi[i][1] * nodes[i][1]; // dy/deta
        J[1][2] += dNdxi[i][2] * nodes[i][1]; // dy/dzeta
        J[2][0] += dNdxi[i][0] * nodes[i][2]; // dz/dxi
        J[2][1] += dNdxi[i][1] * nodes[i][2]; // dz/deta
        J[2][2] += dNdxi[i][2] * nodes[i][2]; // dz/dzeta
    }

    *jacobianDet = det3x3(J);
    *skew = conditionNumber(J);

    if (*jacobianDet <= 0.0) {
        printf("Warning: Jacobian determinant is non-positive (%.4f), indicating an invalid element.\n", *jacobianDet);
    }
}

// Main function to compute volume, Jacobian determinant, and skewness
void hexElementMetrics(double nodes[8][3], double *volume, double *jacobianDet, double *skew) {
    if (!nodes || !volume || !jacobianDet || !skew) {
        printf("Error: Null pointer passed to hexElementMetrics.\n");
        return;
    }

    *volume = computeHexVolume(nodes);
    computeJacobianAndSkew(nodes, jacobianDet, skew);
}

// Example usage
int test() {
    double nodes[8][3] = {
        {0.5, 0.5, 0.0}, // Node 0
        {1.0, 0.0, 0.0}, // Node 1
        {1.0, 1.0, 0.0}, // Node 2
        {0.0, 1.0, 0.0}, // Node 3
        {0.0, 0.0, 1.0}, // Node 4
        {1.0, 0.0, 1.0}, // Node 5
        {1.0, 1.0, 1.0}, // Node 6
        {0.0, 1.0, 1.0}  // Node 7
    };

    double volume, jacobianDet, skew;
    hexElementMetrics(nodes, &volume, &jacobianDet, &skew);

    printf("Volume: %.4f\n", volume);
    printf("Jacobian Determinant: %.4f\n", jacobianDet);
    printf("Skewness: %.4f\n", skew);

    return 0;
}
