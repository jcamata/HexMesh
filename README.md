# HexMesh

## PRESENTATION

This project is focused on the development of a scalable hexahedral mesh generator for large domains based on octree and 27-tree structures. It allows for consideration of topography, bathymetry and coastlines, as well as water bodies and basins for geophysical applications.

The library is mainly developed at NACAD (Universidad Federal do Rio de Janeiro/COPPE, Brazil).

* contact : [José Camata](mailto:camata@nacad.ufrj.br)
* contributors (by order of first commit): L. A. Corrêa, R. Cottereau

It is written in C++/MPI. Additional routines are written in Matlab for the preparation of GTS topography files.
 
## REFERENCES

If you use the library, please cite the following paper:
1. J. Camata, A. Coutinho . Parallel implementation and performance analysis of a linear octree finite element mesh generation scheme, _ Concurrency and Computation: Practice and Experience _ (2013), pp. 826-842. (http://dx.doi.org/10.1002/cpe.2869)

## INSTALLATION

Before using the software, the following libraries need to be installed and available:
1. gts (http://gts.sourceforge.net/)
1. libsc (https://github.com/cburstedde/libsc)
1. hdf5 (https://portal.hdfgroup.org/display/support/Downloads)
1. mesquite (https://software.sandia.gov/mesquite/)
1. CGAL (https://www.cgal.org/), along with its GMP and MPFR dependencies

## COMPILATION

A C++17 compiler is required.

Depending on the OS you are using, modify the paths for GTS_LIB, SC_LIB, HDF5_DIR, MESQUITE_DIR, CGAL_LIB and GLIB_INCLUDE in Make.Linux, Make.mac or Makefile

Compile with (replace OS by Linux or mac)
>> make -f Make.OS

## USE

To prepare the geometry files, write a configuration file with the bounding box and the data sources, and run the Python preprocessor:

>> cd Preproc && python3 preproc.py mauna_loa.input

You need an internet connexion to download the elevation grid and the coastline file (no connexion needed if they are already available on your computer). It writes the topography and bathymetry GTS files directly — no stl2gts step — which should be moved to $(HEXHOME)/input, where $(HEXHOME) is the directory where hexmesh was compiled.

Configuration keys, the coordinate projection, how the bathymetry surface encodes the coastline for the material raytrace, and which elevation databases to use: [`Preproc/README.md`](Preproc/README.md).

The original Matlab workflow (`mainSRTM.m` in ./matlab, run from Matlab, writing STL files that then need `stl2gts`) still works but is superseded — note that its `lonlat2m.m` projection stretches longitude by up to ~27% at mid latitudes.

Before running, edit the input file HexMesh.input, located in $(HEXMESH) (it is read from ./HexMesh.input relative to the working directory). It sets:
* topo / inter: paths to the topography and interface (bathymetry) GTS files, and interfaceNumber (0 or 1 interface admitted for now)
* ref: refinement level of the 27-tree structure (mesh is a cube with 2*3^ref elements in x and y)
* z / zcuts: depth of the model and the depths at which refinement changes from 3 to 1 elements
* movingNodes: whether to move nodes in the smart octree (0/1)
* nmat and the material lines (S/F, with vp, vs, rho): material definitions
* PML, pmlx/pmly/pmlz, nlayersx/y/z, A, npow: PML settings
* meshOpt: mesh optimization flag (not yet implemented)
* CgalUse: use the exact CGAL kernel for intersection computations (1) instead of GTS (0)

To create the mesh, you should run (in a Terminal from the directory $(HEXMESH))

>> mpirun -np <nb_proc> ./hexmesh

where <nb_proc> is an integer specifying the number of processes used to create the mesh (each process creates its own output files). The files referenced by topo and inter in HexMesh.input should be in a repository ./input. By default (this will be made more general later), the depth of the mesh is the larger dimension of the two horizontal dimensions of the topography file; and the depths at which the refinements occur at set in function hexa_tree_cube, in hexa.cpp line 190).

## OUTPUT

Each MPI process writes its own mesh as an HDF5 file mesh_<nb_proc>_<rank>.h5, together with a companion mesh_<nb_proc>_<rank>.h5.xmf XDMF description that can be opened directly in ParaView or VisIt (VTK output is no longer produced; the VTK/VTU writer in hexa_vtk.cpp is disabled due to a known connectivity bug). Each process also writes a Profile_<nb_proc>_<rank>.txt log with timings for each stage. If PML = 1 in HexMesh.input, a material.input.FromHexMesh file is additionally written for use with SEM3D.
