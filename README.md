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

## COMPILATION

Depending on the OS you are using, modify the paths for GTS_LIB, SC_LIB, HDF5_DIR, MESQUITE_DIR and GLIB_INCLUDE in Make.Linux, Make.mac or Makefile

Compile with (replace OS by Linux or mac)
>> make -f Make.OS

## USE

To prepare the geometry files, modify the headers in mainSRTM.m (in particular choose the bounding box in latitude/longitude) and run in Matlab:

>> mainSRTM

You need an internet connexion to download the topography, bathymetry and coastlines files (no connexion needed if they are already available on your computer). The output files are a topography STL file topo.stl and a bathymetry STL file bathy.stl. These files should be transformed to GTS using stl2gts command, and move to directory $(HEXHOME)/input, where $(HEXHOME) is the directory where hexmesh was compiled.

To create the mesh, you should run (in a Terminal from the directory $(HEXMESH))

>> mpirun -np <nb_proc> ./hexmesh <refine_level>

where <nb_proc> is an integer specifying the number of processes used to create the mesh (each process creates its own VTK file), and <refine_level> is an integer specifying the number of level refinements of the 27-tree structure! (For now, the integer must be a power of 3: 1, 3, 9, 27 …). The files topo.gts, bathy.gts and coastline.dat should be in a repository ./input. By default (this will be made more general later), the depth of the mesh is the larger dimension of the two horizontal dimensions of the topography file; and the depths at which the refinements occur at set in function hexa_tree_cube, in hexa.cpp line 190).
