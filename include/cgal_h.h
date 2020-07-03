
#ifndef CGAL_H_
//#define CGAL_H_

//	--- CGAL kernel and common headers
//#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Cartesian_converter.h>
#include <CGAL/Intersections.h>
#include <CGAL/result_of.h>

#include <vector>
#include <tuple>
#include <iostream>
#include <fstream>

// --- Boost/variant and functional
#include <boost/variant/apply_visitor.hpp>
#include <functional>

// --- CGAL polygon and polyhedron
#include <CGAL/Polygon_2.h>
#include <CGAL/Polyhedron_3.h>

// --- CGAL Nef polyhedrons
#include <CGAL/Nef_polyhedron_3.h>
#include <CGAL/Nef_3/SNC_indexed_items.h>

// --- CGAL 3D convex hull and polyhedron decomposition
#include <CGAL/convex_decomposition_3.h>
#include <CGAL/convex_hull_3.h>

// --- CGAL 3D random point generator
#include <CGAL/point_generators_3.h>
/*
 * 		Kernel typedefs
 */

// --- Kernel types
typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef CGAL::Exact_predicates_exact_constructions_kernel ExactKernel;

// --- Converters between exact and inexact kernels
typedef CGAL::Cartesian_converter<Kernel,ExactKernel>	Kernel_to_ExactKernel;
typedef CGAL::Cartesian_converter<ExactKernel,Kernel>	ExactKernel_to_Kernel;

/*
 * 		Kernel-independent typedefs
 */

//typedef CGAL::Bbox_2			Bbox_2;
//typedef CGAL::Bbox_3			Bbox_3;

/*
 * 		EXACT kernel typedefs
 */

//typedef ExactKernel::Point_2	 	ExactPoint_2;
//typedef ExactKernel::Triangle_2		ExactTriangle_2;
//typedef ExactKernel::Segment_2		ExactSegment_2;
//typedef ExactKernel::Vector_2		ExactVector_2;

typedef ExactKernel::Point_3		ExactPoint_3;
typedef ExactKernel::Triangle_3 	ExactTriangle_3;
typedef ExactKernel::Segment_3		ExactSegment_3;
//typedef ExactKernel::Vector_3		ExactVector_3;
//typedef ExactKernel::Plane_3		ExactPlane_3;

//typedef CGAL::Polygon_2<ExactKernel>		ExactPolygon_2;
//typedef CGAL::Polyhedron_3<ExactKernel>		ExactPolyhedron_3;
//typedef CGAL::Tetrahedron_3<ExactKernel>	ExactTetrahedron_3;

//typedef CGAL::Nef_polyhedron_3<ExactKernel, CGAL::SNC_indexed_items>  	Nef_Polyhedron_3;
//typedef Nef_Polyhedron_3::Plane_3			Nef_Plane_3;

//typedef CGAL::Polyhedron_3<ExactKernel,CGAL::Polyhedron_items_3>
//														ExactPolyhedralSurface;

//typedef CGAL::Creator_uniform_3<double, ExactPoint_3>  PointCreator;

/*
 * 		EXACT kernel typedefs - don't even dream to do intersections with these
 */

//typedef Kernel::Point_2	 		Point_2;
//typedef Kernel::Triangle_2		Triangle_2;
//typedef Kernel::Segment_2		Segment_2;
//typedef Kernel::Vector_2		Vector_2;

typedef Kernel::Point_3			Point_3;
//typedef Kernel::Triangle_3  	Triangle_3;
//typedef Kernel::Segment_3		Segment_3;
//typedef Kernel::Vector_3		Vector_3;


//typedef CGAL::Polygon_2<Kernel>			Polygon_2;
//typedef CGAL::Polyhedron_3<Kernel>		Polyhedron_2;
//typedef CGAL::Tetrahedron_3<Kernel>		Tetrahedron_2;

// Intersection typedef (boost::variant)
typedef CGAL::cpp11::result_of<ExactKernel::Intersect_3(ExactTriangle_3, ExactSegment_3)>::type
		Triangle_3_Intersection_Variant;

#endif
