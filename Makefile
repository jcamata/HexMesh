## -*- Makefile -*-
##
## User: camata
## Time: Jul 6, 2015 2:23:03 PM
## Makefile created by Oracle Solaris Studio.
##
## This file is generated automatically.
##

###############################################################################
# File management.  This is where the source, header, and object files are
# defined

#
# source files
srcfiles 	:= $(wildcard src/*.cpp) main.cpp

#
# object files
objects		:= $(patsubst %.cpp, %.o, $(srcfiles))
###############################################################################

GTS_DIR :=/Users/lac/Documents/local
SC_DIR :=/Users/LAC/Codes/local/libsc-master/local
SC_DIR :=/Users/lac/Documents/local
MESQUITE_DIR :=/Users/lac/Documents/local
HDF5_DIR := /opt/local
#CGAL_LIB := -L/opt/local/lib  -Wl,-rpath,/opt/local/lib /opt/local/lib/libmpfr.dylib /opt/local/lib/libgmp.dylib /opt/local/lib/libCGAL_Core.13.0.3.dylib /opt/local/lib/libCGAL.13.0.3.dylib /opt/local/lib/libboost_filesystem-mt.dylib /opt/local/lib/libmpfr.dylib /opt/local/lib/libgmp.dylib

CGAL_LIB := -L/opt/local/lib  -Wl,-rpath,/opt/local/lib /opt/local/lib/libmpfr.dylib /opt/local/lib/libgmp.dylib /opt/local/lib/libboost_filesystem-mt.dylib /opt/local/lib/libmpfr.dylib

GLIB_INCLUDE = -I/opt/local/include/glib-2.0 -I/opt/local/lib/glib-2.0/include

CXX       = mpicxx -lstdc++ -O2 -std=c++17
#CXX       = mpicxx -lstdc++ -g -std=c++17
#CXX       = mpicxx -g -std=c++17 -fsanitize=address
LDFLAGS   =-L$(GTS_DIR)/lib -lgts -L$(SC_DIR)/lib -L$(HDF5_DIR)/lib -lsc -lm -lglib-2.0 -lhdf5_cpp -lhdf5 -L$(MESQUITE_DIR)/lib -lmesquite $(CGAL_LIB)
CXX_FLAGS = -I$(GTS_DIR)/include -I$(SC_DIR)/include $(GLIB_INCLUDE) -I./include -I$(HDF5_DIR)/include -I$(MESQUITE_DIR)/include
# OpenMP for the closeness pass of MeshOptimization and the GLL stability reports
# (optimize_mesh.cpp, stability.cpp; without it they run serially). Apple clang + Homebrew libomp; on Linux/GCC use -fopenmp for both.
OMP_DIR   = /opt/homebrew/opt/libomp
src/optimize_mesh.o src/stability.o: CXX_FLAGS += -Xpreprocessor -fopenmp -I$(OMP_DIR)/include
LDFLAGS  += -L$(OMP_DIR)/lib -lomp -Wl,-rpath,$(OMP_DIR)/lib

#CXX       = mpicxx -O2
#LDFLAGS   = -L/home/lucio/local -lsc -lgts -lm
#CXX_FLAGS = -I/home/camata/local/include -I/usr/include/glib-2.0 -I/usr/lib/x86_64-linux-gnu/glib-2.0/include -Iinclude


## Target: all
all: hexmesh

src/pillow_rebuild.o: include/pillow_geometry.h include/mesh_geom.h include/hexa.h
src/moving_nodes.o src/optimize_mesh.o src/adaptive_refinement.o main.o: include/pillow_geometry.h
$(objects): include/hexa.h
src/refinement_templates.o: include/refinement_templates.h include/pillow_geometry.h
src/conforming_refinement.o: include/conforming_refinement.h include/refinement_templates.h include/pillow_geometry.h
src/adaptive_refinement.o: include/adaptive_refinement.h include/conforming_refinement.h

.PHONY: test-pillow test-refinement-templates test-conforming-refinement
test-conforming-refinement: tests/test_conforming_refinement
	./tests/test_conforming_refinement

tests/test_conforming_refinement: tests/test_conforming_refinement.cpp src/conforming_refinement.o src/refinement_templates.o
	$(CXX) $(CXX_FLAGS) $^ -o $@
test-pillow: tests/test_pillow
	./tests/test_pillow

test-refinement-templates: tests/test_refinement_templates
	./tests/test_refinement_templates

tests/test_refinement_templates: tests/test_refinement_templates.cpp src/refinement_templates.o include/refinement_templates.h
	$(CXX) $(CXX_FLAGS) tests/test_refinement_templates.cpp src/refinement_templates.o -o $@

tests/test_pillow: tests/test_pillow.cpp src/pillow_rebuild.o include/pillow_geometry.h
	$(CXX) $(CXX_FLAGS) tests/test_pillow.cpp src/pillow_rebuild.o -o $@ $(LDFLAGS)


hexmesh: $(objects)
	$(CXX) $(objects) -o hexmesh $(LDFLAGS)

#
# How to compile C++
#
%.o : %.cpp
	@echo "Compiling C++ "$<"..."
	$(CXX) $(CXX_FLAGS) -c $< -o $@



#### Clean target deletes all generated files ####
clean: 
	rm $(objects)


# Enable dependency checking
.KEEP_STATE:
.KEEP_STATE_FILE:.make.state.GNU-amd64-Linux
