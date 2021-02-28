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

GTS_DIR :=/Users/LAC/Codes/local/gts/gts_lib
SC_DIR :=/Users/LAC/Codes/local/libsc-master/local
SC_DIR := /Users/LAC/Codes/local/p4est/libsc/local
MESQUITE_DIR := /Users/LAC/Codes/local
HDF5_DIR := /opt/local
CGAL_LIB := -L/opt/local/lib  -Wl,-rpath,/opt/local/lib /opt/local/lib/libmpfr.dylib /opt/local/lib/libgmp.dylib /opt/local/lib/libCGAL_Core.13.0.3.dylib /opt/local/lib/libCGAL.13.0.3.dylib /opt/local/lib/libboost_filesystem-mt.dylib /opt/local/lib/libmpfr.dylib /opt/local/lib/libgmp.dylib

GLIB_INCLUDE = -I/opt/local/include/glib-2.0 -I/opt/local/lib/glib-2.0/include

CXX       = mpicxx -lstdc++ -g
LDFLAGS   =-L$(GTS_DIR)/lib -lgts -L$(SC_DIR)/lib -L$(HDF5_DIR)/lib -lsc -lm -lglib-2.0 -lhdf5_cpp -lhdf5 -L$(MESQUITE_DIR)/lib -lmesquite $(CGAL_LIB)
CXX_FLAGS = -I$(GTS_DIR)/include -I$(SC_DIR)/include $(GLIB_INCLUDE) -I./include -I$(HDF5_DIR)/include -I$(MESQUITE_DIR)/include

#CXX       = mpicxx -O2
#LDFLAGS   = -L/home/lucio/local -lsc -lgts -lm
#CXX_FLAGS = -I/home/camata/local/include -I/usr/include/glib-2.0 -I/usr/lib/x86_64-linux-gnu/glib-2.0/include -Iinclude


## Target: all
all: hexmesh


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

