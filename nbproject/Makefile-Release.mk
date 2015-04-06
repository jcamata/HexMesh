#
# Generated Makefile - do not edit!
#
# Edit the Makefile in the project folder instead (../Makefile). Each target
# has a -pre and a -post target defined where you can add customized code.
#
# This makefile implements configuration specific macros and targets.


# Environment
MKDIR=mkdir
CP=cp
GREP=grep
NM=nm
CCADMIN=CCadmin
RANLIB=ranlib
CC=gcc
CCC=g++
CXX=g++
FC=gfortran
AS=as

# Macros
CND_PLATFORM=GNU-Linux-x86
CND_DLIB_EXT=so
CND_CONF=Release
CND_DISTDIR=dist
CND_BUILDDIR=build

# Include project Makefile
include Makefile

# Object Directory
OBJECTDIR=${CND_BUILDDIR}/${CND_CONF}/${CND_PLATFORM}

# Object Files
OBJECTFILES= \
	${OBJECTDIR}/_ext/1217637752/hexa.o \
	${OBJECTDIR}/main.o \
	${OBJECTDIR}/src/hexa_face_hanging.o \
	${OBJECTDIR}/src/hexa_mesh.o \
	${OBJECTDIR}/src/hexa_parallel.o \
	${OBJECTDIR}/src/hexa_unv.o \
	${OBJECTDIR}/src/hexa_vtk.o


# C Compiler Flags
CFLAGS=

# CC Compiler Flags
CCFLAGS=
CXXFLAGS=

# Fortran Compiler Flags
FFLAGS=

# Assembler Flags
ASFLAGS=

# Link Libraries and Options
LDLIBSOPTIONS=-L/home/camata/local/lib -lsc

# Build Targets
.build-conf: ${BUILD_SUBPROJECTS}
	"${MAKE}"  -f nbproject/Makefile-${CND_CONF}.mk ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/hexmesh

${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/hexmesh: ${OBJECTFILES}
	${MKDIR} -p ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}
	${LINK.cc} -o ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/hexmesh ${OBJECTFILES} ${LDLIBSOPTIONS}

${OBJECTDIR}/_ext/1217637752/hexa.o: /home/camata/Programming/HexMesh/src/hexa.cpp 
	${MKDIR} -p ${OBJECTDIR}/_ext/1217637752
	${RM} "$@.d"
	$(COMPILE.cc) -O2 -I/home/camata/local/include -Iinclude -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/_ext/1217637752/hexa.o /home/camata/Programming/HexMesh/src/hexa.cpp

${OBJECTDIR}/main.o: main.cpp 
	${MKDIR} -p ${OBJECTDIR}
	${RM} "$@.d"
	$(COMPILE.cc) -O2 -I/home/camata/local/include -Iinclude -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/main.o main.cpp

${OBJECTDIR}/src/hexa_face_hanging.o: src/hexa_face_hanging.cpp 
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} "$@.d"
	$(COMPILE.cc) -O2 -I/home/camata/local/include -Iinclude -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/hexa_face_hanging.o src/hexa_face_hanging.cpp

${OBJECTDIR}/src/hexa_mesh.o: src/hexa_mesh.cpp 
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} "$@.d"
	$(COMPILE.cc) -O2 -I/home/camata/local/include -Iinclude -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/hexa_mesh.o src/hexa_mesh.cpp

${OBJECTDIR}/src/hexa_parallel.o: src/hexa_parallel.cpp 
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} "$@.d"
	$(COMPILE.cc) -O2 -I/home/camata/local/include -Iinclude -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/hexa_parallel.o src/hexa_parallel.cpp

${OBJECTDIR}/src/hexa_unv.o: src/hexa_unv.cpp 
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} "$@.d"
	$(COMPILE.cc) -O2 -I/home/camata/local/include -Iinclude -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/hexa_unv.o src/hexa_unv.cpp

${OBJECTDIR}/src/hexa_vtk.o: src/hexa_vtk.cpp 
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} "$@.d"
	$(COMPILE.cc) -O2 -I/home/camata/local/include -Iinclude -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/hexa_vtk.o src/hexa_vtk.cpp

# Subprojects
.build-subprojects:

# Clean Targets
.clean-conf: ${CLEAN_SUBPROJECTS}
	${RM} -r ${CND_BUILDDIR}/${CND_CONF}
	${RM} ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/hexmesh

# Subprojects
.clean-subprojects:

# Enable dependency checking
.dep.inc: .depcheck-impl

include .dep.inc
