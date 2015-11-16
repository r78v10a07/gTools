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
CND_PLATFORM=GNU-MacOSX
CND_DLIB_EXT=dylib
CND_CONF=Release
CND_DISTDIR=dist
CND_BUILDDIR=build

# Include project Makefile
include Makefile

# Object Directory
OBJECTDIR=${CND_BUILDDIR}/${CND_CONF}/${CND_PLATFORM}

# Object Files
OBJECTFILES= \
	${OBJECTDIR}/src/berror.o \
	${OBJECTDIR}/src/bmemory.o \
	${OBJECTDIR}/src/bstring.o \
	${OBJECTDIR}/src/btree.o \
	${OBJECTDIR}/src/chrentity.o \
	${OBJECTDIR}/src/chromosome.o \
	${OBJECTDIR}/src/gene.o \
	${OBJECTDIR}/src/sam.o


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
LDLIBSOPTIONS=

# Build Targets
.build-conf: ${BUILD_SUBPROJECTS}
	"${MAKE}"  -f nbproject/Makefile-${CND_CONF}.mk ${CND_DISTDIR}/${CND_CONF}/libgtools.a

${CND_DISTDIR}/${CND_CONF}/libgtools.a: ${OBJECTFILES}
	${MKDIR} -p ${CND_DISTDIR}/${CND_CONF}
	${RM} ${CND_DISTDIR}/${CND_CONF}/libgtools.a
	${AR} -rv ${CND_DISTDIR}/${CND_CONF}/libgtools.a ${OBJECTFILES} 
	$(RANLIB) ${CND_DISTDIR}/${CND_CONF}/libgtools.a

${OBJECTDIR}/src/berror.o: src/berror.c 
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} "$@.d"
	$(COMPILE.c) -O2 -Iinclude -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/berror.o src/berror.c

${OBJECTDIR}/src/bmemory.o: src/bmemory.c 
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} "$@.d"
	$(COMPILE.c) -O2 -Iinclude -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/bmemory.o src/bmemory.c

${OBJECTDIR}/src/bstring.o: src/bstring.c 
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} "$@.d"
	$(COMPILE.c) -O2 -Iinclude -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/bstring.o src/bstring.c

${OBJECTDIR}/src/btree.o: src/btree.c 
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} "$@.d"
	$(COMPILE.c) -O2 -Iinclude -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/btree.o src/btree.c

${OBJECTDIR}/src/chrentity.o: src/chrentity.c 
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} "$@.d"
	$(COMPILE.c) -O2 -Iinclude -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/chrentity.o src/chrentity.c

${OBJECTDIR}/src/chromosome.o: src/chromosome.c 
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} "$@.d"
	$(COMPILE.c) -O2 -Iinclude -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/chromosome.o src/chromosome.c

${OBJECTDIR}/src/gene.o: src/gene.c 
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} "$@.d"
	$(COMPILE.c) -O2 -Iinclude -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/gene.o src/gene.c

${OBJECTDIR}/src/sam.o: src/sam.c 
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} "$@.d"
	$(COMPILE.c) -O2 -Iinclude -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/sam.o src/sam.c

# Subprojects
.build-subprojects:

# Clean Targets
.clean-conf: ${CLEAN_SUBPROJECTS}
	${RM} -r ${CND_BUILDDIR}/${CND_CONF}
	${RM} ${CND_DISTDIR}/${CND_CONF}/libgtools.a

# Subprojects
.clean-subprojects:

# Enable dependency checking
.dep.inc: .depcheck-impl

include .dep.inc
