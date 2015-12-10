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
	${OBJECTDIR}/src/reads.o \
	${OBJECTDIR}/src/sam.o

# Test Directory
TESTDIR=${CND_BUILDDIR}/${CND_CONF}/${CND_PLATFORM}/tests

# Test Files
TESTFILES= \
	${TESTDIR}/TestFiles/f1

# Test Object Files
TESTOBJECTFILES= \
	${TESTDIR}/tests/SAMTest.o

# C Compiler Flags
CFLAGS=-g

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

${OBJECTDIR}/src/reads.o: src/reads.c 
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} "$@.d"
	$(COMPILE.c) -O2 -Iinclude -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/reads.o src/reads.c

${OBJECTDIR}/src/sam.o: src/sam.c 
	${MKDIR} -p ${OBJECTDIR}/src
	${RM} "$@.d"
	$(COMPILE.c) -O2 -Iinclude -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/sam.o src/sam.c

# Subprojects
.build-subprojects:

# Build Test Targets
.build-tests-conf: .build-tests-subprojects .build-conf ${TESTFILES}
.build-tests-subprojects:

${TESTDIR}/TestFiles/f1: ${TESTDIR}/tests/SAMTest.o ${OBJECTFILES:%.o=%_nomain.o}
	${MKDIR} -p ${TESTDIR}/TestFiles
	${LINK.c}   -o ${TESTDIR}/TestFiles/f1 $^ ${LDLIBSOPTIONS} -L../../../../../.homebrew/lib -lcunit 


${TESTDIR}/tests/SAMTest.o: tests/SAMTest.c 
	${MKDIR} -p ${TESTDIR}/tests
	${RM} "$@.d"
	$(COMPILE.c) -O2 -Iinclude -Iinclude -I../../../../../.homebrew/include -MMD -MP -MF "$@.d" -o ${TESTDIR}/tests/SAMTest.o tests/SAMTest.c


${OBJECTDIR}/src/berror_nomain.o: ${OBJECTDIR}/src/berror.o src/berror.c 
	${MKDIR} -p ${OBJECTDIR}/src
	@NMOUTPUT=`${NM} ${OBJECTDIR}/src/berror.o`; \
	if (echo "$$NMOUTPUT" | ${GREP} '|main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T _main$$'); \
	then  \
	    ${RM} "$@.d";\
	    $(COMPILE.c) -O2 -Iinclude -Dmain=__nomain -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/berror_nomain.o src/berror.c;\
	else  \
	    ${CP} ${OBJECTDIR}/src/berror.o ${OBJECTDIR}/src/berror_nomain.o;\
	fi

${OBJECTDIR}/src/bmemory_nomain.o: ${OBJECTDIR}/src/bmemory.o src/bmemory.c 
	${MKDIR} -p ${OBJECTDIR}/src
	@NMOUTPUT=`${NM} ${OBJECTDIR}/src/bmemory.o`; \
	if (echo "$$NMOUTPUT" | ${GREP} '|main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T _main$$'); \
	then  \
	    ${RM} "$@.d";\
	    $(COMPILE.c) -O2 -Iinclude -Dmain=__nomain -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/bmemory_nomain.o src/bmemory.c;\
	else  \
	    ${CP} ${OBJECTDIR}/src/bmemory.o ${OBJECTDIR}/src/bmemory_nomain.o;\
	fi

${OBJECTDIR}/src/bstring_nomain.o: ${OBJECTDIR}/src/bstring.o src/bstring.c 
	${MKDIR} -p ${OBJECTDIR}/src
	@NMOUTPUT=`${NM} ${OBJECTDIR}/src/bstring.o`; \
	if (echo "$$NMOUTPUT" | ${GREP} '|main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T _main$$'); \
	then  \
	    ${RM} "$@.d";\
	    $(COMPILE.c) -O2 -Iinclude -Dmain=__nomain -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/bstring_nomain.o src/bstring.c;\
	else  \
	    ${CP} ${OBJECTDIR}/src/bstring.o ${OBJECTDIR}/src/bstring_nomain.o;\
	fi

${OBJECTDIR}/src/btree_nomain.o: ${OBJECTDIR}/src/btree.o src/btree.c 
	${MKDIR} -p ${OBJECTDIR}/src
	@NMOUTPUT=`${NM} ${OBJECTDIR}/src/btree.o`; \
	if (echo "$$NMOUTPUT" | ${GREP} '|main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T _main$$'); \
	then  \
	    ${RM} "$@.d";\
	    $(COMPILE.c) -O2 -Iinclude -Dmain=__nomain -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/btree_nomain.o src/btree.c;\
	else  \
	    ${CP} ${OBJECTDIR}/src/btree.o ${OBJECTDIR}/src/btree_nomain.o;\
	fi

${OBJECTDIR}/src/chrentity_nomain.o: ${OBJECTDIR}/src/chrentity.o src/chrentity.c 
	${MKDIR} -p ${OBJECTDIR}/src
	@NMOUTPUT=`${NM} ${OBJECTDIR}/src/chrentity.o`; \
	if (echo "$$NMOUTPUT" | ${GREP} '|main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T _main$$'); \
	then  \
	    ${RM} "$@.d";\
	    $(COMPILE.c) -O2 -Iinclude -Dmain=__nomain -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/chrentity_nomain.o src/chrentity.c;\
	else  \
	    ${CP} ${OBJECTDIR}/src/chrentity.o ${OBJECTDIR}/src/chrentity_nomain.o;\
	fi

${OBJECTDIR}/src/chromosome_nomain.o: ${OBJECTDIR}/src/chromosome.o src/chromosome.c 
	${MKDIR} -p ${OBJECTDIR}/src
	@NMOUTPUT=`${NM} ${OBJECTDIR}/src/chromosome.o`; \
	if (echo "$$NMOUTPUT" | ${GREP} '|main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T _main$$'); \
	then  \
	    ${RM} "$@.d";\
	    $(COMPILE.c) -O2 -Iinclude -Dmain=__nomain -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/chromosome_nomain.o src/chromosome.c;\
	else  \
	    ${CP} ${OBJECTDIR}/src/chromosome.o ${OBJECTDIR}/src/chromosome_nomain.o;\
	fi

${OBJECTDIR}/src/gene_nomain.o: ${OBJECTDIR}/src/gene.o src/gene.c 
	${MKDIR} -p ${OBJECTDIR}/src
	@NMOUTPUT=`${NM} ${OBJECTDIR}/src/gene.o`; \
	if (echo "$$NMOUTPUT" | ${GREP} '|main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T _main$$'); \
	then  \
	    ${RM} "$@.d";\
	    $(COMPILE.c) -O2 -Iinclude -Dmain=__nomain -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/gene_nomain.o src/gene.c;\
	else  \
	    ${CP} ${OBJECTDIR}/src/gene.o ${OBJECTDIR}/src/gene_nomain.o;\
	fi

${OBJECTDIR}/src/reads_nomain.o: ${OBJECTDIR}/src/reads.o src/reads.c 
	${MKDIR} -p ${OBJECTDIR}/src
	@NMOUTPUT=`${NM} ${OBJECTDIR}/src/reads.o`; \
	if (echo "$$NMOUTPUT" | ${GREP} '|main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T _main$$'); \
	then  \
	    ${RM} "$@.d";\
	    $(COMPILE.c) -O2 -Iinclude -Dmain=__nomain -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/reads_nomain.o src/reads.c;\
	else  \
	    ${CP} ${OBJECTDIR}/src/reads.o ${OBJECTDIR}/src/reads_nomain.o;\
	fi

${OBJECTDIR}/src/sam_nomain.o: ${OBJECTDIR}/src/sam.o src/sam.c 
	${MKDIR} -p ${OBJECTDIR}/src
	@NMOUTPUT=`${NM} ${OBJECTDIR}/src/sam.o`; \
	if (echo "$$NMOUTPUT" | ${GREP} '|main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T main$$') || \
	   (echo "$$NMOUTPUT" | ${GREP} 'T _main$$'); \
	then  \
	    ${RM} "$@.d";\
	    $(COMPILE.c) -O2 -Iinclude -Dmain=__nomain -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/src/sam_nomain.o src/sam.c;\
	else  \
	    ${CP} ${OBJECTDIR}/src/sam.o ${OBJECTDIR}/src/sam_nomain.o;\
	fi

# Run Test Targets
.test-conf:
	@if [ "${TEST}" = "" ]; \
	then  \
	    ${TESTDIR}/TestFiles/f1 || true; \
	else  \
	    ./${TEST} || true; \
	fi

# Clean Targets
.clean-conf: ${CLEAN_SUBPROJECTS}
	${RM} -r ${CND_BUILDDIR}/${CND_CONF}
	${RM} ${CND_DISTDIR}/${CND_CONF}/libgtools.a

# Subprojects
.clean-subprojects:

# Enable dependency checking
.dep.inc: .depcheck-impl

include .dep.inc
