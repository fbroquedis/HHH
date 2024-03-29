%.o: %.f90
	$(F90) $(F90FLAGS1) -c  $< -o $@  $(LIBS)
%.o: %.F90
	$(F90) $(F90FLAGS1) -c $< -o $@  $(LIBS)
%.out: %.f90
	$(F90) $(F90FLAGS2) -o $@ $^  $(LIBS)
%.out: %.F90
	$(F90) $(F90FLAGS2) -o $@ $^  $(LIBS)
!LIBS=	-llapack -lblas
#LIBS=	-ldxml
#LIBS=/usr/lib/liblapack.a /usr/lib/libblas.a 
#LIBS=/usr/lib/liblapack.so.3.1.1 /usr/lib/libblas.so.3.1.1
#LIBS= /home/jdiaz1/FORTRAN90/MUMPS_4.9.2/lib/libzmumps.a -I/home/jdiaz1/FORTRAN90/MUMPS_4.9.2/lib/include/ /home/jdiaz1/FORTRAN90/MUMPS_4.9.2/lib/libmumps_common.a  -L/home/jdiaz1/FORTRAN90/MUMPS_4.9.2/PORD/lib/ -lpord  /usr/lib/libscalapack-openmpi.so.1.8.0  /usr/lib/libblacs-openmpi.so.1.1   /usr/lib/libblacsF77init-openmpi.so.1.1   /usr/lib/libblas/libblas.so.3gf  /usr/lib/libatlas.so.3gf /usr/lib/lapack/liblapack.so.3gf   -lpthread 

#LIBS=/usr/lib/libzmumps-4.9.2.so /usr/lib/libmumps_common-4.9.2.so      /usr/lib/libscalapack-openmpi.so.1.8.0  /usr/lib/libblacs-openmpi.so.1.1   /usr/lib/libblacsF77init-openmpi.so.1.1    /usr/lib/libblas/libblas.so.3gf  /usr/lib/libatlas.so.3gf  /usr/lib/lapack/liblapack.so.3gf   /usr/lib/libmpi_f77.so.0 -lpthread
#LIBS=/usr/lib/libzmumps-4.9.2.so /usr/lib/libmumps_common-4.9.2.so      /usr/lib/libscalapack-openmpi.so.1.8.0  /usr/lib/libblacs-openmpi.so.1.1   /usr/lib/libblacsF77init-openmpi.so.1.1    /usr/lib/libblas/libblas.so.3gf  /usr/lib/libatlas.so.3gf  /usr/lib/lapack/liblapack.so.3gf   /usr/lib/libmpi_f77.so.0 -lpthread
LIBS=/usr/lib/libblas/libblas.so.3gf /usr/lib/lapack/liblapack.so.3gf  -lpthread
#LIBS=/afs/inria.fr/rocq/home/ondes/jdiaz/FORTRAN90/CAGNIARD/flusol/lapack.a  /afs/inria.fr/rocq/home/ondes/jdiaz/FORTRAN90/CAGNIARD/flusol/blas.a


!F90 =/appli_MTS/TMA/Intel_FC/l_fce_c_9.1.043/bin/ifort
F90 = gfortran
!F90FLAGS1 =    -assume byterecl   -Iincludefic/ -Iincludemod/
!F90FLAGS2 =     -assume byterecl -Iincludefic/ -Iincludemod/ 
!F90FLAGS1 =   -O4 -march=k8 -m64 -Iincludefic/ -Iincludemod/ -fbounds-check
!F90FLAGS2 =   -O4   -march=k8 -m64 -Iincludefic/ -Iincludemod/ 
F90FLAGS1 =  -ftree-vectorize -O2 -msse4.2 -ftree-vectorizer-verbose=1 -ftree-vectorizer-verbose=2 -ftree-vectorizer-verbose=3 -ftree-vectorizer-verbose=4 -ftree-vectorizer-verbose=5 -ftree-vectorizer-verbose=6 -ftree-vectorizer-verbose=7 -fopenmp  -Iincludefic/ -Iincludemod/ -Dintel_ -DALLOW_NON_INIT  -I/usr/include/mpich2/ -I.  -fdefault-real-8 -fdefault-double-8 
F90FLAGS2 =  -ftree-vectorize -O2 -msse4.2 -ftree-vectorizer-verbose=1 -ftree-vectorizer-verbose=2 -ftree-vectorizer-verbose=3 -ftree-vectorizer-verbose=4 -ftree-vectorizer-verbose=5 -ftree-vectorizer-verbose=6 -ftree-vectorizer-verbose=7 -fopenmp -Iincludefic/ -Iincludemod/ -Dintel_ -DALLOW_NON_INIT  -I/usr/include/mpich2/ -I. -I/usr/include/  -fdefault-real-8 -fdefault-double-8 
!F90FLAGS1 =  -march='corei7' -ftree-vectorize -O2 -msse4.2 -ftree-vectorizer-verbose=1
!F90FLAGS2 =  -march='corei7' -ftree-vectorize -O2 -msse4.2 -ftree-vectorizer-verbose=1



RM = \rm -f

#F90FLAGS1 =  -ggdb F90FLAGS2 =  -lstdc++

#############
# VARIABLES #
#############
# 1. LIBRAIRIES INDEPENDANTES
# Variables necessaires a la creation de la librairie des Mat
CHEMAT = lib/libmat
LSTMAT = $(foreach DIR,$(CHEMAT),$(wildcard $(DIR)/*.[Ff]90))
OBJMAT = $(LSTMAT:.F90=.o)
LIBMAT = lib/libmat/libmat.a
# Variables necessaires a la creation de la librairie des Mesh
CHEMESH = lib/libmesh
LSTMESH = $(foreach DIR,$(CHEMESH),$(wildcard $(DIR)/*.[Ff]90))
OBJMESH = $(LSTMESH:.F90=.o)
LIBMESH = lib/libmesh/libmesh.a
# Variables necessaires a la creation de la librairie des phiI
CHEBAS = lib/libbasisfunc
LSTBAS = $(foreach DIR,$(CHEBAS),$(wildcard $(DIR)/*.[Ff]90))
OBJBAS = $(LSTBAS:.F90=.o)
LIBBAS = lib/libbasisfunc/libbasisfunc.a
# Variables necessaires a la creation de la librairie des output
CHEOUT = lib/liboutput
LSTOUT = $(foreach DIR,$(CHEOUT),$(wildcard $(DIR)/*.[Ff]90))
OBJOUT = $(LSTOUT:.F90=.o)
LIBOUT = lib/liboutput/liboutput.a
# Variables necessaires a la creation de la librairie des init
CHEINIT = lib/libinit
LSTINIT = $(foreach DIR,$(CHEINIT),$(wildcard $(DIR)/*.[Ff]90))
OBJINIT = $(LSTINIT:.F90=.o)
LIBINIT = lib/libinit/libinit.a
# 2. LIBRAIRIE DES MODULES DE DECLARATION DES VAR GLOBALES
# Variables necessaires a la creation de la librairie des mod
CHEMOD = mod
LSTMOD = $(foreach DIR,$(CHEMOD),$(wildcard $(DIR)/*.[Ff]90)) 
OBJMOD = $(LSTMOD:.f90=.o)
LIBMOD = mod/libmod.a
# 3. LIBRAIRIES DEPENDANT DES MODULES ET DES LIBRAIRIES 1.
# Variables necessaires a la creation de l'executable
CHEBIN = lib/bin
LSTBIN = $(foreach DIR,$(CHEBIN),$(wildcard $(DIR)/*.[Ff]90))
BINBIN = $(addsuffix .out,$(basename $(LSTBIN)))
BIN = $(subst $(CHEBIN),lib,$(BINBIN))

###########
# CHEMINS #
###########
vpath %.F90 $(CHEBIN)
vpath %.m includemod/
vpath %.a lib/libmat  lib/libmesh mod/  lib/liboutput lib/libbasisfunc lib/libinit 

###############
# COMPILATION #
###############
# 1. LIBRAIRIES INDEPENDANTES
srcmat: $(LIBMAT)
$(LIBMAT): $(LIBMAT)($(OBJMAT))
srcmesh: $(LIBMESH)
$(LIBMESH): $(LIBMESH)($(OBJMESH))
srcbas: $(LIBBAS)
$(LIBBAS): $(LIBBAS)($(OBJBAS))
srcout: $(LIBOUT)
$(LIBOUT): $(LIBOUT)($(OBJOUT))
srcinit: $(LIBINIT)
$(LIBINIT): $(LIBINIT)($(OBJINIT))
# Creation des sources 1
src1:  srcmat  srcmesh srcbas srcout srcinit
# 2. LIBRAIRIE DES MODULES DE DECLARATION DES VAR GLOBALES
# Creation de la librairie des mod
srcmod: $(LIBMOD)
$(LIBMOD): $(LIBMOD)($(OBJMOD))
	mv *.mod includemod/.
# 3. LIBRAIRIES DEPENDANT DES MODULES ET DES LIBRAIRIES 1.
# Creation de la librairie libvarpbs.a
# Creation des sources 3
# 5. PROGRAMME PRINCIPAL
# Creation de l'executable du programme principal
all:  $(LIBMOD) $(LIBMAT) $(LIBMESH) $(LIBBAS)  $(LIBOUT) $(LIBINIT) $(BIN)
$(BIN): $(BINBIN)
	mv $(CHEBIN)/$(@F) .
$(BINBIN): -lmod -lmat -lmesh -lbasisfunc -loutput -linit


###############
# SUPPRESSION #
###############
# 1. LIBRAIRIES INDEPENDANTES
# Supression de la librairie des Tab
clean-srcmat:	
	$(RM) lib/libmat/*.a
clean-srcmesh:	
	$(RM) lib/libmesh/*.a
clean-srcbas:	
	$(RM) lib/libbasisfunc/*.a
clean-srcout:	
	$(RM) lib/liboutput/*.a
clean-srcinit:	
	$(RM) lib/libinit/*.a
# Suppression des librairies 1
clean-src1: clean-srcmat  clean-srcmesh  clean-srcbas clean-srcout clean-srcinit
# 2. LIBRAIRIE DES MODULES DE DECLARATION DES VAR GLOBALES
# Suppression de la librairie des mod
clean-srcmod:	
	$(RM) mod/*.a
	$(RM) mod/*.o
	$(RM) includemod/*.mod
# 3. LIBRAIRIES DEPENDANT DES MODULES ET DES LIBRAIRIES 1.
# Suppression de l'executable
clean-out:
	$(RM) *.out

# Suppression de tout
clean: clean-src1 clean-srcmod  clean-out
	$(RM) lib/*/*.o





