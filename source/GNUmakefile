# Fortran compiler and compiler flags
FXX      := gfortran-10
FXXFLAGS := -Wall -std=legacy -O3


# Location of build directory
#BUILD_DIR := /build/


# Location of installation directory
#INSTALL_DIR := /lib/


# Library files to build
#all: libmathhole.a
all: blas_obj libmathhole.a

libmathhole.a: $(blas_obj)
	ar rcs libmathhole.a *.o

# ar rcs libmathhole.a *.o




blas_obj: mathhole_saxpy.o mathhole_snrm2.o mathhole_sgemv.o mathhole_sgemm.o mathhole_strsv.o \
          mathhole_daxpy.o mathhole_dnrm2.o mathhole_dgemv.o mathhole_dgemm.o mathhole_dtrsv.o \
          mathhole_lsame.o mathhole_xerbla.o


# Compile library source files into object files

mathhole_saxpy.o: mathhole_saxpy.f90
	$(FXX) $(FXXFLAGS) -c mathhole_saxpy.f90

mathhole_snrm2.o: mathhole_snrm2.f90
	$(FXX) $(FXXFLAGS) -c mathhole_snrm2.f90

mathhole_sgemv.o: mathhole_sgemv.f90
	$(FXX) $(FXXFLAGS) -c mathhole_sgemv.f90

mathhole_sgemm.o: mathhole_sgemm.f90
	$(FXX) $(FXXFLAGS) -c mathhole_sgemm.f90

mathhole_strsv.o: mathhole_strsv.f90
	$(FXX) $(FXXFLAGS) -c mathhole_strsv.f90

mathhole_daxpy.o: mathhole_daxpy.f90
	$(FXX) $(FXXFLAGS) -c mathhole_daxpy.f90

mathhole_dnrm2.o: mathhole_dnrm2.f90
	$(FXX) $(FXXFLAGS) -c mathhole_dnrm2.f90

mathhole_dgemv.o: mathhole_dgemv.f90
	$(FXX) $(FXXFLAGS) -c mathhole_dgemv.f90

mathhole_dgemm.o: mathhole_dgemm.f90
	$(FXX) $(FXXFLAGS) -c mathhole_dgemm.f90

mathhole_dtrsv.o: mathhole_dtrsv.f90
	$(FXX) $(FXXFLAGS) -c mathhole_dtrsv.f90

mathhole_lsame.o: mathhole_lsame.f90
	$(FXX) $(FXXFLAGS) -c mathhole_lsame.f90

mathhole_xerbla.o: mathhole_xerbla.f90
	$(FXX) $(FXXFLAGS) -c mathhole_xerbla.f90


make install:
	mv libmathhole.a ../lib && rm *.o


# Clean build
make clean:
		rm *.o *.a && clear
