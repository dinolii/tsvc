CC=gcc -std=c99
CXX=g++
FC=gfortran
flags = -mavx2 -O3 -fstrict-aliasing -fivopts
vecflags = -ftree-vectorize 
novecflags = -fno-tree-vectorize
omp_flags=-fopenmp

ifdef VEC_REPORT
vecflags+=-ftree-vectorizer-verbose=5 -fopt-info-vec -fopt-info-vec-missed -fno-inline
endif

ieee_math_flags=
fast_math_flags=-ffast-math

ifdef PRECISE_MATH
$(warning No 'precise' math flags for GNU!)
endif
