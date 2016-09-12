# General variables
NAME:=AlphaPhase
VERSION:=$(shell git rev-parse --short HEAD)
MASTERVERSION:=$(shell git describe --tag)
PROGRAM:=$(NAME)$(MASTERVERSION)

# Compiler
FC:=ifort
#FC:=gfortran

# Options
FFLAGS:=-O3 -m64 -DVERS=""commit-$(VERSION)""

# Set precompilation options by default
DEBUG?=0

# MS Windows
ifeq ($(OS), Windows_NT)
	OSFLAG := "OS_WIN"
	# FFLAGS := $(FFLAGS) /static /fpp  /D $(OSFLAG) -D CLUSTER=0 /Qopenmp /libs:static
	FFLAGS := $(FFLAGS) /fpp  /D $(OSFLAG)

	obj:= .obj
	exe:= .exe

	DEL:= del
else
	OSFLAG := "OS_UNIX"
	# FFLAGS:= $(FFLAGS) -qopenmp -static-intel -fpp -qopenmp-link=static  -D $(OSFLAG) -D CLUSTER=$(CLUSTER)
	FFLAGS:= $(FFLAGS) -fpp -D $(OSFLAG)

	obj:= .o
	exe:=

	DEL:= rm -rf
endif

all: executable

debug: FFLAGS = -DDEBUG=${DEBUG} -g -O0 -check format -check bounds -fpp -traceback -D $(OSFLAG)
debug: executable

OBJS:= Constants$(obj) \
	   ParametersDefinition$(obj) \
	   Random$(obj) \
	   Sorting$(obj) \
	   PedigreeDefinition$(obj) \
	   NRMcode$(obj) \
	   CoreDefinition$(obj) \
	   CoreSubsetDefinition$(obj) \
	   Clustering$(obj) \
	   HaplotypeLibraryDefinition$(obj) \
	   SurrogateDefinition$(obj) \
	   LongRangePhasing$(obj) \
	   MemberManagerDefinition$(obj) \
	   TestResultDefinition$(obj) \
	   InputOutput$(obj) \
	   HaplotypeLibraryPhasing$(obj)

%$(obj):%.f90
	$(FC) $(FFLAGS) -c $<

executable: $(OBJS)
	$(FC) AlphaPhase.f90 $(OBJS) $(FFLAGS) -o $(PROGRAM)$(exe)

clean:
	$(DEL) *$(obj) *.mod *~

veryclean:
	$(DEL) *$(obj)*.mod *~ AlphaImpute$(exe)

.PHONY: make clean veryclean all

