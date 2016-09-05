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
	# FFLAGS:= $(FFLAGS) -openmp -static-intel -fpp -openmp-link=static  -D $(OSFLAG) -D CLUSTER=$(CLUSTER)
	FFLAGS:= $(FFLAGS) -fpp -D $(OSFLAG)

	obj:= .o
	exe:= 

	DEL:= rm -rf
endif

all: executable

debug: FFLAGS = -DDEBUG=${DEBUG} -g -O0 -openmp -check bounds -fpp -traceback -D $(OSFLAG)
debug: executable

ifort  -openmp -static-intel -fpp -openmp-link=static -O3 -m64 -o alphaphase Constants.f90 Parameters.f90 Random.f90 PedigreeDefinition.f90 NRMcode.f90 CoreDefinition.f90 CoreSubsetDefinition.f90 Clustering.f90 HaplotypeLibrary.f90 SurrogateDefinition.f90 Phasing.f90 MemberManagerDefinition.f90 InputOutput.f90 Testing.f90 AlphaPhase.f90

OBJS:=Clustering$(obj) Constants$(obj) CoreDefinition$(obj) CoreSubsetDefinition$(obj) HaplotypeLibrary$(obj) InputOutput$(obj) MemberManagerDefinition$(obj) NRMcode$(obj) Parameters$(obj) PedigreeDefinition$(obj) Phasing$(obj) Random$(obj) SurrogateDefinition$(obj) Testing$(obj)

%$(obj):%.f90
	$(FC) $(FFLAGS) -c $<

executable: $(OBJS)
	$(FC) AlphaPhase.f90 $(OBJS) $(FFLAGS) -o $(PROGRAM)$(exe)
	
clean:
	$(DEL) *$(obj) *.mod *~

veryclean:
	$(DEL) *$(obj)*.mod *~ AlphaImpute$(exe)

.PHONY: make clean veryclean all

