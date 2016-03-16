ifort  -openmp -static-intel -fpp -openmp-link=static -O3 -m64 -o alphaphase Random.f90 DataSubset.f90 HaplotypeLibrary.f90 Phasing.f90 SurrogateDefinition.f90 AlphaPhase.f90
