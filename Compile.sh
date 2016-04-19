ifort -syntax-only AlphaPhase.f90 2> /dev/null
ifort -static-intel -fpp -O3 -m64 -o AlphaPhase1.1 Random.f90 CoreDefinition.f90 DataSubset.f90 SurrogateDefinition.f90 HaplotypeLibrary.f90 Phasing.f90 NRMcode.f90 AlphaPhase.f90
