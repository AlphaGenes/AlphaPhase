rm alphaphase
rm *.mod
ifort -static-intel -O3 -m64 -o alphaphase Constants.f90 Parameters.f90 Random.f90 PedigreeDefinition.f90 NRMcode.f90 CoreDefinition.f90 CoreSubsetDefinition.f90 Clustering.f90 HaplotypeLibraryDefinition.f90 SurrogateDefinition.f90 LongRangePhasing.f90 MemberManagerDefinition.f90 InputOutput.f90 Testing.f90 HaplotypeLibraryPhasing.f90 AlphaPhase.f90
