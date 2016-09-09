rm alphaphase
rm *.mod
ifort -fpp -static-intel -O3 -m64 -o alphaphase -DCOMMIT=" e08d861" Constants.f90 ParametersDefinition.f90 Random.f90 Sorting.f90 PedigreeDefinition.f90 NRMcode.f90 CoreDefinition.f90 CoreSubsetDefinition.f90 Clustering.f90 HaplotypeLibraryDefinition.f90 SurrogateDefinition.f90 LongRangePhasing.f90 MemberManagerDefinition.f90 TestResultDefinition.f90 InputOutput.f90 HaplotypeLibraryPhasing.f90 AlphaPhase.f90
